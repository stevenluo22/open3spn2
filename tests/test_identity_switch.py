"""Tests for the switch forces of the Gibbs identity sampler.

``add_forces(tables, system)`` (or ``tables.add_forces(system)``) adds the DNA forces of the duplex of the identity
tables to an OpenMM system so that the identities of the positions can change during a simulation: the native forces
without the interactions of the tables, and for each kind a switch force with them, in the force group of the kind.
``set_state(tables, context, state)`` (or ``tables.set_state(context, state)``) selects the identities, a name index
per position, in a context of that system. The behavioural requirements tested here:

* S1. For every DNA type, with positions inside the duplex and at its ends, at the wild type and at Watson-Crick
  states set one after the other in one context, the energy of every force group and the forces equal those of the
  native system of the mutant duplex with those identities (``DNA.create_mutant``).
* S2. At states with mixes, the energy of every force group relative to the wild type changes as
  ``tables.state_energies`` of the kinds of that group.
* S3. Setting a state only sets global parameters of the context (in S1, S2 and S5 too), and with the previous state
  given only those of the positions that change. A state of the wrong length or with a name outside the alphabet is
  refused.
* S4. The system holds the native forces and, for each kind with a non-empty table, one switch force: a
  CustomCompoundBondForce in the group of the kind with one bond per interaction of the table. No exclusion removes an
  interaction (OpenMM needs the same exclusions in every CustomNonbondedForce, and its GPU platforms cap the
  exclusions of a CustomHbondForce donor or acceptor): the CustomNonbondedForces keep the native exclusions and the
  CustomHbondForces only exclude donor-acceptor pairs that share a particle. On every platform that runs here, the
  switched system has the energies and forces of the native mutant system on the same platform.
* S5. With periodic boundary conditions, a duplex wrapped across the box has the energies and forces of the native
  mutant system with the same box.
"""
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import openmm
import openmm.unit as unit
import pytest

from open3SPN2 import DNA
from open3SPN2.alchemy import add_forces, alphabet, identity_tables, set_state
from open3SPN2.ff3SPN2 import forces as dna_forces

test_path = Path(__file__).parent
_outdir = Path(tempfile.mkdtemp(prefix='open3spn2_switch_'))

STRUCTURES = {'A': test_path / 'adna', 'B': test_path / 'bdna', 'B_curved': test_path / 'bdna_curv'}
GROUPS = {'Bond': 6, 'Angle': 7, 'Stacking': 8, 'Dihedral': 9, 'BasePair': 10, 'CrossStacking': 11, 'Exclusion': 12,
          'Electrostatics': 13}
SITE_RESIDUES = (12, 13, 16)
END_RESIDUES = (1, 2, 32)
ENERGY_ATOL = 1e-4  # kJ/mol
FORCE_ATOL = 1e-3  # kJ/mol/nm
# Platforms other than Reference differ from it by up to ~1e-3 kJ/mol and ~1e-6 of the largest force, natively too.
PLATFORM_ENERGY_ATOL = 1e-3  # kJ/mol
PLATFORM_FORCE_RTOL = 1e-5  # of the largest force
PLATFORMS = ['CPU', 'CUDA', 'OpenCL', 'HIP']


@lru_cache(maxsize=None)
def _dna(dna_type):
    dna = DNA.fromXYZ(str(STRUCTURES[dna_type] / 'in00_conf.xyz'), dna_type, template_from_X3DNA=False,
                      output_pdb=str(_outdir / f'{dna_type}.pdb'))
    dna.periodic = False
    return dna


def _sites(dna, residues=SITE_RESIDUES):
    """((chain A, residue), (chain B, paired residue)) of each strand-A residue of the antiparallel duplex."""
    atoms = dna.atoms
    first_a = atoms[atoms['chainID'] == 'A']['resSeq'].min()
    last_b = atoms[atoms['chainID'] == 'B']['resSeq'].max()
    return [(('A', r), ('B', int(last_b - (r - first_a)))) for r in residues]


@lru_cache(maxsize=None)
def _tables(dna_type, residues=SITE_RESIDUES):
    dna = _dna(dna_type)
    return identity_tables(dna, _sites(dna, residues), alphabet('wc'))


def _frames(dna_type):
    """The structure and a jitter of it (the excluded volume of the undistorted structure is zero), in nm."""
    x = _dna(dna_type).atoms[['x', 'y', 'z']].to_numpy(float) / 10
    return np.stack([x, x + np.random.default_rng(3).normal(scale=0.03, size=x.shape)])


def _mutant(tables, state):
    mutations = []
    for ((chain_a, residue_a), (chain_b, residue_b)), name in zip(tables.sites, state):
        base_a, base_b = tables.alphabet.identities[tables.alphabet.names[name]]
        mutations += [(chain_a, residue_a, base_a), (chain_b, residue_b, base_b)]
    return tables.dna.create_mutant(mutations)


def _system(n_particles, box=None):
    system = openmm.System()
    for _ in range(n_particles):
        system.addParticle(1.0)
    if box is not None:
        system.setDefaultPeriodicBoxVectors(*np.diag(box))
    return system


def _context(system, platform='Reference'):
    platform = openmm.Platform.getPlatformByName(platform)
    properties = {'Precision': 'double'} if 'Precision' in platform.getPropertyNames() else {}
    return openmm.Context(system, openmm.VerletIntegrator(0.001), platform, properties)


def _native_system(dna, box=None):
    """The OpenMM system of the open3SPN2 DNA forces of `dna`."""
    system = _system(len(dna.atoms), box)
    for force_class in dna_forces.values():
        force_class(dna).addForce(system)
    return system


def _switched_system(tables, box=None):
    system = _system(len(tables.dna.atoms), box)
    tables.add_forces(system)
    return system


def _evaluate(context, coordinates):
    """(energy of each force group in GROUPS (kJ/mol), forces (N, 3) in kJ/mol/nm) at `coordinates` (nm)."""
    context.setPositions(coordinates)
    energies = [context.getState(getEnergy=True, groups={g}).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole) for g in GROUPS.values()]
    forces = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer)
    return np.array(energies), np.asarray(forces)


def _assert_same(context, reference, frames, energy_atol=ENERGY_ATOL, force_rtol=None):
    for x in frames:
        energies, forces = _evaluate(context, x)
        expected_energies, expected_forces = _evaluate(reference, x)
        force_atol = FORCE_ATOL if force_rtol is None else force_rtol * np.abs(expected_forces).max()
        np.testing.assert_allclose(energies, expected_energies, rtol=0, atol=energy_atol)
        np.testing.assert_allclose(forces, expected_forces, rtol=0, atol=force_atol)


class _ParameterContext:
    """Stands for an OpenMM context that only allows setting global parameters: records them and passes them on to
    `context`, if given."""

    def __init__(self, context=None):
        self.context = context
        self.parameters = []

    def setParameter(self, name, value):
        self.parameters.append((name, value))
        if self.context is not None:
            self.context.setParameter(name, value)


@pytest.mark.parametrize('dna_type, residues', [(dna_type, SITE_RESIDUES) for dna_type in STRUCTURES] +
                         [('B_curved', END_RESIDUES)])
def test_switched_system_matches_mutants(dna_type, residues):
    """S1: one context of the switched system, set in turn to the wild type, to states that change every position
    and to random Watson-Crick states, has the group energies and forces of the native system of each mutant."""
    tables = _tables(dna_type, residues)
    frames = _frames(dna_type)
    wild = tables.wild_type
    rng = np.random.default_rng(5)
    states = [wild, (wild + 1) % 4, (wild + 2) % 4] + [rng.integers(0, 4, size=len(wild)) for _ in range(2)] + [wild]

    context = _context(_switched_system(tables))
    previous = None
    for state in states:
        tables.set_state(_ParameterContext(context), state, previous)
        previous = state
        _assert_same(context, _context(_native_system(_mutant(tables, state))), frames)


@pytest.mark.parametrize('dna_type', list(STRUCTURES))
def test_mixes_follow_the_state_energies(dna_type):
    """S2: at states with W, the energy of each force group relative to the wild type equals the change of
    ``state_energies`` of the kinds of that group (none for electrostatics)."""
    tables = _tables(dna_type)
    frames = _frames(dna_type)
    w, at, gc = (tables.alphabet.index(name) for name in ('W', 'AT', 'GC'))
    states = np.array([[w, w, w], [w, at, w], [gc, w, at], [at, gc, w]])
    expected = tables.state_energies(frames, np.vstack([tables.wild_type, states]))
    context = _context(_switched_system(tables))

    for s, state in enumerate(states):
        for f, x in enumerate(frames):
            tables.set_state(_ParameterContext(context), tables.wild_type)
            wild_energies, _ = _evaluate(context, x)
            tables.set_state(_ParameterContext(context), state)
            energies, _ = _evaluate(context, x)
            change = [expected[kind][f, s + 1] - expected[kind][f, 0] if kind in expected else 0.0
                      for kind in GROUPS]
            np.testing.assert_allclose(energies - wild_energies, change, rtol=0, atol=ENERGY_ATOL)


def test_set_state_sets_only_changed_parameters():
    """S3: a state sets one global parameter per position, only those of the changed positions when the previous
    state is given, to the values the real context then holds; malformed states are refused."""
    tables = _tables('B')
    wild = tables.wild_type
    state = wild.copy()
    state[1] = (wild[1] + 1) % 4

    context = _context(_switched_system(tables))
    recorder = _ParameterContext(context)
    set_state(tables, recorder, wild)
    assert len(recorder.parameters) == len(wild)
    recorder = _ParameterContext(context)
    set_state(tables, recorder, state, previous=wild)
    assert len(recorder.parameters) == 1
    (name, value), = recorder.parameters
    assert value == state[1] and context.getParameter(name) == state[1]

    for bad in [wild[:-1], np.append(wild, 0), np.where(np.arange(len(wild)) == 0, len(tables.alphabet), wild),
                np.where(np.arange(len(wild)) == 0, -1, wild)]:
        with pytest.raises(ValueError):
            set_state(tables, _ParameterContext(), bad)


def _exclusions(force):
    return {tuple(force.getExclusionParticles(i)) for i in range(force.getNumExclusions())}


def test_system_structure():
    """S4: the switched system has the native forces plus, for each kind with a non-empty table, one
    CustomCompoundBondForce in the kind's group with a bond per interaction of the table; the CustomNonbondedForces
    keep the native exclusions and every exclusion of a CustomHbondForce joins a donor and an acceptor that share a
    particle."""
    tables = _tables('B_curved')
    native = _native_system(tables.dna)
    system = _system(len(tables.dna.atoms))
    returned = add_forces(tables, system)
    assert set(returned) == set(GROUPS)

    def forces_of(system, group):
        forces = [f for f in system.getForces() if f.getForceGroup() == group]
        compound = sum(f.getNumBonds() for f in forces if isinstance(f, openmm.CustomCompoundBondForce))
        return sorted(type(f).__name__ for f in forces), compound

    switched = {kind: table['term'].nunique() for kind, table in tables.kinds.items() if len(table)}
    assert switched
    assert system.getNumForces() == native.getNumForces() + len(switched)
    for kind, group in GROUPS.items():
        names, bonds = forces_of(system, group)
        native_names, native_bonds = forces_of(native, group)
        expected = sorted(native_names + ['CustomCompoundBondForce']) if kind in switched else native_names
        assert names == expected
        assert bonds == native_bonds + switched.get(kind, 0)

    native_exclusions = {f.getForceGroup(): _exclusions(f) for f in native.getForces()
                         if isinstance(f, openmm.CustomNonbondedForce)}
    for force in system.getForces():
        if isinstance(force, openmm.CustomNonbondedForce):
            assert _exclusions(force) == native_exclusions[force.getForceGroup()]
        if isinstance(force, openmm.CustomHbondForce):
            for donor, acceptor in _exclusions(force):
                shared = set(force.getDonorParameters(donor)[:3]) & set(force.getAcceptorParameters(acceptor)[:3])
                assert shared - {-1}


@pytest.mark.parametrize('platform', PLATFORMS)
def test_platforms(platform):
    """S4: on every platform that can run a context here, the switched system set to states that change every
    position has the group energies and forces of the native mutant system on the same platform."""
    try:
        _context(_system(1), platform)
    except Exception as error:  # absent platform, or no usable device or driver
        pytest.skip(f'{platform} cannot run here: {error}')
    tables = _tables('B_curved')
    context = _context(_switched_system(tables), platform)
    for shift in (1, 2):
        state = (tables.wild_type + shift) % 4
        tables.set_state(_ParameterContext(context), state)
        reference = _context(_native_system(_mutant(tables, state)), platform)
        _assert_same(context, reference, _frames('B_curved'), PLATFORM_ENERGY_ATOL, PLATFORM_FORCE_RTOL)


def test_periodic_box():
    """S5: with periodic boundary conditions, a duplex wrapped across the box boundary, set to states one after the
    other, has the group energies and forces of the native mutant systems with the same box."""
    tables = _tables('B')
    box = np.array([13.7, 12.9, 14.3])
    frames = _frames('B')
    wrapped = (frames - frames.mean(axis=(0, 1)) + box) % box
    states = [(tables.wild_type + 3) % 4, (tables.wild_type + 1) % 4]
    tables.dna.periodic = True
    try:
        context = _context(_switched_system(tables, box))
        references = [_context(_native_system(_mutant(tables, state), box)) for state in states]
    finally:
        tables.dna.periodic = False
    for state, reference in zip(states, references):
        tables.set_state(_ParameterContext(context), state)
        _assert_same(context, reference, wrapped)
