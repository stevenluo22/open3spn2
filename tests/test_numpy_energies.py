"""Tests for the numpy energies of the DNA force classes.

Every force class can score its own interactions in numpy, from the same parameters it gives OpenMM. Its
``table`` is a DataFrame with one row per interaction (particle indices and parameters, read from the OpenMM
force); ``open3SPN2.force.backend.<kind>_energy(table, coordinates)`` returns the energy of each row, and
``force.energy(coordinates)`` evaluates the force's own table with its ``k`` multiplier. The behavioural
requirements tested here:

* N1. For every force class and every reference structure (A, B and curved B DNA), the per-interaction numpy
  energies of a force sum to the energy OpenMM computes for that force alone, for one frame and for a batch of
  frames, including the global ``k`` multiplier. There is one energy per row of ``force.table``, the backend
  function gives the same numbers, and neither path calls OpenMM.
* N2. With periodic boundary conditions, a duplex wrapped across the box boundary has the numpy energy OpenMM
  computes with the same orthorhombic box, which differs from the energy of the same coordinates without a box.
* N3. The box must match the force: a periodic force without a box, a non-periodic force with a box and a
  non-orthorhombic box are refused.
* N4. Evaluating an edited table: when the parameters of a few rows change, only their energies change, and the
  total changes as it does in OpenMM after the same change.
"""
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import openmm
import openmm.unit as unit
import pytest

from open3SPN2 import DNA
from open3SPN2.force import backend, dna as dna_forces
from open3SPN2.parser.xyz import parse_xyz

test_path = Path(__file__).parent
_outdir = Path(tempfile.mkdtemp(prefix='open3spn2_numpy_'))

STRUCTURES = {'A': test_path / 'adna', 'B': test_path / 'bdna', 'B_curved': test_path / 'bdna_curv'}
FORCES = ['Bond', 'Angle', 'Stacking', 'Dihedral']
KERNELS = {'Bond': 'bond', 'Angle': 'angle', 'Stacking': 'stacking', 'Dihedral': 'dihedral'}
K = 0.7
ATOL = 1e-6  # kJ/mol


@lru_cache(maxsize=None)
def _dna(dna_type):
    return DNA.fromXYZ(str(STRUCTURES[dna_type] / 'in00_conf.xyz'), dna_type, template_from_X3DNA=False,
                       output_pdb=str(_outdir / f'{dna_type}.pdb'))


@lru_cache(maxsize=None)
def _frames(dna_type):
    """The structure, three seeded jitters of it and two frames of the LAMMPS trajectory, in nm."""
    x = _dna(dna_type).atoms[['x', 'y', 'z']].to_numpy(float) / 10
    rng = np.random.default_rng(7)
    jitters = [x + rng.normal(scale=0.03, size=x.shape) for _ in range(3)]
    traj = parse_xyz(str(STRUCTURES[dna_type] / 'traj.xyz'))
    steps = sorted(traj['timestep'].unique())
    recorded = [traj[traj['timestep'] == t][['x', 'y', 'z']].to_numpy(float) / 10 for t in (steps[1], steps[-1])]
    return np.stack([x] + jitters + recorded)


def _build(force_name, dna, **kwargs):
    return getattr(dna_forces, force_name)(dna, **kwargs)


def _no_openmm(*args, **kwargs):
    raise AssertionError('numpy energies must not use OpenMM')


def _openmm_energies(force, frames, box=None):
    """Energy of `force` alone, on the Reference platform, for each frame (kJ/mol)."""
    system = openmm.System()
    for _ in range(frames.shape[1]):
        system.addParticle(1.0)
    if box is not None:
        system.setDefaultPeriodicBoxVectors(*np.diag(box))
    force.addForce(system)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName('Reference'))
    energies = []
    for x in frames:
        context.setPositions(x)
        energies.append(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    return np.array(energies)


@pytest.mark.parametrize('dna_type', list(STRUCTURES))
@pytest.mark.parametrize('force_name', FORCES)
def test_numpy_energies_match_openmm(force_name, dna_type, monkeypatch):
    dna = _dna(dna_type)
    dna.periodic = False
    frames = _frames(dna_type)
    force = _build(force_name, dna, k=K)
    table = force.table
    table_energy = getattr(backend, f'{KERNELS[force_name]}_energy')

    with monkeypatch.context() as m:
        m.setattr(openmm, 'Context', _no_openmm)
        energies = force.energy(frames)
        single = force.energy(frames[0])
        functional = K * table_energy(table, frames)
    assert energies.shape == (len(frames), len(table))
    np.testing.assert_allclose(single, energies[0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(functional, energies, rtol=0, atol=1e-12)

    expected = _openmm_energies(_build(force_name, dna, k=K), frames)
    np.testing.assert_allclose(energies.sum(-1), expected, rtol=0, atol=ATOL)


@pytest.mark.parametrize('force_name', FORCES)
def test_periodic_box(force_name):
    dna = _dna('B')
    frames = _frames('B')[:2]
    box = np.array([13.7, 12.9, 14.3])
    wrapped = (frames - frames.mean(axis=(0, 1)) + box) % box  # centred on a box corner
    dna.periodic = True
    try:
        force = _build(force_name, dna)
        reference = _build(force_name, dna)
    finally:
        dna.periodic = False
    energies = force.energy(wrapped, box=box).sum(-1)
    np.testing.assert_allclose(energies, _openmm_energies(reference, wrapped, box=box), rtol=0, atol=ATOL)
    unwrapped = _build(force_name, dna).energy(wrapped).sum(-1)
    assert np.all(np.abs(unwrapped - energies) > 1.0)


def test_box_must_match_the_force():
    dna = _dna('B')
    x = _frames('B')[0]
    dna.periodic = True
    try:
        periodic = _build('Bond', dna)
    finally:
        dna.periodic = False
    plain = _build('Bond', dna)

    with pytest.raises(ValueError):
        periodic.energy(x)
    with pytest.raises(ValueError):
        plain.energy(x, box=np.array([20.0, 20.0, 20.0]))
    with pytest.raises(NotImplementedError):
        periodic.energy(x, box=np.array([[20.0, 0, 0], [5.0, 20.0, 0], [0, 0, 20.0]]))


def _set_parameters(force, i, particles, values):
    particles = [int(p) for p in particles]
    if isinstance(force, openmm.CustomCompoundBondForce):
        force.setBondParameters(i, particles, values)
    elif isinstance(force, openmm.CustomBondForce):
        force.setBondParameters(i, *particles, values)
    elif isinstance(force, openmm.CustomAngleForce):
        force.setAngleParameters(i, *particles, values)
    else:
        force.setTorsionParameters(i, *particles, values)


@pytest.mark.parametrize('force_name', FORCES)
def test_edited_table(force_name):
    dna = _dna('B_curved')
    frames = _frames('B_curved')
    table = _build(force_name, dna).table
    particle_columns = [c for c in ['aai', 'aaj', 'aak', 'aal'] if c in table]
    parameter_columns = [c for c in table if c not in particle_columns]
    changed = np.array([0, len(table) // 2, len(table) - 1])
    edited = table.copy()
    edited.loc[changed, parameter_columns] = edited.loc[changed, parameter_columns] * 1.1 + 0.05

    reference = _build(force_name, dna)
    for i in changed:
        _set_parameters(reference.force, int(i), edited.loc[i, particle_columns],
                        list(edited.loc[i, parameter_columns]))

    table_energy = getattr(backend, f'{KERNELS[force_name]}_energy')
    before, after = table_energy(table, frames), table_energy(edited, frames)
    unchanged = np.setdiff1d(np.arange(len(table)), changed)
    np.testing.assert_array_equal(after[:, unchanged], before[:, unchanged])
    assert np.all(np.abs(after[:, changed] - before[:, changed]).max(axis=1) > 1e-3)
    np.testing.assert_allclose(after.sum(-1), _openmm_energies(reference, frames), rtol=0, atol=ATOL)
