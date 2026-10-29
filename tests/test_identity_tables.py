"""Tests for the identity tables of the Gibbs identity sampler.

An alphabet lists the names a base-pair position can take: identities (one base on each strand) and mixes (the
parameter interpolation between two identities). ``identity_tables(dna, sites, alphabet)`` finds the interactions
whose parameters depend on the identities of the base-pair positions ``sites`` and lists their parameters at every
assignment of names to the positions each interaction depends on: one DataFrame per kind, in the format of the
force's ``parameter_table`` plus the columns ``term`` (the row of the wild-type parameter table), ``position_j`` and
``name_j`` (the positions and the name index of each in this assignment, -1 padded). ``state_energies`` evaluates
the tables for states, arrays with the name index of every position. The behavioural requirements tested here:

* I1. The 'wc' alphabet is GC, AT, TA, CG and the mix W halfway between AT and TA; 'wc4' has no mix. A mix that
  does not interpolate between two identities, reuses a name or has a weight outside (0, 1) is refused.
* I2. For every bonded kind and DNA type, the energy ``state_energies`` gives a Watson-Crick state, relative to the
  wild type, equals the energy change of that force built natively on the mutant duplex (``DNA.create_mutant``),
  for one frame and for a batch of frames; evaluating the tables builds no mutant and does not use OpenMM. The table
  of a kind is empty exactly when the identities do not change its energy.
* I3. A state with W at some positions has the energy of a duplex whose parameters are halfway between those of the
  duplexes with AT and with TA at all of those positions; a circular parameter (an angle such as the dihedral t0) is
  interpolated along the shorter arc.
* I4. Bad inputs are refused: a position whose wild-type base pair is not an identity of the alphabet, an alphabet
  with an identity that is not Watson-Crick (not implemented) and an interaction that would see mixes of different
  weights at once.
"""
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import openmm
import pytest

import pandas

from open3SPN2 import DNA
from open3SPN2.alchemy import Alphabet, alphabet, identity_tables, mix_parameters
from open3SPN2.force import backend, dna as dna_forces
from open3SPN2.force.template import DNAForce

test_path = Path(__file__).parent
_outdir = Path(tempfile.mkdtemp(prefix='open3spn2_identity_'))

STRUCTURES = {'A': test_path / 'adna', 'B': test_path / 'bdna', 'B_curved': test_path / 'bdna_curv'}
BONDED = ['Bond', 'Angle', 'Stacking', 'Dihedral']
# Strand-A residues of the positions: two adjacent ones and one three residues away, sharing some interactions.
SITE_RESIDUES = [12, 13, 16]
ATOL = 1e-6  # kJ/mol


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
def _tables(dna_type):
    dna = _dna(dna_type)
    return identity_tables(dna, _sites(dna), alphabet('wc'))


def _frames(dna_type):
    x = _dna(dna_type).atoms[['x', 'y', 'z']].to_numpy(float) / 10
    rng = np.random.default_rng(11)
    return np.stack([x] + [x + rng.normal(scale=0.03, size=x.shape) for _ in range(2)])


def _mutant(dna, sites, state, names):
    mutations = []
    for ((chain_a, residue_a), (chain_b, residue_b)), name in zip(sites, state):
        base_a, base_b = names.identities[names.names[name]]
        mutations += [(chain_a, residue_a, base_a), (chain_b, residue_b, base_b)]
    return dna.create_mutant(mutations)


def _refuse(*args, **kwargs):
    raise AssertionError('evaluating identity tables must not build mutants, forces or OpenMM contexts')


def _state_energies(tables, coordinates, states, kind, monkeypatch):
    """`tables.state_energies` of `kind`, with DNA topologies, force classes and OpenMM contexts unavailable."""
    with monkeypatch.context() as m:
        m.setattr(DNA, 'create_mutant', _refuse)
        m.setattr(DNA, 'computeTopology', _refuse)
        m.setattr(DNAForce, '__init__', _refuse)
        m.setattr(openmm, 'Context', _refuse)
        return tables.state_energies(coordinates, states)[kind]


def _midpoint_energy(dna, sites, state, names, kind, frames):
    """Energy of the force `kind` with every parameter halfway between those of the duplex with AT and the duplex
    with TA at the W positions of `state`, the dihedral t0 along the shorter arc."""
    w, at, ta = (names.index(name) for name in ('W', 'AT', 'TA'))
    force = getattr(dna_forces, kind)
    low = force(_mutant(dna, sites, np.where(state == w, at, state), names)).parameter_table
    high = force(_mutant(dna, sites, np.where(state == w, ta, state), names)).parameter_table
    particles = [c for c in ('aai', 'aaj', 'aak', 'aal') if c in low]
    parameters = [c for c in low if c not in particles]
    high = high.set_index(particles).loc[pandas.MultiIndex.from_frame(low[particles])].reset_index()
    difference = high[parameters].to_numpy() - low[parameters].to_numpy()
    if kind == 'Dihedral':
        t0 = parameters.index('t0')
        difference[:, t0] = (difference[:, t0] + np.pi) % (2 * np.pi) - np.pi
    middle = low.copy()
    middle[parameters] = low[parameters].to_numpy() + difference / 2
    return getattr(backend, f'{kind.lower()}_energy')(middle, frames).sum(-1)


def test_alphabet():
    """I1: the presets have the expected names and mixes, and inconsistent mixes are refused."""
    wc = alphabet('wc')
    assert wc.names == ['GC', 'AT', 'TA', 'CG', 'W']
    assert len(wc) == 5 and wc.index('W') == 4
    assert wc.identities == {'GC': ('G', 'C'), 'AT': ('A', 'T'), 'TA': ('T', 'A'), 'CG': ('C', 'G')}
    assert wc.mixes == {'W': ('AT', 'TA', 0.5)}
    assert alphabet('wc4').names == ['GC', 'AT', 'TA', 'CG'] and alphabet('wc4').mixes == {}

    identities = dict(wc.identities)
    for mixes in [{'W': ('AT', 'XY', 0.5)}, {'AT': ('GC', 'TA', 0.5)}, {'W': ('AT', 'TA', 0.0)},
                  {'W': ('AT', 'TA', 1.5)}]:
        with pytest.raises(ValueError):
            Alphabet(identities, mixes)
    with pytest.raises(ValueError):
        alphabet('nucleotides')


@pytest.mark.parametrize('dna_type', list(STRUCTURES))
@pytest.mark.parametrize('kind', BONDED)
def test_state_energies_match_mutants(kind, dna_type, monkeypatch):
    """I2: relative to the wild type, the table energy of random Watson-Crick states equals the energy change of the
    force built on the mutant duplex, frame by frame, the tables are evaluated without mutants, forces or OpenMM, and
    only the kinds whose energy changes have rows."""
    dna = _dna(dna_type)
    sites = _sites(dna)
    tables = _tables(dna_type)
    frames = _frames(dna_type)
    rng = np.random.default_rng(3)
    states = np.vstack([tables.wild_type, rng.integers(0, 4, size=(6, len(sites)))])

    energies = _state_energies(tables, frames, states, kind, monkeypatch)
    single = _state_energies(tables, frames[0], states, kind, monkeypatch)
    assert energies.shape == (len(frames), len(states))
    np.testing.assert_allclose(single, energies[0], rtol=0, atol=1e-12)

    force = getattr(dna_forces, kind)
    wild = force(dna).energy(frames).sum(-1)
    expected = np.stack([force(_mutant(dna, sites, s, tables.alphabet)).energy(frames).sum(-1) - wild
                         for s in states], axis=-1)
    np.testing.assert_allclose(energies - energies[:, :1], expected, rtol=0, atol=ATOL)
    assert (len(tables.kinds[kind]) > 0) == (np.abs(expected).max() > ATOL)


@pytest.mark.parametrize('dna_type', list(STRUCTURES))
@pytest.mark.parametrize('kind', BONDED)
def test_mix_is_the_joint_midpoint(kind, dna_type, monkeypatch):
    """I3: relative to the wild type, the table energy of states with W at one or more positions equals the energy
    of the force with the parameters halfway between the duplexes with AT and with TA at those positions."""
    dna = _dna(dna_type)
    sites = _sites(dna)
    tables = _tables(dna_type)
    frames = _frames(dna_type)
    names = tables.alphabet
    w, gc, ta = (names.index(name) for name in ('W', 'GC', 'TA'))
    states = np.array([tables.wild_type, [w, w, w], [w, gc, w], [ta, w, gc], [gc, gc, w]])

    energies = _state_energies(tables, frames, states, kind, monkeypatch)
    expected = np.stack([_midpoint_energy(dna, sites, s, names, kind, frames) for s in states], axis=-1)
    np.testing.assert_allclose(energies - energies[:, :1], expected - expected[:, :1], rtol=0, atol=ATOL)


def test_mix_of_angles_takes_the_shorter_arc():
    """I3: a circular parameter halfway between 170 and -170 degrees is 180 degrees, and a plain one is 0."""
    names = alphabet('wc')
    grid = np.zeros((4, 2))
    grid[names.index('AT')] = np.radians([170.0, 170.0])
    grid[names.index('TA')] = np.radians([-170.0, -170.0])
    parameters = mix_parameters(grid, names, circular=[0])
    angle, plain = parameters[names.index('W')]
    assert np.cos(angle) == pytest.approx(-1.0, abs=1e-12)
    assert plain == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_allclose(parameters[:4], grid, rtol=0, atol=0)


@pytest.mark.parametrize('identity', [('GA', ('G', 'A')), ('Z1', ('A', 'C')), ('XY', ('X', 'Y'))])
def test_non_watson_crick_identities_are_not_implemented(identity):
    """I4: an alphabet with an identity that is not a Watson-Crick pair, whatever its name, is refused."""
    dna = _dna('B_curved')
    name, pair = identity
    with pytest.raises(NotImplementedError):
        identity_tables(dna, _sites(dna), Alphabet({**alphabet('wc4').identities, name: pair}))


def test_bad_inputs():
    """I4: a wild-type base pair outside the alphabet and mixes of different weights in one interaction are
    refused."""
    dna = _dna('B_curved')
    sites = _sites(dna)
    with pytest.raises(ValueError, match='not an identity of the alphabet'):
        identity_tables(dna, sites, Alphabet({'AT': ('A', 'T'), 'TA': ('T', 'A')}))
    unequal = Alphabet(alphabet('wc4').identities, {'W': ('AT', 'TA', 0.5), 'S': ('GC', 'CG', 0.25)})
    with pytest.raises(ValueError, match='mixes of different weights'):
        identity_tables(dna, sites, unequal)
