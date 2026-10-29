"""Identity tables: the parameters of the interactions of a duplex that depend on the identities of a set of
base-pair positions, at every assignment of names to the positions each interaction depends on.

The table of a kind is a DataFrame in the format of the force's ``parameter_table`` (particle and parameter columns)
with one row per (interaction, assignment) and the columns

* ``term``: the row of the interaction in the parameter table of the wild-type force;
* ``position_j``: the j-th position the interaction depends on, -1 padded to the widest interaction of the kind;
* ``name_j``: the index in the alphabet of the name of that position in this assignment, -1 padded.

The rows of an interaction are consecutive, its assignments in ``itertools.product`` order (first position slowest).
A state, an array with the name index of every position, selects one row of each interaction.

Bonded interactions lie on one strand and depend on the bases of that strand within one residue of their particles.
Their parameters are read from a covering design of at most 256 Watson-Crick mutants, in which every window of
positions that share an interaction sees all combinations of identities, so N positions never need 4^N mutants. A
mix has the parameters interpolated jointly between the assignments with all its positions set to its two identities.
"""
import itertools

import numpy as np
import pandas

from ..force import backend, dna as dna_forces
from .alphabet import WATSON_CRICK, alphabet as make_alphabet

FORCES = {'Bond': dna_forces.Bond, 'Angle': dna_forces.Angle, 'Stacking': dna_forces.Stacking,
          'Dihedral': dna_forces.Dihedral}
ENERGIES = {'Bond': backend.bond_energy, 'Angle': backend.angle_energy, 'Stacking': backend.stacking_energy,
            'Dihedral': backend.dihedral_energy}
CIRCULAR = {'Dihedral': ['t0']}
PARTICLES = ['aai', 'aaj', 'aak', 'aal']
COVER = list(WATSON_CRICK.values())
TOLERANCE = 1e-9


class IdentityTables:
    """The identity tables of the positions `sites` of a duplex: `alphabet`, `sites` [((chain, residue) on strand A,
    (chain, residue) on strand B)], `wild_type` (N,) the name indices of the duplex's own identities and `kinds`
    {kind: DataFrame}."""

    def __init__(self, alphabet, sites, wild_type, kinds):
        self.alphabet = alphabet
        self.sites = sites
        self.wild_type = np.asarray(wild_type, dtype=int)
        self.kinds = kinds

    def state_energies(self, coordinates, states, box=None):
        return state_energies(self, coordinates, states, box)


def state_energies(tables, coordinates, states, box=None):
    """{kind: energy (kJ/mol) of the identity-dependent interactions of the kind in each of `states` (S, N), shape
    (..., S)} at `coordinates` (nm)."""
    states = np.atleast_2d(states)
    energies = {}
    for kind, table in tables.kinds.items():
        positions = table[[c for c in table if c.startswith('position_')]].to_numpy(dtype=int)
        names = table[[c for c in table if c.startswith('name_')]].to_numpy(dtype=int)
        selected = ((positions < 0) | (states[:, positions] == names)).all(axis=-1)
        energies[kind] = ENERGIES[kind](table, coordinates, box) @ selected.T
    return energies


def mix_parameters(grid, alphabet, circular=()):
    """(n^k, parameters) the parameters at every assignment of the n names of `alphabet` to k positions, from `grid`
    (4,)*k + (parameters,) the parameters at every combination of Watson-Crick identities (in ``WATSON_CRICK`` order).
    The columns `circular` are angles interpolated along the shorter arc."""
    digit = {name: COVER.index(pair) for name, pair in alphabet.identities.items()}
    low = np.array([digit[name] if name in digit else digit[alphabet.mixes[name][0]] for name in alphabet.names])
    high = np.array([digit[name] if name in digit else digit[alphabet.mixes[name][1]] for name in alphabet.names])
    weight = np.array([alphabet.mixes[name][2] if name in alphabet.mixes else np.nan for name in alphabet.names])

    k = grid.ndim - 1
    assignments = np.array(list(itertools.product(range(len(alphabet)), repeat=k))).reshape(-1, k)
    weights = weight[assignments]
    heaviest, lightest = np.nanmax(weights, axis=1, initial=-np.inf), np.nanmin(weights, axis=1, initial=np.inf)
    mixed = np.isfinite(heaviest)
    if (heaviest[mixed] != lightest[mixed]).any():
        raise ValueError('an interaction would see mixes of different weights at once')
    lo, hi = grid[tuple(low[assignments].T)], grid[tuple(high[assignments].T)]
    difference = hi - lo
    difference[:, circular] = (difference[:, circular] + np.pi) % (2 * np.pi) - np.pi
    return lo + np.where(mixed, heaviest, 0.0)[:, None] * difference


def identity_tables(dna, sites, alphabet=None):
    """The `IdentityTables` of the base-pair positions `sites` [((chain, residue) on strand A, (chain, residue) on
    strand B)] of `dna` in `alphabet` (default 'wc'), for the bonded forces built with ``k = 1``."""
    alphabet = make_alphabet('wc') if alphabet is None else alphabet
    for name, pair in alphabet.identities.items():
        if pair not in COVER:
            raise NotImplementedError(f'identity {name} {pair} is not a Watson-Crick pair')
    sites = [tuple((chain, int(residue)) for chain, residue in pair) for pair in sites]
    position_of = {residue: p for p, pair in enumerate(sites) for residue in pair}

    atoms = dna.atoms
    bases = atoms[atoms['name'].isin(['A', 'C', 'G', 'T'])]
    base_of = dict(zip(zip(bases['chainID'], bases['resSeq'].astype(int)), bases['name']))
    identity_of = {pair: name for name, pair in alphabet.identities.items()}
    wild_type = []
    for pair in sites:
        if not all(residue in base_of for residue in pair):
            raise ValueError(f'position {pair} is not a pair of DNA bases')
        wild = tuple(base_of[residue] for residue in pair)
        if wild not in identity_of:
            raise ValueError(f'the duplex has {wild} at {pair}, which is not an identity of the alphabet')
        wild_type.append(alphabet.index(identity_of[wild]))

    residue_of = dict(zip(atoms.index, zip(atoms['chainID'], atoms['resSeq'].astype(int))))

    def window(particles):
        """The positions whose identity can change an interaction of `particles`: those within one residue of them,
        or only those at them when the former are more than four."""
        residues = {residue_of[i] for i in particles}
        for reach in (1, 0):
            positions = sorted({position_of[(chain, residue + d)] for chain, residue in residues
                                for d in range(-reach, reach + 1) if (chain, residue + d) in position_of})
            if len(positions) <= len(COVER):
                return positions
        raise ValueError(f'an interaction spans {len(positions)} positions; the covering design handles {len(COVER)}')

    wild_tables = {kind: force(dna).parameter_table for kind, force in FORCES.items()}
    windows = {kind: [window(row) for row in table[[c for c in PARTICLES if c in table]].to_numpy()]
               for kind, table in wild_tables.items()}

    # Positions that share a window take different digits of the build number, so every window sees every
    # combination of Watson-Crick identities.
    slot = []
    for p in range(len(sites)):
        taken = {slot[q] for kind_windows in windows.values() for w in kind_windows if p in w for q in w if q < p}
        free = [s for s in range(len(COVER)) if s not in taken]
        if not free:
            raise ValueError(f'position {p} shares interactions with positions in all {len(COVER)} slots')
        slot.append(free[0])
    builds = np.arange(len(COVER) ** (max(slot, default=0) + 1))
    digits = builds[:, None] // len(COVER) ** np.array(slot, dtype=int)[None, :] % len(COVER)

    values = {kind: [] for kind in FORCES}
    for build in digits:
        mutant = dna.create_mutant([(chain, residue, COVER[digit][strand]) for pair, digit in zip(sites, build)
                                    for strand, (chain, residue) in enumerate(pair)])
        for kind, force in FORCES.items():
            wild = wild_tables[kind]
            particles = [c for c in PARTICLES if c in wild]
            table = force(mutant).parameter_table.set_index(particles)
            if len(table) != len(wild) or not table.index.is_unique:
                raise ValueError(f'the {kind} interactions of a mutant differ from those of the wild type')
            table = table.reindex(pandas.MultiIndex.from_frame(wild[particles]))
            if table.isna().any(axis=None):
                raise ValueError(f'the {kind} interactions of a mutant differ from those of the wild type')
            values[kind].append(table.to_numpy())

    kinds = {}
    for kind, wild in wild_tables.items():
        particles = [c for c in PARTICLES if c in wild]
        parameters = [c for c in wild if c not in particles]
        circular = [parameters.index(c) for c in CIRCULAR.get(kind, [])]
        v = np.stack(values[kind])
        terms = []
        for t in np.flatnonzero(np.abs(v - v[0]).max(axis=(0, 2), initial=0) > TOLERANCE):
            window_t = windows[kind][t]
            if not window_t:
                raise ValueError(f'{kind} interaction {t} depends on the identities but is far from every position')
            grid = np.full((len(COVER),) * len(window_t) + (len(parameters),), np.nan)
            grid[tuple(digits[:, window_t].T)] = v[:, t]
            axes = [j for j in range(len(window_t)) if np.abs(grid - np.take(grid, [0], axis=j)).max() > TOLERANCE]
            grid = grid[tuple(slice(None) if j in axes else 0 for j in range(len(window_t)))]
            positions = [window_t[j] for j in axes]
            if not positions or np.abs(grid[tuple(digits[:, positions].T)] - v[:, t]).max() > TOLERANCE:
                raise ValueError(f'{kind} interaction {t} depends on positions outside {window_t}')
            terms.append((t, positions, mix_parameters(grid, alphabet, circular)))

        width = max((len(positions) for _, positions, _ in terms), default=0)
        rows = [len(p) for _, _, p in terms]
        term = np.repeat([t for t, _, _ in terms], rows).astype(int)
        position = np.full((len(term), width), -1)
        name = np.full((len(term), width), -1)
        start = 0
        for (t, positions, _), n in zip(terms, rows):
            position[start:start + n, :len(positions)] = positions
            name[start:start + n, :len(positions)] = list(itertools.product(range(len(alphabet)),
                                                                            repeat=len(positions)))
            start += n
        table = wild.loc[term, particles].reset_index(drop=True)
        table.insert(0, 'term', term)
        for j in range(width):
            table[f'position_{j}'] = position[:, j]
        for j in range(width):
            table[f'name_{j}'] = name[:, j]
        table[parameters] = np.concatenate([p for _, _, p in terms]) if terms else np.empty((0, len(parameters)))
        kinds[kind] = table
    return IdentityTables(alphabet, sites, wild_type, kinds)
