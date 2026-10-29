"""Identity tables: the parameters of the interactions of a duplex that depend on the identities of a set of
base-pair positions, at every assignment of names to the positions each interaction depends on.

The table of a kind is a DataFrame in the format of the force's ``parameter_table`` (particle and parameter columns)
with one row per (interaction, assignment) and the columns

* ``term``: the number of the interaction in the table;
* ``position_j``: the j-th position the interaction depends on, -1 padded to the widest interaction of the kind;
* ``name_j``: the index in the alphabet of the name of that position in this assignment, -1 padded.

The rows of an interaction are consecutive, its assignments in ``itertools.product`` order (first position slowest).
A state, an array with the name index of every position, selects one row of each interaction.

Bonded interactions lie on one strand and depend on the bases of that strand within one residue of their particles.
Their parameters are read from a covering design of at most 256 Watson-Crick mutants, in which every window of
positions that share an interaction sees all combinations of identities, so N positions never need 4^N mutants.

Nonbonded interactions (base pairing, cross stacking and excluded volume) join beads anywhere in the duplex. The
table lists every donor-acceptor or particle pair the force would form for some identities of the positions and that
involves the base of a position, with the parameters the force-field definitions give to the names of its bases. A
pair that does not interact at an assignment has epsilon 0: a mismatch keeps the other parameters of the first
assignment at which it pairs, and a base pair masked from the excluded volume keeps its own.

A mix has the parameters interpolated jointly between the assignments with all its positions set to its two
identities. When a nonbonded pair interacts at only one of them, the mix has the parameters of that one at the
interpolated epsilon.
"""
import itertools

import numpy as np
import openmm.unit as unit
import pandas

from ..force import backend, dna as dna_forces
from .alphabet import WATSON_CRICK, alphabet as make_alphabet

BONDED = {'Bond': dna_forces.Bond, 'Angle': dna_forces.Angle, 'Stacking': dna_forces.Stacking,
          'Dihedral': dna_forces.Dihedral}
ENERGIES = {'Bond': backend.bond_energy, 'Angle': backend.angle_energy, 'Stacking': backend.stacking_energy,
            'Dihedral': backend.dihedral_energy, 'BasePair': backend.base_pair_energy,
            'CrossStacking': backend.cross_stacking_energy, 'Exclusion': backend.exclusion_energy}
CIRCULAR = {'Dihedral': ['t0'], 'BasePair': ['phi0']}
STRENGTH = {'BasePair': 'epsilon', 'CrossStacking': 'epsilon', 'Exclusion': 'epsilon'}
PARTICLES = ['aai', 'aaj', 'aak', 'aal']
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
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


def mix_parameters(grid, alphabet, circular=(), strength=None):
    """(n^k, parameters) the parameters at every assignment of the n names of `alphabet` to k positions, from `grid`
    (4,)*k + (parameters,) the parameters at every combination of Watson-Crick identities (in ``WATSON_CRICK`` order).
    The columns `circular` are angles interpolated along the shorter arc. If `strength` is the column of the strength
    of the interaction, a mix between an assignment where it is 0 and one where it is not has the parameters of the
    latter at the interpolated strength."""
    digit = {name: COVER.index(pair) for name, pair in alphabet.identities.items()}
    low = np.array([digit[name] if name in digit else digit[alphabet.mixes[name][0]] for name in alphabet.names])
    high = np.array([digit[name] if name in digit else digit[alphabet.mixes[name][1]] for name in alphabet.names])
    weight = np.array([alphabet.mixes[name][2] if name in alphabet.mixes else np.nan for name in alphabet.names])

    k = grid.ndim - 1
    assignments = list(itertools.product(range(len(alphabet)), repeat=k))
    assignments = np.array(assignments, dtype=int).reshape(len(assignments), k)
    weights = weight[assignments]
    heaviest, lightest = np.nanmax(weights, axis=1, initial=-np.inf), np.nanmin(weights, axis=1, initial=np.inf)
    mixed = np.isfinite(heaviest)
    if (heaviest[mixed] != lightest[mixed]).any():
        raise ValueError('an interaction would see mixes of different weights at once')
    lo, hi = grid[tuple(low[assignments].T)], grid[tuple(high[assignments].T)]
    if strength is not None:
        off = lo[:, strength] == 0
        lo[off] = hi[off]
        lo[off, strength] = 0.0
        off = hi[:, strength] == 0
        hi[off] = lo[off]
        hi[off, strength] = 0.0
    difference = hi - lo
    difference[:, circular] = (difference[:, circular] + np.pi) % (2 * np.pi) - np.pi
    return lo + np.where(mixed, heaviest, 0.0)[:, None] * difference


def identity_tables(dna, sites, alphabet=None):
    """The `IdentityTables` of the base-pair positions `sites` [((chain, residue) on strand A, (chain, residue) on
    strand B)] of `dna` in `alphabet` (default 'wc'), for the forces built with ``k = 1``."""
    alphabet = make_alphabet('wc') if alphabet is None else alphabet
    for name, pair in alphabet.identities.items():
        if pair not in COVER:
            raise NotImplementedError(f'identity {name} {pair} is not a Watson-Crick pair')
    sites = [tuple((chain, int(residue)) for chain, residue in pair) for pair in sites]

    atoms = dna.atoms
    bases = atoms[atoms['name'].isin(list(COMPLEMENT))]
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

    kinds = {**_bonded_tables(dna, sites, alphabet), **_nonbonded_tables(dna, sites, alphabet)}
    return IdentityTables(alphabet, sites, wild_type, kinds)


def _bonded_tables(dna, sites, alphabet):
    """{kind: identity table} of the bonded forces, read from a covering design of Watson-Crick mutants."""
    position_of = {residue: p for p, pair in enumerate(sites) for residue in pair}
    atoms = dna.atoms
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

    wild_tables = {kind: force(dna).parameter_table for kind, force in BONDED.items()}
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

    values = {kind: [] for kind in BONDED}
    for build in digits:
        mutant = dna.create_mutant([(chain, residue, COVER[digit][strand]) for pair, digit in zip(sites, build)
                                    for strand, (chain, residue) in enumerate(pair)])
        for kind, force in BONDED.items():
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
        v = np.stack(values[kind])
        interactions = []
        for t in np.flatnonzero(np.abs(v - v[0]).max(axis=(0, 2), initial=0) > TOLERANCE):
            window_t = windows[kind][t]
            if not window_t:
                raise ValueError(f'{kind} interaction {t} depends on the identities but is far from every position')
            grid = np.full((len(COVER),) * len(window_t) + (len(parameters),), np.nan)
            grid[tuple(digits[:, window_t].T)] = v[:, t]
            grid, axes = _dependent(grid)
            positions = [window_t[j] for j in axes]
            if not positions or np.abs(grid[tuple(digits[:, positions].T)] - v[:, t]).max() > TOLERANCE:
                raise ValueError(f'{kind} interaction {t} depends on positions outside {window_t}')
            interactions.append((wild.loc[t, particles].to_numpy(dtype=int), positions, grid))
        kinds[kind] = _kind_table(kind, particles, parameters, interactions, alphabet)
    return kinds


def _nonbonded_tables(dna, sites, alphabet):
    """{kind: identity table} of the base pairing, cross stacking and excluded volume forces, from the force-field
    definitions of `dna`."""
    position_of = {residue: (p, strand) for p, pair in enumerate(sites) for strand, residue in enumerate(pair)}
    atoms = dna.atoms
    residue_of = list(zip(atoms['chainID'], atoms['resSeq'].astype(int)))
    seqid = dict(zip(residue_of, dna_forces.chained_sequence_index(atoms)))
    is_dna = atoms['resname'].isin(dna_forces._dnaResidues).to_numpy()
    sugar = {residue_of[i]: i for i in np.flatnonzero(is_dna & (atoms['name'] == 'S').to_numpy())}
    base = {residue_of[i]: i for i in np.flatnonzero(is_dna & atoms['name'].isin(list(COMPLEMENT)).to_numpy())
            if residue_of[i] in sugar}
    name = {residue: atoms['name'].iloc[i] for residue, i in base.items()}

    def interaction(residues, rule, epsilon):
        """(positions, grid (4,)*k + (parameters,)) of an interaction whose parameters `rule(*names)` depend on the
        base names of `residues`, or None if it has no energy at any assignment. `rule` returns None for a pair that
        does not interact; `epsilon` is the column of its strength."""
        positions = sorted({position_of[r][0] for r in residues if r in position_of})
        rows = []
        for digits in itertools.product(range(len(COVER)), repeat=len(positions)):
            names = [COVER[digits[positions.index(position_of[r][0])]][position_of[r][1]] if r in position_of
                     else name[r] for r in residues]
            rows.append(rule(*names))
        interacting = [row for row in rows if row is not None and row[epsilon] != 0]
        if not interacting:
            return None
        grid = np.array([row if row is not None else interacting[0][:epsilon] + [0.0] + interacting[0][epsilon + 1:]
                         for row in rows])
        grid, axes = _dependent(grid.reshape((len(COVER),) * len(positions) + (grid.shape[-1],)))
        return [positions[j] for j in axes], grid

    kinds = {}

    base_pair = dna_forces.BasePair(dna)
    cutoff = next(iter(base_pair.forces.values())).getCutoffDistance().value_in_unit(unit.nanometer)
    pairs = {(row.Base1, row.Base2): row for row in
             dna.pair_definition[dna.pair_definition['DNA'] == dna.DNAtype].itertuples()}

    def pair_rule(d1, a1):
        row = pairs.get((d1, a1))
        if row is None:
            return None
        return [np.radians(row.torsion), row.sigma, np.radians(row.t1), np.radians(row.t2), row.rang, row.epsilon,
                row.alpha, cutoff]

    interactions = []
    for d, a in itertools.permutations(base, 2):
        if (d in position_of or a in position_of) and abs(seqid[a] - seqid[d]) >= base_pair.min_seq_sep:
            found = interaction([d, a], pair_rule, 5)
            if found:
                interactions.append(([base[d], sugar[d], base[a], sugar[a]], *found))
    kinds['BasePair'] = _kind_table('BasePair', ['d1', 'd2', 'a1', 'a2'],
                                    ['phi0', 'sigma', 't01', 't02', 'rng', 'epsilon', 'alpha', 'cutoff'],
                                    interactions, alphabet)

    cross_stacking = dna_forces.CrossStacking(dna)
    force = next(iter(cross_stacking.crossStackingForces.values()))[0]
    cutoff = force.getCutoffDistance().value_in_unit(unit.nanometer)
    crosses = {(row.Base_d1, row.Base_a1, row.Base_a3): row._asdict() for row in
               dna.cross_definition[dna.cross_definition['DNA'] == dna.DNAtype].itertuples()}

    def cross_rule(branch):
        def rule(d1, a1, a3):
            if d1 != COMPLEMENT[a1]:
                return None
            row = crosses[(d1, a1, a3)]
            return [np.radians(row['t03']), np.radians(row[f'T0CS_{branch}']), row[f'rng_cs{branch}'], row['rng_bp'],
                    row[f'eps_cs{branch}'], row[f'alpha_cs{branch}'], row[f'Sigma_{branch}'], cutoff]
        return rule

    # Branch 1 pairs a base and its next base with a base and its previous one; branch 2 the other way round.
    interactions = []
    for branch, (d_shift, a_shift) in ((1, (1, -1)), (2, (-1, 1))):
        for d, a in itertools.permutations(base, 2):
            d3, a3 = (d[0], d[1] + d_shift), (a[0], a[1] + a_shift)
            if d3 not in base or a3 not in base or not {d, a, a3} & position_of.keys():
                continue
            if abs(seqid[a] - seqid[d]) < cross_stacking.min_seq_sep:
                continue
            if {base[d], sugar[d], base[d3]} & {base[a], sugar[a], base[a3]}:
                continue
            found = interaction([d, a, a3], cross_rule(branch), 4)
            if found:
                interactions.append(([base[d], sugar[d], base[a], sugar[a], base[a3]], *found))
    kinds['CrossStacking'] = _kind_table('CrossStacking', ['d1', 'd2', 'a1', 'a2', 'a3'],
                                         ['t03', 't0CS', 'rng_CS', 'rng_BP', 'epsilon', 'alpha', 'sigma', 'cutoff'],
                                         interactions, alphabet)

    exclusion = dna_forces.Exclusion(dna)
    cutoff = exclusion.force.getCutoffDistance().value_in_unit(unit.nanometer)
    fixed = [exclusion.force.getParticleParameters(i) for i in range(exclusion.force.getNumParticles())]
    definition = dna.particle_definition[dna.particle_definition['DNA'] == dna.DNAtype].set_index('name')
    excluded = {tuple(sorted(exclusion.force.getExclusionParticles(i)))
                for i in range(exclusion.force.getNumExclusions())}
    variable = {base[residue]: residue for residue in position_of}

    def exclusion_rule(i, j):
        def particle(index, names):
            if index not in variable:
                return fixed[index][:3]
            bead = next(names)
            return definition.loc[bead, 'epsilon'], definition.loc[bead, 'radius'], dna_forces._base_type[bead]

        def rule(*names):
            names = iter(names)
            (epsilon_i, sigma_i, type_i), (epsilon_j, sigma_j, type_j) = particle(i, names), particle(j, names)
            epsilon = np.sqrt(epsilon_i * epsilon_j) if type_i + type_j != 5 else 0.0
            return [epsilon, 0.5 * (sigma_i + sigma_j), cutoff]
        return rule

    interactions = []
    for i in sorted(variable):
        for j in range(len(atoms)):
            pair = (min(i, j), max(i, j))
            if j == i or (j in variable and j < i) or pair in excluded:
                continue
            if abs(fixed[i][3] - fixed[j][3]) < exclusion.min_seq_sep:
                continue
            residues = [variable[k] for k in pair if k in variable]
            found = interaction(residues, exclusion_rule(*pair), 0)
            if found:
                interactions.append((list(pair), *found))
    kinds['Exclusion'] = _kind_table('Exclusion', ['aai', 'aaj'], ['epsilon', 'sigma', 'cutoff'], interactions,
                                     alphabet)
    return kinds


def _dependent(grid):
    """(grid restricted to the axes along which it varies, those axes) for a grid (4,)*k + (parameters,)."""
    k = grid.ndim - 1
    axes = [j for j in range(k) if np.abs(grid - np.take(grid, [0], axis=j)).max() > TOLERANCE]
    return grid[tuple(slice(None) if j in axes else 0 for j in range(k))], axes


def _kind_table(kind, particles, parameters, interactions, alphabet):
    """The identity table of `interactions` [(particle indices, positions, grid of Watson-Crick parameters)]."""
    circular = [parameters.index(c) for c in CIRCULAR.get(kind, [])]
    strength = parameters.index(STRENGTH[kind]) if kind in STRENGTH else None
    values = [mix_parameters(grid, alphabet, circular, strength) for _, _, grid in interactions]
    rows = [len(v) for v in values]
    width = max((len(positions) for _, positions, _ in interactions), default=0)
    position = np.full((sum(rows), width), -1)
    name = np.full((sum(rows), width), -1)
    start = 0
    for (_, positions, _), n in zip(interactions, rows):
        position[start:start + n, :len(positions)] = positions
        assignments = list(itertools.product(range(len(alphabet)), repeat=len(positions)))
        name[start:start + n, :len(positions)] = np.array(assignments, dtype=int).reshape(n, len(positions))
        start += n
    particle_indices = np.array([p for p, _, _ in interactions], dtype=int).reshape(-1, len(particles))
    table = pandas.DataFrame(np.repeat(particle_indices, rows, axis=0), columns=particles)
    table.insert(0, 'term', np.repeat(np.arange(len(interactions)), rows))
    for j in range(width):
        table[f'position_{j}'] = position[:, j]
    for j in range(width):
        table[f'name_{j}'] = name[:, j]
    table[parameters] = np.concatenate(values) if values else np.empty((0, len(parameters)))
    return table
