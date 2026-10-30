"""Identity tables: the parameters of the interactions of a duplex that depend on the identities of a set of
base-pair positions, at every assignment of names to the positions each interaction depends on.

The table of a kind is a DataFrame in the format of the force's ``parameter_table`` (particle and parameter columns)
with one row per (interaction, assignment) and the columns

* ``term``: the number of the interaction in the table;
* ``position_j``: the j-th position the interaction depends on, -1 padded to the widest interaction of the kind;
* ``name_j``: the index in the alphabet of the name of that position in this assignment, -1 padded.

The rows of an interaction are consecutive, its assignments in ``itertools.product`` order (first position slowest).
A state, an array with the name index of every position, selects one row of each interaction.

There is one function per force. Bonded interactions (``bonded_identity_table``) lie on one strand and depend on the
bases of that strand within one residue of their particles; their parameters are read from the forces built on a
covering design of at most 256 Watson-Crick mutants, in which every window of positions that share an interaction
sees all combinations of identities, so N positions never need 4^N mutants.

Nonbonded interactions (``base_pair_identity_table``, ``cross_stacking_identity_table`` and
``exclusion_identity_table``) join beads anywhere in the duplex. Their table lists every donor-acceptor or particle
pair the force would form for some identities of the positions and that involves the base of a position, with the
parameters the force-field definitions give to the names of its bases. A pair that does not interact at an assignment
has epsilon 0: a mismatch keeps the other parameters of the first assignment at which it pairs, and a base pair
masked from the excluded volume keeps its own.

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
from .switch import add_forces, set_state

BONDED = {'Bond': dna_forces.Bond, 'Angle': dna_forces.Angle, 'Stacking': dna_forces.Stacking,
          'Dihedral': dna_forces.Dihedral}
ENERGIES = {'Bond': backend.bond_energy, 'Angle': backend.angle_energy, 'Stacking': backend.stacking_energy,
            'Dihedral': backend.dihedral_energy, 'BasePair': backend.base_pair_energy,
            'CrossStacking': backend.cross_stacking_energy, 'Exclusion': backend.exclusion_energy}
CIRCULAR = {'Dihedral': ['t0']}
PARTICLES = ['aai', 'aaj', 'aak', 'aal']
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
COVER = list(WATSON_CRICK.values())
TOLERANCE = 1e-9


class IdentityTables:
    """The identity tables of the positions `sites` of the duplex `dna`: `alphabet`, `sites` [((chain, residue) on
    strand A, (chain, residue) on strand B)], `wild_type` (N,) the name indices of the duplex's own identities and
    `kinds` {kind: DataFrame}."""

    def __init__(self, dna, alphabet, sites, wild_type, kinds):
        self.dna = dna
        self.alphabet = alphabet
        self.sites = sites
        self.wild_type = np.asarray(wild_type, dtype=int)
        self.kinds = kinds

    def state_energies(self, coordinates, states, box=None):
        return state_energies(self, coordinates, states, box)

    def add_forces(self, system):
        return add_forces(self, system)

    def set_state(self, context, state, previous=None):
        set_state(self, context, state, previous)


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
    if any(pair not in COVER for pair in alphabet.identities.values()):
        raise NotImplementedError('identities that are not Watson-Crick pairs')
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

    bases = base_table(dna, sites)
    identity_of = {pair: name for name, pair in alphabet.identities.items()}
    wild_type = []
    for pair in sites:
        wild = tuple(bases.at[residue, 'name'] for residue in pair)
        if wild not in identity_of:
            raise ValueError(f'the duplex has {wild} at {pair}, which is not an identity of the alphabet')
        wild_type.append(alphabet.index(identity_of[wild]))

    design = covering_design(dna, sites)
    kinds = {kind: bonded_identity_table(force(dna), sites, alphabet, design) for kind, force in BONDED.items()}
    kinds['BasePair'] = base_pair_identity_table(dna_forces.BasePair(dna), sites, alphabet)
    kinds['CrossStacking'] = cross_stacking_identity_table(dna_forces.CrossStacking(dna), sites, alphabet)
    kinds['Exclusion'] = exclusion_identity_table(dna_forces.Exclusion(dna), sites, alphabet)
    return IdentityTables(dna, alphabet, sites, wild_type, kinds)


def base_table(dna, sites=()):
    """DataFrame of the DNA bases of `dna`, indexed by (chainID, resSeq): the particle indices of the ``base`` and of
    its ``sugar``, the sequence index ``seqid`` of the force masks, the base ``name`` and, for the bases of `sites`,
    their ``position`` and ``strand`` (0 for strand A, 1 for strand B; -1 for the other bases)."""
    atoms = dna.atoms.assign(seqid=dna_forces.chained_sequence_index(dna.atoms), particle=dna.atoms.index,
                             resSeq=dna.atoms['resSeq'].astype(int))
    atoms = atoms[atoms['resname'].isin(dna_forces._dnaResidues)]
    bases = atoms[atoms['name'].isin(list(COMPLEMENT))].set_index(['chainID', 'resSeq'])
    sugars = atoms[atoms['name'] == 'S'].set_index(['chainID', 'resSeq'])
    table = pandas.DataFrame({'base': bases['particle'], 'sugar': sugars['particle'].reindex(bases.index),
                              'seqid': bases['seqid'], 'name': bases['name'], 'position': -1, 'strand': -1})
    table = table.dropna().astype({'sugar': int})
    for p, pair in enumerate(sites):
        for strand, (chain, residue) in enumerate(pair):
            if (chain, int(residue)) not in table.index:
                raise ValueError(f'position {pair} is not a pair of DNA bases')
            table.loc[(chain, int(residue)), ['position', 'strand']] = p, strand
    return table


def covering_design(dna, sites):
    """(digits (builds, N), mutants): Watson-Crick mutants of `dna` in which every set of positions `sites` that share
    a bonded interaction takes all combinations of identities; digits[b, p] is the identity (in ``WATSON_CRICK``
    order) of position p in mutant b."""
    windows = [w for force in BONDED.values() for w in _windows(dna, sites, force(dna).parameter_table)]
    # Positions that share a window take different digits of the build number.
    slot = []
    for p in range(len(sites)):
        taken = {slot[q] for w in windows if p in w for q in w if q < p}
        free = [s for s in range(len(COVER)) if s not in taken]
        if not free:
            raise ValueError(f'position {p} shares interactions with positions in all {len(COVER)} slots')
        slot.append(free[0])
    builds = np.arange(len(COVER) ** (max(slot, default=0) + 1))
    digits = builds[:, None] // len(COVER) ** np.array(slot, dtype=int)[None, :] % len(COVER)
    mutants = [dna.create_mutant([(chain, residue, COVER[digit][strand]) for pair, digit in zip(sites, build)
                                  for strand, (chain, residue) in enumerate(pair)]) for build in digits]
    return digits, mutants


def bonded_identity_table(force, sites, alphabet, design=None):
    """The identity table of a bonded force (Bond, Angle, Stacking or Dihedral) built on the wild-type duplex, read
    from the same force built on the mutants of `design` (default ``covering_design(force.dna, sites)``)."""
    digits, mutants = covering_design(force.dna, sites) if design is None else design
    wild = force.parameter_table
    particles = [c for c in PARTICLES if c in wild]
    parameters = [c for c in wild if c not in particles]
    values = []
    for mutant in mutants:
        table = type(force)(mutant).parameter_table.set_index(particles)
        if len(table) != len(wild) or not table.index.is_unique:
            raise ValueError(f'the interactions of {type(force).__name__} on a mutant differ from the wild type')
        table = table.reindex(pandas.MultiIndex.from_frame(wild[particles]))
        if table.isna().any(axis=None):
            raise ValueError(f'the interactions of {type(force).__name__} on a mutant differ from the wild type')
        values.append(table[parameters].to_numpy())
    values = np.stack(values)

    windows = _windows(force.dna, sites, wild)
    interactions = []
    for t in np.flatnonzero(np.abs(values - values[0]).max(axis=(0, 2), initial=0) > TOLERANCE):
        if not windows[t]:
            raise ValueError(f'interaction {t} depends on the identities but is far from every position')
        grid = np.full((len(COVER),) * len(windows[t]) + (len(parameters),), np.nan)
        grid[tuple(digits[:, windows[t]].T)] = values[:, t]
        grid, axes = _dependent(grid)
        positions = [windows[t][j] for j in axes]
        if not positions or np.abs(grid[tuple(digits[:, positions].T)] - values[:, t]).max() > TOLERANCE:
            raise ValueError(f'interaction {t} depends on positions outside {windows[t]}')
        interactions.append((wild.loc[t, particles].to_numpy(dtype=int), positions, grid.reshape(-1, len(parameters))))
    return _identity_table(particles, parameters, interactions, alphabet,
                           circular=CIRCULAR.get(type(force).__name__, []))


def base_pair_identity_table(force, sites, alphabet):
    """The identity table of a BasePair force: every ordered pair of bases (donor d1 with sugar d2, acceptor a1 with
    sugar a2) far enough in sequence, of which at least one is the base of a position."""
    dna = force.dna
    bases = base_table(dna, sites).to_dict('index')
    cutoff = next(iter(force.forces.values())).getCutoffDistance().value_in_unit(unit.nanometer)
    definitions = dna.pair_definition[dna.pair_definition['DNA'] == dna.DNAtype]
    pairs = {(row.Base1, row.Base2): row for row in definitions.itertuples()}

    parameters = ['phi0', 'sigma', 't01', 't02', 'rng', 'epsilon', 'alpha', 'cutoff']
    interactions = []
    for d, a in itertools.permutations(bases, 2):
        if bases[d]['position'] < 0 and bases[a]['position'] < 0:
            continue
        if abs(bases[a]['seqid'] - bases[d]['seqid']) < force.min_seq_sep:
            continue
        positions, names = _assignments(bases, [d, a])
        rows = []
        for name_d, name_a in names:
            row = pairs.get((name_d, name_a))
            rows.append(None if row is None else [np.radians(row.torsion), row.sigma, np.radians(row.t1),
                                                  np.radians(row.t2), row.rang, row.epsilon, row.alpha, cutoff])
        particles = [bases[d]['base'], bases[d]['sugar'], bases[a]['base'], bases[a]['sugar']]
        interactions.append((particles, positions, rows))
    return _identity_table(['d1', 'd2', 'a1', 'a2'], parameters, interactions, alphabet, circular=['phi0'],
                           strength='epsilon')


def cross_stacking_identity_table(force, sites, alphabet):
    """The identity table of a CrossStacking force: every pair of a donor base d1 (sugar d2) and an acceptor base a1
    (sugar a2) with the neighbour a3 of the acceptor, far enough in sequence, of which d1, a1 or a3 is the base of a
    position. Branch 1 pairs a donor followed by a base with an acceptor preceded by a base; branch 2 the opposite."""
    dna = force.dna
    bases = base_table(dna, sites).to_dict('index')
    cutoff = next(iter(force.crossStackingForces.values()))[0].getCutoffDistance().value_in_unit(unit.nanometer)
    definitions = dna.cross_definition[dna.cross_definition['DNA'] == dna.DNAtype]
    crosses = {(row['Base_d1'], row['Base_a1'], row['Base_a3']): row for _, row in definitions.iterrows()}

    parameters = ['t03', 't0CS', 'rng_CS', 'rng_BP', 'epsilon', 'alpha', 'sigma', 'cutoff']
    interactions = []
    for branch, d_shift, a_shift in ((1, 1, -1), (2, -1, 1)):
        for d, a in itertools.permutations(bases, 2):
            d3, a3 = (d[0], d[1] + d_shift), (a[0], a[1] + a_shift)
            if d3 not in bases or a3 not in bases:
                continue
            if all(bases[r]['position'] < 0 for r in (d, a, a3)):
                continue
            if abs(bases[a]['seqid'] - bases[d]['seqid']) < force.min_seq_sep:
                continue
            donor = {bases[d]['base'], bases[d]['sugar'], bases[d3]['base']}
            if donor & {bases[a]['base'], bases[a]['sugar'], bases[a3]['base']}:
                continue
            positions, names = _assignments(bases, [d, a, a3])
            rows = []
            for name_d, name_a, name_a3 in names:
                if name_d != COMPLEMENT[name_a]:
                    rows.append(None)
                    continue
                row = crosses[(name_d, name_a, name_a3)]
                rows.append([np.radians(row['t03']), np.radians(row[f'T0CS_{branch}']), row[f'rng_cs{branch}'],
                             row['rng_bp'], row[f'eps_cs{branch}'], row[f'alpha_cs{branch}'], row[f'Sigma_{branch}'],
                             cutoff])
            particles = [bases[d]['base'], bases[d]['sugar'], bases[a]['base'], bases[a]['sugar'],
                         bases[a3]['base']]
            interactions.append((particles, positions, rows))
    return _identity_table(['d1', 'd2', 'a1', 'a2', 'a3'], parameters, interactions, alphabet, strength='epsilon')


def exclusion_identity_table(force, sites, alphabet):
    """The identity table of an Exclusion force: every pair of particles far enough in sequence of which at least one
    is the base of a position, with the combined epsilon (0 for a base pair) and sigma."""
    dna = force.dna
    bases = base_table(dna, sites).to_dict('index')
    cutoff = force.force.getCutoffDistance().value_in_unit(unit.nanometer)
    names = [force.force.getPerParticleParameterName(i) for i in range(force.force.getNumPerParticleParameters())]
    particles = pandas.DataFrame([force.force.getParticleParameters(i) for i in range(force.force.getNumParticles())],
                                 columns=names)
    definitions = dna.particle_definition[dna.particle_definition['DNA'] == dna.DNAtype].set_index('name')
    bead = {name: (definitions.at[name, 'epsilon'], definitions.at[name, 'radius'], dna_forces._base_type[name])
            for name in COMPLEMENT}
    excluded = {tuple(sorted(force.force.getExclusionParticles(i))) for i in range(force.force.getNumExclusions())}
    residue_of = {base['base']: residue for residue, base in bases.items() if base['position'] >= 0}
    seqid = particles['seqid'].to_numpy()
    fixed = particles[['epsilon', 'sigma', 'bt']].to_numpy()

    parameters = ['epsilon', 'sigma', 'cutoff']
    interactions = []
    for i, residue in residue_of.items():
        for j in range(len(particles)):
            if j == i or (j in residue_of and j < i) or tuple(sorted((i, j))) in excluded:
                continue
            if abs(seqid[i] - seqid[j]) < force.min_seq_sep:
                continue
            residues = [residue] if j not in residue_of else [residue, residue_of[j]]
            positions, assigned = _assignments(bases, residues)
            rows = []
            for names_ij in assigned:
                epsilon_i, sigma_i, type_i = bead[names_ij[0]]
                epsilon_j, sigma_j, type_j = bead[names_ij[1]] if j in residue_of else fixed[j]
                epsilon = np.sqrt(epsilon_i * epsilon_j) if type_i + type_j != 5 else 0.0
                rows.append([epsilon, 0.5 * (sigma_i + sigma_j), cutoff])
            interactions.append(([min(i, j), max(i, j)], positions, rows))
    return _identity_table(['aai', 'aaj'], parameters, interactions, alphabet, strength='epsilon')


def _windows(dna, sites, table):
    """The positions whose identity can change each bonded interaction of `table`: those within one residue of its
    particles, or only those at them when the former are more than four."""
    atoms = dna.atoms
    residue_of = dict(zip(atoms.index, zip(atoms['chainID'], atoms['resSeq'].astype(int))))
    position_of = {(chain, int(residue)): p for p, pair in enumerate(sites) for chain, residue in pair}
    windows = []
    for particles in table[[c for c in PARTICLES if c in table]].to_numpy():
        residues = {residue_of[i] for i in particles}
        for reach in (1, 0):
            positions = sorted({position_of[(chain, residue + d)] for chain, residue in residues
                                for d in range(-reach, reach + 1) if (chain, residue + d) in position_of})
            if len(positions) <= len(COVER):
                break
        else:
            raise ValueError(f'an interaction spans {len(positions)} positions; at most {len(COVER)} are handled')
        windows.append(positions)
    return windows


def _assignments(bases, residues):
    """(positions, names): the sorted positions of the bases `residues` of `bases` (a `base_table` as a dict of
    rows), and the names of those bases at every Watson-Crick assignment of identities to the positions
    (``itertools.product`` order)."""
    rows = [bases[r] for r in residues]
    positions = sorted({row['position'] for row in rows if row['position'] >= 0})
    names = [tuple(identities[positions.index(row['position'])][row['strand']] if row['position'] >= 0
                   else row['name'] for row in rows)
             for identities in itertools.product(COVER, repeat=len(positions))]
    return positions, names


def _dependent(grid):
    """(grid restricted to the axes along which it varies, those axes) for a grid (4,)*k + (parameters,)."""
    k = grid.ndim - 1
    axes = [j for j in range(k) if np.abs(grid - np.take(grid, [0], axis=j)).max() > TOLERANCE]
    return grid[tuple(slice(None) if j in axes else 0 for j in range(k))], axes


def _identity_table(particles, parameters, interactions, alphabet, circular=(), strength=None):
    """The identity table of `interactions` [(particle indices, positions, rows)], where rows holds the parameters at
    every Watson-Crick assignment of the positions, or None where the interaction is off. Interactions whose
    `strength` column is 0 at every assignment are left out; positions that do not change the parameters are
    dropped."""
    circular = [parameters.index(c) for c in circular]
    strength = None if strength is None else parameters.index(strength)
    kept = []
    for particle_indices, positions, rows in interactions:
        if strength is not None:
            on = [row for row in rows if row is not None and row[strength] != 0]
            if not on:
                continue
            off = list(on[0])
            off[strength] = 0.0
            rows = [off if row is None else row for row in rows]
        grid = np.asarray(rows, dtype=float).reshape((len(COVER),) * len(positions) + (len(parameters),))
        grid, axes = _dependent(grid)
        values = mix_parameters(grid, alphabet, circular, strength)
        kept.append((particle_indices, [positions[j] for j in axes], values))

    rows = [len(values) for _, _, values in kept]
    width = max((len(positions) for _, positions, _ in kept), default=0)
    position = np.full((sum(rows), width), -1)
    name = np.full((sum(rows), width), -1)
    start = 0
    for (_, positions, _), n in zip(kept, rows):
        position[start:start + n, :len(positions)] = positions
        assignments = list(itertools.product(range(len(alphabet)), repeat=len(positions)))
        name[start:start + n, :len(positions)] = np.array(assignments, dtype=int).reshape(n, len(positions))
        start += n
    particle_indices = np.array([p for p, _, _ in kept], dtype=int).reshape(-1, len(particles))
    table = pandas.DataFrame(np.repeat(particle_indices, rows, axis=0), columns=particles)
    table.insert(0, 'term', np.repeat(np.arange(len(kept)), rows))
    for j in range(width):
        table[f'position_{j}'] = position[:, j]
    for j in range(width):
        table[f'name_{j}'] = name[:, j]
    table[parameters] = np.concatenate([v for _, _, v in kept]) if kept else np.empty((0, len(parameters)))
    return table
