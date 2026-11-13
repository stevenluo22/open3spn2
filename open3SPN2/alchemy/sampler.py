"""Sampler of the identities of the positions at fixed coordinates.

At fixed coordinates the identity-dependent energy of a state s (the name index of every position) is a sum of
factors: factor f depends on the names at the sorted positions ``sets[f]`` (padded with -1) and holds the energy of
each of the n^k assignments of names to them, in ``itertools.product`` order (first position slowest), at
``ptr[f]:ptr[f + 1]``. A model is the distribution p(s) ~ exp(-(U_s + V_s) / kT) of one frame, where the offsets
V_s = sum_p offsets[p, s_p] bias the occupancy of the names; its factors hold the log-weights -(U_f + V_f) / kT, with
the offsets in single-position factors.

Every operation is built on ``log_weights(model, states)``: a conditional enumerates the assignments of a block with
the rest of the state fixed, summing only the factors that touch the block. Each function taking a model is also a
method of `Model`.
"""
import functools
import itertools

import numpy as np

from .tables import ENERGIES, GEOMETRIES

PARTICLES = ('aai', 'aaj', 'aak', 'aal', 'd1', 'd2', 'a1', 'a2', 'a3')
# {kind: (distance, bound)}: a row has no energy where the geometry `distance` of its interaction exceeds `bound`.
SCREENS = {'BasePair': ('r', 'cutoff'), 'CrossStacking': ('r_pair', 'cutoff'), 'Exclusion': ('r', 'sigma')}

MAX_ENUMERATION = 3125


class Factors:
    """The factors of the identity-dependent energy at a batch of frames: `sets` (F, w) sorted positions padded with
    -1, `ptr` (F + 1,) and `energies` (frames, ptr[-1]) in kJ/mol, for `n_names` names at `n_positions` positions."""

    def __init__(self, sets, ptr, energies, n_names, n_positions):
        self.sets = np.asarray(sets, dtype=int).reshape(len(ptr) - 1, -1)
        self.ptr = np.asarray(ptr, dtype=int)
        self.energies = np.atleast_2d(np.asarray(energies, dtype=float))
        self.n_names = int(n_names)
        self.n_positions = int(n_positions)


class Model:
    """The log-weights of the states of one frame as factors: `sets` (F, w), `ptr` (F + 1,) and `factor_log_weights`
    (ptr[-1],), for `n_names` names at `n_positions` positions."""

    def __init__(self, sets, ptr, factor_log_weights, n_names, n_positions):
        self.sets = np.asarray(sets, dtype=int).reshape(len(ptr) - 1, -1)
        self.ptr = np.asarray(ptr, dtype=int)
        self.factor_log_weights = np.asarray(factor_log_weights, dtype=float)
        self.n_names = int(n_names)
        self.n_positions = int(n_positions)
        valid = self.sets >= 0
        exponent = valid.sum(axis=1, keepdims=True) - 1 - np.arange(self.sets.shape[1])
        self.powers = np.where(valid, self.n_names ** np.maximum(exponent, 0), 0)

    def log_weights(self, states):
        return log_weights(self, states)

    def conditional(self, state, block):
        return conditional(self, state, block)

    def window_conditionals(self, state, size):
        return window_conditionals(self, state, size)

    def joint(self, rng):
        return joint(self, rng)

    def sweep(self, state, rng, block_size=2):
        return sweep(self, state, rng, block_size)


class FactorLayout:
    """Where the rows of identity tables go in their factors: `sets` (F, w) sorted positions padded with -1, `ptr`
    (F + 1,), and per kind in `kinds` {kind: (terms, cells, cell_ptr, row_column)}:

    * ``terms``: the particles of every interaction, in term order, as a dict of column arrays;
    * ``cells``: the term and the parameters of the distinct rows (interaction and parameters) whose energy can be
      non-zero, i.e. with epsilon != 0, as a dict of column arrays;
    * ``cell_ptr`` (C + 1,), ``row_column``: the factor entries of the rows of cell c are
      ``row_column[cell_ptr[c]:cell_ptr[c + 1]]``.
    """

    def __init__(self, sets, ptr, kinds, n_names, n_positions):
        self.sets = sets
        self.ptr = ptr
        self.kinds = kinds
        self.n_names = n_names
        self.n_positions = n_positions


def factor_layout(tables):
    """The `FactorLayout` of the identity tables `tables` (cached as ``tables.layout``): the interactions that depend
    on the same set of positions share a factor."""
    n = len(tables.alphabet)
    term_sets = {}
    for kind, table in tables.kinds.items():
        first = np.flatnonzero(np.diff(table['term'].to_numpy(), prepend=-1))
        positions = table[[c for c in table if c.startswith('position_')]].to_numpy(dtype=int)[first]
        term_sets[kind] = first, [tuple(p[p >= 0]) for p in positions]
    sets = sorted({s for _, kind_sets in term_sets.values() for s in kind_sets}, key=lambda s: (len(s), s))
    index = {s: f for f, s in enumerate(sets)}
    ptr = np.concatenate([[0], np.cumsum([n ** len(s) for s in sets], dtype=int)])

    kinds = {}
    for kind, table in tables.kinds.items():
        if not len(table):
            continue
        first, kind_sets = term_sets[kind]
        term = table['term'].to_numpy()
        column = ptr[np.array([index[s] for s in kind_sets], dtype=int)[term]] + np.arange(len(table)) - first[term]
        key = [c for c in table if not c.startswith(('position_', 'name_'))]
        row_cell = table.groupby(key, sort=False).ngroup().to_numpy()
        cells = table.drop_duplicates(key)
        live = cells['epsilon'].to_numpy() != 0 if 'epsilon' in cells else np.ones(len(cells), dtype=bool)
        live_rows = np.flatnonzero(live[row_cell])
        live_cell = (np.cumsum(live) - 1)[row_cell[live_rows]]
        order = np.argsort(live_cell, kind='stable')
        cell_ptr = np.concatenate([[0], np.cumsum(np.bincount(live_cell, minlength=live.sum()))])
        particles = [c for c in key if c in PARTICLES]
        kinds[kind] = ({c: table[c].to_numpy()[first] for c in particles},
                       {c: cells[c].to_numpy()[live] for c in key if c not in particles}, cell_ptr,
                       column[live_rows[order]])

    padded = np.full((len(sets), max((len(s) for s in sets), default=0)), -1)
    for f, s in enumerate(sets):
        padded[f, :len(s)] = s
    return FactorLayout(padded, ptr, kinds, n, len(tables.sites))


def factor_tables(tables, coordinates, box=None):
    """The `Factors` of the identity tables `tables` at `coordinates` (N, 3) or (frames, N, 3) in nm: the energies of
    the interactions that depend on the same set of positions summed into one factor.

    The geometry of every interaction is computed once, the energy once per distinct parameter row with epsilon != 0,
    and, for the nonbonded kinds, only for the rows whose particles are within the range of the potential (``SCREENS``)
    at some frame; every other entry is exactly 0."""
    layout = tables.layout
    frames = np.asarray(coordinates, dtype=float).reshape(-1, len(tables.dna.atoms), 3)
    columns, weights = [], []
    for kind, (terms, cells, cell_ptr, row_column) in layout.kinds.items():
        geometry = GEOMETRIES[kind](terms, frames, box)
        cell_term = cells['term']
        if kind in SCREENS:
            distance, bound = SCREENS[kind]
            near = np.flatnonzero((geometry[distance][:, cell_term] <= cells[bound]).any(axis=0))
        else:
            near = np.arange(len(cell_term))
        energies = ENERGIES[kind]({c: v[near] for c, v in cells.items()}, None,
                                  geometry={k: np.take(g, cell_term[near], axis=-1) for k, g in geometry.items()})
        counts = cell_ptr[near + 1] - cell_ptr[near]
        rows = np.repeat(cell_ptr[near] - np.cumsum(counts) + counts, counts) + np.arange(counts.sum())
        columns.append(row_column[rows])
        weights.append(np.repeat(energies, counts, axis=1))
    columns = np.concatenate(columns) if columns else np.zeros(0, dtype=int)
    weights = np.concatenate(weights, axis=1) if weights else np.zeros((len(frames), 0))
    energies = np.stack([np.bincount(columns, w, minlength=layout.ptr[-1]) for w in weights])
    return Factors(layout.sets.copy(), layout.ptr.copy(), energies, layout.n_names, layout.n_positions)


def model(factors, kT, offsets=None, frame=0):
    """The `Model` of frame `frame` of `factors` at temperature `kT` (kJ/mol), with the offsets (N, n) in kJ/mol
    added to the single-position factors (created where there are none)."""
    if kT <= 0:
        raise ValueError('kT must be positive')
    n, N = factors.n_names, factors.n_positions
    sets, ptr, energies = factors.sets, factors.ptr, factors.energies[frame].copy()
    if offsets is not None:
        offsets = np.asarray(offsets, dtype=float)
        if offsets.shape != (N, n):
            raise ValueError(f'offsets must have shape {(N, n)}')
        single = {int(s[0]): f for f, s in enumerate(sets) if (s >= 0).sum() == 1}
        missing = [p for p in range(N) if p not in single]
        width = max(sets.shape[1], 1)
        sets = np.vstack([np.pad(sets, ((0, 0), (0, width - sets.shape[1])), constant_values=-1),
                          np.pad(np.array(missing, dtype=int).reshape(-1, 1), ((0, 0), (0, width - 1)),
                                 constant_values=-1)])
        single.update({p: len(ptr) - 1 + j for j, p in enumerate(missing)})
        ptr = np.concatenate([ptr, ptr[-1] + n * np.arange(1, len(missing) + 1)])
        energies = np.concatenate([energies, np.zeros(n * len(missing))])
        for p in range(N):
            energies[ptr[single[p]]:ptr[single[p]] + n] += offsets[p]
    return Model(sets, ptr, -energies / kT, n, N)


def log_weights(model, states):
    """(S,) the log-weights -(U_s + V_s) / kT of the states (S, N) or of one state (N,)."""
    return _log_weights(model, np.atleast_2d(states), np.arange(len(model.sets)))


def conditional(model, state, block):
    """(n^len(block),) the probabilities of the assignments of names to the positions `block` (in that order,
    itertools.product order) given the names of the other positions in `state`."""
    block = [int(p) for p in block]
    if not block or len(set(block)) != len(block) or not all(0 <= p < model.n_positions for p in block):
        raise ValueError(f'a block is distinct positions in [0, {model.n_positions})')
    candidates = np.repeat(_checked(model, state)[None], model.n_names ** len(block), axis=0)
    candidates[:, block] = _assignments(model.n_names, len(block))
    touching = np.flatnonzero(np.isin(model.sets, block).any(axis=1))
    return _probabilities(_log_weights(model, candidates, touching))


def window_conditionals(model, state, size):
    """(N - size + 1,) + (n,) * size: entry [q, a_0, ..., a_size-1] is the probability that positions q..q+size-1
    hold a_0..a_size-1 given the rest of `state`."""
    if not 1 <= size <= model.n_positions:
        raise ValueError(f'size must be in [1, {model.n_positions}]')
    shape = (model.n_names,) * size
    return np.stack([conditional(model, state, range(q, q + size)).reshape(shape)
                     for q in range(model.n_positions - size + 1)])


def joint(model, rng):
    """A state drawn exactly from p(s) by enumerating every state (at most 3125)."""
    if model.n_names ** model.n_positions > MAX_ENUMERATION:
        raise ValueError(f'{model.n_names}^{model.n_positions} states; at most {MAX_ENUMERATION} are enumerated')
    states = _assignments(model.n_names, model.n_positions)
    p = _probabilities(log_weights(model, states))
    return states[rng.choice(len(p), p=p)].copy()


def sweep(model, state, rng, block_size=2):
    """The state after redrawing every position once from its conditional, in consecutive blocks of `block_size`
    positions whose partition is shifted by a random offset."""
    if block_size < 1:
        raise ValueError('block_size must be at least 1')
    state = _checked(model, state).copy()
    N = model.n_positions
    shift = int(rng.integers(block_size))
    edges = sorted({0, N, *range(shift, N, block_size)})
    for start, stop in zip(edges[:-1], edges[1:]):
        p = conditional(model, state, range(start, stop))
        state[start:stop] = _assignments(model.n_names, stop - start)[rng.choice(len(p), p=p)]
    return state


def adapt_offsets(offsets, singles, kT, gamma, target=None):
    """One self-adjusted mixture sampling update of `offsets` (N, n) from `singles` (N, n), the distributions of the
    names at each position: offsets + gamma kT (singles / target - 1), which raises the offsets of the names that
    occur more than their occupancy `target` (n,) (default uniform)."""
    singles = np.asarray(singles, dtype=float)
    target = 1.0 / singles.shape[-1] if target is None else np.asarray(target, dtype=float)
    return np.asarray(offsets, dtype=float) + gamma * kT * (singles / target - 1)


def w_target(alphabet, w):
    """(n,) the occupancy that gives the mixes of `alphabet` a total of `w` and the identities 1 - w, evenly."""
    if not 0 < w < 1 or not alphabet.mixes:
        raise ValueError('w must be in (0, 1) and the alphabet needs at least one mix')
    return np.array([w / len(alphabet.mixes) if name in alphabet.mixes else (1 - w) / len(alphabet.identities)
                     for name in alphabet.names])


def _log_weights(model, states, factors):
    """(S,) the sum of the log-weights of `factors` at the states (S, N)."""
    padded_states = np.concatenate([states, np.zeros((len(states), 1), dtype=int)], axis=1)
    sets = np.where(model.sets[factors] >= 0, model.sets[factors], model.n_positions)
    codes = (padded_states[:, sets] * model.powers[factors]).sum(axis=-1)
    return model.factor_log_weights[model.ptr[factors] + codes].sum(axis=-1)


def _checked(model, state):
    state = np.asarray(state, dtype=int)
    if state.shape != (model.n_positions,):
        raise ValueError(f'a state is {model.n_positions} name indices')
    return state


@functools.lru_cache(maxsize=None)
def _assignments(n, k):
    """(n^k, k) every assignment of n names to k positions, the first slowest (read-only)."""
    assignments = np.array(list(itertools.product(range(n), repeat=k)), dtype=int).reshape(n ** k, k)
    assignments.setflags(write=False)
    return assignments


def _probabilities(log_weights):
    weights = np.exp(log_weights - log_weights.max())
    return weights / weights.sum()
