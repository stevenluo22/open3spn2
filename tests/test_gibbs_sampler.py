"""Tests for the sampler and estimator of the Gibbs identity sampler.

At fixed coordinates the identity-dependent energy of a state is a sum of factors: ``factor_tables(tables,
coordinates)`` sums the interactions of the identity tables that depend on the same set of positions into one table
per set, with the energy at every assignment of names to its positions. ``model(factors, kT, offsets, frame)`` is the
distribution p(s) ~ exp(-(U_s + V_s) / kT) of the states of one frame, where the offsets V add offsets[p, s_p] for
every position. Its ``log_weights`` are the exponent, ``conditional`` and ``window_conditionals`` the distributions of
blocks of positions given the rest of a state, ``joint`` draws a state exactly and ``sweep`` redraws every position
once by blocked Gibbs. ``adapt_offsets`` is the self-adjusted mixture sampling update of the offsets towards a target
occupancy of the names (``w_target``). The estimator turns the window distributions of a run into sequence free
energies: ``markov_log_probability`` is ln p(s) of the order-(k-1) Markov chain with those windows, and
``free_energies`` removes the offsets. The behavioural requirements tested here:

* G1. For the identity tables of a duplex, the factors are one table per distinct set of positions of the
  interactions, with n^k entries per frame; a model built from those arrays alone gives, at every frame of a batch
  (with interactions in range at some frames only) and for any states (mixes included), the log-weights
  -(U + V) / kT, with U the sum over kinds of ``tables.state_energies``.
* G2. ``conditional`` of any block of positions (adjacent or not) and ``window_conditionals`` of windows of 1, 2 and 3
  consecutive positions equal the distributions obtained by enumerating every state.
* G3. The states drawn by ``joint`` and by ``sweep`` with blocks of 1 and 2 positions have the single-position and
  adjacent-pair marginals of p(s) (total variation 0.02 and 0.03 over 20000 draws), for a model whose pair
  marginals are far from independent; consecutive ``joint`` draws are independent. One ``sweep`` from a given state
  follows the Gibbs kernel: the blocks of the partition, shifted by an offset uniform in [0, block size), each redrawn
  in turn from its conditional (total variation 0.035 over 20000 sweeps). ``joint`` refuses models of more than 3125
  states.
* G4. ``w_target`` gives each identity (1 - w) / identities and each mix w / mixes. ``adapt_offsets`` adds
  gamma kT (singles / target - 1), with a uniform target by default. Iterated with the exact occupancies it brings
  them to a non-uniform target; driven by sweeps and the single-position conditionals with a decreasing step, the
  offsets reach occupancies within 0.05 of the target in 3000 sweeps.
* G5. With the exact window marginals of k positions (k = 1 to 4), ``markov_log_probability`` is ln p(s) of every
  state when every factor spans at most k consecutive positions (and not otherwise). ``free_energies`` of the exact ln p of a model with
  offsets gives U_s - U_wild, and refuses a wild type that is not among the states.
"""
import itertools
from pathlib import Path

import numpy as np
import pytest

from open3SPN2 import DNA
from open3SPN2.alchemy import (Alphabet, Factors, alphabet, adapt_offsets, factor_tables, free_energies,
                               identity_tables, markov_log_probability, model, sequence_offsets, w_target)
from open3SPN2.alchemy.alphabet import WATSON_CRICK

test_path = Path(__file__).parent
KT = 0.0083144626 * 300  # kJ/mol
N_NAMES = 5
N_POSITIONS = 4
CHAIN = [(0,), (1,), (2,), (3,), (0, 1), (1, 2), (2, 3)]
PAIRS = [(0, 1), (1, 2), (2, 3)]
SINGLES = [(0,), (1,), (2,), (3,)]
SPAN_3 = CHAIN + [(0, 1, 2), (1, 3)]
RICH = CHAIN + [(0, 3), (0, 2, 3)]


def _factors(sets, seed=0, scale=1.5, n=N_NAMES, N=N_POSITIONS):
    """Factors of one frame with random energies (in units of `scale` kT) for the position sets `sets`."""
    width = max(len(s) for s in sets)
    padded = np.full((len(sets), width), -1)
    for f, s in enumerate(sets):
        padded[f, :len(s)] = s
    ptr = np.concatenate([[0], np.cumsum([n ** len(s) for s in sets])])
    energies = np.random.default_rng(seed).normal(scale=scale * KT, size=(1, ptr[-1]))
    return Factors(padded, ptr, energies, n, N)


def _all_states(n=N_NAMES, N=N_POSITIONS):
    return np.array(list(itertools.product(range(n), repeat=N)))


def _exact_log_p(factors, offsets=None, kT=KT):
    """ln p(s) of every state (itertools.product order), summing the factor entries directly."""
    states = _all_states(factors.n_names, factors.n_positions)
    energy = np.zeros(len(states))
    for f in range(len(factors.sets)):
        positions = factors.sets[f][factors.sets[f] >= 0]
        code = np.zeros(len(states), dtype=int)
        for p in positions:
            code = code * factors.n_names + states[:, p]
        energy += factors.energies[0, factors.ptr[f] + code]
    if offsets is not None:
        energy += offsets[np.arange(factors.n_positions), states].sum(axis=1)
    log_weights = -energy / kT
    return log_weights - np.logaddexp.reduce(log_weights)


def _enumerated_conditional(factors, state, block):
    """(n,) * len(block) the distribution of the names of `block` (in that order) given the rest of `state`."""
    p = np.exp(_exact_log_p(factors))
    states = _all_states(factors.n_names, factors.n_positions)
    rest = [q for q in range(factors.n_positions) if q not in block]
    matching = (states[:, rest] == state[rest]).all(axis=1)
    expected = np.zeros((factors.n_names,) * len(block))
    for s, weight in zip(states[matching], p[matching]):
        expected[tuple(s[list(block)])] += weight
    return expected / expected.sum()


def _marginal(p, keep, n=N_NAMES, N=N_POSITIONS):
    """Marginal of the positions `keep` (sorted) of the joint probabilities `p` (product order), flattened."""
    joint = p.reshape((n,) * N)
    return joint.sum(axis=tuple(a for a in range(N) if a not in keep)).reshape(-1)


@pytest.fixture(scope='module')
def duplex_tables(tmp_path_factory):
    dna = DNA.fromXYZ(str(test_path / 'bdna_curv' / 'in00_conf.xyz'), 'B_curved', template_from_X3DNA=False,
                      output_pdb=str(tmp_path_factory.mktemp('sampler') / 'B_curved.pdb'))
    dna.periodic = False
    sites = [(('A', r), ('B', 65 - r)) for r in (12, 13, 16)]
    return identity_tables(dna, sites, alphabet('wc'))


def test_log_weights_match_state_energies(duplex_tables):
    """G1: the factors are one table of n^k entries per distinct set of positions of the identity tables, and a model
    of those bare arrays gives random states (with W), at each frame of a batch whose middle frame, compressed, brings
    interactions in range that the others do not, the log-weights -(U + V)/kT with U from ``tables.state_energies``."""
    tables = duplex_tables
    n, N = len(tables.alphabet), len(tables.sites)
    x = tables.dna.atoms[['x', 'y', 'z']].to_numpy(float) / 10
    compressed = x.mean(axis=0) + 0.95 * (x - x.mean(axis=0))
    frames = np.stack([x, compressed, x + np.random.default_rng(2).normal(scale=0.03, size=x.shape)])
    rng = np.random.default_rng(4)
    states = rng.integers(0, len(tables.alphabet), size=(40, len(tables.sites)))
    offsets = rng.normal(size=(len(tables.sites), len(tables.alphabet)))
    energies = sum(tables.state_energies(frames, states).values())
    factors = factor_tables(tables, frames)

    interaction_sets = set()
    for table in tables.kinds.values():
        positions = table.drop_duplicates('term')[[c for c in table if c.startswith('position_')]].to_numpy(dtype=int)
        interaction_sets |= {tuple(p[p >= 0]) for p in positions}
    sets = [tuple(s[s >= 0]) for s in factors.sets]
    assert len(sets) == len(set(sets)) and set(sets) == interaction_sets
    np.testing.assert_array_equal(np.diff(factors.ptr), [n ** len(s) for s in sets])
    assert factors.energies.shape == (len(frames), factors.ptr[-1])

    bare = Factors(factors.sets.copy(), factors.ptr.copy(), factors.energies.copy(), n, N)
    for frame in range(len(frames)):
        log_weights = model(bare, KT, offsets, frame).log_weights(states)
        expected = -(energies[frame] + offsets[np.arange(len(tables.sites)), states].sum(axis=1)) / KT
        np.testing.assert_allclose(log_weights, expected, rtol=0, atol=1e-8)


@pytest.mark.parametrize('block', [(1,), (0, 1), (0, 2), (3, 1), (0, 1, 3)])
def test_conditional_matches_enumeration(block):
    """G2: the conditional of a block (in the order given) given the rest of a state is the enumerated one."""
    factors = _factors(RICH)
    state = np.array([2, 4, 0, 3])
    expected = _enumerated_conditional(factors, state, block)
    conditional = model(factors, KT).conditional(state, block)
    np.testing.assert_allclose(conditional.reshape(expected.shape), expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize('size', [1, 2, 3])
def test_window_conditionals_match_enumeration(size):
    """G2: the conditional of every window of `size` consecutive positions given the rest of a state is the
    enumerated one."""
    factors = _factors(RICH, seed=1)
    state = np.array([1, 0, 4, 2])
    windows = model(factors, KT).window_conditionals(state, size)
    assert windows.shape == (N_POSITIONS - size + 1,) + (N_NAMES,) * size
    for q in range(N_POSITIONS - size + 1):
        expected = _enumerated_conditional(factors, state, tuple(range(q, q + size)))
        np.testing.assert_allclose(windows[q], expected, rtol=0, atol=1e-12)


def _total_variation(a, b):
    return 0.5 * np.abs(a - b).sum()


@pytest.mark.parametrize('move', ['joint', 'sweep1', 'sweep2'])
def test_draws_have_the_exact_marginals(move):
    """G3: 20000 draws of `move` have the single and adjacent-pair marginals of p(s) of a strongly coupled model."""
    factors = _factors(PAIRS, seed=0, scale=2.0)
    p = np.exp(_exact_log_p(factors))
    pairs = [_marginal(p, (q, q + 1)) for q in range(N_POSITIONS - 1)]
    singles = [_marginal(p, (q,)) for q in range(N_POSITIONS)]
    for q, pair in enumerate(pairs):
        assert _total_variation(pair, np.outer(singles[q], singles[q + 1]).reshape(-1)) > 0.1

    m = model(factors, KT)
    rng = np.random.default_rng(8)
    draws = []
    state = np.zeros(N_POSITIONS, dtype=int)
    for i in range(20100):
        if move == 'joint':
            state = m.joint(rng)
        else:
            state = m.sweep(state, rng, block_size=int(move[-1]))
        if move == 'joint' or i >= 100:
            draws.append(state)
    draws = np.array(draws[-20000:])

    for q in range(N_POSITIONS):
        counts = np.bincount(draws[:, q], minlength=N_NAMES) / len(draws)
        assert _total_variation(counts, singles[q]) < 0.02
    for q in range(N_POSITIONS - 1):
        counts = np.bincount(draws[:, q] * N_NAMES + draws[:, q + 1], minlength=N_NAMES ** 2) / len(draws)
        assert _total_variation(counts, pairs[q]) < 0.03
    if move == 'joint':
        for q in range(N_POSITIONS):
            indicator = (draws[:, q] == np.argmax(singles[q])).astype(float)
            assert abs(np.corrcoef(indicator[:-1], indicator[1:])[0, 1]) < 0.03


def _sweep_kernel(factors, state, blocks):
    """(n^N,) the distribution of the state after redrawing `blocks` of `state` in turn from their conditionals."""
    n, N = factors.n_names, factors.n_positions
    states = _all_states(n, N)
    distribution = (states == state).all(axis=1).astype(float)
    for block in blocks:
        moved = np.zeros_like(distribution)
        for s in np.flatnonzero(distribution):
            for assignment, weight in zip(itertools.product(range(n), repeat=len(block)),
                                          _enumerated_conditional(factors, states[s], block).reshape(-1)):
                target = states[s].copy()
                target[list(block)] = assignment
                moved[np.flatnonzero((states == target).all(axis=1))[0]] += distribution[s] * weight
        distribution = moved
    return distribution


@pytest.mark.parametrize('block_size', [1, 2])
def test_sweep_follows_the_gibbs_kernel(block_size):
    """G3: the states after one sweep from a fixed state follow the blocked Gibbs kernel averaged over the shifts of
    the partition, which differs from p(s) and, for blocks of 2, from the kernel of an unshifted partition."""
    n, N = 3, 3
    factors = _factors([(0, 1), (1, 2), (0, 2)], seed=2, scale=2.0, n=n, N=N)
    start = np.array([2, 0, 1])
    partitions = {1: [[(0,), (1,), (2,)]], 2: [[(0, 1), (2,)], [(0,), (1, 2)]]}[block_size]
    kernels = [_sweep_kernel(factors, start, blocks) for blocks in partitions]
    expected = np.mean(kernels, axis=0)
    assert _total_variation(expected, np.exp(_exact_log_p(factors))) > 0.1
    if block_size == 2:
        assert _total_variation(expected, kernels[0]) > 0.05

    m = model(factors, KT)
    rng = np.random.default_rng(12)
    codes = [int(np.ravel_multi_index(tuple(m.sweep(start, rng, block_size)), (n,) * N)) for _ in range(20000)]
    assert _total_variation(np.bincount(codes, minlength=n ** N) / len(codes), expected) < 0.035


def test_joint_refuses_large_models():
    """G3: enumerating more than 3125 states is refused."""
    factors = _factors([(p,) for p in range(6)], N=6)
    with pytest.raises(ValueError):
        model(factors, KT).joint(np.random.default_rng(0))


def test_w_target():
    """G4: the target occupancy splits w among the mixes and 1 - w among the identities; w outside (0, 1) and an
    alphabet without mixes are refused."""
    np.testing.assert_allclose(w_target(alphabet('wc'), 0.1), [0.225, 0.225, 0.225, 0.225, 0.1])
    two_mixes = Alphabet(WATSON_CRICK, {'W': ('AT', 'TA', 0.5), 'S': ('GC', 'CG', 0.5)})
    np.testing.assert_allclose(w_target(two_mixes, 0.3), [0.175, 0.175, 0.175, 0.175, 0.15, 0.15])
    three = Alphabet({'AT': ('A', 'T'), 'TA': ('T', 'A'), 'GC': ('G', 'C')}, {'W': ('AT', 'TA', 0.5)})
    np.testing.assert_allclose(w_target(three, 0.4), [0.2, 0.2, 0.2, 0.4])
    for w in (0.0, 1.0, 1.5):
        with pytest.raises(ValueError):
            w_target(alphabet('wc'), w)
    with pytest.raises(ValueError):
        w_target(alphabet('wc4'), 0.1)


def _occupancy(factors, offsets):
    p = np.exp(_exact_log_p(factors, offsets))
    N = factors.n_positions
    return np.array([_marginal(p, (q,), N=N) for q in range(N)])


def test_offsets_reach_the_target():
    """G4: the update with the default (uniform) target is gamma kT (singles n - 1); iterating adapt_offsets with the
    exact occupancies converges to a non-uniform target; driven by sweeps and the
    single-position conditionals with the step 0.1 / max(1, t / 100), it gets within 0.05 of it in 3000 sweeps (the
    noise of the stochastic updates is 0.01-0.03 at that length)."""
    singles = np.random.default_rng(3).dirichlet(np.ones(N_NAMES), size=N_POSITIONS)
    offsets = np.random.default_rng(4).normal(size=(N_POSITIONS, N_NAMES))
    np.testing.assert_allclose(adapt_offsets(offsets, singles, KT, 0.3),
                               offsets + 0.3 * KT * (singles * N_NAMES - 1), rtol=0, atol=1e-12)

    factors = _factors(CHAIN, seed=5, scale=2.0)
    target = w_target(alphabet('wc'), 0.1)

    offsets = np.zeros((N_POSITIONS, N_NAMES))
    for _ in range(3000):
        offsets = adapt_offsets(offsets, _occupancy(factors, offsets), KT, 0.1, target)
    np.testing.assert_allclose(_occupancy(factors, offsets), np.tile(target, (N_POSITIONS, 1)), atol=1e-3)

    factors = _factors([(0,), (1,), (0, 1)], seed=1, scale=2.0, N=2)
    offsets = np.zeros((2, N_NAMES))
    rng = np.random.default_rng(101)
    state = np.zeros(2, dtype=int)
    for t in range(1, 3001):
        m = model(factors, KT, offsets)
        state = m.sweep(state, rng)
        offsets = adapt_offsets(offsets, m.window_conditionals(state, 1), KT, 0.1 / max(1.0, t / 100), target)
    assert np.abs(_occupancy(factors, offsets) - target).max() < 0.05


@pytest.mark.parametrize('sets, k, exact', [(SINGLES, 1, True), (CHAIN, 1, False), (CHAIN, 2, True), (SPAN_3, 3, True),
                                            (SPAN_3, 2, False), (RICH, 4, True)])
def test_markov_estimator(sets, k, exact):
    """G5: with exact window marginals of k positions, the Markov chain gives ln p(s) of every state exactly when
    every factor spans at most k consecutive positions."""
    factors = _factors(sets, seed=7)
    log_p = _exact_log_p(factors)
    p = np.exp(log_p)
    windows = np.array([_marginal(p, tuple(range(q, q + k))).reshape((N_NAMES,) * k)
                        for q in range(N_POSITIONS - k + 1)])
    estimate = markov_log_probability(windows, _all_states())
    if exact:
        np.testing.assert_allclose(estimate, log_p, rtol=0, atol=1e-9)
    else:
        assert np.abs(estimate - log_p).max() > 1e-3


def test_free_energies_remove_the_offsets():
    """G5: the free energies of the exact ln p of a model with offsets are the energy differences to the wild type;
    the offsets of a state are the sum of its single offsets; a wild type outside the states is refused."""
    factors = _factors(CHAIN, seed=9)
    offsets = np.random.default_rng(10).normal(scale=KT, size=(N_POSITIONS, N_NAMES))
    states = _all_states()
    np.testing.assert_allclose(sequence_offsets(states, offsets),
                               offsets[np.arange(N_POSITIONS), states].sum(axis=1), rtol=0, atol=1e-12)

    log_p = _exact_log_p(factors, offsets)
    energies = -KT * (_exact_log_p(factors) - _exact_log_p(factors)[0])
    wild = np.array([1, 2, 3, 0])
    w = np.flatnonzero((states == wild).all(axis=1))[0]
    np.testing.assert_allclose(free_energies(log_p, states, wild, offsets, KT), energies - energies[w], rtol=0,
                               atol=1e-9)
    with pytest.raises(ValueError):
        free_energies(log_p[1:], states[1:], states[0], offsets, KT)
