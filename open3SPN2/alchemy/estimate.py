"""Sequence free energies from the window distributions of a run, with the offsets of the run removed.

The distributions of the names of every window of k consecutive positions, averaged over a run (for instance of
``window_conditionals``), estimate the window marginals of p(s) ~ exp(-(U_s + V_s) / kT). The order-(k - 1) Markov
chain with those marginals gives ln p(s) of any sequence, exactly when every factor of the energy spans at most k
consecutive positions, and F_s - F_wild = -kT (ln p(s) - ln p(wild)) - (V_s - V_wild).
"""
import numpy as np


def markov_log_probability(windows, states):
    """(S,) ln p(s) of the states (S, N) under the order-(k - 1) Markov chain with window marginals `windows`
    (N - k + 1,) + (n,) * k:

        ln P_0(s_0..s_k-1) + sum_{q >= 1} [ln P_q(s_q..s_q+k-1) - ln sum_b P_q(s_q..s_q+k-2, b)]
    """
    windows, states = np.asarray(windows, dtype=float), np.atleast_2d(states)
    k, count = windows.ndim - 1, windows.shape[0]
    with np.errstate(divide='ignore'):
        log_windows, log_overlaps = np.log(windows), np.log(windows.sum(axis=-1))
    q = np.arange(count)[None]
    names = [states[:, j:j + count] for j in range(k)]
    joint = log_windows[(q, *names)]
    overlap = log_overlaps[(q[:, 1:], *(name[:, 1:] for name in names[:-1]))]
    with np.errstate(invalid='ignore'):
        steps = np.where(np.isneginf(joint[:, 1:]), -np.inf, joint[:, 1:] - overlap)
    return joint[:, 0] + steps.sum(axis=1)


def sequence_offsets(states, offsets):
    """(S,) V_s = sum_p offsets[p, s_p] of the states (S, N), with `offsets` (N, n)."""
    states = np.atleast_2d(states)
    return np.asarray(offsets, dtype=float)[np.arange(states.shape[1])[None], states].sum(axis=1)


def free_energies(log_p, states, wild, offsets, kT):
    """(S,) F_s - F_wild = -kT (log_p - log_p[wild]) - (V_s - V_wild) of the states (S, N) from their log-probabilities
    `log_p` in a run with `offsets` (N, n) at temperature `kT`; `wild` (N,) must be one of the states."""
    states, log_p = np.atleast_2d(states), np.asarray(log_p, dtype=float)
    match = np.flatnonzero((states == np.asarray(wild)).all(axis=1))
    if not len(match):
        raise ValueError('the wild type must be one of the states')
    v = sequence_offsets(states, offsets)
    return -kT * (log_p - log_p[match[0]]) - (v - v[match[0]])
