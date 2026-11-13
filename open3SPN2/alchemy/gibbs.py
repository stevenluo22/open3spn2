"""The Gibbs identity sampler: molecular dynamics of a duplex alternated with redraws of the identities of its
positions.

A move runs `steps_per_move` steps of the simulation in the current state, builds the model of the identity tables at
the positions the dynamics ended at with the offsets in force, records its window conditionals at the current state,
redraws the state from the model (``sweep`` or ``joint``) and sets it in the context. With a step `gamma` the offsets
are adapted after every move by self-adjusted mixture sampling from the single-position conditionals, with the step
gamma / max(1, t / 100) at the t-th move of the run.
"""
import numpy as np
import openmm.unit as unit

from .sampler import MAX_ENUMERATION, adapt_offsets, factor_tables, joint, model, sweep, window_conditionals
from .switch import set_state

MOVES = ('sweep', 'joint')


class GibbsSampler:
    """The identities of the positions of the identity tables `tables` in `simulation`, an OpenMM Simulation of a
    system built with ``tables.add_forces``.

    The states are drawn at `kT` in kJ/mol (default that of the integrator) with `offsets` (N, n) in kJ/mol (default
    zero), by `move` 'sweep' (blocked Gibbs in blocks of `block_size` positions) or 'joint' (exact, at most 3125
    states), with a generator seeded with `seed`. The runs record the window conditionals of the sizes `windows`
    (default 1 and 2, up to the number of positions). The context is set to `state` (default the wild type), which
    the sampler then keeps in ``state``.
    """

    def __init__(self, simulation, tables, kT=None, offsets=None, move='sweep', block_size=2, windows=None,
                 seed=None, state=None):
        N, n = len(tables.wild_type), len(tables.alphabet)
        if kT is None:
            if not hasattr(simulation.integrator, 'getTemperature'):
                raise ValueError('kT is needed with an integrator that has no temperature')
            kT = (unit.MOLAR_GAS_CONSTANT_R * simulation.integrator.getTemperature()).value_in_unit(
                unit.kilojoule_per_mole)
        if move not in MOVES:
            raise ValueError(f'move must be one of {MOVES}')
        if move == 'joint' and n ** N > MAX_ENUMERATION:
            raise ValueError(f'{n}^{N} states; the joint move enumerates at most {MAX_ENUMERATION}')
        windows = tuple(k for k in (1, 2) if k <= N) if windows is None else tuple(windows)
        if not all(1 <= k <= N for k in windows):
            raise ValueError(f'window sizes must be in [1, {N}]')
        offsets = np.zeros((N, n)) if offsets is None else np.array(offsets, dtype=float)
        if offsets.shape != (N, n):
            raise ValueError(f'offsets must have shape {(N, n)}')
        self.simulation = simulation
        self.tables = tables
        self.kT = float(kT)
        self.offsets = offsets
        self.move = move
        self.block_size = block_size
        self.windows = windows
        self.rng = np.random.default_rng(seed)
        self.state = tables.wild_type.copy() if state is None else np.array(state, dtype=int)
        set_state(tables, simulation.context, self.state)

    def run(self, n_moves, steps_per_move=25, gamma=None, target=None):
        return run(self, n_moves, steps_per_move, gamma, target)


def run(sampler, n_moves, steps_per_move=25, gamma=None, target=None):
    """Runs `n_moves` moves of `sampler`, adapting its offsets towards the occupancy `target` (n,) (default uniform)
    when `gamma` is given, and returns per move

    * ``state`` (n_moves, N): the state the dynamics ran in;
    * ``offsets`` (n_moves, N, n): the offsets in force;
    * ``window{k}`` (n_moves, N - k + 1) + (n,) * k: the window conditionals at that state, for k in
      ``sampler.windows``.
    """
    if n_moves < 1:
        raise ValueError('n_moves must be at least 1')
    tables, context = sampler.tables, sampler.simulation.context
    periodic = context.getSystem().usesPeriodicBoundaryConditions()
    records = {'state': [], 'offsets': [], **{f'window{k}': [] for k in sampler.windows}}
    for t in range(1, n_moves + 1):
        sampler.simulation.step(steps_per_move)
        positions = context.getState(getPositions=True)
        x = positions.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        box = positions.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer) if periodic else None
        current = model(factor_tables(tables, x, box), sampler.kT, sampler.offsets)
        windows = {k: window_conditionals(current, sampler.state, k) for k in sampler.windows}
        records['state'].append(sampler.state)
        records['offsets'].append(sampler.offsets)
        for k, window in windows.items():
            records[f'window{k}'].append(window)

        if sampler.move == 'sweep':
            new = sweep(current, sampler.state, sampler.rng, sampler.block_size)
        else:
            new = joint(current, sampler.rng)
        set_state(tables, context, new, sampler.state)
        if gamma is not None:
            singles = windows[1] if 1 in windows else window_conditionals(current, sampler.state, 1)
            sampler.offsets = adapt_offsets(sampler.offsets, singles, sampler.kT, gamma / max(1.0, t / 100), target)
        sampler.state = new
    return {key: np.array(values) for key, values in records.items()}
