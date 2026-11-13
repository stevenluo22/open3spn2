"""Tests for the driver of the Gibbs identity sampler.

``GibbsSampler(simulation, tables, ...)`` samples the identities of the positions of the identity tables in an OpenMM
Simulation of a system built with ``tables.add_forces``. ``run(sampler, n_moves, steps_per_move, gamma, target)`` (or
``sampler.run(...)``) alternates dynamics and redraws: each move runs the dynamics in the current state, builds the
model of the positions the dynamics ended at with the offsets in force, records its window conditionals at the
current state, redraws the state and sets it in the context; with `gamma` it adapts the offsets. The behavioural
requirements tested here:

* D1. A run of m moves advances the simulation by m times the steps per move and returns per move the state the
  dynamics ran in (the sampler's state before the move, the previous move's draw after the first), the offsets in
  force and the window conditionals of the requested sizes (by default 1 and 2), each a distribution. After every
  move the context holds the sampler's state.
* D2. The records are the window conditionals, at the recorded state, of the model of the positions the dynamics
  ended at, at kT (by default R T of the integrator) with the offsets in force; with periodic boundary conditions,
  with the box of the context.
* D3. At fixed positions (moves without steps), the states drawn by the 'joint' and 'sweep' moves follow p(s) of the
  model with the offsets: the single-position marginals of 1000 moves are within total variation 0.06 of the exact
  ones, for offsets that change them by far more. The context holds the sampler's state as the positions revisit
  earlier names.
* D4. With a step gamma, the offsets after the t-th move of a run are ``adapt_offsets`` of those in force with the
  single-position conditionals of the move, the step gamma / max(1, t / 100) and the target; without it they stay.
* D5. Runs with the same seeds of the integrator and of the sampler are identical.
* D6. Unknown moves, window sizes beyond the positions, offsets of the wrong shape, a missing kT with an integrator
  without temperature and runs of no moves are refused.
"""
import itertools
from pathlib import Path

import numpy as np
import openmm
import openmm.unit as unit
import pytest

import open3SPN2
from open3SPN2 import DNA
from open3SPN2.alchemy import GibbsSampler, adapt_offsets, alphabet, factor_tables, identity_tables, model, \
    w_target, window_conditionals

test_path = Path(__file__).parent
TEMPERATURE = 300  # K
KT = (unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE * unit.kelvin).value_in_unit(unit.kilojoule_per_mole)
BOX = np.array([13.7, 12.9, 14.3])  # nm


@pytest.fixture(scope='module')
def tables(tmp_path_factory):
    dna = DNA.fromXYZ(str(test_path / 'bdna_curv' / 'in00_conf.xyz'), 'B_curved', template_from_X3DNA=False,
                      output_pdb=str(tmp_path_factory.mktemp('gibbs_run') / 'B_curved.pdb'))
    dna.periodic = False
    sites = [(('A', r), ('B', 65 - r)) for r in (12, 13, 16)]
    return identity_tables(dna, sites, alphabet('wc'))


def _simulation(tables, seed=1, integrator=None, box=None):
    """A Reference-platform Simulation of the switched system of `tables`, with a Langevin integrator at TEMPERATURE
    seeded with `seed` unless `integrator` is given, periodic with the orthorhombic `box` (nm) if given."""
    system = open3SPN2.System(tables.dna)
    system.clearForces()
    tables.dna.periodic = box is not None
    try:
        tables.add_forces(system)
    finally:
        tables.dna.periodic = False
    if box is not None:
        system.setDefaultPeriodicBoxVectors(*np.diag(box))
    if integrator is None:
        integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE * unit.kelvin, 2 / unit.picosecond,
                                                     2 * unit.femtosecond)
        integrator.setRandomNumberSeed(seed)
    simulation = openmm.app.Simulation(system.top, system._wrapped_system, integrator,
                                       openmm.Platform.getPlatformByName('Reference'))
    x = np.asarray(system.coord.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    if box is not None:
        x = (x - x.mean(axis=0) + box) % box
    simulation.context.setPositions(x)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin, seed)
    return simulation


def _positions(simulation):
    state = simulation.context.getState(getPositions=True)
    return (state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
            state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))


def _context_state(simulation, n_positions):
    return np.array([simulation.context.getParameter(f'state_{p}') for p in range(n_positions)])


def _offsets(tables, scale, seed=4):
    return np.random.default_rng(seed).normal(scale=scale, size=(len(tables.wild_type), len(tables.alphabet)))


def test_run_records_and_context(tables):
    """D1: the run advances the dynamics by moves times steps, returns per move the state the dynamics ran in, the
    offsets in force and window distributions of the requested sizes, and leaves the context in the sampler's
    state after every move."""
    N, n = len(tables.wild_type), len(tables.alphabet)
    simulation = _simulation(tables)
    start = (tables.wild_type + 1) % 4
    sampler = GibbsSampler(simulation, tables, windows=(1, 2, 3), seed=0, state=start)
    np.testing.assert_array_equal(_context_state(simulation, N), start)

    states = [start]
    for steps in (3, 0, 2):
        before = sampler.state.copy()
        records = sampler.run(4, steps_per_move=steps)
        assert simulation.currentStep == simulation.context.getStepCount()
        np.testing.assert_array_equal(records['state'][0], before)
        np.testing.assert_array_equal(_context_state(simulation, N), sampler.state)
        states += list(records['state']) + [sampler.state]
        assert records['state'].shape == (4, N)
        assert records['offsets'].shape == (4, N, n)
        for k in (1, 2, 3):
            window = records[f'window{k}']
            assert window.shape == (4, N - k + 1) + (n,) * k
            assert (window >= 0).all()
            np.testing.assert_allclose(window.sum(axis=tuple(range(2, 2 + k))), 1.0, rtol=0, atol=1e-12)
    assert simulation.currentStep == 4 * (3 + 0 + 2)
    assert len({tuple(s) for s in states}) > 1


@pytest.mark.parametrize('periodic', [False, True])
def test_records_are_the_model_conditionals(tables, periodic):
    """D2: each record is the window conditionals, at the recorded state, of the model of the positions the
    dynamics ended at, at R T of the integrator (or the kT given) and with the offsets in force; with a duplex
    wrapped across a periodic box, with the box of the context."""
    box = BOX if periodic else None
    for kT in (None, 2 * KT):
        simulation = _simulation(tables, box=box)
        offsets = _offsets(tables, scale=KT)
        sampler = GibbsSampler(simulation, tables, kT=kT, offsets=offsets, seed=1)
        assert sampler.kT == pytest.approx(KT if kT is None else kT, rel=1e-12)
        assert sampler.windows == (1, 2)
        for _ in range(3):
            state = sampler.state.copy()
            records = sampler.run(1, steps_per_move=2)
            x, box_vectors = _positions(simulation)
            factors = factor_tables(tables, x, box_vectors if periodic else None)
            if periodic:
                np.testing.assert_allclose(box_vectors, np.diag(BOX), rtol=0, atol=1e-12)
                assert np.abs(factor_tables(tables, x).energies - factors.energies).max() > 1.0
            expected = model(factors, sampler.kT, offsets)
            np.testing.assert_array_equal(records['state'][0], state)
            np.testing.assert_array_equal(records['offsets'][0], offsets)
            for k in (1, 2):
                np.testing.assert_allclose(records[f'window{k}'][0], window_conditionals(expected, state, k),
                                           rtol=0, atol=1e-12)


def _singles(states, n, weights=None):
    """(N, n) the occupancy of the names at each position in `states` (S, N), weighted by `weights` (S,)."""
    weights = np.full(len(states), 1 / len(states)) if weights is None else weights
    return np.array([np.bincount(states[:, q], weights=weights, minlength=n) for q in range(states.shape[1])])


def _exact_singles(model_):
    """(N, n) the single-position marginals of p(s) of `model_`."""
    every = np.array(list(itertools.product(range(model_.n_names), repeat=model_.n_positions)))
    log_weights = model_.log_weights(every)
    p = np.exp(log_weights - log_weights.max())
    return _singles(every, model_.n_names, p / p.sum())


@pytest.mark.parametrize('move', ['joint', 'sweep'])
def test_draws_follow_the_model(tables, move):
    """D3: without dynamics the moves redraw the state from p(s) of the model of the fixed positions with the
    offsets in force."""
    simulation = _simulation(tables)
    factors = factor_tables(tables, _positions(simulation)[0])
    flattening = KT * np.log(window_conditionals(model(factors, KT), tables.wild_type, 1))
    offsets = flattening + _offsets(tables, scale=KT, seed=1)
    sampler = GibbsSampler(simulation, tables, offsets=offsets, move=move, seed=2)
    exact = _exact_singles(model(factors, KT, offsets))
    states = []
    for _ in range(100):
        states.append(sampler.run(10, steps_per_move=0)['state'])
        np.testing.assert_array_equal(_context_state(simulation, len(tables.wild_type)), sampler.state)
    empirical = _singles(np.concatenate(states), len(tables.alphabet))

    for other in (model(factors, KT), model(factors, 2 * KT, offsets)):
        assert 0.5 * np.abs(_exact_singles(other) - exact).sum(axis=1).max() > 0.2
    assert 0.5 * np.abs(empirical - exact).sum(axis=1).max() < 0.06


def test_offsets_adapt_by_sams(tables):
    """D4: with gamma, the offsets after move t of a run are the SAMS update of those in force with the move's
    single-position conditionals, the step gamma / max(1, t / 100) and the target; without gamma they stay."""
    simulation = _simulation(tables)
    target = w_target(tables.alphabet, 0.2)
    sampler = GibbsSampler(simulation, tables, seed=3)
    gamma, moves = 0.5, 120
    records = sampler.run(moves, steps_per_move=0, gamma=gamma, target=target)
    offsets = list(records['offsets']) + [sampler.offsets]
    np.testing.assert_array_equal(offsets[0], 0.0)
    for t in range(1, moves + 1):
        expected = adapt_offsets(offsets[t - 1], records['window1'][t - 1], KT, gamma / max(1.0, t / 100), target)
        np.testing.assert_allclose(offsets[t], expected, rtol=0, atol=1e-12)

    adapted = sampler.offsets.copy()
    records = sampler.run(3, steps_per_move=0)
    np.testing.assert_array_equal(records['offsets'], [adapted] * 3)
    np.testing.assert_array_equal(sampler.offsets, adapted)


def test_same_seeds_give_the_same_run(tables):
    """D5: two runs with the same integrator and sampler seeds give the same states, records and positions; another
    sampler seed changes the states."""
    def run(seed):
        simulation = _simulation(tables, seed=7)
        records = GibbsSampler(simulation, tables, seed=seed).run(15, steps_per_move=2, gamma=0.1)
        return records, _positions(simulation)[0]

    (first, x_first), (second, x_second), (other, _) = run(5), run(5), run(6)
    for key in first:
        np.testing.assert_array_equal(first[key], second[key])
    np.testing.assert_array_equal(x_first, x_second)
    assert (first['state'] != other['state']).any()


def test_refusals(tables):
    """D6: malformed samplers and runs are refused."""
    simulation = _simulation(tables)
    N, n = len(tables.wild_type), len(tables.alphabet)
    with pytest.raises(ValueError):
        GibbsSampler(simulation, tables, move='metropolis')
    with pytest.raises(ValueError):
        GibbsSampler(simulation, tables, windows=(1, N + 1))
    with pytest.raises(ValueError):
        GibbsSampler(simulation, tables, offsets=np.zeros((N, n - 1)))
    with pytest.raises(ValueError):
        GibbsSampler(_simulation(tables, integrator=openmm.VerletIntegrator(0.001)), tables)
    with pytest.raises(ValueError):
        GibbsSampler(simulation, tables).run(0)
