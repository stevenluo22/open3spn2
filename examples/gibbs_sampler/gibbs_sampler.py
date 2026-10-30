"""Free energies of the base pairs at three adjacent positions of a duplex with the Gibbs identity sampler.

The identities of positions 13-15 of strand A (and their partners on strand B) change during the simulation: every
move runs 25 steps of dynamics and redraws the identities from their distribution at the current coordinates. The
offsets are adapted by self-adjusted mixture sampling during equilibration so that every name is visited, then fixed
for production. The averaged conditionals of pairs of adjacent positions give the free energy of every Watson-Crick
sequence of the three positions relative to the wild type, exactly when the identity-dependent energy couples only
adjacent positions (the nearest-neighbour model).

    python gibbs_sampler.py [--adapt-ns 0.5] [--prod-ns 2] [--platform OpenCL]
"""
import argparse
import itertools

import numpy as np
import openmm
import openmm.unit as unit
import pandas

import open3SPN2
from open3SPN2.alchemy import GibbsSampler, alphabet, free_energies, identity_tables, markov_log_probability

SEQUENCE = 'ACGTCACTTCCTCTGTACGT'
POSITIONS = (13, 14, 15)
STEPS_PER_MOVE = 25
TIMESTEP = 2 * unit.femtosecond

parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
parser.add_argument('--adapt-ns', type=float, default=0.5)
parser.add_argument('--prod-ns', type=float, default=2.0)
parser.add_argument('--platform', default=None, help='OpenMM platform (default the fastest)')
args = parser.parse_args()
moves_in = lambda ns: round(ns * unit.nanosecond / (STEPS_PER_MOVE * TIMESTEP))

dna = open3SPN2.DNA.fromSequence(SEQUENCE, dna_type='B_curved', output_pdb='initial.pdb')
sites = [(('A', i), ('B', len(SEQUENCE) + 1 - i)) for i in POSITIONS]
tables = identity_tables(dna, sites, alphabet('wc'))

system = open3SPN2.System(dna)
system.clearForces()
tables.add_forces(system)
integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 2 / unit.picosecond, TIMESTEP)
platform = None if args.platform is None else openmm.Platform.getPlatformByName(args.platform)
simulation = openmm.app.Simulation(system.top, system._wrapped_system, integrator, platform)
simulation.context.setPositions(system.coord.positions)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

sampler = GibbsSampler(simulation, tables, windows=(2,))
sampler.run(moves_in(args.adapt_ns), STEPS_PER_MOVE, gamma=0.1)
records = sampler.run(moves_in(args.prod_ns), STEPS_PER_MOVE)

names = tables.alphabet.names
identities = [names.index(name) for name in tables.alphabet.identities]
sequences = np.array(list(itertools.product(identities, repeat=len(POSITIONS))))
log_p = markov_log_probability(records['window2'].mean(axis=0), sequences)
energies = free_energies(log_p, sequences, tables.wild_type, sampler.offsets, sampler.kT)
table = pandas.DataFrame({'sequence': ['-'.join(names[s] for s in sequence) for sequence in sequences],
                          'F - F_wild (kJ/mol)': energies}).sort_values('F - F_wild (kJ/mol)')
print(table.to_string(index=False, float_format='%.2f'))
