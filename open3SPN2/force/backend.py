"""Numpy implementation of the 3SPN2 potentials.

Every interaction has two functions:

* the potential, e.g. ``angle(theta, t0, k)``, the mathematical expression evaluated on arrays of geometry and
  parameters;
* the table evaluator, e.g. ``angle_energy(table, coordinates, box=None)``, which computes the geometry of every
  row of an interaction table and returns the energy of each row.

An interaction table is a DataFrame with one row per interaction: the particle indices in the columns ``aai``,
``aaj``, ``aak``, ``aal`` and one column per parameter, named as in the OpenMM force (the ``table`` of a force
class). Coordinates are (N, 3) or (frames, N, 3); a periodic ``box`` is orthorhombic. Units are those of the
OpenMM forces: nm, rad and kJ/mol. The global ``k`` multiplier of a force is not included.
"""
import numpy as np

from ..geometry import distances, angles, dihedrals


def gate(dt):
    """Angular modulation of the 3SPN2 attractions: 1 for |dt| <= pi/2, sin(dt)^2 for pi/2 < |dt| <= pi, 0 beyond."""
    inside = np.abs(dt) <= np.pi
    return np.where(inside, np.maximum(np.sin(dt) ** 2, np.abs(dt) <= np.pi / 2), 0.0)


def morse(dr, epsilon, alpha):
    """(repulsive, attractive) parts of the 3SPN2 Morse well at dr = r - sigma: the repulsion acts for dr <= 0 and
    the attraction, shifted to -epsilon at the minimum, for dr >= 0."""
    well = epsilon * (1 - np.exp(-alpha * dr)) ** 2
    return np.where(dr <= 0, well, 0.0), np.where(dr >= 0, well, 0.0) - epsilon


def bond(r, r0, Kb2, Kb3, Kb4):
    """Quartic bond potential."""
    dr = r - r0
    return Kb2 * dr ** 2 + Kb3 * dr ** 3 + Kb4 * dr ** 4


def bond_energy(table, coordinates, box=None):
    r = distances(coordinates, table[['aai', 'aaj']].to_numpy(), box)
    return bond(r, table['r0'].to_numpy(), table['Kb2'].to_numpy(), table['Kb3'].to_numpy(),
                table['Kb4'].to_numpy())


def angle(theta, t0, k):
    """Harmonic angle potential."""
    return 0.5 * k * (theta - t0) ** 2


def angle_energy(table, coordinates, box=None):
    theta = angles(coordinates, table[['aai', 'aaj', 'aak']].to_numpy(), box)
    return angle(theta, table['t0'].to_numpy(), table['k'].to_numpy())


def stacking(r, theta, epsilon, sigma, t0, alpha, rng):
    """Base stacking: a Morse well in the distance r of the stacked bases, whose attraction is modulated by the
    angle theta (sugar, base, next base)."""
    repulsion, attraction = morse(r - sigma, epsilon, alpha)
    return repulsion + gate(rng * (theta - t0)) * attraction


def stacking_energy(table, coordinates, box=None):
    """The rows of a stacking table are (sugar ``aai``, base ``aaj``, stacked base ``aak``) triples."""
    r = distances(coordinates, table[['aaj', 'aak']].to_numpy(), box)
    theta = angles(coordinates, table[['aai', 'aaj', 'aak']].to_numpy(), box)
    return stacking(r, theta, table['epsilon'].to_numpy(), table['sigma'].to_numpy(), table['t0'].to_numpy(),
                    table['alpha'].to_numpy(), table['rng'].to_numpy())


def dihedral(theta, K_periodic, K_gaussian, sigma, t0):
    """Periodic dihedral potential with a Gaussian well at t0."""
    dt = theta - t0
    dt_periodic = dt - np.floor((dt + np.pi) / (2 * np.pi)) * (2 * np.pi)
    return K_periodic * (1 - np.cos(dt)) - K_gaussian * np.exp(-dt_periodic ** 2 / 2 / sigma ** 2)


def dihedral_energy(table, coordinates, box=None):
    theta = dihedrals(coordinates, table[['aai', 'aaj', 'aak', 'aal']].to_numpy(), box)
    return dihedral(theta, table['K_periodic'].to_numpy(), table['K_gaussian'].to_numpy(),
                    table['sigma'].to_numpy(), table['t0'].to_numpy())
