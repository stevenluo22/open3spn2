"""Numpy implementation of the 3SPN2 potentials.

Every interaction has two functions:

* the potential, e.g. ``angle(theta, t0, k)``, the mathematical expression evaluated on arrays of geometry and
  parameters;
* the table evaluator, e.g. ``angle_energy(table, coordinates, box=None)``, which computes the geometry of every
  row of an interaction table and returns the energy of each row.

An interaction table is a DataFrame with one row per interaction (the ``parameter_table`` of a force class): the
particle indices in the columns ``aai``, ``aaj``, ``aak``, ``aal``, or ``d1``, ``d2``, ``a1``, ``a2``, ``a3`` for the
donor-acceptor forces, and one column per parameter, named as in the OpenMM force. Parameters that OpenMM combines
from two particles are stored combined, and the cutoff and global parameters are columns too. Coordinates are
(N, 3) or (frames, N, 3); a periodic ``box`` is orthorhombic. Units are those of the OpenMM forces: nm, rad and
kJ/mol. The global ``k`` multiplier of a force is not included.
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


def base_pair(r, theta1, theta2, phi, phi0, sigma, t01, t02, rng, epsilon, alpha):
    """Base pairing: a Morse well in the distance r of the paired bases, whose attraction is modulated by the
    angles theta1 (sugar, base, paired base) and theta2 (paired sugar, paired base, base) and by the dihedral phi
    (sugar, base, paired base, paired sugar)."""
    repulsion, attraction = morse(r - sigma, epsilon, alpha)
    return repulsion + 0.5 * (1 + np.cos(phi - phi0)) * gate(rng * (theta1 - t01)) * gate(rng * (theta2 - t02)) \
        * attraction


def base_pair_energy(table, coordinates, box=None):
    """The rows of a base pair table are (base ``d1``, sugar ``d2``) donors paired with (base ``a1``, sugar ``a2``)
    acceptors; pairs whose bases are farther apart than ``cutoff`` have no energy."""
    r = distances(coordinates, table[['d1', 'a1']].to_numpy(), box)
    theta1 = angles(coordinates, table[['d2', 'd1', 'a1']].to_numpy(), box)
    theta2 = angles(coordinates, table[['a2', 'a1', 'd1']].to_numpy(), box)
    phi = dihedrals(coordinates, table[['d2', 'd1', 'a1', 'a2']].to_numpy(), box)
    energy = base_pair(r, theta1, theta2, phi, table['phi0'].to_numpy(), table['sigma'].to_numpy(),
                       table['t01'].to_numpy(), table['t02'].to_numpy(), table['rng'].to_numpy(),
                       table['epsilon'].to_numpy(), table['alpha'].to_numpy())
    return np.where(r < table['cutoff'].to_numpy(), energy, 0.0)


def cross_stacking(r, theta1, theta2, phi, theta_cs, t03, t0CS, rng_CS, rng_BP, epsilon, alpha, sigma):
    """Cross stacking: the attractive part of a Morse well in the distance r between a base and the neighbour of
    its paired base, modulated by the angle theta_cs (sugar, base, neighbour) and by the angle theta3 between the
    base-sugar bonds of the pair, computed from the base pair angles theta1, theta2 and dihedral phi."""
    _, attraction = morse(r - sigma, epsilon, alpha)
    cos_theta3 = np.clip(np.sin(theta1) * np.sin(theta2) * np.cos(phi) - np.cos(theta1) * np.cos(theta2), -0.99, 0.99)
    theta3 = np.arccos(cos_theta3)
    return gate(rng_BP * (theta3 - t03)) * gate(rng_CS * (theta_cs - t0CS)) * attraction / 2


def cross_stacking_energy(table, coordinates, box=None):
    """The rows of a cross stacking table are (base ``d1``, sugar ``d2``) donors paired with (base ``a1``, sugar
    ``a2``, neighbour ``a3``) acceptors; pairs whose bases ``d1``, ``a1`` are farther apart than ``cutoff`` have no
    energy."""
    r = distances(coordinates, table[['d1', 'a3']].to_numpy(), box)
    theta1 = angles(coordinates, table[['d2', 'd1', 'a1']].to_numpy(), box)
    theta2 = angles(coordinates, table[['d1', 'a1', 'a2']].to_numpy(), box)
    phi = dihedrals(coordinates, table[['d2', 'd1', 'a1', 'a2']].to_numpy(), box)
    theta_cs = angles(coordinates, table[['d2', 'd1', 'a3']].to_numpy(), box)
    energy = cross_stacking(r, theta1, theta2, phi, theta_cs, table['t03'].to_numpy(), table['t0CS'].to_numpy(),
                            table['rng_CS'].to_numpy(), table['rng_BP'].to_numpy(), table['epsilon'].to_numpy(),
                            table['alpha'].to_numpy(), table['sigma'].to_numpy())
    cutoff = distances(coordinates, table[['d1', 'a1']].to_numpy(), box) < table['cutoff'].to_numpy()
    return np.where(cutoff, energy, 0.0)


def exclusion(r, epsilon, sigma):
    """Excluded volume: a Lennard-Jones repulsion shifted to zero at, and truncated beyond, r = sigma."""
    return np.where(r <= sigma, epsilon * ((sigma / r) ** 12 - 2 * (sigma / r) ** 6) + epsilon, 0.0)


def exclusion_energy(table, coordinates, box=None):
    r = distances(coordinates, table[['aai', 'aaj']].to_numpy(), box)
    energy = exclusion(r, table['epsilon'].to_numpy(), table['sigma'].to_numpy())
    return np.where(r < table['cutoff'].to_numpy(), energy, 0.0)


def electrostatics(r, qq, dh_length, denominator):
    """Debye-Hueckel screened Coulomb interaction of two charges with product qq."""
    return qq * np.exp(-r / dh_length) / denominator / r


def electrostatics_energy(table, coordinates, box=None):
    r = distances(coordinates, table[['aai', 'aaj']].to_numpy(), box)
    energy = electrostatics(r, table['qq'].to_numpy(), table['dh_length'].to_numpy(),
                            table['denominator'].to_numpy())
    return np.where(r < table['cutoff'].to_numpy(), energy, 0.0)
