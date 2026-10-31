"""Numpy implementation of the 3SPN2 potentials.

Every interaction has two functions:

* the potential, e.g. ``angle(theta, t0, k)``, the mathematical expression evaluated on arrays of geometry and
  parameters;
* the geometry, e.g. ``angle_geometry(table, coordinates, box=None)``, which depends only on the particles of every
  row of an interaction table: a dict of arrays (..., rows) named as the arguments of the potential;
* the table evaluator, e.g. ``angle_energy(table, coordinates, box=None, geometry=None)``, which returns the energy of
  each row at the coordinates, or at the given `geometry` of the rows (the coordinates are then not used).

An interaction table is a DataFrame with one row per interaction (the ``parameter_table`` of a force class), or a dict
of column arrays: the particle indices in the columns ``aai``, ``aaj``, ``aak``, ``aal``, or ``d1``, ``d2``, ``a1``,
``a2``, ``a3`` for the donor-acceptor forces, and one column per parameter, named as in the OpenMM force. Parameters
that OpenMM combines from two particles are stored combined, and the cutoff and global parameters are columns too.
Coordinates are (N, 3) or (frames, N, 3); a periodic ``box`` is orthorhombic. Units are those of the OpenMM forces: nm,
rad and kJ/mol. The global ``k`` multiplier of a force is not included.
"""
import numpy as np

from ..geometry import distances, angles, dihedrals


def _columns(table, *names):
    """(rows, len(names)) the columns `names` of `table` side by side."""
    return np.stack([np.asarray(table[name]) for name in names], axis=-1)


def gate(dt):
    """Angular modulation of the 3SPN2 attractions: 1 for abs(dt) <= pi/2, sin(dt)^2 for pi/2 < abs(dt) <= pi, 0 beyond."""
    magnitude = np.abs(dt)
    out = (magnitude <= np.pi / 2).astype(float)
    edge = (magnitude > np.pi / 2) & (magnitude <= np.pi)
    out[edge] = np.sin(dt[edge]) ** 2
    return out


def morse(dr, epsilon, alpha):
    """(repulsive, attractive) parts of the 3SPN2 Morse well at dr = r - sigma: the repulsion acts for dr <= 0 and
    the attraction, shifted to -epsilon at the minimum, for dr >= 0."""
    well = epsilon * (1 - np.exp(-alpha * dr)) ** 2
    return np.where(dr <= 0, well, 0.0), np.where(dr >= 0, well, 0.0) - epsilon


def bond(r, r0, Kb2, Kb3, Kb4):
    """Quartic bond potential."""
    dr = r - r0
    return Kb2 * dr ** 2 + Kb3 * dr ** 3 + Kb4 * dr ** 4


def bond_geometry(table, coordinates, box=None):
    return {'r': distances(coordinates, _columns(table, 'aai', 'aaj'), box)}


def bond_energy(table, coordinates, box=None, geometry=None):
    geometry = bond_geometry(table, coordinates, box) if geometry is None else geometry
    return bond(geometry['r'], np.asarray(table['r0']), np.asarray(table['Kb2']), np.asarray(table['Kb3']),
                np.asarray(table['Kb4']))


def angle(theta, t0, k):
    """Harmonic angle potential."""
    return 0.5 * k * (theta - t0) ** 2


def angle_geometry(table, coordinates, box=None):
    return {'theta': angles(coordinates, _columns(table, 'aai', 'aaj', 'aak'), box)}


def angle_energy(table, coordinates, box=None, geometry=None):
    geometry = angle_geometry(table, coordinates, box) if geometry is None else geometry
    return angle(geometry['theta'], np.asarray(table['t0']), np.asarray(table['k']))


def stacking(r, theta, epsilon, sigma, t0, alpha, rng):
    """Base stacking: a Morse well in the distance r of the stacked bases, whose attraction is modulated by the
    angle theta (sugar, base, next base)."""
    repulsion, attraction = morse(r - sigma, epsilon, alpha)
    return repulsion + gate(rng * (theta - t0)) * attraction


def stacking_geometry(table, coordinates, box=None):
    """The rows of a stacking table are (sugar ``aai``, base ``aaj``, stacked base ``aak``) triples."""
    return {'r': distances(coordinates, _columns(table, 'aaj', 'aak'), box),
            'theta': angles(coordinates, _columns(table, 'aai', 'aaj', 'aak'), box)}


def stacking_energy(table, coordinates, box=None, geometry=None):
    geometry = stacking_geometry(table, coordinates, box) if geometry is None else geometry
    return stacking(geometry['r'], geometry['theta'], np.asarray(table['epsilon']), np.asarray(table['sigma']),
                    np.asarray(table['t0']), np.asarray(table['alpha']), np.asarray(table['rng']))


def dihedral(theta, K_periodic, K_gaussian, sigma, t0):
    """Periodic dihedral potential with a Gaussian well at t0."""
    dt = theta - t0
    dt_periodic = dt - np.floor((dt + np.pi) / (2 * np.pi)) * (2 * np.pi)
    return K_periodic * (1 - np.cos(dt)) - K_gaussian * np.exp(-dt_periodic ** 2 / 2 / sigma ** 2)


def dihedral_geometry(table, coordinates, box=None):
    return {'theta': dihedrals(coordinates, _columns(table, 'aai', 'aaj', 'aak', 'aal'), box)}


def dihedral_energy(table, coordinates, box=None, geometry=None):
    geometry = dihedral_geometry(table, coordinates, box) if geometry is None else geometry
    return dihedral(geometry['theta'], np.asarray(table['K_periodic']), np.asarray(table['K_gaussian']),
                    np.asarray(table['sigma']), np.asarray(table['t0']))


def base_pair(r, theta1, theta2, phi, phi0, sigma, t01, t02, rng, epsilon, alpha):
    """Base pairing: a Morse well in the distance r of the paired bases, whose attraction is modulated by the
    angles theta1 (sugar, base, paired base) and theta2 (paired sugar, paired base, base) and by the dihedral phi
    (sugar, base, paired base, paired sugar)."""
    repulsion, attraction = morse(r - sigma, epsilon, alpha)
    return repulsion + 0.5 * (1 + np.cos(phi - phi0)) * gate(rng * (theta1 - t01)) * gate(rng * (theta2 - t02)) \
        * attraction


def base_pair_geometry(table, coordinates, box=None):
    """The rows of a base pair table are (base ``d1``, sugar ``d2``) donors paired with (base ``a1``, sugar ``a2``)
    acceptors."""
    return {'r': distances(coordinates, _columns(table, 'd1', 'a1'), box),
            'theta1': angles(coordinates, _columns(table, 'd2', 'd1', 'a1'), box),
            'theta2': angles(coordinates, _columns(table, 'a2', 'a1', 'd1'), box),
            'phi': dihedrals(coordinates, _columns(table, 'd2', 'd1', 'a1', 'a2'), box)}


def base_pair_energy(table, coordinates, box=None, geometry=None):
    """Pairs whose bases are farther apart than ``cutoff`` have no energy."""
    geometry = base_pair_geometry(table, coordinates, box) if geometry is None else geometry
    energy = base_pair(geometry['r'], geometry['theta1'], geometry['theta2'], geometry['phi'],
                       np.asarray(table['phi0']), np.asarray(table['sigma']), np.asarray(table['t01']),
                       np.asarray(table['t02']), np.asarray(table['rng']), np.asarray(table['epsilon']),
                       np.asarray(table['alpha']))
    return np.where(geometry['r'] < np.asarray(table['cutoff']), energy, 0.0)


def cross_stacking(r, theta3, theta_cs, t03, t0CS, rng_CS, rng_BP, epsilon, alpha, sigma):
    """Cross stacking: the attractive part of a Morse well in the distance r between a base and the neighbour of
    its paired base, modulated by the angle theta_cs (sugar, base, neighbour) and by the angle theta3 between the
    base-sugar bonds of the pair."""
    _, attraction = morse(r - sigma, epsilon, alpha)
    return gate(rng_BP * (theta3 - t03)) * gate(rng_CS * (theta_cs - t0CS)) * attraction / 2


def cross_stacking_geometry(table, coordinates, box=None):
    """The rows of a cross stacking table are (base ``d1``, sugar ``d2``) donors paired with (base ``a1``, sugar
    ``a2``, neighbour ``a3``) acceptors. theta3 is computed, as in OpenMM, from the angles theta1 (d2, d1, a1) and
    theta2 (d1, a1, a2) and the dihedral phi (d2, d1, a1, a2) of the pair, its cosine clipped to +-0.99; ``r_pair`` is
    the distance of the paired bases ``d1``, ``a1``."""
    theta1 = angles(coordinates, _columns(table, 'd2', 'd1', 'a1'), box)
    theta2 = angles(coordinates, _columns(table, 'd1', 'a1', 'a2'), box)
    phi = dihedrals(coordinates, _columns(table, 'd2', 'd1', 'a1', 'a2'), box)
    cos_theta3 = np.sin(theta1) * np.sin(theta2) * np.cos(phi) - np.cos(theta1) * np.cos(theta2)
    return {'r': distances(coordinates, _columns(table, 'd1', 'a3'), box),
            'theta3': np.arccos(np.clip(cos_theta3, -0.99, 0.99)),
            'theta_cs': angles(coordinates, _columns(table, 'd2', 'd1', 'a3'), box),
            'r_pair': distances(coordinates, _columns(table, 'd1', 'a1'), box)}


def cross_stacking_energy(table, coordinates, box=None, geometry=None):
    """Pairs whose bases ``d1``, ``a1`` are farther apart than ``cutoff`` have no energy."""
    geometry = cross_stacking_geometry(table, coordinates, box) if geometry is None else geometry
    energy = cross_stacking(geometry['r'], geometry['theta3'], geometry['theta_cs'], np.asarray(table['t03']),
                            np.asarray(table['t0CS']), np.asarray(table['rng_CS']), np.asarray(table['rng_BP']),
                            np.asarray(table['epsilon']), np.asarray(table['alpha']), np.asarray(table['sigma']))
    return np.where(geometry['r_pair'] < np.asarray(table['cutoff']), energy, 0.0)


def exclusion(r, epsilon, sigma):
    """Excluded volume: a Lennard-Jones repulsion shifted to zero at, and truncated beyond, r = sigma."""
    return np.where(r <= sigma, epsilon * ((sigma / r) ** 12 - 2 * (sigma / r) ** 6) + epsilon, 0.0)


def exclusion_geometry(table, coordinates, box=None):
    return {'r': distances(coordinates, _columns(table, 'aai', 'aaj'), box)}


def exclusion_energy(table, coordinates, box=None, geometry=None):
    geometry = exclusion_geometry(table, coordinates, box) if geometry is None else geometry
    energy = exclusion(geometry['r'], np.asarray(table['epsilon']), np.asarray(table['sigma']))
    return np.where(geometry['r'] < np.asarray(table['cutoff']), energy, 0.0)


def electrostatics(r, qq, dh_length, denominator):
    """Debye-Hueckel screened Coulomb interaction of two charges with product qq."""
    return qq * np.exp(-r / dh_length) / denominator / r


def electrostatics_geometry(table, coordinates, box=None):
    return {'r': distances(coordinates, _columns(table, 'aai', 'aaj'), box)}


def electrostatics_energy(table, coordinates, box=None, geometry=None):
    geometry = electrostatics_geometry(table, coordinates, box) if geometry is None else geometry
    energy = electrostatics(geometry['r'], np.asarray(table['qq']), np.asarray(table['dh_length']),
                            np.asarray(table['denominator']))
    return np.where(geometry['r'] < np.asarray(table['cutoff']), energy, 0.0)
