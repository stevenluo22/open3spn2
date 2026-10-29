import numpy as np
import openmm.unit as unit

class Transform:
    def __init__(self, dx=0, dy=0, dz=0, rx=0, ry=0, rz=0):
        """
        Initialize the Transform object with displacements and rotations.
        
        Parameters:
        - dx, dy, dz: Displacements in the x, y, z directions respectively (Angstroms).
        - rx, ry, rz: Rotations around the x, y, z axes respectively (Radians).
        """
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rx = rx
        self.ry = ry
        self.rz = rz
        
        # Angle needed to rotate around the z axis to achieve a single rotation around the y axis
        self.phi = np.arctan2(self.rx, self.ry)
        
        # Magnitude of the combined x and y rotation
        self.yx_ang = np.sqrt(self.rx**2 + self.ry**2)

    @property
    def full_rotation(self):
        """
        Compute the full rotation transformation matrix that transitions from one state to another.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        rotz1 = self.rotz(0.5 * self.rz - self.phi)
        rotyx = self.roty(self.yx_ang)
        rotz2 = self.rotz(0.5 * self.rz + self.phi)
        return rotz1 @ rotyx @ rotz2

    @property
    def half_rotation(self):
        """
        Compute the rotation matrix for half the rotation.
        
        This aims to rotate two objects equally relative to a reference frame.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        rotz1 = self.rotz(0.5 * self.rz - self.phi)
        rotyx = self.roty(0.5 * self.yx_ang)
        rotz2 = self.rotz(self.phi)
        return rotz1 @ rotyx @ rotz2

    @property
    def full_displacement(self):
        """
        Compute the displacement vector for transitioning from one state to another.
        
        Returns:
        - A 1x3 displacement vector.
        """
        displacement = np.array([self.dx, self.dy, self.dz])
        return np.dot(displacement, self.half_rotation.T)

    @property
    def half_displacement(self):
        """
        Compute the midpoint of the displacement vector.
        
        Returns:
        - A 1x3 displacement vector.
        """
        return 0.5 * self.full_displacement
    
    @staticmethod
    def rotz(theta):
        """
        Compute the rotation matrix around the Z-axis.
        
        Parameters:
        - theta: Rotation angle in radians.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    @staticmethod
    def roty(theta):
        """
        Compute the rotation matrix around the Y-axis.
        
        Parameters:
        - theta: Rotation angle in radians.
        
        Returns:
        - A 3x3 rotation matrix.
        """
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def _nanometers(x):
    return np.asarray(x.value_in_unit(unit.nanometer) if unit.is_quantity(x) else x, dtype=float)


def minimum_image(v, box=None):
    """Displacement vectors `v` (..., 3) (nm) mapped to their nearest periodic image.

    `box` is an orthorhombic box, as its three edge lengths (3,) or as diagonal box vectors (3, 3), in nm; with
    no box the vectors are returned unchanged."""
    if box is None:
        return v
    box = _nanometers(box)
    if box.shape == (3, 3):
        if np.any(box[~np.eye(3, dtype=bool)] != 0):
            raise NotImplementedError('Only orthorhombic periodic boxes are supported')
        box = np.diag(box)
    if box.shape != (3,):
        raise ValueError(f'A box has shape (3,) or (3, 3), not {box.shape}')
    return v - box * np.floor(v / box + 0.5)


def _vectors(x, start, end, box):
    return minimum_image(x[..., end, :] - x[..., start, :], box)


def _dot(a, b):
    return (a * b).sum(-1)


def distances(x, pairs, box=None):
    """Distance (nm) of each pair (P, 2) of particles at positions `x` ((N, 3) or (frames, N, 3), nm)."""
    x, pairs = _nanometers(x), np.asarray(pairs)
    v = _vectors(x, pairs[:, 0], pairs[:, 1], box)
    return np.sqrt(_dot(v, v))


def angles(x, triples, box=None):
    """Angle (rad) at the middle particle of each triple (A, 3) of particles."""
    x, triples = _nanometers(x), np.asarray(triples)
    v1 = _vectors(x, triples[:, 1], triples[:, 0], box)
    v2 = _vectors(x, triples[:, 1], triples[:, 2], box)
    return np.arccos(np.clip(_dot(v1, v2) / np.sqrt(_dot(v1, v1) * _dot(v2, v2)), -1.0, 1.0))


def dihedrals(x, quadruplets, box=None):
    """Dihedral angle (rad, in [-pi, pi], OpenMM sign convention) of each quadruplet (D, 4) of particles."""
    x, quadruplets = _nanometers(x), np.asarray(quadruplets)
    b0, b1, b2 = (_vectors(x, quadruplets[:, i], quadruplets[:, i + 1], box) for i in range(3))
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    m1 = np.cross(n1, b1 / np.sqrt(_dot(b1, b1))[..., None])
    return -np.arctan2(_dot(m1, n2), _dot(n1, n2))
