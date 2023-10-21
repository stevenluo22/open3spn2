from open3SPN2 import *

# Unit testing
def test_DNA_from_pdb():
    """ Test correct DNA initialization from PDB"""
    mol = DNA.fromPDB("Tests/1svc/1svc.pdb", template_from_X3DNA=False)


def test_DNA_from_gro():
    """ Test correct DNA initialization from gromacs files"""
    pass


def test_DNA_from_seq():
    """ Test correct DNA initialization from sequence files"""
    #return True #Needs X3DNA
    seq = 'ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
    dna = DNA.fromSequence(seq, dna_type='B_curved')

    # Compute the topology for the DNA structure.
    # Since the dna was generated from the sequence using X3DNA,
    # it is not necesary to recompute the geometry.

    dna.computeTopology(template_from_X3DNA=False)

    # Create the system.
    # To set periodic boundary conditions (periodicBox=[50,50,50]).
    # The periodic box size is in nanometers.
    dna.periodic = False
    s = System(dna, periodicBox=None)

    # Add 3SPN2 forces
    s.add3SPN2forces(verbose=True)

    # Initialize Molecular Dynamics simulations
    s.initializeMD(temperature=300 * openmm.unit.kelvin, platform_name='OpenCL')
    simulation = s.simulation

    # Set initial positions
    simulation.context.setPositions(s.coord.getPositions())

    energy_unit = openmm.unit.kilojoule_per_mole
    # Total energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print('TotalEnergy', round(energy, 6), energy_unit.get_symbol())

    # Detailed energy
    energies = {}
    for force_name, force in s.forces.items():
        group = force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2 ** group)
        energies[force_name] = state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in s.forces.keys():
        print(force_name, round(energies[force_name], 6), energy_unit.get_symbol())


def test_DNA_from_xyz():
    """Tests the correct parsing from an xyz file"""
    mol = DNA.fromXYZ('Tests/adna/in00_conf.xyz', template_from_X3DNA=False)
    assert mol.atoms.at[8, 'name'] == 'P'
    assert round(mol.atoms.at[188, 'y'], 6) == -8.779343


# Functional testing

def parse_xyz(filename=''):
    columns = ['N', 'timestep', 'id', 'name', 'x', 'y', 'z']
    data = []
    with open(filename, 'r') as traj_file:
        atom = pandas.Series(index=columns)
        atom['id'] = None
        for line in traj_file:
            s = line.split()
            if len(s) == 1:
                atom['N'] = int(s[0])
                if atom['id'] > -1:
                    assert atom['id'] == atoms
                atoms = int(s[0])
            elif len(s) == 3:
                atom['timestep'] = int(s[2])
                atom['id'] = 0
            elif len(s) > 3:
                atom['name'] = int(s[0])
                atom['x'], atom['y'], atom['z'] = [float(a) for a in s[1:4]]
                data += [atom.copy()]
                atom['id'] += 1
    xyz_data = pandas.concat(data, axis=1).T
    for i in ['N', 'timestep', 'id', 'name']:
        xyz_data[i] = xyz_data[i].astype(int)
    return xyz_data


def parse_log(filename=''):
    columns = ''
    log_data = []
    with open(filename, 'r') as log_file:
        start = False
        for line in log_file:
            if line[:4] == 'Step':
                columns = line.split()
                start = True
                continue
            if start:
                try:
                    log_data += [[float(a) for a in line.split()]]
                except ValueError:
                    break
    log_data = pandas.DataFrame(log_data, columns=columns)

    try:
        for i in ['Step', 'nbp']:
            log_data[i] = log_data[i].astype(int)
    except KeyError:
        for i in ['Step', 'v_nbp']:
            log_data[i] = log_data[i].astype(int)
    return log_data


def test_parse_xyz():
    """Tests the example trajectory parsing"""
    xyz_data = parse_xyz('Tests/adna/traj.xyz')
    assert xyz_data.at[1, 'name'] == 7
    assert xyz_data.at[1, 'x'] == 4.34621


def test_parse_log():
    """Tests the example log parsing"""
    log_data = parse_log('Tests/adna/sim.log')
    assert log_data.at[1, 'Step'] == 2000
    assert log_data.at[1, 'eexcl'] == 0.45734636