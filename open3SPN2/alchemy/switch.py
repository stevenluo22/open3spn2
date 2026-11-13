"""Switch forces: the identity-dependent interactions of a duplex in an OpenMM system, following the state of the
positions set in the context.

``add_forces`` adds to a system the open3SPN2 DNA forces of the duplex without the interactions of the identity
tables (``native_forces``) and, for every kind with a non-empty table, a switch force with them in the force group of
the kind (``switch_force``). A switch force is a CustomCompoundBondForce with one bond per interaction (``term``) of
the table, whose parameters are the tabulated values at the row

    first + sum_j n^(k-1-j) * state_{position_j}

of the table: ``first`` is the first row of the term, k the number of its positions and n the size of the alphabet
(the rows of a term are in ``itertools.product`` order). The global parameter ``state_p`` holds the name index of
position p, so ``set_state`` changes the identities without rebuilding the context.

The native forces lose the switched interactions without new exclusions, since OpenMM needs the same exclusions in
every CustomNonbondedForce and its GPU platforms cap the exclusions of a CustomHbondForce donor or acceptor: the
bonded interactions keep their place with zero strength, the bases of the positions have zero excluded-volume epsilon,
and the base-pair and cross-stacking forces are built without the donors and acceptors of those bases.
"""
import numpy as np
import openmm
import pandas

from ..ff3SPN2 import forces as dna_forces

STRENGTH = {'Bond': ['Kb2', 'Kb3', 'Kb4'], 'Angle': ['k'], 'Stacking': ['epsilon'],
            'Dihedral': ['K_periodic', 'K_gaussian']}

# kind -> (particle columns of the table, as p1, p2, ...; energy in the parameters of the table)
FORMULAS = {
    'Bond': (['aai', 'aaj'], 'Kb2*dr^2+Kb3*dr^3+Kb4*dr^4; dr=distance(p1,p2)-r0'),
    'Angle': (['aai', 'aaj', 'aak'], '0.5*k*(angle(p1,p2,p3)-t0)^2'),
    'Stacking': (['aai', 'aaj', 'aak'], '''rep+f2*attr;
                 rep=epsilon*(1-exp(-alpha*(dr)))^2*step(-dr);
                 attr=epsilon*(1-exp(-alpha*(dr)))^2*step(dr)-epsilon;
                 dr=distance(p2,p3)-sigma;
                 f2=max(f*pair2,pair1);
                 pair1=step(dt+pi/2)*step(pi/2-dt);
                 pair2=step(dt+pi)*step(pi-dt);
                 f=1-cos(dt)^2;
                 dt=rng*(angle(p1,p2,p3)-t0)'''),
    'Dihedral': (['aai', 'aaj', 'aak', 'aal'], '''K_periodic*(1-cs)-K_gaussian*exp(-dt_periodic^2/2/sigma^2);
                 cs=cos(dt);
                 dt_periodic=dt-floor((dt+pi)/(2*pi))*(2*pi);
                 dt=dihedral(p1,p2,p3,p4)-t0'''),
    'BasePair': (['d1', 'd2', 'a1', 'a2'], '''step(cutoff-distance(p1,p3))*(rep+1/2*(1+cos(dphi))*fdt1*fdt2*attr);
                 rep=epsilon*(1-exp(-alpha*dr))^2*(1-step(dr));
                 attr=epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                 fdt1=max(f1*pair0t1,pair1t1);
                 fdt2=max(f2*pair0t2,pair1t2);
                 pair1t1=step(pi/2+dt1)*step(pi/2-dt1);
                 pair1t2=step(pi/2+dt2)*step(pi/2-dt2);
                 pair0t1=step(pi+dt1)*step(pi-dt1);
                 pair0t2=step(pi+dt2)*step(pi-dt2);
                 f1=1-cos(dt1)^2;
                 f2=1-cos(dt2)^2;
                 dphi=dihedral(p2,p1,p3,p4)-phi0;
                 dr=distance(p1,p3)-sigma;
                 dt1=rng*(angle(p2,p1,p3)-t01);
                 dt2=rng*(angle(p4,p3,p1)-t02)'''),
    'CrossStacking': (['d1', 'd2', 'a1', 'a2', 'a3'], '''step(cutoff-distance(p1,p3))*fdt3*fdtCS*attr/2;
                      attr=epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                      fdt3=max(f1*pair0t3,pair1t3);
                      fdtCS=max(f2*pair0tCS,pair1tCS);
                      pair0t3=step(pi+dt3)*step(pi-dt3);
                      pair0tCS=step(pi+dtCS)*step(pi-dtCS);
                      pair1t3=step(pi/2+dt3)*step(pi/2-dt3);
                      pair1tCS=step(pi/2+dtCS)*step(pi/2-dtCS);
                      f1=1-cos(dt3)^2;
                      f2=1-cos(dtCS)^2;
                      dr=distance(p1,p5)-sigma;
                      dt3=rng_BP*(t3-t03);
                      dtCS=rng_CS*(tCS-t0CS);
                      tCS=angle(p2,p1,p5);
                      t3=acos(cost3lim);
                      cost3lim=min(max(cost3,-0.99),0.99);
                      cost3=sin(t1)*sin(t2)*cos(phi)-cos(t1)*cos(t2);
                      t1=angle(p2,p1,p3);
                      t2=angle(p1,p3,p4);
                      phi=dihedral(p2,p1,p3,p4)'''),
    'Exclusion': (['aai', 'aaj'], '''step(cutoff-r)*(epsilon*((sigma/r)^12-2*(sigma/r)^6)+epsilon)*step(sigma-r);
                  r=distance(p1,p2)'''),
}


def add_forces(tables, system):
    """Adds to `system` the open3SPN2 DNA forces of ``tables.dna`` whose identity-dependent interactions follow the
    state set with `set_state` (initially the wild type) and returns {name: native force}."""
    forces = native_forces(tables)
    for name, force in forces.items():
        force.addForce(system)
        if len(tables.kinds.get(name, ())):
            switch = switch_force(tables, name)
            switch.setForceGroup(force.force_group)
            system.addForce(switch)
    return forces


def set_state(tables, context, state, previous=None):
    """Selects in `context` the names `state` (the name index of every position) for the switch forces of `tables`;
    with `previous`, the state the context holds, only the positions that change are set."""
    state = np.asarray(state)
    if state.shape != tables.wild_type.shape or not np.isin(state, np.arange(len(tables.alphabet))).all():
        raise ValueError(f'a state is {len(tables.wild_type)} name indices of {tables.alphabet}')
    changed = range(len(state)) if previous is None else np.flatnonzero(state != np.asarray(previous))
    for p in changed:
        context.setParameter(f'state_{p}', float(state[p]))


def native_forces(tables):
    """{name: force} the open3SPN2 DNA forces of ``tables.dna`` without the interactions of the identity tables."""
    dna = tables.dna
    residues = {residue for pair in tables.sites for residue in pair}
    atoms = dna.atoms
    is_variable = atoms['name'].isin(['A', 'T', 'G', 'C']) & pandas.Series(
        [(chain, int(residue)) in residues for chain, residue in zip(atoms['chainID'], atoms['resSeq'])],
        index=atoms.index)
    variable = atoms.index[is_variable].tolist()

    forces = {}
    for name, force_class in dna_forces.items():
        if name in ('BasePair', 'CrossStacking'):
            forces[name] = force_class(dna, skip_bases=variable if name in tables.kinds else ())
        else:
            forces[name] = force_class(dna)
    for kind, strength in STRENGTH.items():
        if kind in tables.kinds:
            _remove_bonded(forces[kind], tables.kinds[kind], FORMULAS[kind][0], strength)
    if 'Exclusion' in tables.kinds:
        exclusion = forces['Exclusion']
        for i in variable:
            epsilon, sigma, bt, seqid = exclusion.force.getParticleParameters(i)
            exclusion.force.setParticleParameters(i, [0.0, sigma, bt, seqid])
    return forces


def switch_force(tables, kind):
    """The CustomCompoundBondForce of the interactions of `kind` in `tables`, one bond per term, following the
    global parameters ``state_p`` (default the wild type); its force group is left at 0."""
    table = tables.kinds[kind]
    particles, formula = FORMULAS[kind]
    position_columns = [c for c in table if c.startswith('position_')]
    name_columns = [c for c in table if c.startswith('name_')]
    parameters = [c for c in table if c not in ['term'] + particles + position_columns + name_columns]

    first = np.flatnonzero(np.diff(table['term'].to_numpy(), prepend=-1))
    positions = table[position_columns].to_numpy(dtype=int)[first]
    strides = np.zeros((len(first), len(tables.wild_type)), dtype=int)
    for t, term_positions in enumerate(positions):
        term_positions = term_positions[term_positions >= 0]
        strides[t, term_positions] = len(tables.alphabet) ** np.arange(len(term_positions))[::-1]
    used = np.flatnonzero(strides.any(axis=0))

    lookups = ''.join(f'; {c}={c}_table(row)' for c in parameters)
    row = 'first' + ''.join(f'+stride_{p}*state_{p}' for p in used)
    force = openmm.CustomCompoundBondForce(len(particles), f'{formula}{lookups}; row={row}')
    force.addPerBondParameter('first')
    for p in used:
        force.addPerBondParameter(f'stride_{p}')
    for p, name in enumerate(tables.wild_type):
        force.addGlobalParameter(f'state_{p}', float(name))
    force.addGlobalParameter('pi', np.pi)
    for c in parameters:
        force.addTabulatedFunction(f'{c}_table', openmm.Discrete1DFunction(table[c].to_numpy(dtype=float).tolist()))
    bond_particles = table[particles].to_numpy(dtype=int)[first]
    for t, row_index in enumerate(first):
        force.addBond(bond_particles[t].tolist(), [float(row_index)] + strides[t, used].astype(float).tolist())
    force.setUsesPeriodicBoundaryConditions(tables.dna.periodic)
    force.setName(f'{kind} switch')
    return force


def _remove_bonded(force, table, particles, strength):
    """Sets to zero the `strength` parameters of the interactions of the bonded `force` that are in `table`."""
    native = force.parameter_table
    parameters = [c for c in native if c not in particles]
    switched = pandas.MultiIndex.from_frame(native[particles]).isin(pandas.MultiIndex.from_frame(table[particles]))
    openmm_force = force.force
    for i in np.flatnonzero(switched):
        indices = [int(p) for p in native.loc[i, particles]]
        values = native.loc[i, parameters].copy()
        values[strength] = 0.0
        values = values.tolist()
        if isinstance(openmm_force, openmm.CustomCompoundBondForce):
            openmm_force.setBondParameters(int(i), indices, values)
        elif isinstance(openmm_force, openmm.CustomBondForce):
            openmm_force.setBondParameters(int(i), *indices, values)
        elif isinstance(openmm_force, openmm.CustomAngleForce):
            openmm_force.setAngleParameters(int(i), *indices, values)
        else:
            openmm_force.setTorsionParameters(int(i), *indices, values)
    del force.parameter_table
