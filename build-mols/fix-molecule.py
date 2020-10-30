import mbuild
import foyer
import pyomo.environ as pyo
import numpy as np

# Load/apply force field and save initial structure
ff = foyer.forcefields.load_OPLSAA()
mol = mbuild.load("C1CCCCC1CCOC", smiles=True)
mol_ff = ff.apply(mol)
mol_ff.save("unconstrained.pdb", overwrite=True)

# Extract coordinates and bonds
xyz = mol_ff.coordinates
idxs = list(range(len(mol_ff.atoms)))
atom_is = []
atom_js = []
reqs = []
for bond in mol_ff.bonds:
    atom_is.append(bond.atom1.idx) 
    atom_js.append(bond.atom2.idx) 
    reqs.append(bond.type.req)
    

# Modify one of the bond lengths to make the changes larger
reqs[2] = 2.0

#######################################
################ Pyomo ################

# Create pyomo model
m = pyo.ConcreteModel()

# Create list of atom indicies
m.ATOMS = pyo.Set(initialize=idxs)

# Create target values
m.x_start = pyo.Param(m.ATOMS, initialize={idx: x for idx, x in zip(idxs, xyz[:,0])},within=pyo.Reals)
# , mutable = True
m.y_start = pyo.Param(m.ATOMS, initialize={idx: x for idx, x in zip(idxs, xyz[:,1])},within=pyo.Reals)
m.z_start = pyo.Param(m.ATOMS, initialize={idx: x for idx, x in zip(idxs, xyz[:,2])},within=pyo.Reals)

# Create bond lengths
I = atom_is
J = atom_js
L = reqs

# Construct an empty dictionary
B = {}

# Loop over bonds
for n in range(len(I)):
    B[(I[n],J[n])] = L[n]

# Create bond length set
# m.BONDS = pyo.Set(B.keys())    

# Create bond length parameter
m.bond_lengths = pyo.Param(m.ATOMS, m.ATOMS, initialize=B)

# Create position variables
m.x = pyo.Var(m.ATOMS, initialize={idx: x for idx, x in zip(idxs, xyz[:,0])},within=pyo.Reals)
m.y = pyo.Var(m.ATOMS, initialize={idx: x for idx, x in zip(idxs, xyz[:,1])},within=pyo.Reals)
m.z = pyo.Var(m.ATOMS, initialize={idx: x for idx, x in zip(idxs, xyz[:,2])},within=pyo.Reals)

# Add bond length constraint
# TODO: define this over the set BONDS
def calc_bond_length(m, i, j):
    if (i, j) in B.keys():
        return (m.x[i] - m.x[j])**2 + (m.y[i] - m.y[j])**2 + (m.z[i] - m.z[j])**2 == m.bond_lengths[(i,j)]**2
    else:
        return pyo.Constraint.Skip

m.eq_calc_bound_length = pyo.Constraint(m.ATOMS, m.ATOMS, rule=calc_bond_length)


m.obj = pyo.Objective(expr = sum((m.x[i] - m.x_start[i])**2
                                 + (m.y[i] - m.y_start[i])**2 
                                 + (m.z[i] - m.z_start[i])**2 
                                 for i in m.ATOMS))

results = pyo.SolverFactory('ipopt').solve(m,tee=True)

############## End Pyomo ##############
#######################################


# Save the results
new_xyz = np.zeros(xyz.shape)
for i in m.ATOMS:
    new_xyz[i, 0] = m.x[i].value
    new_xyz[i, 1] = m.y[i].value
    new_xyz[i, 2] = m.z[i].value
mol_ff.coordinates = new_xyz

mol_ff.save("constrained.pdb", overwrite=True)
