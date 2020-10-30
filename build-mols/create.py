import mbuild
import foyer
import pandas as pd


ff = foyer.forcefields.load_OPLSAA()

mol = mbuild.load("CC", smiles=True)
mol_ff = ff.apply(mol)

xyz = mol_ff.coordinates

atom_is = []
atom_js = []
reqs = []
for bond in mol_ff.bonds:
    atom_is.append(bond.atom1.idx) 
    atom_js.append(bond.atom2.idx) 
    reqs.append(bond.type.req) 

df_atoms = pd.DataFrame(columns=["id", "x", "y", "z"])

df_atoms["id"] = range(xyz.shape[0])
df_atoms["x"] = xyz[:,0]
df_atoms["y"] = xyz[:,1]
df_atoms["z"] = xyz[:,2]

df_bonds = pd.DataFrame(columns=["id_i", "id_j", "req"])
df_bonds["id_i"] = atom_is
df_bonds["id_j"] = atom_js
df_bonds["req"] = reqs

df_atoms.to_csv("positions.csv")
df_bonds.to_csv("constraints.csv")
