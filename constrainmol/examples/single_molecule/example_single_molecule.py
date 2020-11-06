import mbuild
import foyer
from constrainmol import ConstrainedMolecule


ff = foyer.forcefields.load_OPLSAA()
mol = mbuild.load("c1ccccc1O", smiles=True)
mol_ff = ff.apply(mol)
mol_ff.save("unconstrained.pdb", overwrite=True)

constrain_mol = ConstrainedMolecule(mol_ff)
constrain_mol.solve()

# Before we update the coordinates let's see
# how much they changed
diff = mol_ff.coordinates - constrain_mol.xyz
print("Final (x,y,z) - initial (x,y,z)")
print(diff)
mol_ff.coordinates = constrain_mol.xyz
mol_ff.save("constrained.pdb", overwrite=True)

