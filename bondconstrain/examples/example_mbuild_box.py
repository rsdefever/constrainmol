import mbuild
import foyer
import bondconstrain


ff = foyer.forcefields.load_OPLSAA()
dee = mbuild.load("CCOCC", smiles=True)
dee_ff = ff.apply(dee)
box = mbuild.fill_box(dee, 500, density=600)
box.save("unconstrained.pdb", overwrite=True)

constrain_mol = bondconstrain.ConstrainedMolecule(dee_ff)
for mol in box.children:
    constrain_mol.update_xyz(mol.xyz * 10)  # nm to angstrom
    constrain_mol.solve()
    mol.xyz = constrain_mol.xyz / 10.0  # angstrom to nm

box.save("constrained.pdb", overwrite=True)
