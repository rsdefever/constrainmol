import parmed
import numpy as np

from constrainmol import ConstrainedMolecule

# Load from gro/topology
system = parmed.load_file("system.top", xyz="system.gro")
system.save("unconstrained.pdb", overwrite=True)

constrained_coordinates = np.zeros(system.coordinates.shape)

unique_res = system.split()
for (res, resids) in unique_res:
    constrain_mol = ConstrainedMolecule(res)
    for resid in resids:
        constrain_mol.update_xyz(system[resid, :].coordinates)
        constrain_mol.solve()
        
        atom_ids = [atom.idx for atom in system.residues[resid].atoms] 
        constrained_coordinates[atom_ids] = constrain_mol.structure.coordinates

system.coordinates = constrained_coordinates
system.save("constrained.pdb", overwrite=True)

