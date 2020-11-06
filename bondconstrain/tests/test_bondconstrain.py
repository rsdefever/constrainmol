import pytest
import parmed
import pyomo.environ as pyo
import numpy as np

from copy import deepcopy
from bondconstrain.tests.base_test import BaseTest
from bondconstrain import ConstrainedMolecule


class TestConstrainedMolecule(BaseTest):
    def test_invalid_init(self):
        with pytest.raises(TypeError):
            ConstrainedMolecule("structure")

    def test_no_bonds(self):
        with pytest.raises(ValueError):
            ConstrainedMolecule(parmed.Structure())

    def test_copy_coordinates(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        assert np.allclose(constrain_mol.structure.coordinates, ethane_ua.coordinates)

    def test_copy_bonds(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        assert len(constrain_mol.structure.bonds) == len(ethane_ua.bonds)

    def test_create_model(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        assert isinstance(constrain_mol.model, pyo.ConcreteModel)

    def test_model_solved(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        assert constrain_mol.model_solved is False

    def test_model_xyz(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        assert np.allclose(constrain_mol.model.x_start[0].value, ethane_ua.coordinates[0, 0])
        assert np.allclose(constrain_mol.model.x_start[1].value, ethane_ua.coordinates[1, 0])
        assert np.allclose(constrain_mol.model.y_start[0].value, ethane_ua.coordinates[0, 1])
        assert np.allclose(constrain_mol.model.y_start[1].value, ethane_ua.coordinates[1, 1])
        assert np.allclose(constrain_mol.model.z_start[1].value, ethane_ua.coordinates[1, 2])
        assert np.allclose(constrain_mol.model.z_start[0].value, ethane_ua.coordinates[0, 2])
        assert np.allclose(constrain_mol.model.x[0].value, ethane_ua.coordinates[0, 0])
        assert np.allclose(constrain_mol.model.x[1].value, ethane_ua.coordinates[1, 0])
        assert np.allclose(constrain_mol.model.y[0].value, ethane_ua.coordinates[0, 1])
        assert np.allclose(constrain_mol.model.y[1].value, ethane_ua.coordinates[1, 1])
        assert np.allclose(constrain_mol.model.z[0].value, ethane_ua.coordinates[0, 2])
        assert np.allclose(constrain_mol.model.z[1].value, ethane_ua.coordinates[1, 2])

    def test_model_constraints(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        assert np.allclose(constrain_mol.model.bond_lengths[(0, 1)], ethane_ua.bonds[0].type.req)

    def test_solved_model(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        constrain_mol.solve()
        assert np.allclose(constrain_mol.model.x_start[0].value, ethane_ua.coordinates[0, 0])
        assert np.allclose(constrain_mol.model.x_start[1].value, ethane_ua.coordinates[1, 0])
        assert np.allclose(constrain_mol.model.y_start[0].value, ethane_ua.coordinates[0, 1])
        assert np.allclose(constrain_mol.model.y_start[1].value, ethane_ua.coordinates[1, 1])
        assert np.allclose(constrain_mol.model.z_start[0].value, ethane_ua.coordinates[0, 2])
        assert np.allclose(constrain_mol.model.z_start[1].value, ethane_ua.coordinates[1, 2])
        assert np.allclose(constrain_mol.model.x[0].value, ethane_ua.coordinates[0, 0])
        assert np.allclose(constrain_mol.model.x[1].value, ethane_ua.coordinates[1, 0])
        assert np.allclose(constrain_mol.model.y[0].value, ethane_ua.coordinates[0, 1])
        assert np.allclose(constrain_mol.model.y[1].value, ethane_ua.coordinates[1, 1])
        assert np.allclose(constrain_mol.model.z[0].value, ethane_ua.coordinates[0, 2])
        assert np.allclose(constrain_mol.model.z[1].value, ethane_ua.coordinates[1, 2])
        assert constrain_mol.model_solved is True

    def test_invalid_update(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        with pytest.raises(TypeError):
            constrain_mol.update_xyz("string")
        with pytest.raises(ValueError):
            constrain_mol.update_xyz([[1.0, 1.0, 1.0]])

    def test_update_xyz(self, ethane_ua):
        constrain_mol = ConstrainedMolecule(ethane_ua)
        ethane_ua.coordinates[0, 0] = -1.0
        ethane_ua.coordinates[0, 1] = -0.3
        ethane_ua.coordinates[0, 2] = 0.2
        constrain_mol.update_xyz(ethane_ua.coordinates)
        assert np.allclose(constrain_mol.structure.coordinates, ethane_ua.coordinates)
        assert np.allclose(constrain_mol.model.x_start[0].value, ethane_ua.coordinates[0, 0])
        assert np.allclose(constrain_mol.model.x_start[1].value, ethane_ua.coordinates[1, 0])
        assert np.allclose(constrain_mol.model.y_start[0].value, ethane_ua.coordinates[0, 1])
        assert np.allclose(constrain_mol.model.y_start[1].value, ethane_ua.coordinates[1, 1])
        assert np.allclose(constrain_mol.model.z_start[0].value, ethane_ua.coordinates[0, 2])
        assert np.allclose(constrain_mol.model.z_start[1].value, ethane_ua.coordinates[1, 2])
        assert np.allclose(constrain_mol.model.x[0].value, ethane_ua.coordinates[0, 0])
        assert np.allclose(constrain_mol.model.x[1].value, ethane_ua.coordinates[1, 0])
        assert np.allclose(constrain_mol.model.y[0].value, ethane_ua.coordinates[0, 1])
        assert np.allclose(constrain_mol.model.y[1].value, ethane_ua.coordinates[1, 1])
        assert np.allclose(constrain_mol.model.z[0].value, ethane_ua.coordinates[0, 2])
        assert np.allclose(constrain_mol.model.z[1].value, ethane_ua.coordinates[1, 2])

    def test_resolve_model(self, propane_ua):
        constrain_mol = ConstrainedMolecule(propane_ua)
        constrain_mol.solve()
        assert constrain_mol.model_solved is True
        propane_solved = deepcopy(constrain_mol.structure)
        constrain_mol.update_xyz(propane_ua.coordinates)
        assert constrain_mol.model_solved is False
        constrain_mol.solve()
        assert constrain_mol.model_solved is True
        assert np.allclose(propane_solved.coordinates, constrain_mol.structure.coordinates)

    def test_dimethylether(self, dimehtylether_oplsaa):
        constrain_mol = ConstrainedMolecule(dimehtylether_oplsaa)
        constrain_mol.solve()
        optimized = constrain_mol.structure
        xyz = optimized.coordinates
        for bond in optimized.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx
            dist = np.sqrt(np.sum((xyz[idx2] - xyz[idx1])**2))
            assert np.allclose(dist, bond.type.req)

    def test_benzene(self, benzene_oplsaa):
        constrain_mol = ConstrainedMolecule(benzene_oplsaa)
        constrain_mol.solve()
        optimized = constrain_mol.structure
        xyz = optimized.coordinates
        for bond in optimized.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx
            dist = np.sqrt(np.sum((xyz[idx2] - xyz[idx1])**2))
            assert np.allclose(dist, bond.type.req)

    def test_box_diethylether(self, diethylether_box):
        (dee_ff, box) = diethylether_box
        constrain_mol = ConstrainedMolecule(dee_ff)
        for mol in box.children:
            constrain_mol.update_xyz(mol.xyz * 10)  # nm to angstrom
            constrain_mol.solve()
            optimized = constrain_mol.structure
            xyz = optimized.coordinates
            for bond in optimized.bonds:
                idx1 = bond.atom1.idx
                idx2 = bond.atom2.idx
                dist = np.sqrt(np.sum((xyz[idx2] - xyz[idx1])**2))
                assert np.allclose(dist, bond.type.req)

    def test_xyz_getter(self, benzene_oplsaa):
        mol = ConstrainedMolecule(benzene_oplsaa)
        assert np.allclose(mol.xyz, benzene_oplsaa.coordinates)

    def test_xyz_setter(self, benzene_oplsaa):
        mol = ConstrainedMolecule(benzene_oplsaa)
        new_xyz = np.zeros(benzene_oplsaa.coordinates.shape)
        with pytest.raises(ValueError):
            mol.xyz = new_xyz
