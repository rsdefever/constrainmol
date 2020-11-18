import pytest
import mbuild
import foyer
import parmed
import numpy as np


class BaseTest:

    @pytest.fixture
    def ethane_ua(self):
        ethane = parmed.Structure()
        a1 = parmed.topologyobjects.Atom(name="C", atomic_number=6)
        a2 = parmed.topologyobjects.Atom(name="C", atomic_number=6)
        ethane.add_atom(a1, resname="RES", resnum=1)
        ethane.add_atom(a2, resname="RES", resnum=1)
        bond_type = parmed.topologyobjects.BondType(1.0, 1.0)
        bond = parmed.topologyobjects.Bond(a1, a2, type=bond_type)
        ethane.bonds.append(bond)
        ethane.bond_types.append(bond_type)
        ethane.coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        return ethane

    @pytest.fixture
    def propane_ua(self):
        propane = parmed.Structure()
        a1 = parmed.topologyobjects.Atom(name="C", atomic_number=6)
        a2 = parmed.topologyobjects.Atom(name="C", atomic_number=6)
        a3 = parmed.topologyobjects.Atom(name="C", atomic_number=6)
        propane.add_atom(a1, resname="RES", resnum=1)
        propane.add_atom(a2, resname="RES", resnum=1)
        propane.add_atom(a3, resname="RES", resnum=1)
        bond_type = parmed.topologyobjects.BondType(1.5, 1.5)
        b1 = parmed.topologyobjects.Bond(a1, a2, type=bond_type)
        b2 = parmed.topologyobjects.Bond(a2, a3, type=bond_type)
        angle_type = parmed.topologyobjects.AngleType(1.0, 120)
        ang1 = parmed.topologyobjects.Angle(a1, a2, a3, type=angle_type)
        propane.bonds.append(b1)
        propane.bonds.append(b2)
        propane.angles.append(ang1)
        propane.bond_types.append(bond_type)
        propane.angle_types.append(angle_type)
        propane.coordinates = np.array(
            [
                [0.0, 0.1, 0.2],
                [1.0, 0.3, 0.1],
                [2.0, 0.4, 0.2]
            ]
        )
        return propane

    @pytest.fixture
    def water_spce(self):
        ff = foyer.forcefields.load_OPLSAA()
        water = mbuild.load("O", smiles=True)
        water = ff.apply(water)
        return water

    @pytest.fixture
    def dimehtylether_oplsaa(self):
        ff = foyer.forcefields.load_OPLSAA()
        dme = mbuild.load("COC", smiles=True)
        dme = ff.apply(dme)
        return dme

    @pytest.fixture
    def benzene_oplsaa(self):
        ff = foyer.forcefields.load_OPLSAA()
        benzene = mbuild.load("c1ccccc1", smiles=True)
        benzene = ff.apply(benzene)
        return benzene

    @pytest.fixture
    def diethylether_box(self):
        ff = foyer.forcefields.load_OPLSAA()
        dee = mbuild.load("CCOCC", smiles=True)
        dee_ff = ff.apply(dee)
        box = mbuild.fill_box(dee, 50, density=600)
        return dee_ff, box
