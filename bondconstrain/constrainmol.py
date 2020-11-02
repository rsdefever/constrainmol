import pyomo.environ as pyo
import parmed
import numpy as np
from copy import deepcopy


class ConstrainMol(object):

    def __init__(self, structure):
        """Initialize the ConstrainedMolecule from a parmed.Structure

        Parameters
        ----------
        structure: parmed.Structure
            parmed structure with coordinates and bond lengths

        Notes
        -----
        The pyomo model is created on initialization
        """

        if not isinstance(structure, parmed.Structure):
            raise TypeError(
                f"structure: {structure} must be a parmed.Structure"
            )
        if len(structure.bonds) == 0:
            raise ValueError(
                f"structure: {structure} contains no bonds"
            )

        # Extract coordinates and bonds
        xyz = structure.coordinates
        constraints = {}
        for bond in structure.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx
            constraints[(idx1, idx2)] = bond.type.req

        self.structure = deepcopy(structure)
        self.model = self._create_model(xyz, constraints)

    def _create_model(self, xyz, constraints):
        """Create a pyomo model to make xyz satisfy bond constraints

        Parameters
        ----------
        xyz: np.ndarray, shape=(n_atoms, 3)
            initial xyz coordinates
        constraints: dict, keys=(idx1, idx2), values=bond length
            dictionary with bond length constraints

        Returns
        -------
        pyomo.ConcreteModel
            the pyomo model
        """

        # Create pyomo model
        m = pyo.ConcreteModel()

        # Create list of atom indicies
        ids = range(xyz.shape[0])
        m.atom_ids = pyo.Set(initialize=ids)

        # Create target values
        m.x_start = pyo.Param(
            m.atom_ids,
            initialize={idx: x for idx, x in zip(ids, xyz[:, 0])},
            within=pyo.Reals,
            mutable=True,
        )
        m.y_start = pyo.Param(
            m.atom_ids,
            initialize={idx: y for idx, y in zip(ids, xyz[:, 1])},
            within=pyo.Reals,
            mutable=True,
        )
        m.z_start = pyo.Param(
            m.atom_ids,
            initialize={idx: z for idx, z in zip(ids, xyz[:, 2])},
            within=pyo.Reals,
            mutable=True,
        )

        # Create position variables
        m.x = pyo.Var(
            m.atom_ids,
            initialize={idx: x for idx, x in zip(ids, xyz[:, 0])},
            within=pyo.Reals
        )
        m.y = pyo.Var(
            m.atom_ids,
            initialize={idx: y for idx, y in zip(ids, xyz[:, 1])},
            within=pyo.Reals
        )
        m.z = pyo.Var(
            m.atom_ids,
            initialize={idx: z for idx, z in zip(ids, xyz[:, 2])},
            within=pyo.Reals
        )

        # Create bond length parameter
        m.bond_lengths = pyo.Param(
            m.atom_ids, m.atom_ids, initialize=constraints
        )

        # Add bond length constraints
        m.eq_calc_bound_length = pyo.Constraint(
            m.atom_ids, m.atom_ids, rule=self._calc_bond_length
        )

        m.obj = pyo.Objective(
            expr=sum(
                (m.x[i] - m.x_start[i])**2 +
                (m.y[i] - m.y_start[i])**2 +
                (m.z[i] - m.z_start[i])**2
                for i in m.atom_ids
            )
        )

        return m

    @staticmethod
    def _calc_bond_length(m, i, j):
        if (i, j) in m.bond_lengths.keys():
            return (
                    (m.x[i] - m.x[j])**2 +
                    (m.y[i] - m.y[j])**2 +
                    (m.z[i] - m.z[j])**2 ==
                    m.bond_lengths[(i, j)]**2
            )
        else:
            return pyo.Constraint.Skip

    def solve(self):
        """Solve the pyomo model to find the constrained coordinates

        Updates the structure.coordinates if solve is successful
        """
        result = pyo.SolverFactory('ipopt').solve(self.model, tee=True)
        success = (
                str(result["Solver"][0]["Termination condition"]) == 'optimal'
        )
        if not success:
            raise ValueError("Optimal solution not found.")
        constrained_xyz = np.zeros((len(self.model.atom_ids), 3))
        for idx in self.model.atom_ids:
            constrained_xyz[idx, 0] = self.model.x[idx].value
            constrained_xyz[idx, 1] = self.model.y[idx].value
            constrained_xyz[idx, 2] = self.model.z[idx].value

        self.structure.coordinates = constrained_xyz

    def update_xyz(self, xyz):
        """
        Update the initial unconstrained coordinates

        Parameters
        ----------
        xyz: np.ndarray, shape=(n_atoms, 3)
            the new coordinates

        Returns
        -------

        """
        try:
            xyz = np.array(xyz, dtype=np.float64)
        except ValueError:
            raise TypeError(f"New xyz: {xyz} must be array-like")

        if xyz.shape != self.structure.coordinates.shape:
            raise ValueError(
                f"New xyz shape: {xyz.shape} does not match the coordinates "
                f"used to instantiate the class: "
                f"{self.structure.coordinates.shape}."
            )
        self.structure.coordinates = xyz
        for idx in range(xyz.shape[0]):
            self.model.x_start[idx] = xyz[idx, 0]
            self.model.y_start[idx] = xyz[idx, 1]
            self.model.z_start[idx] = xyz[idx, 2]
            self.model.x[idx] = xyz[idx, 0]
            self.model.y[idx] = xyz[idx, 1]
            self.model.z[idx] = xyz[idx, 2]
