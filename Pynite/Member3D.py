from __future__ import annotations  # Allows more recent type hints features
from typing import TYPE_CHECKING, Literal, Union, List
from math import isclose

from numpy import array, zeros, add, subtract, matmul, insert, dot, cross, divide, count_nonzero, concatenate
from numpy import linspace, vstack, hstack, allclose, radians, sin, cos
from numpy.linalg import inv, pinv, norm

import Pynite.FixedEndReactions
from Pynite.BeamSegZ import BeamSegZ
from Pynite.BeamSegY import BeamSegY

if TYPE_CHECKING:

    from typing import Dict, List, Tuple, Optional, Any, Literal

    from numpy import float64
    from numpy.typing import NDArray

    from Pynite.FEModel3D import FEModel3D
    from Pynite.Node3D import Node3D
    from Pynite.Material import Material
    from Pynite.Section import Section
    from Pynite.LoadCombo import LoadCombo


class Member3D():
    """
    A class representing a 3D frame element in a finite element model.

    Most users will not need to interface with this class directly. Rather, the physical member class, which inherits from this class and stitches together a seires of colinear `Member3D` objects will be more useful.
    """

    # '__plt' is used to store the 'pyplot' from matplotlib once it gets imported. Setting it to 'None' for now allows us to defer importing it until it's actually needed.
    __plt = None

    def __init__(self, model: FEModel3D, name: str, i_node: Node3D,
                 j_node: Node3D, material_name: str, section_name: str,
                 rotation: float = 0.0, tension_only: bool = False,
                 comp_only: bool = False) -> None:
        """
        Initializes a new member.

        :param model: The finite element model this member belongs to
        :type model: FEModel3D
        :param name: A unique name for the member
        :type name: str
        :param i_node: The element's i-node (start node)
        :type i_node: Node3D
        :param j_node: The element's j-node (end node)
        :type j_node: Node3D
        :param material_name: The name of the material to use for this member
        :type material_name: str
        :param section_name: The name of the section to use for this member
        :type section_name: str
        :param rotation: Member rotation (degrees) about its local x-axis, defaults to 0.0
        :type rotation: float, optional
        :param tension_only: Indicates if the member is tension-only, defaults to False
        :type tension_only: bool, optional
        :param comp_only: Indicates if the member is compression-only, defaults to False
        :type comp_only: bool, optional
        """

        self.name: str = name      # A unique name for the member given by the user
        self.ID: int | None = None        # Unique index number for the member assigned by the program
        self.i_node: Node3D = i_node  # The element's i-node
        self.j_node: Node3D = j_node  # The element's j-node

        try:
            self.material: Material = model.materials[material_name]  # The element's material
        except KeyError:
            raise NameError(f"No material named '{material_name}'")

        try:
            self.section: Section = model.sections[section_name]  # The element's section
        except KeyError:
            raise NameError(f"No section names '{section_name}'")

        # Variables used to track nonlinear material member end forces
        self._fxi: float = 0
        self._myi: float = 0
        self._mzi: float = 0
        self._fxj: float = 0
        self._myj: float = 0
        self._mzj: float = 0

        # Variable used to track plastic load reveral
        self.i_reversal: bool = False
        self.j_reversal: bool = False

        self.rotation: float = rotation  # Member rotation (degrees) about its local x-axis
        self.PtLoads: List[Tuple] = []  # A list of point loads & moments applied to the element (Direction, P, x, case='Case 1') or (Direction, M, x, case='Case 1')
        self.DistLoads: List[Tuple] = []       # A list of linear distributed loads applied to the element (Direction, w1, w2, x1, x2, case='Case 1', self_weight=False)
        self.SegmentsZ: List[BeamSegZ] = []       # A list of mathematically continuous beam segments for z-bending
        self.SegmentsY: List[BeamSegY] = []       # A list of mathematically continuous beam segments for y-bending
        self.SegmentsX: List[BeamSegZ] = []       # A list of mathematically continuous beam segments for torsion
        self.Releases: List[bool] = [False, False, False, False, False, False, False, False, False, False, False, False]
        self.tension_only: bool = tension_only  # Indicates whether the member is tension-only
        self.comp_only: bool = comp_only  # Indicates whether the member is compression-only

        # Members need to track whether they are active or not for any given load combination. They may become inactive for a load combination during a tension/compression-only analysis. This dictionary will be used when the model is solved.
        self.active: Dict[str, bool] = {}  # Key = load combo name, Value = True or False

        # The 'Member3D' object will store results for one load combination at a time. To reduce repetative calculations the '_solved_combo' variable will be used to track whether the member needs to be resegmented before running calculations for any given load combination.
        self._solved_combo: LoadCombo | None = None  # The current solved load combination

        # Members need a link to the model they belong to
        self.model: FEModel3D = model

# %%
    def L(self) -> float:
        """
        Returns the length of the member.

        :return: The length of the member
        :rtype: float
        """

        # Return the distance between the two nodes
        return self.i_node.distance(self.j_node)

# %%
    def _partition_D(self) -> Tuple[List[int], List[int]]:
        """
        Builds lists of unreleased and released degree of freedom indices for the member.

        Returns
        -------
        R1_indices : list
            A list of the indices for the unreleased DOFs
        R2_indices : list
            A list of the indices for the released DOFs
        """

        R1_indices = []
        R2_indices = []
        for i in range(12):
            if self.Releases[i] == False:
                R1_indices.append(i)
            else:
                R2_indices.append(i)

        return R1_indices, R2_indices

# %%
    def k(self) -> NDArray[Any]:
        """
        Returns the condensed (and expanded) local stiffness matrix for the member.

        :return: The condensed local stiffness matrix
        :rtype: ndarray
        """

        # Partition the local stiffness matrix as 4 submatrices in
        # preparation for static condensation
        k11, k12, k21, k22 = self._partition(self._k_unc())

        # Calculate the condensed local stiffness matrix
        k_Condensed = subtract(k11, matmul(matmul(k12, inv(k22)), k21))

        # Expand the condensed local stiffness matrix
        i = 0
        for DOF in self.Releases:

            if DOF == True:
                k_Condensed = insert(k_Condensed, i, 0, axis=0)
                k_Condensed = insert(k_Condensed, i, 0, axis=1)

            i += 1

        # Return the local stiffness matrix, with end releases applied
        return k_Condensed

# %%
    def _k_unc(self) -> NDArray[float64]:
        """
        Returns the uncondensed local stiffness matrix for the member.

        :return: The uncondensed local stiffness matrix
        :rtype: NDArray[float64]
        """

        # Get the properties needed to form the local stiffness matrix
        E = self.material.E
        G = self.material.G
        Iy = self.section.Iy
        Iz = self.section.Iz
        J = self.section.J
        A = self.section.A
        L = self.L()

        # Create the uncondensed local stiffness matrix
        k = array([[A*E/L,  0,             0,             0,      0,            0,            -A*E/L, 0,             0,             0,      0,            0           ],
                   [0,      12*E*Iz/L**3,  0,             0,      0,            6*E*Iz/L**2,  0,      -12*E*Iz/L**3, 0,             0,      0,            6*E*Iz/L**2 ],
                   [0,      0,             12*E*Iy/L**3,  0,      -6*E*Iy/L**2, 0,            0,      0,             -12*E*Iy/L**3, 0,      -6*E*Iy/L**2, 0           ],
                   [0,      0,             0,             G*J/L,  0,            0,            0,      0,             0,             -G*J/L, 0,            0           ],
                   [0,      0,             -6*E*Iy/L**2,  0,      4*E*Iy/L,     0,            0,      0,             6*E*Iy/L**2,   0,      2*E*Iy/L,     0           ],
                   [0,      6*E*Iz/L**2,   0,             0,      0,            4*E*Iz/L,     0,      -6*E*Iz/L**2,  0,             0,      0,            2*E*Iz/L    ],
                   [-A*E/L, 0,             0,             0,      0,            0,            A*E/L,  0,             0,             0,      0,            0           ],
                   [0,      -12*E*Iz/L**3, 0,             0,      0,            -6*E*Iz/L**2, 0,      12*E*Iz/L**3,  0,             0,      0,            -6*E*Iz/L**2],
                   [0,      0,             -12*E*Iy/L**3, 0,      6*E*Iy/L**2,  0,            0,      0,             12*E*Iy/L**3,  0,      6*E*Iy/L**2,  0           ],
                   [0,      0,             0,             -G*J/L, 0,            0,            0,      0,             0,             G*J/L,  0,            0           ],
                   [0,      0,             -6*E*Iy/L**2,  0,      2*E*Iy/L,     0,            0,      0,             6*E*Iy/L**2,   0,      4*E*Iy/L,     0           ],
                   [0,      6*E*Iz/L**2,   0,             0,      0,            2*E*Iz/L,     0,      -6*E*Iz/L**2,  0,             0,      0,            4*E*Iz/L    ]])

        # Return the uncondensed local stiffness matrix
        return k

    def kg(self, P: float = 0) -> NDArray[float64]:
        """
        Returns the condensed (expanded) local geometric stiffness matrix for the member.

        Parameters
        ----------
        P : number, optional
            The axial force acting on the member (compression = +, tension = -)

        :return: The condensed local geometric stiffness matrix
        :rtype: NDArray[float64]
        """

        # Get the properties needed to form the local geometric stiffness matrix
        Ip = self.section.Iy + self.section.Iz
        A = self.section.A
        L = self.L()

        # Create the uncondensed local geometric stiffness matrix
        kg = array([[1,  0,    0,     0,     0,         0,         -1, 0,     0,    0,     0,         0        ],
                    [0,  6/5,  0,     0,     0,         L/10,      0,  -6/5,  0,    0,     0,         L/10     ],
                    [0,  0,    6/5,   0,     -L/10,     0,         0,  0,     -6/5, 0,     -L/10,     0        ],
                    [0,  0,    0,     Ip/A,  0,         0,         0,  0,     0,    -Ip/A, 0,         0        ],
                    [0,  0,    -L/10, 0,     2*L**2/15, 0,         0,  0,     L/10, 0,     -L**2/30,  0        ],
                    [0,  L/10, 0,     0,     0,         2*L**2/15, 0,  -L/10, 0,    0,     0,         -L**2/30 ],
                    [-1, 0,    0,     0,     0,         0,         1,  0,     0,    0,     0,         0        ],
                    [0,  -6/5, 0,     0,     0,         -L/10,     0,  6/5,   0,    0,     0,         -L/10    ],
                    [0,  0,    -6/5,  0,     L/10,      0,         0,  0,     6/5,  0,     L/10,      0        ],
                    [0,  0,    0,     -Ip/A, 0,         0,         0,  0,     0,    Ip/A,  0,         0        ],
                    [0,  0,    -L/10, 0,     -L**2/30,  0,         0,  0,     L/10, 0,     2*L**2/15, 0        ],
                    [0,  L/10, 0,     0,     0,         -L**2/30,  0,  -L/10, 0,    0,     0,         2*L**2/15]])

        kg = kg*P/L

        # Partition the geometric stiffness matrix as 4 submatrices in
        # preparation for static condensation
        kg11, kg12, kg21, kg22 = self._partition(kg)

        # Calculate the condensed local geometric stiffness matrix
        # Note that a matrix of zeros cannot be inverted, so if P is 0 an error will occur
        if isclose(P, 0.0):
            kg_Condensed = zeros(kg11.shape)
        else:
            kg_Condensed = subtract(kg11, matmul(matmul(kg12, inv(kg22)), kg21))

        # Expand the condensed local geometric stiffness matrix
        i = 0
        for DOF in self.Releases:

            if DOF == True:
                kg_Condensed = insert(kg_Condensed, i, 0, axis=0)
                kg_Condensed = insert(kg_Condensed, i, 0, axis=1)

            i += 1

        # Return the local geomtric stiffness matrix, with end releases applied
        return kg_Condensed

    def km(self, combo_name: str = 'Combo 1') -> NDArray[float64]:
        """Returns the local plastic reduction matrix for the element.

        :param combo_name: The name of the load combination to get the plastic reduction matrix for. Defaults to 'Combo 1'.
        :type combo_name: str, optional
        :return: The plastic reduction matrix for the element
        :rtype: NDArray[float64]
        """

        # List the degrees of freedom associated with axial and bending stiffnesses
        # dofs = [0, 3, 4, 6, 9, 10]

        # Get the elastic local stiffness matrix (for only axial and bending)
        # Note that using the entire stiffness matrix with all terms would lead to an uninvertible term later on
        ke = self.k()  # [dofs][:, dofs]

        # Get the member's axial force
        P = self._fxi

        # Get the geometric local stiffness matrix (for only axial and bending)
        kg = self.kg(P)  # [dofs][:, dofs]

        # Get the total elastic local stiffness matrix
        ke = add(ke, kg)

        # Get the gradient to the failure surface at at each end of the element
        if self.section is None:
            raise Exception(f'Nonlinear material analysis requires member sections to be defined. A section definition is missing for element {self.name}.')
        else:

            # Check for load reversal at the i-node
            if self.i_reversal == True:
                # Gi is a null vector if load reversal is occuring
                Gi = zeros((6, 1))
            else:
                Gi = self.section.G(self._fxi, self._myi, self._mzi)

            # Check for load reversal at the j-node
            if self.j_reversal == True:
                # Gj is a null vector if load reversal is occuring
                Gj = zeros((6, 1))
            else:
                Gj = self.section.G(self._fxj, self._myj, self._mzj)

        # Combine the gradients at the i and j-nodes
        zeros_array = zeros((6, 1))
        Gi = vstack((Gi, zeros_array))
        Gj = vstack((zeros_array, Gj))
        G = hstack((Gi, Gj))

        if self.name == 'M1a':
            M1a_Phi_i = self.section.Phi(self._fxi, self._myi, self._mzi)
            print(f'M1a Phi_i: {M1a_Phi_i}')
            # print(f'M1a Gi: {Gi}')

        # Calculate the plastic reduction matrix for each end of the element
        # TODO: Note that `ke` below already accounts for P-Delta effects and any member end releases which should spill into `km`. I believe end releases will resolve themselves because of this. We'll see how this tests when we get to testing. If it causes problems when end releases are applied we may need to adjust our calculation of G when end releases are present.
        # Check that G is not a zero matrix, which indicates no plastic behavior
        if allclose(G, 0, atol=1e-14):
            return zeros((12, 12))
        else:
            # Solve for `km` using a psuedo-inverse (pinv). The psuedo-inverse takes into account that we may have rows of zeros that make the matrix otherwise uninvertable.
            return -ke @ G @ pinv(G.T @ ke @ G) @ G.T @ ke

    def _m_unc(self, mass_combo_name: str, mass_direction: str = 'Y', gravity: float = 1.0) -> NDArray[float64]:
        """
        Returns the uncondensed mass matrix for the member in local coordinates.

        :param mass_combo_name: Load combination name to define mass via forces.
        :type mass_combo_name: str
        :param mass_direction: Direction for load-to-mass conversion ('X', 'Y', or 'Z'). Any loads applied in this direction (positive or negative) will be converted to mass. Default is 'Y'.
        :type mass_direction: str, optional
        :param gravity: The acceleration due to gravity. Defaults to 1.0. In most cases you'll want to change this to be in units consistent with your model.
        :type gravity: float
        :return: The uncondensed local mass matrix
        :rtype: NDArray[float64]
        """

        # Calculate the non-material load-based mass from the load combination
        load_mass, x = self._calc_load_mass(mass_combo_name, mass_direction, gravity)

        # Calculate the lumped mass from the loads, and the consitent mass from the self-weight loads
        lumped_mass = self.lumped_m(load_mass, x)
        material_mass = self.consistent_m(mass_combo_name, gravity)

        return lumped_mass + material_mass

    def consistent_m(self, mass_combo_name, gravity: float = 1.0) -> NDArray[float64]:

        # Get the section properties needed to form the local mass matrix
        J = self.section.J
        L = self.L()
        A = self.section.A

        # Initialize the material mass to zero
        material_mass = 0.0

        # Get the mass load combination
        mass_combo = self.model.load_combos[mass_combo_name]

        # Step through each member distributed load (looking for self-weight loads)
        for dist_load in self.DistLoads:

            # Extract values from this distributed load
            load_dir, w1, w2, x1, x2, case, self_weight = dist_load

            # Check if this is a self-weight load and if it's part of the mass combo
            if self_weight and case in mass_combo.factors.keys():

                # Find the load factor the user has specified for this load
                factor = mass_combo.factors[case]

                # Calculate the factored mass
                rho = self.material.rho
                material_mass += factor*rho*L*A/gravity

        # Consistent mass matrix for 3D beam element
        #   [dxi     dyi     dzi      rxi      ryi      rzi      dxj  dyj     dzj    rxj      ryj      rzj   ]
        m_coeff = array([
            [140,    0,      0,       0,       0,       0,       70,  0,      0,     0,       0,       0     ],
            [0,      156,    0,       0,       0,       22*L,    0,   54,     0,     0,       0,       -13*L ],
            [0,      0,      156,     0,       -22*L,   0,       0,   0,      54,    0,       13*L,    0     ],
            [0,      0,      0,       140*J/A, 0,       0,       0,   0,      0,     70*J/A,  0,       0     ],
            [0,      0,      -22*L,   0,       4*L**2,  0,       0,   0,      -13*L, 0,       -3*L**2, 0     ],
            [0,      22*L,   0,       0,       0,       4*L**2,  0,   13*L,   0,     0,       0,       -3*L**2],
            [70,     0,      0,       0,       0,       0,       140, 0,      0,     0,       0,       0     ],
            [0,      54,     0,       0,       0,       13*L,    0,   156,    0,     0,       0,       -22*L ],
            [0,      0,      54,      0,       -13*L,   0,       0,   0,      156,   0,       22*L,    0     ],
            [0,      0,      0,       70*J/A,  0,       0,       0,   0,      0,     140*J/A, 0,       0     ],
            [0,      0,      13*L,    0,       -3*L**2, 0,       0,   0,      22*L,  0,       4*L**2,  0     ],
            [0,      -13*L,  0,       0,       0,       -3*L**2, 0,   -22*L,  0,     0,       0,       4*L**2]
        ])

        # # Calculate sum of translational coefficients
        # trans_coeff_sum = 0
        # for i in [0, 1, 2, 6, 7, 8]:  # Translational DOFs
        #     for j in range(12):
        #         trans_coeff_sum += m_coeff[i, j]

        # print(f"DEBUG: Total mass = {total_mass:.2f} kg")
        # print(f"DEBUG: Sum of translational coefficients = {trans_coeff_sum}")
        # print(f"DEBUG: Expected sum for 2 nodes = 420")
        # print(f"DEBUG: Ratio = {trans_coeff_sum/420:.3f}")

        m = m_coeff*(abs(material_mass)/420)

        return m

    def lumped_m(self, load_mass, x) -> NDArray[float64]:
        """Lumps the load's mass to each end node based on its position, `x`.

        :param load_mass: The member's total massp from loads.
        :type load_mass: float
        :param x: The location of the load mass along the member.
        :type x: float
        :return: The member's lumped mass matrix for the load mass.
        :rtype: NDArray[float64]
        """

        # Get the properties needed to form the local mass matrix
        L = self.L()
        m = zeros((12, 12))

        # Distribute half mass to each node's translational DOFs
        # TODO: Distribute the mass based on distance from each load instead
        i_node_mass = load_mass*(1 - x)/L
        j_node_mass = load_mass*x/L
        m[0, 0] = i_node_mass   # FX i-node
        m[1, 1] = i_node_mass   # FY i-node  
        m[2, 2] = i_node_mass   # FZ i-node
        m[6, 6] = j_node_mass   # FX j-node
        m[7, 7] = j_node_mass   # FY j-node
        m[8, 8] = j_node_mass   # FZ j-node

        # Add rotational inertia for the mass
        # Calculate the rotational inertia: I = m*x²
        i_node_rot_inertia = i_node_mass*x**2
        j_node_rot_inertia = j_node_mass*(L - x)**2

        # Apply rotational inertia about the major and minor axes
        # There is no rotational inertia about the torsional axis (zero lever arm)
        # Any numerical instability from zeros in the torsional direction will be fixed by the global solver
        m[4, 4] = i_node_rot_inertia    # RZ i-node
        m[5, 5] = i_node_rot_inertia    # RZ i-node

        m[10, 10] = j_node_rot_inertia    # RZ j-node
        m[11, 11] = j_node_rot_inertia  # RZ j-node

        return m

    def _calc_load_mass(self, mass_combo_name: str, mass_direction: str = 'Y', gravity: float = 1.0) -> float:
        """ Calculates the total mass from a load combination in a specific local direction.

        :param mass_combo_name: Name of the load combination to use for mass calculation.
        :type mass_combo_name: str
        :param mass_direction: Local direction component to extract: 'X', 'Y', or 'Z'. Defaults to 'Y'.
        :type mass_direction: str, optional
        :param gravity: The acceleration due to gravity used for load to mass conversion. Defaults to 1.0.
        :type gravity: float, optional
        :raises NameError: Occurs if `mass_combo_name` is invalid.
        :return: Total mass in the specified local direction (absolute value) and it's location along the member
        :rtype: List[float, float]
        """

        # Get the mass load combination
        try:
            mass_combo = self.model.load_combos[mass_combo_name]
        except KeyError:
            raise NameError(f"No load combination named '{mass_combo_name}'")

        # Initialize the total force to zero
        sum_force = 0.0
        sum_force_x = 0.0

        # Get the transformation matrix once for efficiency
        T_local = self.T()[:3, :3]  # 3x3 rotation matrix

        # Define vector for the mass direction
        if mass_direction == 'X':
            m_vector = array([1.0, 0.0, 0.0])
        elif mass_direction == 'Y':
            m_vector = array([0.0, 1.0, 0.0])
        elif mass_direction == 'Z':
            m_vector = array([0.0, 0.0, 1.0])

        # Sum forces from point loads
        # Step through each point load in the member
        for pt_load in self.PtLoads:

            # Step through each load case and load factor in the mass load combo
            for case, factor in mass_combo.factors.items():

                # Retrive the load's components for clearer reference below
                load_dir, P, x, load_case = pt_load

                # Check if this point load is used in this load combo
                if load_case == case:

                    # Convert the point load to a global vector
                    if load_dir == 'Fx':
                        P_global = T_local.T() @ array([P, 0, 0]) @ T_local
                        P_global = T_local.T() @ array([P, 0, 0]) @ T_local
                    elif load_dir == 'Fy':
                        P_global = T_local.T() @ array([0, P, 0]) @ T_local
                        P_global = T_local.T() @ array([0, P, 0]) @ T_local
                    elif load_dir == 'Fz':
                        P_global = T_local.T() @ array([0, 0, P]) @ T_local
                        P_global = T_local.T() @ array([0, 0, P]) @ T_local
                    elif load_dir == 'FX':
                        P_global = array([P, 0, 0])
                        P_global = array([P, 0, 0])
                    elif load_dir == 'FY':
                        P_global = array([0, P, 0])
                        P_global = array([0, P, 0])
                    elif load_dir == 'FZ':
                        P_global = array([0, 0, P])
                        P_global = array([0, 0, P])
                    else:
                        # Assume zero for any other load directions
                        P_global = array([0, 0, 0])

                    # Calculate the load component acting in the mass direction
                    P_m_comp = dot(P_global, m_vector)

                    # Sum the total for the load component
                    sum_force += factor*P_m_comp
                    sum_force_x += factor*P_m_comp*x

        # Sum forces from distributed loads
        for dist_load in self.DistLoads:

            # Retrive the load's components for clearer reference below
            load_dir, w1, w2, x1, x2, load_case, self_weight = dist_load

            # Don't add masses for self-weight loads. They will be added elsewhere using a consistent mass matrix, rather than a lumped mass matrix
            if not self_weight:

                # Step through each load case and factor in the mass combo
                for case, factor in mass_combo.factors.items():

                    # Check if this load's case is in the mass combo
                    if load_case == case:

                        # Calculate the length of the distributed load
                        length_loaded = x2 - x1

                        # Convert the distributed load to a global vector
                        if load_dir == 'Fx':
                            w1_global = T_local.T() @ array([w1, 0, 0]) @ T_local
                            w2_global = T_local.T() @ array([w2, 0, 0]) @ T_local
                        elif load_dir == 'Fy':
                            w1_global = T_local.T() @ array([0, w1, 0]) @ T_local
                            w2_global = T_local.T() @ array([0, w2, 0]) @ T_local
                        elif load_dir == 'Fz':
                            w1_global = T_local.T() @ array([0, 0, w1]) @ T_local
                            w2_global = T_local.T() @ array([0, 0, w2]) @ T_local
                        elif load_dir == 'FX':
                            w1_global = array([w1, 0, 0])
                            w2_global = array([w2, 0, 0])
                        elif load_dir == 'FY':
                            w1_global = array([0, w1, 0])
                            w2_global = array([0, w2, 0])
                        elif load_dir == 'FZ':
                            w1_global = array([0, 0, w1])
                            w2_global = array([0, 0, w2])
                        else:
                            # Assume zero for any other load directions
                            w1_global = array([0, 0, 0])
                            w2_global = array([0, 0, 0])

                        # Calculate the load component acting in the mass direction
                        w1_m_comp = dot(w1_global, m_vector)
                        w2_m_comp = dot(w2_global, m_vector)

                        # Sum the average for the load component
                        avg_load = (w1_m_comp + w2_m_comp)/2
                        sum_force += factor*avg_load*length_loaded

                        # Identify the point through which the load acts
                        sum_force_x += factor*avg_load*(w2 - w1)*(w1_m_comp*(x2 - x1)**2/2
                                                                  + 0.5*(w2_m_comp - w1_m_comp)*(x2 - x1)**2*(2/3))

        # Identify the load's center of gravity
        if sum_force_x != 0.0:
            total_x = sum_force_x/sum_force
        else:
            total_x = self.L()/2

        return [abs(sum_force/gravity), total_x]  # Mass is always positive

    def m(self, mass_combo_name: str, mass_direction: str = 'Y', gravity: float = 1.0) -> NDArray[Any]:
        """
        Returns the condensed (and expanded) local mass matrix for the member. Condensing the matrix is used to account for member end releases.

        :param mass_combo_name: Name of the load combination used to define mass
        :param mass_direction: Direction for load-to-mass conversion ('X', 'Y', or 'Z'). Any loads applied in this direction (positive or negative) will be converted to mass. Default is 'Y'.
        :type mass_direction: str, optional
        :param gravity: Acceleration due to gravity. Default is 1.0.
        :type gravity: float, optional
        :return: The condensed local mass matrix
        :rtype: ndarray
        """

        # Check if there are any member end releases
        if True not in self.Releases:

            # If no releases, return the full uncondensed mass matrix
            return self._m_unc(mass_combo_name, mass_direction, gravity)

        # Partition the local mass matrix as 4 submatrices in
        # preparation for static condensation (same as stiffness matrix)
        R1_indices, R2_indices = self._partition_D()

        # Get the uncondensed mass matrix
        m_unc = self._m_unc(mass_combo_name, mass_direction, gravity)

        # Partition the mass matrix
        m11 = m_unc[R1_indices, :][:, R1_indices]
        m12 = m_unc[R1_indices, :][:, R2_indices] 
        m21 = m_unc[R2_indices, :][:, R1_indices]
        m22 = m_unc[R2_indices, :][:, R2_indices]

        # Calculate the condensed local mass matrix
        m_condensed = m11 - m12 @ pinv(m22) @ m21

        # Expand the condensed local mass matrix back to full size
        m_expanded = zeros((12, 12))
        for i, idx_i in enumerate(R1_indices):
            for j, idx_j in enumerate(R1_indices):
                m_expanded[idx_i, idx_j] = m_condensed[i, j]

        return m_expanded

    def M(self, mass_combo_name: str | None = None, mass_direction: str = 'Y', gravity: float = 1.0) -> NDArray[Any]:
        """Returns the member's global mass matrix.

        The mass is created from loads in the mass combination.

        :param mass_combo_name: Load combination name for force-based mass calculation, defaults to ""
        :type mass_combo_name: str, optional
        :param mass_direction: Direction for load-to-mass conversion ('X', 'Y', or 'Z'). Any loads applied in this direction (positive or negative) will be converted to mass. Default is 'Y'.
        :type mass_direction: str, optional
        :return: Global mass matrix of shape (12, 12)
        :rtype: numpy.ndarray
        """

        # Get the member's local mass matrix
        m = self.m(mass_combo_name, mass_direction, gravity)

        # Get the member's transformation matrix
        T = self.T()

        # Calculate the global mass matrix
        # M_global = T^T * m_local * T
        M_global = T.T @ m @ T

        return M_global

    def lamb(self, model_Delta_D: NDArray[float64], combo_name: str = 'Combo 1', push_combo: str = 'Push', step_num: int = 1) -> NDArray[float64]:
        """
        Returns the `lambda` vector used in pushover analysis.

        `lambda` is a vector representing the magnitude of the plastic deformations in the member.

        :param model_Delta_D: The change in the global displacement vector calculated from the latest load step
        :type model_Delta_D: ndarray
        :param combo_name: The load combination name, defaults to 'Combo 1'
        :type combo_name: str, optional
        :param push_combo: The pushover load combination name, defaults to 'Push'
        :type push_combo: str, optional
        :param step_num: The current load step number, defaults to 1
        :type step_num: int, optional
        :return: The lambda vector
        :rtype: NDArray[float64]
        """

        # Obtain the change in the member's end displacements from the calculated displacement change vector
        Delta_D = array([model_Delta_D[self.i_node.ID*6 + 0],
                         model_Delta_D[self.i_node.ID*6 + 1],
                         model_Delta_D[self.i_node.ID*6 + 2],
                         model_Delta_D[self.i_node.ID*6 + 3],
                         model_Delta_D[self.i_node.ID*6 + 4],
                         model_Delta_D[self.i_node.ID*6 + 5],
                         model_Delta_D[self.j_node.ID*6 + 0],
                         model_Delta_D[self.j_node.ID*6 + 1],
                         model_Delta_D[self.j_node.ID*6 + 2],
                         model_Delta_D[self.j_node.ID*6 + 3],
                         model_Delta_D[self.j_node.ID*6 + 4],
                         model_Delta_D[self.j_node.ID*6 + 5]]).reshape(12, 1)

        # Convert the global changes in displacement to local coordinates
        Delta_d = self.T() @ Delta_D

        # Get the elastic local stiffness matrix (includeing goemetric stiffness)
        d_total = self.d(combo_name)  # Total displacements acting on the member at the current load stp
        delta_dx_total = d_total[6, 0] - d_total[0, 0]  # Change in displacement across the lenght of the member
        P = self.section.A*self.material.E/self.L()*delta_dx_total  # Axial load acting on the member at the current load step
        ke = self.k() + self.kg(P)  # Elastic stiffness (including geometric stiffness)

        # Get the gradient to the failure surface at at each end of the element
        if self.section is None:
            raise Exception(f'Nonlinear material analysis requires member sections to be defined. A section definition is missing for element {self.name}.')
        else:
            Gi = self.section.G(self._fxi, self._myi, self._mzi)
            Gj = self.section.G(self._fxj, self._myj, self._mzj)

        # Combine the gradients for the i and j-nodes
        zeros_array = zeros((6, 1))
        Gi = vstack((Gi, zeros_array))
        Gj = vstack((zeros_array, Gj))
        G = hstack((Gi, Gj))

        # Check if all terms in [G] are zero
        if allclose(G, 0, atol=1e-14):
            # No plasticity is occuring, so `lambda` is a 2x1 zero matrix
            return array([[0], [0]])
        else:
            # Note: `pinv` accounts for rows of zeros that would normally make the matrix singular
            return pinv(G.T @ ke @ G) @ G.T @ ke @ Delta_d

    def fer(self, combo_name: str = 'Combo 1') -> NDArray[float64]:
        """
        Returns the condensed (and expanded) local fixed end reaction vector for the member for the given load combination.

        :param combo_name: The name of the load combination to construct the fixed end reaction vector for
        :type combo_name: str
        :return: The fixed end reaction vector for the member
        :rtype: NDArray[float64]
        """

        # Get the lists of unreleased and released degree of freedom indices
        R1_indices, R2_indices = self._partition_D()

        # Partition the local stiffness matrix and local fixed end reaction vector
        k11, k12, k21, k22 = self._partition(self._k_unc())
        fer1, fer2 = self._partition(self._fer_unc(combo_name))

        # Calculate the condensed fixed end reaction vector
        ferCondensed = subtract(fer1, matmul(matmul(k12, inv(k22)), fer2))

        # Expand the condensed fixed end reaction vector
        i = 0
        for DOF in self.Releases:

            if DOF == True:
                ferCondensed = insert(ferCondensed, i, 0, axis=0)

            i += 1

        # Return the fixed end reaction vector        
        return ferCondensed

    def _fer_unc(self, combo_name:str = 'Combo 1') -> NDArray[float64]:
        """
        Returns the member's local fixed end reaction vector, ignoring the effects of end releases.
        Needed to apply the slope-deflection equation properly.
        """

        # Initialize the fixed end reaction vector
        fer = zeros((12, 1))

        # Get the requested load combination
        combo = self.model.load_combos[combo_name]

        # Loop through each load case and factor in the load combination
        for case, factor in combo.factors.items():

            # Sum the fixed end reactions for the point loads & moments
            for ptLoad in self.PtLoads:

                # Check if the current point load corresponds to the current load case
                if ptLoad[3] == case:

                    if ptLoad[0] == 'Fx':
                        fer = add(fer, Pynite.FixedEndReactions.FER_AxialPtLoad(factor*ptLoad[1], ptLoad[2], self.L()))
                    elif ptLoad[0] == 'Fy':
                        fer = add(fer, Pynite.FixedEndReactions.FER_PtLoad(factor*ptLoad[1], ptLoad[2], self.L(), 'Fy'))
                    elif ptLoad[0] == 'Fz':
                        fer = add(fer, Pynite.FixedEndReactions.FER_PtLoad(factor*ptLoad[1], ptLoad[2], self.L(), 'Fz'))
                    elif ptLoad[0] == 'Mx':
                        fer = add(fer, Pynite.FixedEndReactions.FER_Torque(factor*ptLoad[1], ptLoad[2], self.L()))
                    elif ptLoad[0] == 'My':
                        fer = add(fer, Pynite.FixedEndReactions.FER_Moment(factor*ptLoad[1], ptLoad[2], self.L(), 'My'))
                    elif ptLoad[0] == 'Mz':     
                        fer = add(fer, Pynite.FixedEndReactions.FER_Moment(factor*ptLoad[1], ptLoad[2], self.L(), 'Mz'))
                    elif ptLoad[0] == 'FX' or ptLoad[0] == 'FY' or ptLoad[0] == 'FZ':
                        FX, FY, FZ = 0, 0, 0
                        if ptLoad[0] == 'FX': FX = 1
                        if ptLoad[0] == 'FY': FY = 1
                        if ptLoad[0] == 'FZ': FZ = 1
                        f = self.T()[:3, :][:, :3] @ array([FX*ptLoad[1], FY*ptLoad[1], FZ*ptLoad[1]])
                        fer = add(fer, Pynite.FixedEndReactions.FER_AxialPtLoad(factor*f[0], ptLoad[2], self.L()))
                        fer = add(fer, Pynite.FixedEndReactions.FER_PtLoad(factor*f[1], ptLoad[2], self.L(), 'Fy'))
                        fer = add(fer, Pynite.FixedEndReactions.FER_PtLoad(factor*f[2], ptLoad[2], self.L(), 'Fz'))
                    elif ptLoad[0] == 'MX' or ptLoad[0] == 'MY' or ptLoad[0] == 'MZ':
                        MX, MY, MZ = 0, 0, 0
                        if ptLoad[0] == 'MX': MX = 1
                        if ptLoad[0] == 'MY': MY = 1
                        if ptLoad[0] == 'MZ': MZ = 1
                        f = self.T()[:3, :][:, :3] @ array([MX*ptLoad[1], MY*ptLoad[1], MZ*ptLoad[1]])
                        fer = add(fer, Pynite.FixedEndReactions.FER_Torque(factor*f[0], ptLoad[2], self.L()))
                        fer = add(fer, Pynite.FixedEndReactions.FER_Moment(factor*f[1], ptLoad[2], self.L(), 'My'))
                        fer = add(fer, Pynite.FixedEndReactions.FER_Moment(factor*f[2], ptLoad[2], self.L(), 'Mz'))
                    else:
                        raise Exception('Invalid member point load direction specified.')

            # Sum the fixed end reactions for the distributed loads
            for distLoad in self.DistLoads:

                # Check if the current distributed load corresponds to the current load case
                if distLoad[5] == case:

                    if distLoad[0] == 'Fx':
                        fer = add(fer, Pynite.FixedEndReactions.FER_AxialLinLoad(factor*distLoad[1], factor*distLoad[2], distLoad[3], distLoad[4], self.L()))
                    elif distLoad[0] == 'Fy' or distLoad[0] == 'Fz':
                        fer = add(fer, Pynite.FixedEndReactions.FER_LinLoad(factor*distLoad[1], factor*distLoad[2], distLoad[3], distLoad[4], self.L(), distLoad[0]))
                    elif distLoad[0] == 'FX' or distLoad[0] == 'FY' or distLoad[0] == 'FZ':
                        FX, FY, FZ = 0, 0, 0
                        if distLoad[0] == 'FX': FX = 1
                        if distLoad[0] == 'FY': FY = 1
                        if distLoad[0] == 'FZ': FZ = 1
                        w1 = self.T()[:3, :][:, :3] @ array([FX*distLoad[1], FY*distLoad[1], FZ*distLoad[1]])
                        w2 = self.T()[:3, :][:, :3] @ array([FX*distLoad[2], FY*distLoad[2], FZ*distLoad[2]])
                        fer = add(fer, Pynite.FixedEndReactions.FER_AxialLinLoad(factor*w1[0], factor*w2[0], distLoad[3], distLoad[4], self.L()))
                        fer = add(fer, Pynite.FixedEndReactions.FER_LinLoad(factor*w1[1], factor*w2[1], distLoad[3], distLoad[4], self.L(), 'Fy'))
                        fer = add(fer, Pynite.FixedEndReactions.FER_LinLoad(factor*w1[2], factor*w2[2], distLoad[3], distLoad[4], self.L(), 'Fz'))

        # Return the fixed end reaction vector, uncondensed
        return fer

    def _partition(self, unp_matrix: NDArray[float64])-> tuple[Any, Any] | tuple[Any, Any, Any, Any]:
        """
        Partitions a matrix into sub-matrices based on unreleased and released degree of freedom indices.
        """

        # Create auxiliary lists of released/unreleased DOFs
        R1_indices, R2_indices = self._partition_D()

        # Partition the matrix by slicing
        if unp_matrix.shape[1] == 1:
            m1 = unp_matrix[R1_indices, :]
            m2 = unp_matrix[R2_indices, :]
            return m1, m2
        else:
            m11 = unp_matrix[R1_indices, :][:, R1_indices]
            m12 = unp_matrix[R1_indices, :][:, R2_indices]
            m21 = unp_matrix[R2_indices, :][:, R1_indices]
            m22 = unp_matrix[R2_indices, :][:, R2_indices]
            return  m11, m12, m21, m22

    def f(self, combo_name: str='Combo 1', push_combo: str = None, step_num: int = None) -> NDArray[float64]:
        """Returns the member's local end force vector for the given load combination.

        :param combo_name: The load combination to get the local end for vector for. Defaults to 'Combo 1'.
        :type combo_name: str, optional
        :return: The member's local end force vector for the given load combination.
        :rtype: array
        """

        # Calculate and return the member's local end force vector
        if self.model.solution == 'P-Delta':

            # Back-calculate the axial force on the member from the axial strain
            P = (self.d(combo_name)[6, 0] - self.d(combo_name)[0, 0])*self.section.A*self.material.E/self.L()

            return add(matmul(add(self.k(), self.kg(P)), self.d(combo_name)), self.fer(combo_name))

        # Check for a pushover analysis
        elif push_combo is not None and step_num is not None:

            # Calculate the axial force on the member from the latest elasto-plastic member end forces
            P = self._fxi

            # Calculate the total stiffness matrix
            kt = self.k() + self.kg(P) + self.km(combo_name)

            # Calculate the fixed end reaction vector for this load step
            fer = self.fer(combo_name) + self.fer(push_combo)*step_num

            # Retern the new member end forces
            return kt @ self.d(combo_name) + fer

        else:

            return self.k() @ self.d(combo_name) + self.fer(combo_name)

    def d(self, combo_name='Combo 1') -> NDArray[float64]:
        """
        Returns the member's local displacement vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to construct the displacement vector for (not the load combination itself).
        """

        # Calculate and return the local displacement vector
        return self.T() @ self.D(combo_name)

    # Transformation matrix
    def T(self) -> NDArray[float64]:
        """
        Returns the transformation matrix for the member.
        """

        # Get the global coordinates for the two ends
        Xi = self.i_node.X
        Xj = self.j_node.X

        Yi = self.i_node.Y
        Yj = self.j_node.Y

        Zi = self.i_node.Z
        Zj = self.j_node.Z

        # Calculate the length of the member
        L = self.L()

        # Calculate the direction cosines for the local x-axis
        x = [(Xj - Xi)/L, (Yj - Yi)/L, (Zj - Zi)/L]

        # Calculate the remaining direction cosines.
        # For now, the local z-axis will be kept parallel to the global XZ plane in all cases. It will be adjusted later if a rotation has been applied to the member.
        # Vertical members
        if isclose(Xi, Xj) and isclose(Zi, Zj):

            # For vertical members, keep the local y-axis in the XY plane to make 2D problems easier to solve in the XY plane
            if Yj > Yi:
                y = [-1, 0, 0]
                z = [0, 0, 1]
            else:
                y = [1, 0, 0]
                z = [0, 0, 1]

        # Horizontal members
        elif isclose(Yi, Yj):

            # Find a vector in the direction of the local z-axis by taking the cross-product
            # of the local x-axis and the local y-axis. This vector will be perpendicular to
            # both the local x-axis and the local y-axis.
            y = [0, 1, 0]
            z = cross(x, y)

            # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
            z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)

        # Members neither vertical or horizontal
        else:

            # Find the projection of x on the global XZ plane
            proj = [Xj - Xi, 0, Zj - Zi]

            # Find a vector in the direction of the local z-axis by taking the cross-product
            # of the local x-axis and its projection on a plane parallel to the XZ plane. This
            # produces a vector perpendicular to both the local x-axis and its projection. This
            # vector will always be horizontal since it's parallel to the XZ plane. The order
            # in which the vectors are 'crossed' has been selected to ensure the y-axis always
            # has an upward component (i.e. the top of the beam is always on top).
            if Yj > Yi:
                z = cross(proj, x)
            else:
                z = cross(x, proj)

            # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
            z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)

            # Find the direction cosines for the local y-axis
            y = cross(z, x)
            y = divide(y, (y[0]**2 + y[1]**2 + y[2]**2)**0.5)

        # Check if the member is rotated
        if self.rotation != 0.0:

            # Convert the rotation angle to radians
            theta = radians(self.rotation)

            # Shorthand cosine and sine functions
            c = cos(theta)
            s = sin(theta)

            # Convert `x` to a numpy array
            u = array(x)

            # Rotate y using the Rodrigues formula
            v = array(y)
            y = v * c + cross(u, v) * s + u * (dot(u, v)) * (1 - c)
            y /= norm(y)

            # Rotate z using the Rodrigues formula
            v = array(z)
            z = v * c + cross(u, v) * s + u * (dot(u, v)) * (1 - c)
            z /= norm(z)

        # Create the direction cosines matrix
        dirCos = array([x, y, z])

        # Build the transformation matrix
        transMatrix = zeros((12, 12))
        transMatrix[0:3, 0:3] = dirCos
        transMatrix[3:6, 3:6] = dirCos
        transMatrix[6:9, 6:9] = dirCos
        transMatrix[9:12, 9:12] = dirCos

        return transMatrix

    # Member global stiffness matrix
    def K(self) -> NDArray[float64]:
        """Returns the global elastic stiffness matrix for the member.

        :return: The global elastic stiffness matrix for the member.
        :rtype: array
        """

        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(inv(self.T()), self.k()), self.T())

    def Kg(self, P: float=0.0):
        """Returns the global geometric stiffness matrix for the member. Used for P-Delta analysis.

        :param P: Member axial load. Defaults to 0.
        :type P: float, optional
        :return: The global geometric stiffness matrix for the member.
        :rtype: array
        """

        # Calculate and return the geometric stiffness matrix in global coordinates
        return matmul(matmul(inv(self.T()), self.kg(P)), self.T())

    def Km(self, combo_name: str) -> NDArray[float64]:
        """Returns the global plastic reduction matrix for the member. Used to modify member behavior for plastic hinges at the ends.

        :param combo_name: The name of the load combination to get the plastic reduction matrix for.
        :type combo_name: str
        :return: The global plastic reduction matrix for the member.
        :rtype: array
        """

        # Calculate and return the plastic reduction matrix in global coordinates
        return matmul(matmul(inv(self.T()), self.km(combo_name)), self.T())

    def F(self, combo_name: str='Combo 1') -> NDArray[float64]:
        """
        Returns the member's global end force vector for the given load combination.
        """

        # Calculate and return the global force vector
        return matmul(inv(self.T()), self.f(combo_name))

    def FER(self, combo_name: str = 'Combo 1') -> NDArray[float64]:
        """
        Returns the global fixed end reaction vector

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the fixed end reaction vector for (not the load combination itself).
        """

        # Calculate and return the fixed end reaction vector
        return matmul(inv(self.T()), self.fer(combo_name))

    def D(self, combo_name: str = 'Combo 1') -> NDArray[float64]:
        """
        Returns the member's global displacement vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to construct the global
            displacement vector for (not the load combination itelf).
        """

        # Initialize the displacement vector
        D = zeros((12, 1))

        # TODO: I'm not sure this next block is the best way to handle inactive members - need to review
        # Read in the global displacements from the nodes
        # Apply axial displacements only if the member is active
        if self.active[combo_name] == True:
            D[0, 0] = self.i_node.DX[combo_name]
            D[6, 0] = self.j_node.DX[combo_name]

        # Apply the remaining displacements
        D[1, 0] = self.i_node.DY[combo_name]
        D[2, 0] = self.i_node.DZ[combo_name]
        D[3, 0] = self.i_node.RX[combo_name]
        D[4, 0] = self.i_node.RY[combo_name]
        D[5, 0] = self.i_node.RZ[combo_name]
        D[7, 0] = self.j_node.DY[combo_name]
        D[8, 0] = self.j_node.DZ[combo_name]
        D[9, 0] = self.j_node.RX[combo_name]
        D[10, 0] = self.j_node.RY[combo_name]
        D[11, 0] = self.j_node.RZ[combo_name]

        # Return the global displacement vector
        return D

    def shear(self, Direction: Literal['Fy', 'Fz'], x: float, combo_name: str = 'Combo 1') -> float:
        """
        Returns the shear at a point along the member's length.

        Parameters
        ----------
        Direction : string
            The direction in which to find the shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis.
                'Fz' = Shear acting on the local z-axis.
        x : number
            The location at which to find the shear.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        # Only calculate results if the member is currently active
        if self.active[combo_name]:

            # Segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Check which direction is of interest
            if Direction == 'Fy':

                # Check which segment 'x' falls on
                for segment in self.SegmentsZ:
                    if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                        return segment.Shear(x - segment.x1)

                if isclose(x, self.L()):  
                    lastIndex = len(self.SegmentsZ) - 1
                    return self.SegmentsZ[lastIndex].Shear(x - self.SegmentsZ[lastIndex].x1)

            elif Direction == 'Fz':

                for segment in self.SegmentsY:

                    if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):

                        return segment.Shear(x - segment.x1)

                if isclose(x, self.L()):

                    lastIndex = len(self.SegmentsY) - 1
                    return self.SegmentsY[lastIndex].Shear(x - self.SegmentsY[lastIndex].x1)

        else:

            return 0

    def max_shear(self, Direction: Literal['Fy', 'Fz'], combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the maximum shear force in the member for the specified direction
        and load combination(s).

        Parameters
        ----------
        Direction : Literal['Fy', 'Fz']
            The direction in which to find the maximum shear:
            - 'Fy' : Shear acting along the local y-axis.
            - 'Fz' : Shear acting along the local z-axis.
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the maximum shear
            across those combos will be returned.

        Returns
        -------
        float
            The maximum shear value in the specified direction for the given
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.

        Notes
        -----
        - Unlike moments, shear does not depend on P-Delta effects.
        - Automatically segments the member if it has not yet been segmented
        for a load combination.
        """
        # Normalize input: single name vs list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]  # treat as a single combo name
        else:
            # Filter combos that contain *any* of the given tags
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        Vmax_global = None  # will store the global maximum across combos

        for combo_name in combo_names:
            # Skip inactive combos
            if not self.active.get(combo_name, False):
                continue

            # If member not yet segmented for this combo, do it
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Select the correct segment list
            segments = self.SegmentsZ if Direction == 'Fy' else self.SegmentsY

            if not segments:
                continue

            # Get the maximum shear for this combo
            Vmax = max(segment.max_shear() for segment in segments)

            # Update global maximum
            if Vmax_global is None or Vmax > Vmax_global:
                Vmax_global = Vmax

        # Return 0 if no valid combos were found
        return Vmax_global if Vmax_global is not None else 0

    def min_shear(self, Direction: Literal['Fy', 'Fz'], combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the minimum shear force in the member for the specified direction
        and load combination(s).

        Parameters
        ----------
        Direction : Literal['Fy', 'Fz']
            The direction in which to find the minimum shear:
            - 'Fy' : Shear acting along the local y-axis.
            - 'Fz' : Shear acting along the local z-axis.
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the minimum shear
            across those combos will be returned.

        Returns
        -------
        float
            The minimum shear value in the specified direction for the given
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.

        Notes
        -----
        - Unlike moments, shear does not depend on P-Delta effects.
        - Automatically segments the member if it has not yet been segmented
        for a load combination.
        """
        # Normalize input: single name vs list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]  # treat as a single combo name
        else:
            # Filter combos that contain *any* of the given tags
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        Vmin_global = None  # will store the global minimum across combos

        for combo_name in combo_names:
            # Skip inactive combos
            if not self.active.get(combo_name, False):
                continue

            # If member not yet segmented for this combo, do it
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Select the correct segment list
            segments = self.SegmentsZ if Direction == 'Fy' else self.SegmentsY

            if not segments:
                continue

            # Get the minimum shear for this combo
            Vmin = min(segment.min_shear() for segment in segments)

            # Update global minimum
            if Vmin_global is None or Vmin < Vmin_global:
                Vmin_global = Vmin

        # Return 0 if no valid combos were found
        return Vmin_global if Vmin_global is not None else 0

    def plot_shear(self, Direction: Literal['Fy', 'Fz'], combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the shear diagram for the member

        Parameters
        ----------
        Direction : string
            The direction in which to find the moment. Must be one of the following:
                'Fy' = Shear acting on the local y-axis.
                'Fz' = Shear acting on the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        x, V = self.shear_array(Direction, n_points, combo_name)

        Member3D.__plt.plot(x, V)
        Member3D.__plt.ylabel('Shear')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.name + '\n' + combo_name)
        Member3D.__plt.show()    

    def shear_array(self, Direction: Literal['Fy', 'Fz'], n_points: int, combo_name='Combo 1', x_array=None) -> NDArray[float64]:
        """
        Returns the array of the shear in the member for the given direction

        Parameters
        ----------
        Direction : string
            The direction to plot the shear for. Must be one of the following:
                'Fy' = Shear acting on the local y-axis.
                'Fz' = Shear acting on the local z-axis.
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated.
            Values must be provided in local member coordinates (between 0 and L) and be in ascending order
        """

        # Segment the member into segments with mathematically continuous loads if not already done
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        L = self.L()
        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        # Check which axis is of interest
        if Direction == 'Fz':
            return self._extract_vector_results(self.SegmentsY, x_array, 'shear')

        elif Direction == 'Fy':
            return self._extract_vector_results(self.SegmentsZ, x_array, 'shear')

        else:
            raise ValueError(f"Direction must be 'Fy' or 'Fz'. {Direction} was given.")

    def moment(self, Direction: Literal['My', 'Mz'], x: float, combo_name: str = 'Combo 1') -> float:
        """
        Returns the moment at a point along the member's length

        Parameters
        ----------
        Direction : string
            The direction in which to find the moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = moment about the local z-axis.
        x : number
            The location at which to find the moment.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        # Only calculate results if the member is currently active
        if self.active[combo_name]:

            # Segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Determine if a P-Delta analysis has been run
            if self.model.solution == 'P-Delta' or self.model.solution == 'Pushover':
                # Include P-little-delta effects in the moment results
                P_delta = True
            else:
                # Do not include P-little delta effects in the moment results
                P_delta = False

            # Check which axis is of interest
            if Direction == 'My':

                # Check which segment 'x' falls on
                for segment in self.SegmentsY:

                    if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):

                        return segment.moment(x - segment.x1, P_delta)

                if isclose(x, self.L()):

                    return self.SegmentsY[-1].moment(x - self.SegmentsY[-1].x1, P_delta)

            elif Direction == 'Mz':

                for segment in self.SegmentsZ:

                    if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):

                        return segment.moment(x - segment.x1, P_delta)

                if isclose(x, self.L()):

                    return self.SegmentsZ[-1].moment(x - self.SegmentsZ[-1].x1, P_delta)

            else:
                raise ValueError(f"Direction must be 'My' or 'Mz'. {Direction} was given.")

        else:

            return 0

    def max_moment(self, Direction: Literal['My', 'Mz'], combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the maximum bending moment in the member for the specified direction
        and load combination(s).

        Parameters
        ----------
        Direction : Literal['My', 'Mz']
            The direction in which to find the maximum moment:
            - 'My' : Moment about the local y-axis.
            - 'Mz' : Moment about the local z-axis.
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the maximum moment
            across those combos will be returned.

        Returns
        -------
        float
            The maximum moment value in the specified direction for the given
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.

        Notes
        -----
        - Includes P-Delta effects if the model's solution is 'P-Delta' or 'Pushover'.
        - Automatically segments the member if it has not yet been segmented
        for a load combination.
        """
        # Normalize input: single name vs list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]  # treat as a single combo name
        else:
            # Filter combos that contain *any* of the given tags
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        Mmax_global = None  # will store the global maximum across combos

        for combo_name in combo_names:
            # Skip inactive combos
            if not self.active.get(combo_name, False):
                continue

            # Determine if P-Delta (or Pushover) effects should be included
            P_delta = self.model.solution in ('P-Delta', 'Pushover')

            # If member not yet segmented for this combo, do it
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Select the correct segment list
            segments = self.SegmentsY if Direction == 'My' else self.SegmentsZ

            if not segments:
                continue

            # Get the maximum moment for this combo
            Mmax = max(segment.max_moment() for segment in segments)

            # Update global maximum
            if Mmax_global is None or Mmax > Mmax_global:
                Mmax_global = Mmax

        # Return 0 if no valid combos were found
        return Mmax_global if Mmax_global is not None else 0

    def min_moment(self, Direction: Literal['My', 'Mz'], combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the minimum bending moment in the member for the specified direction
        and load combination(s).

        Parameters
        ----------
        Direction : Literal['My', 'Mz']
            The direction in which to find the minimum moment:
            - 'My' : Moment about the local y-axis.
            - 'Mz' : Moment about the local z-axis.
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the minimum moment
            across those combos will be returned.

        Returns
        -------
        float
            The minimum moment value in the specified direction for the given
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.

        Notes
        -----
        - Includes P-Delta effects if the model's solution is 'P-Delta' or 'Pushover'.
        - Automatically segments the member if it has not yet been segmented
        for a load combination.
        """
        # Normalize input: single name vs list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]  # treat as a single combo name
        else:
            # Filter combos that contain *any* of the given tags
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        Mmin_global = None  # will store the global minimum across combos

        for combo_name in combo_names:
            # Skip inactive combos
            if not self.active.get(combo_name, False):
                continue

            # Determine if P-Delta (or Pushover) effects should be included
            P_delta = self.model.solution in ('P-Delta', 'Pushover')

            # If member not yet segmented for this combo, do it
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Select the correct segment list
            segments = self.SegmentsY if Direction == 'My' else self.SegmentsZ

            if not segments:
                continue

            # Get the minimum moment for this combo
            Mmin = min(segment.min_moment(P_delta) for segment in segments)

            # Update global minimum
            if Mmin_global is None or Mmin < Mmin_global:
                Mmin_global = Mmin

        # Return 0 if no valid combos were found
        return Mmin_global if Mmin_global is not None else 0


    def plot_moment(self, Direction: Literal['My', 'Mz'], combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the moment diagram for the member

        Parameters
        ----------
        Direction : string
            The direction in which to plot the moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = moment about the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        # Generate the moment diagram coordinates
        x, M = self.moment_array(Direction, n_points, combo_name)

        Member3D.__plt.plot(x, M)
        Member3D.__plt.ylabel('Moment')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.name + '\n' + combo_name)
        Member3D.__plt.show()

    def moment_array(self, Direction: Literal['My', 'Mz'], n_points: int, combo_name: str = 'Combo 1', x_array: Optional[NDArray[float64]] = None) -> NDArray[float64]:
        """
        Returns the array of the moment in the member for the given direction

        Parameters
        ----------
        Direction : string
            The direction in which to find the moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = moment about the local z-axis.
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated.
            Values must be provided in local member coordinates (between 0 and L) and be in ascending order.
        """
        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        # Determine if a P-Delta analysis has been run
        if self.model.solution == 'P-Delta' or self.model.solution == 'Pushover':
            # Include P-little-delta effects in the moment results
            P_delta = True
        else:
            # Do not include P-little delta effects in the moment results
            P_delta = False

        L = self.L()

        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        if P_delta:
            # P-delta analysis is not vectorised yet, do it element-wise
            y_arr = array([self.moment(Direction, x, combo_name) for x in x_array])
            return array([x_array, y_arr])

        else:
            # Check which axis is of interest
            if Direction == 'My':
                return self._extract_vector_results(self.SegmentsY, x_array, 'moment', P_delta)

            elif Direction == 'Mz':
                return self._extract_vector_results(self.SegmentsZ, x_array, 'moment', P_delta)

            else:
                raise ValueError(f"Direction must be 'My' or 'Mz'. {Direction} was given.")

    def torque(self, x: float, combo_name: str = 'Combo 1') -> float:
        """
        Returns the torsional moment at a point along the member's length

        Parameters
        ----------
        x : number
            The location at which to find the torque
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        # Only calculate results if the member is currently active
        if self.active[combo_name]:

            # Segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Check which segment 'x' falls on
            for segment in self.SegmentsX:
                if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                    return segment.Torsion()

                if isclose(x, self.L()):  
                    lastIndex = len(self.SegmentsX) - 1
                    return self.SegmentsX[lastIndex].Torsion()

        else:

            return 0

    def max_torque(self, combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the maximum torsional moment in the member across the
        specified load combination(s).

        Parameters
        ----------
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the maximum
            torsional moment across those combos will be returned.

        Returns
        -------
        float
            The maximum torsional moment in the member for the specified
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.
        """
        # Normalize input: single combo name or list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]
        else:
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        # Track global maximum torsion across all combos
        Tmax_global = None

        # Loop through each candidate combo
        for combo_name in combo_names:
            # Skip if inactive for this combo
            if not self.active.get(combo_name, False):
                continue

            # Re-segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Skip if no segments exist
            if not self.SegmentsX:
                continue

            # Get maximum torsion across all segments for this combo
            Tmax = max(segment.MaxTorsion() for segment in self.SegmentsX)

            # Update global maximum
            if Tmax_global is None or Tmax > Tmax_global:
                Tmax_global = Tmax

        # Return global maximum or 0 if nothing found
        return Tmax_global if Tmax_global is not None else 0


    def min_torque(self, combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the minimum torsional moment in the member across the
        specified load combination(s).

        Parameters
        ----------
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the minimum
            torsional moment across those combos will be returned.

        Returns
        -------
        float
            The minimum torsional moment in the member for the specified
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.
        """
        # Normalize input: single combo name or list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]
        else:
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        # Track global minimum torsion across all combos
        Tmin_global = None

        # Loop through each candidate combo
        for combo_name in combo_names:
            # Skip if inactive for this combo
            if not self.active.get(combo_name, False):
                continue

            # Re-segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Skip if no segments exist
            if not self.SegmentsX:
                continue

            # Get minimum torsion across all segments for this combo
            Tmin = min(segment.MinTorsion() for segment in self.SegmentsX)

            # Update global minimum
            if Tmin_global is None or Tmin < Tmin_global:
                Tmin_global = Tmin

        # Return global minimum or 0 if nothing found
        return Tmin_global if Tmin_global is not None else 0

    def plot_torque(self, combo_name='Combo 1', n_points=20):
        """
        Plots the torque diagram for the member.

        Paramters
        ---------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        x, T = self.torque_array(n_points, combo_name)

        Member3D.__plt.plot(x, T)
        Member3D.__plt.ylabel('Torsional Moment (Warping Torsion Not Included)') # Torsion results are for pure torsion. Torsional warping has not been considered
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.name + '\n' + combo_name)
        Member3D.__plt.show()

    def torque_array(self, n_points, combo_name='Combo 1', x_array = None) -> NDArray[float64]:
        """
        Returns the array of the torque in the member

        Parameters
        ----------
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated.
            Values must be provided in local member coordinates (between 0 and L) and be in ascending order.
        """
        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        L = self.L()

        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        return self._extract_vector_results(self.SegmentsZ, x_array, 'torque')

    def axial(self, x: float, combo_name: str = 'Combo 1') -> float:
        """
        Returns the axial force at a point along the member's length.

        Parameters
        ----------
        x : number
            The location at which to find the axial force.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        # Only calculate results if the member is currently active
        if self.active[combo_name]:

            # Segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Check which segment 'x' falls on
            for segment in self.SegmentsZ:
                if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                    return segment.axial(x - segment.x1)

                if isclose(x, self.L()):  
                    lastIndex = len(self.SegmentsZ) - 1
                    return self.SegmentsZ[lastIndex].axial(x - self.SegmentsZ[lastIndex].x1)

        else:

            return 0

    def max_axial(self, combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the maximum axial force in the member across the specified
        load combination(s).

        Parameters
        ----------
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the maximum axial
            force across those combos will be returned.

        Returns
        -------
        float
            The maximum axial force in the member for the specified
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.
        """
        # Normalize input: if a string is passed, treat it as a single combo name
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]
        else:
            # If tags are provided, gather all combos that match ANY of the tags
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        # Track the global maximum axial force across all combos
        Pmax_global = None

        # Loop through each candidate combo
        for combo_name in combo_names:
            # Skip if the member is inactive for this combo
            if not self.active.get(combo_name, False):
                continue

            # Re-segment the member if this combo hasn’t been solved yet
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Skip if no segments are available
            if not self.SegmentsZ:
                continue

            # Get the maximum axial force across all segments for this combo
            Pmax = max(segment.max_axial() for segment in self.SegmentsZ)

            # Update the global maximum
            if Pmax_global is None or Pmax > Pmax_global:
                Pmax_global = Pmax

        # Return the global maximum, or 0 if nothing was found
        return Pmax_global if Pmax_global is not None else 0

    def min_axial(self, combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the minimum axial force in the member across the specified
        load combination(s).

        Parameters
        ----------
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the minimum axial
            force across those combos will be returned.

        Returns
        -------
        float
            The minimum axial force in the member for the specified
            combination(s). Returns 0 if the member is inactive for all
            specified combinations or if no segments are available.
        """
        # Normalize input: if a string is passed, treat it as a single combo name
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]
        else:
            # If tags are provided, gather all combos that match ANY of the tags
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        # Track the global minimum axial force across all combos
        Pmin_global = None

        # Loop through each candidate combo
        for combo_name in combo_names:
            # Skip if the member is inactive for this combo
            if not self.active.get(combo_name, False):
                continue

            # Re-segment the member if this combo hasn’t been solved yet
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Skip if no segments are available
            if not self.SegmentsZ:
                continue

            # Get the minimum axial force across all segments for this combo
            Pmin = min(segment.min_axial() for segment in self.SegmentsZ)

            # Update the global minimum
            if Pmin_global is None or Pmin < Pmin_global:
                Pmin_global = Pmin

        # Return the global minimum, or 0 if nothing was found
        return Pmin_global if Pmin_global is not None else 0

    def plot_axial(self, combo_name: str = 'Combo 1', n_points=20) -> None:
        """
        Plots the axial force diagram for the member.
        
        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]
        
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x, P = self.axial_array(n_points, combo_name)

        Member3D.__plt.plot(x, P)
        Member3D.__plt.ylabel('Axial Force')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.name + '\n' + combo_name)
        Member3D.__plt.show()    

    def axial_array(self, n_points: int, combo_name: str = 'Combo 1', x_array: Optional[NDArray[float64]] = None) -> NDArray[float64]:
        """
        Returns the array of the axial force in the member for the given direction
        
        Parameters
        ----------
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated.
            Values must be provided in local member coordinates (between 0 and L) and be in ascending order.
        """

        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        L = self.L()
        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array<0) or any(x_array>L):
                raise ValueError(f"All x values must be in the range 0 to {L}")
            
        return self._extract_vector_results(self.SegmentsZ, x_array, 'axial')

                        
    def deflection(self, Direction: Literal['dx', 'dy', 'dz'], x: float, combo_name: str = 'Combo 1') -> float:
        """
        Returns the deflection at a point along the member's length.
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the deflection. Must be one of the following:
                'dx' = Deflection in the local x-axis.
                'dy' = Deflection in the local y-axis.
                'dz' = Deflection in the local z-axis.
        x : number
            The location at which to find the deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """
        
        # Only calculate results if the member is currently active
        if self.active[combo_name]:

            # Segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]
            
            if self.model.solution == 'P-Delta' or self.model.solution == 'Pushover':
                P_delta = True
            else:
                P_delta = False

            # Check which axis is of interest
            if Direction == 'dx':
                
                # Check which segment 'x' falls on
                for segment in self.SegmentsZ:
                    
                    if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                        return segment.axial_deflection(x - segment.x1)
                    
                if isclose(x, self.L()):
                    
                    lastIndex = len(self.SegmentsZ) - 1
                    return self.SegmentsZ[lastIndex].axial_deflection(x - self.SegmentsZ[lastIndex].x1)

            elif Direction == 'dy':
                
                # Check which segment 'x' falls on
                for segment in self.SegmentsZ:
                    
                    if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                        
                        return segment.deflection(x - segment.x1, P_delta)
                    
                if isclose(x, self.L()):
                    
                    lastIndex = len(self.SegmentsZ) - 1
                    return self.SegmentsZ[lastIndex].deflection(x - self.SegmentsZ[lastIndex].x1, P_delta)
                    
            elif Direction == 'dz':
                
                for segment in self.SegmentsY:
                    
                    if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                        
                        return segment.deflection(x - segment.x1)
                    
                if isclose(x, self.L()):
                    
                    lastIndex = len(self.SegmentsY) - 1
                    return self.SegmentsY[lastIndex].deflection(x - self.SegmentsY[lastIndex].x1) 
        
        else:

            return 0

    def max_deflection(self, Direction: Literal['dx', 'dy', 'dz'], combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the maximum deflection in the member across the specified
        load combination(s).

        Parameters
        ----------
        Direction : {'dx', 'dy', 'dz'}
            The direction in which to find the deflection:
                - 'dx' = Deflection along the local x-axis
                - 'dy' = Deflection along the local y-axis
                - 'dz' = Deflection along the local z-axis
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the maximum
            deflection across those combos will be returned.

        Returns
        -------
        float
            The maximum deflection in the member for the specified
            combination(s). Returns 0 if the member is inactive for all
            specified combinations.
        """
        # Normalize input: single name or list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]
        else:
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        # Track global maximum deflection across all combos
        dmax_global = None

        # Loop through each candidate combo
        for combo_name in combo_names:
            # Skip if member inactive for this combo
            if not self.active.get(combo_name, False):
                continue

            # Re-segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Initialize the maximum deflection for this combo (at start)
            dmax = self.deflection(Direction, 0, combo_name)

            # Sample deflections at 100 points along the member length
            for i in range(100):
                d = self.deflection(Direction, self.L() * i / 99, combo_name)
                if d > dmax:
                    dmax = d

            # Update global maximum
            if dmax_global is None or dmax > dmax_global:
                dmax_global = dmax

        # Return global maximum or 0 if nothing found
        return dmax_global if dmax_global is not None else 0


    def min_deflection(self, Direction: Literal['dx', 'dy', 'dz'], combo_tags: Union[str, List[str]] = 'Combo 1') -> float:
        """
        Returns the minimum deflection in the member across the specified
        load combination(s).

        Parameters
        ----------
        Direction : {'dx', 'dy', 'dz'}
            The direction in which to find the deflection:
                - 'dx' = Deflection along the local x-axis
                - 'dy' = Deflection along the local y-axis
                - 'dz' = Deflection along the local z-axis
        combo_tags : str or list of str, optional
            Either:
            - A single load combination name (default: 'Combo 1').
            - A list of tags. In this case, all load combinations that contain
            **any** of the given tags will be evaluated, and the minimum
            deflection across those combos will be returned.

        Returns
        -------
        float
            The minimum deflection in the member for the specified
            combination(s). Returns 0 if the member is inactive for all
            specified combinations.
        """
        # Normalize input: single name or list of tags
        if isinstance(combo_tags, str):
            combo_names = [combo_tags]
        else:
            combo_names = [
                name for name, combo in self.model.load_combos.items()
                if any(tag in combo.combo_tags for tag in combo_tags)
            ]

        # Track global minimum deflection across all combos
        dmin_global = None

        # Loop through each candidate combo
        for combo_name in combo_names:
            # Skip if member inactive for this combo
            if not self.active.get(combo_name, False):
                continue

            # Re-segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]

            # Initialize the minimum deflection for this combo (at start)
            dmin = self.deflection(Direction, 0, combo_name)

            # Sample deflections at 100 points along the member length
            for i in range(100):
                d = self.deflection(Direction, self.L() * i / 99, combo_name)
                if d < dmin:
                    dmin = d

            # Update global minimum
            if dmin_global is None or dmin < dmin_global:
                dmin_global = dmin

        # Return global minimum or 0 if nothing found
        return dmin_global if dmin_global is not None else 0
              
    def plot_deflection(self, Direction: Literal['dx', 'dy', 'dz'], combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the deflection diagram for the member
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the deflection. Must be one of the following:
                'dx' = Deflection in the local x-axis.
                'dy' = Deflection in the local y-axis.
                'dz' = Deflection in the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        n_points: int
            The number of points used to generate the plot
        """
        
        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]
                
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt
        
        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x, d = self.deflection_array(Direction, n_points, combo_name)

        Member3D.__plt.plot(x, d)
        Member3D.__plt.ylabel('Deflection')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.name + '\n' + combo_name)
        Member3D.__plt.show()

    def deflection_array(self, Direction: Literal['dx', 'dy', 'dz'], n_points: int, combo_name: str = 'Combo 1', x_array: Optional[NDArray[float64]] = None) -> NDArray[float64]:
        """
        Returns the array of the deflection in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the deflection. Must be one of the following:
                'dx' = Deflection in the local x-axis.
                'dy' = Deflection in the local y-axis.
                'dz' = Deflection in the local z-axis.
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated.
            Values must be provided in local member coordinates (between 0 and L) and be in ascending order.
        """
        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        # Determine if a P-Delta analysis has been run
        if self.model.solution == 'P-Delta' or self.model.solution == 'Pushover':
            # Include P-little-delta effects in the moment results
            P_delta = True
        else:
            # Do not include P-little delta effects in the moment results
            P_delta = False

        L = self.L()

        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array<0) or any(x_array>L):
                raise ValueError(f"All x values must be in the range 0 to {L}")
                
        if P_delta:
            # P-delta analysis is not vectorised yet, do it element-wise
            y_arr = array([self.deflection(Direction, x, combo_name) for x in x_array])
            return array([x_array, y_arr])

        else:
            # Check which axis is of interest
            if Direction == 'dz':
                return self._extract_vector_results(self.SegmentsY, x_array, 'deflection', P_delta)

            elif Direction == 'dy':
                return self._extract_vector_results(self.SegmentsZ, x_array, 'deflection', P_delta)
                                
            elif Direction == 'dx':
                return self._extract_vector_results(self.SegmentsZ, x_array, 'axial_deflection', P_delta)
            
            else:
                raise ValueError(f"Direction must be 'dx', 'dy' or 'dz'. {Direction} was given.")

    def rel_deflection(self, Direction: Literal['dx', 'dy', 'dz'], x: float, combo_name: str = 'Combo 1') -> float:
        """
        Returns the relative deflection at a point along the member's length
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the relative deflection. Must be one of the following:
                'dx' = Deflection in the local x-axis
                'dy' = Deflection in the local y-axis
                'dz' = Deflection in the local z-axis
        x : number
            The location at which to find the relative deflection
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        """

        # Only calculate results if the member is currently active
        if self.active[combo_name]:

            # Segment the member if necessary
            if self._solved_combo is None or combo_name != self._solved_combo.name:
                self._segment_member(combo_name)
                self._solved_combo = self.model.load_combos[combo_name]
            
            d = self.d(combo_name)
            dyi = d[1,0]
            dyj = d[7,0]
            dzi = d[2,0]
            dzj = d[8,0]
            L = self.L()
        
            # Check which axis is of interest
            if Direction == 'dy':
                
                # Check which segment 'x' falls on
                for segment in self.SegmentsZ:
                    
                    if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                        
                        return (segment.deflection(x - segment.x1)) - (dyi + (dyj-dyi)/L*x)
                    
                if isclose(x, self.L()):
                    
                    lastIndex = len(self.SegmentsZ) - 1
                    return (self.SegmentsZ[lastIndex].deflection(x - self.SegmentsZ[lastIndex].x1))-dyj
                    
            elif Direction == 'dz':
                
                for segment in self.SegmentsY:
                    
                    if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                        
                        return (segment.deflection(x - segment.x1)) - (dzi + (dzj-dzi)/L*x)
                    
                if isclose(x, self.L()):
                    
                    lastIndex = len(self.SegmentsY) - 1
                    return (self.SegmentsY[lastIndex].deflection(x - self.SegmentsY[lastIndex].x1)) - dzj
        
        else:

            return 0

    def plot_rel_deflection(self, Direction: Literal['dx', 'dy', 'dz'], combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the deflection diagram for the member
        
        Parameters
        ----------
        Direction : string
            The direction in which to plot the deflection. Must be one of the following:
                'dx' = Deflection in the local x-axis.
                'dy' = Deflection in the local y-axis.
                'dz' = Deflection in the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        """
        
        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]
                
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt
        
        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x, d_relative = self.rel_deflection_array(Direction, n_points, combo_name)

        Member3D.__plt.plot(x, d_relative)
        Member3D.__plt.ylabel('Relative Deflection')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.name + '\n' + combo_name)
        Member3D.__plt.show()

    def rel_deflection_array(self, Direction: Literal['dx', 'dy', 'dz'], n_points: int, combo_name: str = 'Combo 1', x_array: Optional[NDArray[float64]] = None) -> NDArray[float64]:
        """
        Returns the array of the relative deflection in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the deflection. Must be one of the following:
                'dx' = Deflection in the local x-axis.
                'dy' = Deflection in the local y-axis.
                'dz' = Deflection in the local z-axis.
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated.
            Values must be provided in local member coordinates (between 0 and L) and be in ascending order.
        """
        # Segment the member if necessary
        if self._solved_combo is None or combo_name != self._solved_combo.name:
            self._segment_member(combo_name)
            self._solved_combo = self.model.load_combos[combo_name]

        d = self.d(combo_name)
        dyi = d[1, 0]
        dyj = d[7, 0]
        dzi = d[2, 0]
        dzj = d[8, 0]

        L = self.L()
        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        # Check which axis is of interest
        if Direction == 'dy':
            deflections = self._extract_vector_results(self.SegmentsZ, x_array, 'deflection')[1]
            return vstack((x_array, deflections - (dyi + (dyj-dyi)/L*x_array)))

        elif Direction == 'dz':
            deflections = self._extract_vector_results(self.SegmentsY, x_array, 'deflection')[1]
            return vstack((x_array, deflections - (dzi + (dzj-dzi)/L*x_array)))

    def _segment_member(self, combo_name='Combo 1'):
        """
        Divides the element up into mathematically continuous segments along each axis
        """

        # Get the member's length and stiffness properties
        L = self.L()
        E = self.material.E
        A = self.section.A
        Iz = self.section.Iz
        Iy = self.section.Iy
        SegmentsZ = self.SegmentsZ
        SegmentsY = self.SegmentsY
        SegmentsX = self.SegmentsX

        # Get the load combination to segment the member for
        combo = self.model.load_combos[combo_name]

        # Create a list of discontinuity locations
        disconts = [0, L]  # Member ends

        for load in self.PtLoads:
            disconts.append(load[2])  # Point load locations

        for load in self.DistLoads:
            disconts.append(load[3])  # Distributed load start locations
            disconts.append(load[4])  # Distributed load end locations

        # Sort the list and eliminate duplicate values
        disconts = sorted(set(disconts))

        # Clear out old data from any previous analyses
        SegmentsZ.clear()
        SegmentsY.clear()
        SegmentsX.clear()

        # Create a list of mathematically continuous segments for each direction
        for index in range(len(disconts) - 1):

            # z-direction segments (bending about local z-axis)
            newSeg = BeamSegZ()            # Create the new segment
            newSeg.x1 = disconts[index]    # Segment start location
            newSeg.x2 = disconts[index+1]  # Segment end location
            newSeg.EI = E*Iz               # Segment flexural stiffness
            newSeg.EA = E*A                # Segment axial stiffness
            SegmentsZ.append(newSeg)       # Add the segment to the list

            # y-direction segments (bending about local y-axis)
            newSeg = BeamSegY()            # Create the new segment
            newSeg.x1 = disconts[index]    # Segment start location
            newSeg.x2 = disconts[index+1]  # Segment end location
            newSeg.EI = E*Iy               # Segment flexural stiffness
            newSeg.EA = E*A                # Segment axial stiffness
            SegmentsY.append(newSeg)       # Add the segment to the list

            # x-direction segments (for torsional moment)
            newSeg = BeamSegZ()            # Create the new segment
            newSeg.x1 = disconts[index]    # Segment start location
            newSeg.x2 = disconts[index + 1]  # Segment end location
            newSeg.EA = E*A                # Segment axial stiffness
            SegmentsX.append(newSeg)       # Add the segment to the list

        # Get the element local end forces, local fixed end reactions, and local displacements
        f = self.f(combo_name)           # Member local end force vector
        fer = self._fer_unc(combo_name)  # Member local fixed end reaction vector
        d = self.d(combo_name)           # Member local displacement vector

        # Get the local deflections and calculate the slope at the start of the member
        # Note 1: The slope may not be available directly from the local displacement vector if member end releases have been used, so slope-deflection has been applied to solve for it.
        # Note 2: The traditional slope-deflection equations assume a sign convention opposite of what Pynite uses for moments about the local y-axis, so a negative value has been applied to those values specifically.
        m1z = f[5, 0]        # local z-axis moment at start of member
        m2z = f[11, 0]       # local z-axis moment at end of member
        m1y = -f[4, 0]       # local y-axis moment at start of member
        m2y = -f[10, 0]      # local y-axis moment at end of member
        fem1z = fer[5, 0]    # local z-axis fixed end moment at start of member
        fem2z = fer[11, 0]   # local z-axis fixed end moment at end of member
        fem1y = -fer[4, 0]   # local y-axis fixed end moment at start of member
        fem2y = -fer[10, 0]  # local y-axis fixed end moment at end of member
        delta1y = d[1, 0]    # local y displacement at start of member
        delta2y = d[7, 0]    # local y displacement at end of member
        delta1z = d[2, 0]    # local z displacement at start of member
        delta2z = d[8, 0]    # local z displacement at end of member
        SegmentsZ[0].delta1 = delta1y
        SegmentsY[0].delta1 = delta1z
        SegmentsZ[0].theta1 = 1/3*((m1z - fem1z)*L/(E*Iz) - (m2z - fem2z)*L/(2*E*Iz) + 3*(delta2y - delta1y)/L)
        SegmentsY[0].theta1 = -1/3*((m1y - fem1y)*L/(E*Iy) - (m2y - fem2y)*L/(2*E*Iy) + 3*(delta2z - delta1z)/L)

        # Add the axial deflection at the start of the member
        SegmentsZ[0].delta_x1 = d[0, 0]
        SegmentsY[0].delta_x1 = d[0, 0]
        SegmentsX[0].delta_x1 = d[0, 0]

        # Add loads to each segment
        for i in range(len(SegmentsZ)):

            # Get the starting point of the segment
            x = SegmentsZ[i].x1

            # Initialize the distributed loads on the segment to zero
            SegmentsZ[i].w1 = 0
            SegmentsZ[i].w2 = 0
            SegmentsZ[i].p1 = 0
            SegmentsZ[i].p2 = 0
            SegmentsY[i].w1 = 0
            SegmentsY[i].w2 = 0
            SegmentsY[i].p1 = 0
            SegmentsY[i].p2 = 0

            # Initialize the slope and displacement at the start of the segment
            if i > 0:  # The first segment has already been initialized
                SegmentsZ[i].theta1 = SegmentsZ[i-1].slope(SegmentsZ[i-1].Length())
                SegmentsZ[i].delta1 = SegmentsZ[i-1].deflection(SegmentsZ[i-1].Length())
                SegmentsZ[i].delta_x1 = SegmentsZ[i-1].axial_deflection(SegmentsZ[i-1].Length())
                SegmentsY[i].theta1 = SegmentsY[i-1].slope(SegmentsY[i-1].Length())
                SegmentsY[i].delta1 = SegmentsY[i-1].deflection(SegmentsY[i-1].Length())
                SegmentsY[i].delta_x1 = SegmentsY[i-1].axial_deflection(SegmentsY[i-1].Length())

            # Add the effects of the beam end forces to the segment
            SegmentsZ[i].P1 = f[0, 0]
            SegmentsZ[i].V1 = f[1, 0]
            SegmentsZ[i].M1 = f[5, 0] - f[1, 0]*x
            SegmentsY[i].P1 = f[0, 0]
            SegmentsY[i].V1 = f[2, 0]
            SegmentsY[i].M1 = f[4, 0] + f[2, 0]*x
            SegmentsX[i].T1 = f[3, 0]

            # Step through each load case in the specified load combination
            for case, factor in combo.factors.items():

                # Add effects of point loads occuring prior to this segment
                for ptLoad in self.PtLoads:

                    if round(ptLoad[2], 10) <= round(x, 10) and case == ptLoad[3]:

                        if ptLoad[0] == 'Fx':
                            SegmentsZ[i].P1 += factor*ptLoad[1]
                        elif ptLoad[0] == 'Fy':
                            SegmentsZ[i].V1 += factor*ptLoad[1]
                            SegmentsZ[i].M1 -= factor*ptLoad[1]*(x - ptLoad[2])
                        elif ptLoad[0] == 'Fz':
                            SegmentsY[i].V1 += factor*ptLoad[1]
                            SegmentsY[i].M1 += factor*ptLoad[1]*(x - ptLoad[2])
                        elif ptLoad[0] == 'Mx':
                            SegmentsX[i].T1 += factor*ptLoad[1]    
                        elif ptLoad[0] == 'My':
                            SegmentsY[i].M1 += factor*ptLoad[1]
                        elif ptLoad[0] == 'Mz':
                            SegmentsZ[i].M1 += factor*ptLoad[1]
                        elif ptLoad[0] == 'FX' or ptLoad[0] == 'FY' or ptLoad[0] == 'FZ':
                            FX, FY, FZ = 0, 0, 0
                            if ptLoad[0] == 'FX':
                                FX = 1
                            if ptLoad[0] == 'FY':
                                FY = 1
                            if ptLoad[0] == 'FZ':
                                FZ = 1
                            force = self.T()[:3, :][:, :3] @ array([FX*ptLoad[1], FY*ptLoad[1], FZ*ptLoad[1]])
                            SegmentsZ[i].P1 += factor*force[0]
                            SegmentsZ[i].V1 += factor*force[1]
                            SegmentsZ[i].M1 -= factor*force[1]*(x - ptLoad[2])
                            SegmentsY[i].V1 += factor*force[2]
                            SegmentsY[i].M1 += factor*force[2]*(x - ptLoad[2])
                        elif ptLoad[0] == 'MX' or ptLoad[0] == 'MY' or ptLoad[0] == 'MZ':
                            MX, MY, MZ = 0, 0, 0
                            if ptLoad[0] == 'MX':
                                MX = 1
                            if ptLoad[0] == 'MY':
                                MY = 1
                            if ptLoad[0] == 'MZ':
                                MZ = 1
                            force = self.T()[:3, :][:, :3] @ array([MX*ptLoad[1], MY*ptLoad[1], MZ*ptLoad[1]])
                            SegmentsX[i].T1 += factor*force[0]
                            SegmentsY[i].M1 += factor*force[1]
                            SegmentsZ[i].M1 += factor*force[2]

                # Add distributed loads to the segment
                for distLoad in self.DistLoads:

                    if case == distLoad[5]:

                        # Get the parameters for the distributed load
                        Direction = distLoad[0]
                        w1 = factor*distLoad[1]
                        w2 = factor*distLoad[2]
                        x1 = distLoad[3]
                        x2 = distLoad[4]

                        # Determine if the load affects the segment
                        if round(x1, 10) <= round(x, 10):

                            if Direction == 'Fx':

                                # Determine if the load is on this segment
                                if round(x2, 10) > round(x, 10):

                                    # Break up the load and place it on the segment
                                    # Note that 'w1' and 'w2' are really the axial loads 'p1' and 'p2' here
                                    SegmentsZ[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                                    SegmentsZ[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1
                                    SegmentsY[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                                    SegmentsY[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1

                                    # Calculate the magnitude of the load at the start of the segment
                                    w2 = w1+(w2-w1)/(x2-x1)*(x-x1)
                                    x2 = x

                                # Calculate the axial force at the start of the segment
                                SegmentsZ[i].P1 += (w1 + w2)/2*(x2 - x1)
                                SegmentsY[i].P1 += (w1 + w2)/2*(x2 - x1)

                            elif Direction == 'Fy':

                                # Determine if the load is on this segment
                                if round(x2, 10) > round(x, 10):

                                    # Break up the load and place it on the segment
                                    SegmentsZ[i].w1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                                    SegmentsZ[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1

                                    # Calculate the magnitude of the load at the start of the segment
                                    # This will be used as the 'x2' value for the load prior to the start of the segment
                                    w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                                    x2 = x

                                # Calculate the shear and moment at the start of the segment due to the load
                                SegmentsZ[i].V1 += (w1 + w2)/2*(x2 - x1)
                                SegmentsZ[i].M1 -= (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6

                            elif Direction == 'Fz':

                                # Determine if the load is on this segment
                                if round(x2, 10) > round(x, 10):

                                    # Break up the load and place it on the segment
                                    SegmentsY[i].w1 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x1 - x1) + w1
                                    SegmentsY[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1

                                    # Calculate the magnitude of the load at the start of the segment
                                    w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                                    x2 = x

                                # Calculate the shear and moment at the start of the segment due to the load
                                SegmentsY[i].V1 += (w1 + w2)/2*(x2 - x1)
                                SegmentsY[i].M1 += (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6

                            elif Direction == 'FX' or Direction == 'FY' or Direction == 'FZ':

                                FX, FY, FZ = 0, 0, 0
                                if Direction == 'FX': FX = 1
                                if Direction == 'FY': FY = 1
                                if Direction == 'FZ': FZ = 1
                                T = self.T()[:3, :][:, :3]
                                f1 = T @ array([FX*w1, FY*w1, FZ*w1])
                                f2 = T @ array([FX*w2, FY*w2, FZ*w2])

                                # Determine if the load is on this segment
                                if round(x2, 10) > round(x, 10):

                                    # Break up the load and place it on the segment
                                    SegmentsZ[i].p1 += (f2[0] - f1[0])/(x2 - x1)*(x - x1) + f1[0]
                                    SegmentsZ[i].p2 += (f2[0] - f1[0])/(x2 - x1)*(SegmentsZ[i].x2 - x1) + f1[0]
                                    SegmentsY[i].p1 += (f2[0] - f1[0])/(x2 - x1)*(x - x1) + f1[0]
                                    SegmentsY[i].p2 += (f2[0] - f1[0])/(x2 - x1)*(SegmentsY[i].x2 - x1) + f1[0]

                                    SegmentsZ[i].w1 += (f2[1] - f1[1])/(x2 - x1)*(x - x1) + f1[1]
                                    SegmentsZ[i].w2 += (f2[1] - f1[1])/(x2 - x1)*(SegmentsZ[i].x2 - x1) + f1[1]

                                    SegmentsY[i].w1 += (f2[2] - f1[2])/(x2 - x1)*(SegmentsY[i].x1 - x1) + f1[2]
                                    SegmentsY[i].w2 += (f2[2] - f1[2])/(x2 - x1)*(SegmentsY[i].x2 - x1) + f1[2]

                                    # Calculate the magnitude of the load at the start of the segment
                                    w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                                    f2 = T @ array([FX*w2, FY*w2, FZ*w2])
                                    x2 = x

                                # Calculate the axial force, shear and moment at the start of the segment
                                SegmentsZ[i].P1 += (f1[0] + f2[0])/2*(x2 - x1)
                                SegmentsY[i].P1 += (f1[0] + f2[0])/2*(x2 - x1)

                                SegmentsZ[i].V1 += (f1[1] + f2[1])/2*(x2 - x1)
                                SegmentsZ[i].M1 -= (x1 - x2)*(2*f1[1]*x1 - 3*f1[1]*x + f1[1]*x2 + f2[1]*x1 - 3*f2[1]*x + 2*f2[1]*x2)/6

                                SegmentsY[i].V1 += (f1[2] + f2[2])/2*(x2 - x1)
                                SegmentsY[i].M1 += (x1 - x2)*(2*f1[2]*x1 - 3*f1[2]*x + f1[2]*x2 + f2[2]*x1 - 3*f2[2]*x + 2*f2[2]*x2)/6

    def _extract_vector_results(self, segments: List, x_array: NDArray[float64], result_name: Literal['moment', 'shear', 'axial', 'torque', 'deflection', 'axial_deflection'], P_delta: bool = False) -> NDArray[float64]:
        """
        Extracts result values at specified locations along a structural member using efficient, 
        vectorized evaluation of piecewise segment functions.

        Assumes:
            - `segments` are ordered from left to right (increasing x).
            - `x_array` is sorted in ascending order.
            - Each `segment` provides a method for evaluating the requested result type.

        Parameters
        ----------
        segments : List
            List of segment objects. Each segment represents a portion of a structural member
            and must have properties `x1`, `x2`, and appropriate result methods (e.g. `moment()`, `Shear()`, etc.).

        x_array : NDArray[float64]
            1D NumPy array of x-coordinates (global positions along the member) at which to evaluate results.

        result_name : Literal
            The type of result to extract. Must be one of:
            'moment', 'shear', 'axial', 'torque', 'deflection', or 'axial_deflection'.

        P_delta : bool, optional
            Whether to include second-order (P-Delta) effects for applicable result types 
            ('moment' and 'deflection'). Default is False.

        Returns
        -------
        NDArray[float64]
            A 2D NumPy array of shape (2, N), where:
                - Row 0 contains the x-locations (same order as x_array, filtered by segment coverage).
                - Row 1 contains the corresponding result values (same order).

        Raises
        ------
        ValueError
            If `result_name` is not a supported option.
        """

        # Prepare storage for results
        segment_results = []  # Stores y-values per segment
        x_results = []        # Stores x-values that belong to each segment

        # Dispatch table maps the result name to the correct evaluation function
        method_map = {
            "moment": lambda s, x: s.moment(x, P_delta),
            "shear": lambda s, x: s.Shear(x),
            "axial": lambda s, x: s.axial(x),
            "torque": lambda s, x: s.Torsion(x),
            "deflection": lambda s, x: s.deflection(x, P_delta),
            "axial_deflection": lambda s, x: s.axial_deflection(x),
        }

        # Lookup the result computation method
        compute_result = method_map.get(result_name)
        if compute_result is None:
            raise ValueError(f"Unsupported result_name: {result_name}")

        idx = 0  # Tracks current position in x_array
        n = x_array.size

        for i, segment in enumerate(segments):

            # For the last segment, include points up to and including x2
            if i == len(segments) - 1:
                filter = (x_array[idx:] >= segment.x1) & (x_array[idx:] <= segment.x2)
            else:
                # For intermediate segments, include points up to but not including x2
                filter = (x_array[idx:] >= segment.x1) & (x_array[idx:] < segment.x2)

            # Count how many x-values fall within this segment
            count = count_nonzero(filter)
            if count == 0:
                continue  # No points to evaluate in this segment

            # Extract the relevant x-values for this segment
            segment_x = x_array[idx:][filter]
            local_x = segment_x - segment.x1  # Convert global x to local segment coordinates

            # Evaluate the selected result at each local x
            segment_y = compute_result(segment, local_x)

            # Store for final output
            x_results.append(segment_x)
            segment_results.append(segment_y)

            # Advance the index so we don't re-check already-processed x values
            idx += count
            if idx >= n:
                break  # All x-values have been processed

        # Assemble full x and y arrays
        all_x = concatenate(x_results)
        all_y = concatenate(segment_results)

        return vstack((all_x, all_y))

# %%
