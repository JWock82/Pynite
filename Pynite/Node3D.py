# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:04:56 2017

@author: D. Craig Brinck, SE
"""
from __future__ import annotations  # Allows more recent type hints features
# %%
from numpy import array, zeros

from typing import List, Tuple, Dict, Optional,TYPE_CHECKING
if TYPE_CHECKING:

    from typing import Dict, List, Tuple, Optional, Any, Literal

    from numpy import float64
    from numpy.typing import NDArray

    from Pynite.FEModel3D import FEModel3D
    from Pynite.LoadCombo import LoadCombo


class Node3D():
    """
    A class representing a node in a 3D finite element model.
    """

    def __init__(self, model: FEModel3D, name: str, X: float, Y: float, Z: float):

        self.name = name                 # A unique name for the node assigned by the user
        self.ID: Optional[int] = None    # A unique index number for the node assigned by the program

        self.X = X          # Global X coordinate
        self.Y = Y          # Global Y coordinate
        self.Z = Z          # Global Z coordinate

        self.NodeLoads: List[Tuple[str, float, str]] = []  # A list of loads applied to the node (Direction, P, case) or (Direction, M, case)

        # Initialize the dictionaries of calculated node displacements
        self.DX: Dict[str, float] = {}
        self.DY: Dict[str, float] = {}
        self.DZ: Dict[str, float] = {}
        self.RX: Dict[str, float] = {}
        self.RY: Dict[str, float] = {}
        self.RZ: Dict[str, float] = {}

        # Initialize the dictionaries of calculated node reactions
        self.RxnFX: Dict[str, float] = {}
        self.RxnFY: Dict[str, float] = {}
        self.RxnFZ: Dict[str, float] = {}
        self.RxnMX: Dict[str, float] = {}
        self.RxnMY: Dict[str, float] = {}
        self.RxnMZ: Dict[str, float] = {}

        # Initialize all support conditions to `False`
        self.support_DX: bool = False
        self.support_DY: bool = False
        self.support_DZ: bool = False
        self.support_RX: bool = False
        self.support_RY: bool = False
        self.support_RZ: bool = False

        # Inititialize all support springs
        self.spring_DX: List[float | str | bool | None] = [None, None, None]  # [stiffness, direction, active]
        self.spring_DY: List[float | str | bool | None] = [None, None, None]
        self.spring_DZ: List[float | str | bool | None] = [None, None, None]
        self.spring_RX: List[float | str | bool | None] = [None, None, None]
        self.spring_RY: List[float | str | bool | None] = [None, None, None]
        self.spring_RZ: List[float | str | bool | None] = [None, None, None]

        # Initialize all enforced displacements to `None`
        self.EnforcedDX: float | None = None
        self.EnforcedDY: float | None = None
        self.EnforcedDZ: float | None = None
        self.EnforcedRX: float | None = None
        self.EnforcedRY: float | None = None
        self.EnforcedRZ: float | None = None

        # Initialize the color contour value for the node. This will be used for contour smoothing.
        self.contour: List[float] = []

        # The 'Node3D' object will store results for one load combination at a time. 

        # Adding a link to the model that Nodes belong to 
        self.model: FEModel3D = model

    def distance(self, other: 'Node3D') -> float:
        """
        Returns the distance to another node.

        Parameters
        ----------
        other : Node3D
            A node object to compare coordinates with.
        """
        return ((self.X - other.X)**2 + (self.Y - other.Y)**2 + (self.Z - other.Z)**2)**0.5

    def M(self, mass_combo_name: str | None = None, mass_direction: str = 'Y', gravity: float = 1.0, characteristic_length: float | None = None) -> NDArray:
        """Returns the node's mass matrix (6x6 diagonal). In member-based models, nodes provide translational mass only to prevent double-counting of rotational inertia. For member-less models, rotational inertia can be added by providing a characteristic_length.

        :param mass_combo_name: Load combination name for force-based mass calculation, defaults to ''.
        :type mass_combo_name: str, optional
        :param mass_direction: Direction for load-to-mass conversion ('X', 'Y', or 'Z'). Any loads applied in this direction (positive or negative) will be converted to mass. Default is 'Y'.
        :type mass_direction: str, optional
        :param gravity: The acceleration due to gravity. Defaults to 1.0. In most cases you'll want to change this to be in units consistent with your model.
        :type gravity: float
        :param characteristic_length: If provided, adds rotational inertia scaled by I = m * (0.01L²).
                                      Use None for translational-only mass, defaults to None
        :type characteristic_length: float, optional
        :return: Node mass matrix of shape (6, 6)
        :rtype: numpy.ndarray
        """

        # Initialize the nodal mass matrix
        m = zeros((6, 6))

        # Get mass from nodal loads in the specified load combination
        if mass_combo_name not in self.model.load_combos:
            raise ValueError(f'Load combination {mass_combo_name} not found in the model.')
        else:
            total_mass = self._calc_mass(mass_combo_name, mass_direction, gravity)

        # Check if there is any mass at this node
        if total_mass > 0:

            # Create lumped mass matrix for the node
            # Mass is distributed to translational DOFs only (standard practice)
            m[0, 0] = total_mass  # FX
            m[1, 1] = total_mass  # FY
            m[2, 2] = total_mass  # FZ

            # Smart rotational inertia - only for free rotational DOFs
            # Use a small value based on mass and a characteristic length
            # Scale by characteristic length squared (I = m * r²)
            if characteristic_length is not None:
                rotational_inertia = total_mass * (characteristic_length ** 2) * 0.01  # 1% scaling

                if not self.support_RX:
                    m[3, 3] = rotational_inertia
                if not self.support_RY:
                    m[4, 4] = rotational_inertia
                if not self.support_RZ:
                    m[5, 5] = rotational_inertia

        return m

    def _calc_mass(self, mass_combo_name: str, mass_direction: str = 'Y', gravity: float = 1.0) -> float:
        """Calculates the total mass from nodal loads in a load combination.

        For modal analysis, nodal "masses" are typically specified as force loads that get converted to mass using F = m*a (where a is typically gravity).

        :param mass_combo_name: The name of the load combination to use for mass calculation.
        :type mass_combo_name: str
        :param mass_direction: Direction for load-to-mass conversion ('X', 'Y', or 'Z'). Any loads applied in this direction (positive or negative) will be converted to mass. Default is 'Y'.
        :type mass_direction: str, optional
        :param gravity: The acceleration due to gravity. Defaults to 1.0. In most cases you'll want to change this to be in units consistent with your model.
        :type gravity: float
        :raises NameError: _description_
        :return: The total mass calculated from forces in the specified direction
        :rtype: float
        """

        # Initialize the force summation to zero
        total_force = 0.0

        # Get the mass combo from the model
        try:
            mass_combo = self.model.load_combos[mass_combo_name]
        except KeyError:
            raise NameError(f"No load combination named '{mass_combo_name}'")

        # Step through each nodal load acting on this node
        for load in self.NodeLoads:

            # Read in the nodal load parameters
            load_direction, load_value, load_case = load

            # Apply load factors to load cases in the mass combo
            if load_case in mass_combo.factors:

                factor = mass_combo.factors[load_case]
                load_magnitude = factor * load_value

                # Sum forces in the specified direction
                if mass_direction == 'X' and load_direction == 'FX':
                    total_force += load_magnitude
                elif mass_direction == 'Y' and load_direction == 'FY':
                    total_force += load_magnitude
                elif mass_direction == 'Z' and load_direction == 'FZ':
                    total_force += load_magnitude

        # Convert force to mass
        total_mass = total_force / gravity if gravity != 0.0 else 0.0

        # Debug output (you can remove this later)
        # if abs(total_mass) > 1e-10:
        #     print(f'{self.name}: total_force: {total_force:.6f}, total_mass: {total_mass:.6f}, combo_factors: {combo.factors}')

        return abs(total_mass)   # Mass is always positive
