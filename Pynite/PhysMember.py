from __future__ import annotations # Allows more recent type hints features
from typing import Dict, List, Literal, Tuple, TYPE_CHECKING
from Pynite.Member3D import Member3D

if TYPE_CHECKING:

    from Pynite.Node3D import Node3D
    from Pynite.FEModel3D import FEModel3D
    import numpy.typing as npt
    from numpy import float64
    from numpy.typing import NDArray

from numpy import array, dot, linspace, hstack, empty
from numpy.linalg import norm
from math import isclose, acos

class PhysMember(Member3D):
    """
    A physical member.

    Physical members can detect internal nodes and subdivide themselves into sub-members at those
    nodes.
    """

    # '__plt' is used to store the 'pyplot' from matplotlib once it gets imported. Setting it to 'None' for now allows us to defer importing it until it's actually needed.
    __plt = None

    def __init__(self, model: FEModel3D, name: str, i_node: Node3D, j_node: Node3D, material_name: str, section_name: str, rotation: float = 0.0,
                 tension_only: bool = False, comp_only: bool = False) -> None:

        super().__init__(model, name, i_node, j_node, material_name, section_name, rotation, tension_only, comp_only)
        self.sub_members: Dict[str, Member3D] = {}

    def descritize(self) -> None:
        """
        Subdivides the physical member into sub-members at each node along the physical member
        """

        # Clear out any old sub_members
        self.sub_members = {}

        # Start a new list of nodes along the member
        int_nodes: List[Tuple[Node3D, float]] = []

        # Create a vector from the i-node to the j-node
        Xi, Yi, Zi = self.i_node.X, self.i_node.Y, self.i_node.Z
        Xj, Yj, Zj = self.j_node.X, self.j_node.Y, self.j_node.Z
        v_ij = array([Xj - Xi, Yj - Yi, Zj - Zi])
        L = norm(v_ij)

        # Guard against zero-length members
        if L == 0:
            return

        u = v_ij / L  # Unit vector along member axis

        # Add the i-node and j-node to the list
        int_nodes.append((self.i_node, 0.0))
        int_nodes.append((self.j_node, L))

        # Absolute/per-length tolerance for colinearity (consistent with previous strictness)
        # Using a small multiple of length maintains scale-invariance like the prior approach.
        tol = 1e-12 * (1.0 + L)

        # Fast axis-aligned bounding box (AABB) prefilter to cull distant nodes quickly
        # Use a slightly larger tolerance than `tol` to avoid false negatives near corners
        bb_tol = 1e-9 * (1.0 + L)
        xmin = min(Xi, Xj) - bb_tol
        xmax = max(Xi, Xj) + bb_tol
        ymin = min(Yi, Yj) - bb_tol
        ymax = max(Yi, Yj) + bb_tol
        zmin = min(Zi, Zj) - bb_tol
        zmax = max(Zi, Zj) + bb_tol

        # Step through each node in the model
        for node in self.model.nodes.values():

            # Skip the end nodes
            if node is self.i_node or node is self.j_node:
                continue

            # Bounding-box reject (very fast)
            if (
                node.X < xmin or node.X > xmax or
                node.Y < ymin or node.Y > ymax or
                node.Z < zmin or node.Z > zmax
            ):
                continue

            # Vector from i-node to this node
            dx = node.X - Xi
            dy = node.Y - Yi
            dz = node.Z - Zi

            # Parametric location along the member (projection onto axis)
            t = dx * u[0] + dy * u[1] + dz * u[2]

            # Quickly reject nodes outside the i-j segment
            if t <= 0.0 or t >= L:
                continue

            # Perpendicular distance from the member line
            px = dx - t * u[0]
            py = dy - t * u[1]
            pz = dz - t * u[2]
            d_perp = (px**2 + py**2 + pz**2) ** 0.5

            # Consider node on member if perpendicular distance is negligible
            if d_perp <= tol:
                int_nodes.append((node, t))

        # Create a list of sorted intermediate nodes by distance from the i-node
        int_nodes = sorted(int_nodes, key=lambda x: x[1])

        # Break up the member into sub-members at each intermediate node
        for i in range(len(int_nodes) - 1):

            # Generate the sub-member's name (physical member name + a, b, c, etc.)
            name = self.name + chr(i+97)

            # Find the i and j nodes for the sub-member, and their positions along the physical
            # member's local x-axis
            i_node = int_nodes[i][0]
            j_node = int_nodes[i+1][0]
            xi = int_nodes[i][1]
            xj = int_nodes[i+1][1]

            # Create a new sub-member
            new_sub_member = Member3D(self.model, name, i_node, j_node, self.material.name, 
                                      self.section.name, self.rotation, self.tension_only, 
                                      self.comp_only)

            # Flag the sub-member as active
            for combo_name in self.model.load_combos.keys():
                new_sub_member.active[combo_name] = True

            # Apply end releases if applicable
            if i == 0:
                new_sub_member.Releases[0:6] = self.Releases[0:6]
            if i == len(int_nodes) - 2:
                new_sub_member.Releases[6:12] = self.Releases[6:12]

            # Add distributed to the sub-member
            for dist_load in self.DistLoads:

                # Find the start and end points of the distributed load in the physical member's
                # local coordinate system
                x1_load = dist_load[3]
                x2_load = dist_load[4]

                # Determine if the distributed load should be applied to this segment
                if x1_load <= xj and x2_load > xi: 

                    direction = dist_load[0]
                    w1 = dist_load[1]
                    w2 = dist_load[2]
                    case = dist_load[5]
                    self_weight = dist_load[6]

                    # Equation describing the load as a function of x
                    # Linear interpolation of distributed load
                    w = lambda x: (w2 - w1)/(x2_load - x1_load)*(x - x1_load) + w1

                    # Chop up the distributed load for the sub-member
                    if x1_load > xi:
                        x1 = x1_load - xi
                    else:
                        x1 = 0
                        w1 = w(xi)

                    if x2_load < xj:
                        x2 = x2_load - xi
                    else:
                        x2 = xj - xi
                        w2 = w(xj)

                    # Add the load to the sub-member
                    new_sub_member.DistLoads.append([direction, w1, w2, x1, x2, case, self_weight])

            # Add point loads to the sub-member
            for pt_load in self.PtLoads:

                direction = pt_load[0]
                P = pt_load[1]
                x = pt_load[2]
                case = pt_load[3]

                # Determine if the point load should be applied to this segment
                if x >= xi and x < xj or (isclose(x, xj) and isclose(xj, self.L())):

                    x = x - xi

                    # Add the load to the sub-member
                    new_sub_member.PtLoads.append([direction, P, x, case])

            # Add the new sub-member to the sub-member dictionary for this physical member
            self.sub_members[name] = new_sub_member

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
            The name of the load combination to get the results for (not the combination itself).
        """

        member, x_mod = self.find_member(x)
        return member.shear(Direction, x_mod, combo_name)

    def max_shear(self, Direction: Literal['Fy', 'Fz'], combo_name: str = 'Combo 1') -> float:
        """
        Returns the maximum shear in the member for the given direction

        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis
                'Fz' = Shear acting on the local z-axis
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        """

        Vmax = None
        for member in self.sub_members.values():
            V = member.max_shear(Direction, combo_name)
            if Vmax is None or V > Vmax:
                Vmax = V
        return Vmax

    def min_shear(self, Direction: Literal['Fy', 'Fz'], combo_name: str = 'Combo 1') -> float:
        """
        Returns the minimum shear in the member for the given direction

        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis
                'Fz' = Shear acting on the local z-axis
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        Vmin = None
        for member in self.sub_members.values():
            V = member.min_shear(Direction, combo_name)
            if Vmin is None or V < Vmin:
                Vmin = V
        return Vmin

    def plot_shear(self, Direction: Literal['Fy', 'Fz'], combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the shear diagram for the member

        Parameters
        ----------
        Direction : string
            The direction in which to plot the shear force. Must be one of the following:
                'Fy' = Shear in the local y-axis.
                'Fz' = Shear in the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Import 'pyplot' if not already done
        if PhysMember.__plt is None:
            from matplotlib import pyplot as plt
            PhysMember.__plt = plt

        fig, ax = PhysMember.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        # Generate the shear diagram
        V_array = self.shear_array(Direction, n_points, combo_name)
        x = V_array[0]
        V = V_array[1]

        PhysMember.__plt.plot(x, V)
        PhysMember.__plt.ylabel('Shear')
        PhysMember.__plt.xlabel('Location')
        PhysMember.__plt.title('Member ' + self.name + '\n' + combo_name)
        PhysMember.__plt.show()

    def shear_array(self, Direction: Literal['Fy', 'Fz'], n_points: int, combo_name='Combo 1', x_array=None) -> NDArray[float64]:
        """
        Returns the array of the shear in the physical member for the given direction

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
            A custom array of x values that may be provided by the user, otherwise an array is generated. Values must be provided in local member coordinates (between 0 and L) and be in ascending order
        """

        # `v_array2` will be used to store the shear values for the overall member
        v_array2 = empty((2, 1))

        # Create an array of locations along the physical member to obtain results at
        L = self.L()
        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        # Step through each submember in the physical member
        x_o = 0
        for i, submember in enumerate(self.sub_members.values()):

            # Segment the submember into segments with mathematically continuous loads if not already done
            if submember._solved_combo is None or combo_name != submember._solved_combo.name:
                submember._segment_member(combo_name)
                submember._solved_combo = self.model.load_combos[combo_name]

            # Check if this is the last submember
            if i == len(self.sub_members.values()) - 1:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array <= x_o + submember.L())

            # Not the last submember
            else:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array < x_o + submember.L())

            x_subm_array = x_array[filter] - x_o

            # Check which axis is of interest
            if Direction == 'Fz':
                v_array = self._extract_vector_results(submember.SegmentsY, x_subm_array, 'shear')
            elif Direction == 'Fy':
                v_array = self._extract_vector_results(submember.SegmentsZ, x_subm_array, 'shear')
            else:
                raise ValueError(f"Direction must be 'Fy' or 'Fz'. {Direction} was given.")

            # Adjust from the submember's coordinate system to the physical member's coordinate system
            v_array[0] = [x_o + x for x in v_array[0]]

            # Add the submember shear values to the overall member shear values in `v_array2`
            if i != 0:
                v_array2 = hstack((v_array2, v_array))
            else:
                v_array2 = v_array

            # Get the starting position of the next submember
            x_o += submember.L()

        # Return the results
        return v_array2

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

        member, x_mod = self.find_member(x)
        return member.moment(Direction, x_mod, combo_name)

    def max_moment(self, Direction: Literal['My', 'Mz'], combo_name: str = 'Combo 1') -> float:
        """
        Returns the maximum moment in the member for the given direction.
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = Moment about the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        """

        Mmax = None
        for member in self.sub_members.values():
            M = member.max_moment(Direction, combo_name)
            if Mmax is None or M > Mmax:
                Mmax = M
        return Mmax

    def min_moment(self, Direction: Literal['My', 'Mz'], combo_name: str = 'Combo 1') -> float:
        """
        Returns the minimum moment in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = Moment about the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """
        
        Mmin = None
        for member in self.sub_members.values():
            M = member.min_moment(Direction, combo_name)
            if Mmin is None or M < Mmin:
                Mmin = M
        return Mmin

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

        # Import 'pyplot' if not already done
        if PhysMember.__plt is None:
            from matplotlib import pyplot as plt
            PhysMember.__plt = plt

        fig, ax = PhysMember.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        # Generate the moment diagram
        M_array = self.moment_array(Direction, n_points, combo_name)
        x = M_array[0]
        M = M_array[1]

        PhysMember.__plt.plot(x, M)
        PhysMember.__plt.ylabel('Moment')
        PhysMember.__plt.xlabel('Location')
        PhysMember.__plt.title('Member ' + self.name + '\n' + combo_name)
        PhysMember.__plt.show()

    def moment_array(self, Direction: Literal['My', 'Mz'], n_points: int, combo_name='Combo 1', x_array=None) -> NDArray[float64]:
        """
        Returns the array of the moment in the physical member for the given direction

        Parameters
        ----------
        Direction : string
            The direction to plot the moment for. Must be one of the following:
                'My' = Moment acting about the local y-axis (usually the weak-axis).
                'Mz' = Moment acting about the local z-axis (usually the strong-axis).
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated.
            Values must be provided in local member coordinates (between 0 and L) and be in ascending order
        """

        # `m_array2` will be used to store the moment values for the overall member
        m_array2 = empty((2, 1))

        # Create an array of locations along the physical member to obtain results at
        L = self.L()
        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        # Step through each submember in the physical member
        x_o = 0
        for i, submember in enumerate(self.sub_members.values()):

            # Segment the submember into segments with mathematically continuous loads if not already done
            if submember._solved_combo is None or combo_name != submember._solved_combo.name:
                submember._segment_member(combo_name)
                submember._solved_combo = self.model.load_combos[combo_name]

            # Check if this is the last submember
            if i == len(self.sub_members.values()) - 1:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array <= x_o + submember.L())

            # Not the last submember
            else:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array < x_o + submember.L())

            x_subm_array = x_array[filter] - x_o

            # Check if P-Delta analysis was run
            if self.model.solution == 'P-Delta':
                PDelta = True
            else:
                PDelta = False

            # Check which axis is of interest
            if Direction == 'My':
                m_array = self._extract_vector_results(submember.SegmentsY, x_subm_array, 'moment', PDelta)
            elif Direction == 'Mz':
                m_array = self._extract_vector_results(submember.SegmentsZ, x_subm_array, 'moment', PDelta)
            else:
                raise ValueError(f"Direction must be 'My' or 'Mz'. {Direction} was given.")

            # Adjust from the submember's coordinate system to the physical member's coordinate system
            m_array[0] = [x_o + x for x in m_array[0]]

            # Add the submember moment values to the overall member shear values in `m_array2`
            if i != 0:
                m_array2 = hstack((m_array2, m_array))
            else:
                m_array2 = m_array

            # Get the starting position of the next submember
            x_o += submember.L()

        # Return the results
        return m_array2

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
        
        member, x_mod = self.find_member(x)
        return member.torque(x_mod, combo_name)

    def max_torque(self, combo_name: str = 'Combo 1') -> float:
        
        Tmax = None
        for member in self.sub_members.values():
            T = member.max_torque(combo_name)
            if Tmax is None or T > Tmax:
                Tmax = T
        return Tmax

    def min_torque(self, combo_name: str = 'Combo 1') -> float:
        """
        Returns the minimum torsional moment in the member.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        Tmin = None
        for member in self.sub_members.values():
            T = member.min_torque(combo_name)
            if Tmin is None or T < Tmin:
                Tmin = T
        return Tmin

    def plot_torque(self, combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the torque diagram for the member

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Import 'pyplot' if not already done
        if PhysMember.__plt is None:
            from matplotlib import pyplot as plt
            PhysMember.__plt = plt

        fig, ax = PhysMember.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        # Generate the torque diagram
        T_array = self.torque_array(n_points, combo_name)
        x = T_array[0]
        T = T_array[1]

        PhysMember.__plt.plot(x, T)
        PhysMember.__plt.ylabel('Torque')
        PhysMember.__plt.xlabel('Location')
        PhysMember.__plt.title('Member ' + self.name + '\n' + combo_name)
        PhysMember.__plt.show()

    def torque_array(self, n_points: int, combo_name='Combo 1', x_array=None) -> NDArray[float64]:
        """
        Returns the array of the torque in the physical member.

        Parameters
        ----------
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated. Values must be provided in local member coordinates (between 0 and L) and be in ascending order
        """

        # `t_array2` will be used to store the torque values for the overall member
        t_array2 = empty((2, 1))

        # Create an array of locations along the physical member to obtain results at
        L = self.L()
        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        # Step through each submember in the physical member
        x_o = 0
        for i, submember in enumerate(self.sub_members.values()):

            # Segment the submember into segments with mathematically continuous loads if not already done
            if submember._solved_combo is None or combo_name != submember._solved_combo.name:
                submember._segment_member(combo_name)
                submember._solved_combo = self.model.load_combos[combo_name]

            # Check if this is the last submember
            if i == len(self.sub_members.values()) - 1:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array <= x_o + submember.L())

            # Not the last submember
            else:

                # Find any points from `x_array` that lie along this submember
                # x_subm_array = [x - x_o for x in x_array if x >= x_o and x < x_o + submember.L()]
                filter = (x_array >= x_o) & (x_array < x_o + submember.L())

            x_subm_array = x_array[filter] - x_o

            # Check which axis is of interest
            t_array = self._extract_vector_results(submember.SegmentsX, x_subm_array, 'torque')

            # Adjust from the submember's coordinate system to the physical member's coordinate system
            t_array[0] = [x_o + x for x in t_array[0]]

            # Add the submember torque values to the overall member shear values in `t_array2`
            if i != 0:
                t_array2 = hstack((t_array2, t_array))
            else:
                t_array2 = t_array

            # Get the starting position of the next submember
            x_o += submember.L()

        # Return the results
        return t_array2

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

        member, x_mod = self.find_member(x)
        return member.axial(x_mod, combo_name)

    def max_axial(self, combo_name: str = 'Combo 1') -> float:
        
        Pmax = None
        for member in self.sub_members.values():
            P = member.max_axial(combo_name)
            if Pmax is None or P > Pmax:
                Pmax = P
        return Pmax

    def min_axial(self, combo_name: str = 'Combo 1') -> float:

        Pmin = None
        for member in self.sub_members.values():
            P = member.min_axial(combo_name)
            if Pmin is None or P < Pmin:
                Pmin = P
        return Pmin

    def plot_axial(self, combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the axial force diagram for the member

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Import 'pyplot' if not already done
        if PhysMember.__plt is None:
            from matplotlib import pyplot as plt
            PhysMember.__plt = plt

        fig, ax = PhysMember.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        # Generate the axial force array
        P_array = self.axial_array(n_points, combo_name)
        x = P_array[0]
        P = P_array[1]

        PhysMember.__plt.plot(x, P)
        PhysMember.__plt.ylabel('Axial Force')
        PhysMember.__plt.xlabel('Location')
        PhysMember.__plt.title('Member ' + self.name + '\n' + combo_name)
        PhysMember.__plt.show()

    def axial_array(self, n_points: int, combo_name='Combo 1', x_array=None) -> NDArray[float64]:
        """
        Returns the array of the axial force in the physical member.

        Parameters
        ----------
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated. Values must be provided in local member coordinates (between 0 and L) and be in ascending order
        """

        # `a_array2` will be used to store the axial force values for the overall member
        a_array2 = empty((2, 1))

        # Create an array of locations along the physical member to obtain results at
        L = self.L()
        if x_array is None:
            x_array = linspace(0, L, n_points)
        else:
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        # Step through each submember in the physical member
        x_o = 0
        for i, submember in enumerate(self.sub_members.values()):

            # Segment the submember into segments with mathematically continuous loads if not already done
            if submember._solved_combo is None or combo_name != submember._solved_combo.name:
                submember._segment_member(combo_name)
                submember._solved_combo = self.model.load_combos[combo_name]

            # Check if this is the last submember
            if i == len(self.sub_members.values()) - 1:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array <= x_o + submember.L())

            # Not the last submember
            else:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array < x_o + submember.L())

            x_subm_array = x_array[filter] - x_o

            # Check which axis is of interest
            a_array = self._extract_vector_results(submember.SegmentsZ, x_subm_array, 'axial')

            # Adjust from the submember's coordinate system to the physical member's coordinate system
            a_array[0] = [x_o + x for x in a_array[0]]

            # Add the submember axial values to the overall member shear values in `a_array2`
            if i != 0:
                a_array2 = hstack((a_array2, a_array))
            else:
                a_array2 = a_array

            # Get the starting position of the next submember
            x_o += submember.L()
        
        # Return the results
        return a_array2

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

        member, x_mod = self.find_member(x)
        return member.deflection(Direction, x_mod, combo_name)

    def max_deflection(self, Direction: Literal['dx', 'dy', 'dz'], combo_name: str = 'Combo 1') -> float:
        """
        Returns the maximum deflection in the member.
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to find the maximum deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        dmax = None
        for member in self.sub_members.values():
            d = member.max_deflection(Direction, combo_name)
            if dmax is None or d > dmax:
                dmax = d
        return dmax

    def min_deflection(self, Direction: Literal['dx', 'dy', 'dz'], combo_name: str = 'Combo 1') -> float:
        """
        Returns the minimum deflection in the member.
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to find the minimum deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        dmin = None
        for member in self.sub_members.values():
            d = member.min_deflection(Direction, combo_name)
            if dmin is None or d < dmin:
                dmin = d
        return dmin

    def rel_deflection(self, Direction: Literal['dx', 'dy', 'dz'], x: float, combo_name: str = 'Combo 1') -> float:
        """
        Returns the relative deflection at a point along the member's length
        
        Parameters
        ----------

        Direction : string
            The direction in which to find the relative deflection. Must be one of the following:
                'dy' = Deflection in the local y-axis
                'dz' = Deflection in the local x-axis
        x : number
            The location at which to find the relative deflection
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        """
        
        member, x_mod = self.find_member(x)
        return member.rel_deflection(Direction, x_mod, combo_name)

    def plot_deflection(self, Direction: Literal['dx', 'dy', 'dz'], combo_name: str = 'Combo 1', n_points: int = 20) -> None:
        """
        Plots the deflection diagram for the member

        Parameters
        ----------
        Direction : string
            The direction in which to plot the deflection. Must be one of the following:
                'dy' = Deflection in the local y-axis.
                'dz' = Deflection in the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        n_points: int
            The number of points used to generate the plot
        """

        # Import 'pyplot' if not already done
        if PhysMember.__plt is None:
            from matplotlib import pyplot as plt
            PhysMember.__plt = plt

        fig, ax = PhysMember.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()

        d_array = self.deflection_array(Direction, n_points, combo_name)
        x = d_array[0]
        d = d_array[1]

        PhysMember.__plt.plot(x, d)
        PhysMember.__plt.ylabel('Deflection')
        PhysMember.__plt.xlabel('Location')
        PhysMember.__plt.title('Member ' + self.name + '\n' + combo_name)
        PhysMember.__plt.show()

    def deflection_array(self, Direction: Literal['dx', 'dy', 'dz'], n_points: int, combo_name='Combo 1', x_array=None) -> NDArray[float64]:
        """
        Returns the array of the deflection in the physical member for the given direction

        Parameters
        ----------
        Direction : string
            The direction to plot the deflection for. Must be one of the following:
                'dx' = Deflection in the local x-direction (axial deflection)
                'dy' = Deflection in the local y-direction (usually the strong-axis).
                'dz' = Deflection in the local z-direction (usually the weak-axis).
        n_points: int
            The number of points in the array to generate over the full length of the member.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        x_array : array = None
            A custom array of x values that may be provided by the user, otherwise an array is generated. Values must be provided in local member coordinates (between 0 and L) and be in ascending order
        """

        # `d_array2` will be used to store the deflection values for the overall member
        d_array2 = empty((2, 1))

        # Create an array of locations along the physical member to obtain results at
        L = self.L()
        if x_array is None:
            # Create an array of evenly spaced points
            x_array = linspace(0, L, n_points)
        else:
            # Ensure the requested points are within the member
            if any(x_array < 0) or any(x_array > L):
                raise ValueError(f"All x values must be in the range 0 to {L}")

        # Step through each submember in the physical member
        x_o = 0
        for i, submember in enumerate(self.sub_members.values()):

            # Segment the submember into segments with mathematically continuous loads if not already done
            if submember._solved_combo is None or combo_name != submember._solved_combo.name:
                submember._segment_member(combo_name)
                submember._solved_combo = self.model.load_combos[combo_name]

            # Check if this is the last submember
            if i == len(self.sub_members.values()) - 1:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array <= x_o + submember.L())

            # Not the last submember
            else:

                # Find any points from `x_array` that lie along this submember
                filter = (x_array >= x_o) & (x_array < x_o + submember.L())

            x_subm_array = x_array[filter] - x_o

            # Check which axis is of interest
            if Direction == 'dx':
                d_array = self._extract_vector_results(submember.SegmentsZ, x_subm_array, 'axial_deflection')
            elif Direction == 'dy':
                d_array = self._extract_vector_results(submember.SegmentsZ, x_subm_array, 'deflection')
            elif Direction == 'dz':
                d_array = self._extract_vector_results(submember.SegmentsY, x_subm_array, 'deflection')
            else:
                raise ValueError(f"Direction must be 'dy' or 'dz'. {Direction} was given.")

            # Adjust from the submember's coordinate system to the physical member's coordinate system
            d_array[0] = [x_o + x for x in d_array[0]]

            # Add the submember deflection values to the overall member shear values in `d_array2`
            if i != 0:
                d_array2 = hstack((d_array2, d_array))
            else:
                d_array2 = d_array

            # Get the starting position of the next submember
            x_o += submember.L()

        # Return the results
        return d_array2

    def find_member(self, x: float) -> Tuple[Member3D, float]:
        """
        Returns the sub-member that the physical member's local point 'x' lies on, and 'x' modified for that sub-member's local coordinate system.
        """

        # Initialize a summation of sub-member lengths
        L = 0

        # Step through each sub-member (in order from start to end)
        for i, member in enumerate(self.sub_members.values()):

            # Sum the sub-member's length
            L += member.L()

            # Check if 'x' lies on this sub-member
            if x < L or (isclose(x, L) and i == len(self.sub_members.values()) - 1):

                # Return the sub-member, and a modified value for 'x' relative to the sub-member's
                # i-node
                return member, x - (L - member.L())

                # Exit the 'for' loop
                break
        else:
            raise ValueError(f"Location x={x} does not lie on this member")
