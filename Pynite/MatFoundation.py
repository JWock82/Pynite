"""Mat foundation mesh and loading helpers.

This module defines the ``MatFoundation`` class, a thin convenience wrapper around ``Pynite.Mesh.RectangleMesh`` for modeling mat (raft) foundations. It adds utilities to:

- register rectangular cutouts/openings in the mat
- register concentrated loads at X-Z coordinates on the mat surface
- generate soil springs at every mesh node using a Winkler foundation model
  with tributary areas and a modulus of subgrade reaction ``ks``

Typical workflow:

1. Create a ``MatFoundation`` tied to an existing ``FEModel3D`` instance
2. Add optional cutouts and point loads
3. Call ``generate()`` to build the mesh, place the loads, and define springs

Notes
- The mat is meshed in the ``XZ`` plane by default so node vertical springs are defined about the model's global Y direction (``DY``).
- Point load directions must match the host model's accepted node load keys (for example: ``FX``, ``FY``, ``FZ``, ``MX``, ``MY``, ``MZ``).
"""

import numpy as np

from Pynite.Mesh import RectangleMesh


class MatFoundation(RectangleMesh):
    """Mat foundation modeled as a rectangular shell mesh on Winkler springs.

    This class specializes ``RectangleMesh`` for mats by tracking point loads defined in the mat's local X-Z plane and by assigning vertical soil springs to each generated node based on tributary area and a supplied modulus of subgrade reaction ``ks``.
    """

    def __init__(self, name, mesh_size, length_X, length_Z, thickness, material_name, model, ks, origin=[0, 0, 0], x_control=None, y_control=None):
        """Initialize a new mat foundation mesh.

        :param str name: Unique name for the mat mesh within the model.
        :param float mesh_size: Target element size used by the mesher.
        :param float length_X: Mat length in the global X direction.
        :param float length_Z: Mat length in the global Z direction.
        :param float thickness: Plate/shell thickness for the mat elements.
        :param str material_name: Name of an existing material in the model.
        :param model: The host ``FEModel3D`` to which this mesh belongs.
        :param float ks: Modulus of subgrade reaction for soil springs (units consistent with model, e.g., force/length^3).
        :param list origin: Global coordinates of the mesh origin ``[X, Y, Z]``. Defaults to ``[0, 0, 0]``.
        :param list x_control: Optional control points along local x (global X) used to force node locations in the mesh.
        :param list y_control: Optional control points along local y (global Z) used to force node locations in the mesh.
        :rtype: None
        """

        super().__init__(mesh_size, length_X, length_Z, thickness, material_name, model, 1, 1, origin, 'XZ', x_control, y_control)

        self.name = name
        self.ks = ks
        self.pt_loads = []  # [XZ_coord, direction, magnitude, case]

    def add_rect_opening(self, name, X_min, Z_min, X_max, Z_max):
        """Add a rectangular opening to the mat by corner coordinates.

        The opening is defined in the model's global X-Z plane using minimum and maximum corner coordinates. Internally this delegates to the ``RectangleMesh.add_rect_opening`` method that uses width/height.

        :param str name: Unique opening name within this mesh.
        :param float X_min: Lower-left X coordinate of the opening.
        :param float Z_min: Lower-left Z coordinate of the opening.
        :param float X_max: Upper-right X coordinate of the opening.
        :param float Z_max: Upper-right Z coordinate of the opening.
        :return: ``None``
        :rtype: None
        """

        super().add_rect_opening(name, X_min, Z_min, X_max - X_min, Z_max - Z_min)

    def add_mat_pt_load(self, XZ_coord, direction, magnitude, case='Case 1'):
        """Register a concentrated load at an X-Z coordinate on the mat.

        Loads are stored until ``generate()`` is called, at which time they are applied to the node that coincides with the specified coordinate. The given X and Z are also inserted as mesh control points to ensure a node exists at the target location.

        :param list XZ_coord: Two floats ``[X, Z]`` in the mat plane.
        :param str direction: Load direction key understood by the model's ``add_node_load`` API (eg., ``FX``, ``FY``, ``FZ``, ``MX``, ``MY``, ``MZ``).
        :param float magnitude: Load magnitude (units consistent with the model).
        :param str case: Load case name to assign the load to. Defaults to ``'Case 1'``.
        :return: ``None``
        :rtype: None
        """

        self.x_control.append(XZ_coord[0])
        self.y_control.append(XZ_coord[1])
        self.pt_loads.append([XZ_coord, direction, magnitude, case])

    def generate(self):
        """Generate the mesh, apply point loads, and define soil springs.

        Actions performed:

        - Calls the base ``RectangleMesh.generate()`` to create nodes/elements.
        - Applies any registered point loads to nodes that coincide with the specified X-Z coordinates (within floating-point tolerance).
        - Computes a tributary area for each node as one quarter of the area of each attached plate and assigns a vertical spring in ``DY`` with stiffness ``ks * tributary_area`` to represent subgrade reaction.

        This method mutates the attached model by adding node loads and defining support springs.

        :return: ``None``
        :rtype: None
        """

        # Generate the mesh
        super().generate()

        # Add point loads to the model
        for node in self.nodes.values():
            for pt_load in self.pt_loads:
                if np.isclose(node.X, pt_load[0][0]) and np.isclose(node.Z, pt_load[0][1]):
                    self.model.add_node_load(node.name, pt_load[1], pt_load[2], pt_load[3])

        # Step through each node in the mat
        for node in self.nodes.values():

            # Initialize the tributary area to the node to zero
            trib = 0

            # Step through each plate in the model
            for plate in self.elements.values():

                # Determine if the plate is attached to the node
                if node.name in [plate.i_node.name, plate.j_node.name, plate.m_node.name, plate.n_node.name]:

                    # Add 1/4 the plate's area to the tributary area to the node
                    trib += abs(plate.j_node.X - plate.i_node.X)*abs(plate.m_node.Z - plate.j_node.Z)/4

            # Add a soil spring to the node
            self.model.def_support_spring(node.name, 'DY', self.ks*trib, '-')

    def soil_pressure(self, x, y, combo_name: str = 'Combo 1'):
        """Return the soil contact pressure at a point in the mat.

        The mat is modeled on Winkler springs, so contact pressure at a point is taken as the surrounding nodal vertical reactions divided by their tributary areas, bilinearly interpolated within the plate element that contains the point.

        This requires analysis results to have been computed so that nodal reactions (``RxnFY``) are available for the specified load combination.

        :param float x: Global X of the query point (mat plane).
        :param float y: Global Z of the query point (mat plane).
        :param str combo_name: Load combination to read reactions from.
        :return: Soil contact pressure at ``(x, y)``. ``None`` if the point is outside the mat mesh or falls within an opening.
        :rtype: float | 0
        :raises RuntimeError: If analysis results are unavailable for ``combo_name``.
        """

        # Helper: compute area of a rectangular element
        def plate_area(plate):
            length_x = abs(plate.j_node.X - plate.i_node.X)
            length_y = abs(plate.n_node.Z - plate.i_node.Z)
            return length_x*length_y

        # Helper: tributary area for a node (sum 1/4 of connected plate areas)
        def tributary_area(node):
            area = 0.0
            for plate in self.elements.values():
                if node.name in (plate.i_node.name, plate.j_node.name, plate.m_node.name, plate.n_node.name):
                    area += plate_area(plate) / 4.0
            return area

        # Find the plate (rectangular) whose X-Z bounds contain the point
        containing_plate = None
        for plate in self.elements.values():
            X_vals = [plate.i_node.X, plate.j_node.X, plate.m_node.X, plate.n_node.X]
            Z_vals = [plate.i_node.Z, plate.j_node.Z, plate.m_node.Z, plate.n_node.Z]
            Xmin, Xmax = min(X_vals), max(X_vals)
            Zmin, Zmax = min(Z_vals), max(Z_vals)
            # Include boundary points with small tolerance
            if (x >= Xmin - 1e-9) and (x <= Xmax + 1e-9) and (y >= Zmin - 1e-9) and (y <= Zmax + 1e-9):
                containing_plate = plate
                break

        # If no plate contains the point, we are outside the mesh or in an opening
        if containing_plate is None:
            return 0

        # Ensure reactions exist for this combo
        if self.model.solution is None:
            raise RuntimeError(f'No nodal reactions found for combination {combo_name!r}. Run analysis first.')

        # Compute nodal pressures at each node on the containing plate = RxnFY / tributary area
        Rxn_i = containing_plate.i_node.RxnFY[combo_name]
        Rxn_j = containing_plate.j_node.RxnFY[combo_name]
        Rxn_m = containing_plate.m_node.RxnFY[combo_name]
        Rxn_n = containing_plate.n_node.RxnFY[combo_name]

        trib_i = tributary_area(containing_plate.i_node)
        trib_j = tributary_area(containing_plate.j_node)
        trib_m = tributary_area(containing_plate.m_node)
        trib_n = tributary_area(containing_plate.n_node)

        p_i = Rxn_i/trib_i
        p_j = Rxn_j/trib_j
        p_m = Rxn_m/trib_m
        p_n = Rxn_n/trib_n

        # Use bilinear interpolation for the point lying between the nodes
        # xi and eta measure distance from the i node in the local x and y direction in terms of the element's length in each direction
        xi_denom = (containing_plate.j_node.X - containing_plate.i_node.X)
        eta_denom = (containing_plate.n_node.Z - containing_plate.i_node.Z)

        xi = (x - containing_plate.i_node.X) / xi_denom
        eta = (y - containing_plate.i_node.Z) / eta_denom

        # Clamp to [0, 1] to be robust to small numerical tolerances
        xi = min(max(xi, 0.0), 1.0)
        eta = min(max(eta, 0.0), 1.0)

        # Shape functions for a bilinear quad aligned with axes
        N_i = (1 - xi) * (1 - eta)
        N_j = xi * (1 - eta)
        N_m = xi * eta
        N_n = (1 - xi) * eta

        # Calculate and return the interpolated pressure
        p = N_i*p_i + N_j*p_j + N_m*p_m + N_n*p_n
        return p
