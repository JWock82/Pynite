from __future__ import annotations # Allows more recent type hints features
from typing import TYPE_CHECKING
from math import pi, sin, cos, ceil, isclose

from Pynite.Node3D import Node3D
from Pynite.Quad3D import Quad3D
from Pynite.Plate3D import Plate3D

if TYPE_CHECKING:
    from typing import List, Union, Dict
    from Pynite.FEModel3D import FEModel3D

#%%
class Mesh():
    """
    A parent class for meshes to inherit from.
    """

    def __init__(self, thickness: float, material_name: str, model: FEModel3D, kx_mod: float = 1, ky_mod: float = 1, start_node: str = 'N1', start_element: str = 'Q1') -> None:
        """Base mesh container storing common settings and membership.

        :param thickness: Element thickness used by generated elements.
        :type thickness: float
        :param material_name: Name of the element material.
        :type material_name: str
        :param model: Owning FEModel3D; nodes/elements are inserted here on generate().
        :type model: FEModel3D
        :param kx_mod: In-plane stiffness modifier along local x for generated elements. Defaults to 1.
        :type kx_mod: float, optional
        :param ky_mod: In-plane stiffness modifier along local y for generated elements. Defaults to 1.
        :type ky_mod: float, optional
        :param start_node: First node name to use when numbering; e.g., ``'N1'``. Defaults to ``'N1'``.
        :type start_node: str, optional
        :param start_element: First element name to use when numbering; e.g., ``'Q1'`` or ``'R1'``. Defaults to ``'Q1'``.
        :type start_element: str, optional
        """

        self.thickness = thickness          # Thickness
        self.material_name = material_name            # The name of the element material
        self.model = model                  # Meshes need a link to the model they belong to
        self.kx_mod = kx_mod                # Local x stiffness modification factor for elements in the mesh
        self.ky_mod = ky_mod                # Local y stiffness modification factor for elements in the mesh
        self.start_node = start_node        # The name of the first node in the mesh
        self.last_node = None               # The name of the last node in the mesh
        self.start_element = start_element  # The name of the first element in the mesh
        self.last_element = None            # The name of the last element in the mesh
        self.nodes: Dict[str, Node3D] = {}  # A dictionary containing the nodes in the mesh
        self.elements: Dict[str, Union[Quad3D, Plate3D]] = {}  # A dictionary containing the elements in the mesh
        self.element_type = 'Quad'          # The type of element used in the mesh
        self.is_generated = False          # A flag indicating whether the mesh has been generated
    
    def _rename_duplicates(self) -> None:
        """Renames any nodes or elements in the mesh that are already in the model
        """

        # Initialize lists to track node and element name changes
        revised_nodes: Dict[str, Node3D] = {}
        revised_elements: Dict[str, Union[Quad3D, Plate3D]] = {}

        # Step through each node in the mesh
        for node in self.nodes.values():

            # Check if this node name is already being used in the model's `Nodes` dictionary
            if node.name in self.model.nodes.keys():

                # Come up with a new node name
                node.name = self.model.unique_name(self.model.nodes, 'N')
            
            # Save the node to the model
            self.model.nodes[node.name] = node

            # Add this node to the mesh's new/replacement `nodes` dictionary
            revised_nodes[node.name] = node
            
        # Step through each element in the mesh
        for element in self.elements.values():

            # Check which type of element this is
            if element.type == 'Rect':

                # Check if this element name is already being used in the model's `Plates` dictionary
                if element.name in self.model.plates.keys():

                    # Come up with a new element name
                    element.name = self.model.unique_name(self.model.plates, 'R')

                # Save the element to the model
                self.model.plates[element.name] = element

            elif element.type == 'Quad':

                # Check if this element name is already being used in the model's `Quads` dictionary
                if element.name in self.model.quads.keys():

                    # Come up with a new element name
                    element.name = self.model.unique_name(self.model.quads, prefix='Q')

                # Save the element to the model
                self.model.quads[element.name] = element
            
            # Add this element to the mesh's new/replacement `elements` dictionary
            revised_elements[element.name] = element
        
        # Replace the old dictionaries of nodes and elements with the revised dictionaries
        self.nodes = revised_nodes
        self.elements = revised_elements
            
    def generate(self) -> None:
        """
        A placeholder to be overwritten by subclasses inheriting from this class
        """
        pass

    def max_shear(self, direction: str = 'Qx', combo_tags: str | list[str] = 'Combo 1') -> float:
        """Returns the maximum shear in the mesh.

        Checks corner and center shears in all elements contained in the mesh. The mesh
        must belong to a solved model prior to calling this method.

        :param direction: Shear component to evaluate. Use local components ``'Qx'`` or
            ``'Qy'``, or global components ``'QX'`` or ``'QY'``. Defaults to ``'Qx'``.
        :type direction: str
        :param combo_tags: Either a single load combination name (``str``) or a list of
            tags (``list[str]``). If a list is provided, any load combination with at
            least one of these tags will be considered. Defaults to ``'Combo 1'``.
        :type combo_tags: str | list[str]
        :return: The maximum shear found for the requested component (0.0 if none).
        :rtype: float
        """

        # Determine if the shear is requested in local or global axes
        if direction in ['QX', 'QY']:
            local = False
        else:
            local = True

        # Map direction string to index in element.shear() results
        if direction.upper() == 'QX':
            i = 0
        elif direction.upper() == 'QY':
            i = 1
        else:
            raise Exception(
                "Invalid direction specified for mesh shear results. "
                "Valid values are 'Qx', 'Qy', 'QX', or 'QY'."
            )

        # Initialize the maximum value
        Q_max = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine evaluation points depending on element type
            if element.type == 'Rect':
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn = 0, element.height()
            elif element.type == 'Quad':
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1
            else:
                continue  # Skip unsupported element types

            # Step through each load combination
            for load_combo in self.model.load_combos.values():

                # Decide whether this combo should be checked
                if isinstance(combo_tags, str):
                    include = (load_combo.name == combo_tags)
                elif isinstance(combo_tags, list):
                    include = any(tag in load_combo.tags for tag in combo_tags)
                else:
                    include = False

                if not include:
                    continue

                # Evaluate corner and center shears
                Q_element = max([
                    element.shear(xi, yi, local, load_combo.name)[i, 0],
                    element.shear(xj, yj, local, load_combo.name)[i, 0],
                    element.shear(xm, ym, local, load_combo.name)[i, 0],
                    element.shear(xn, yn, local, load_combo.name)[i, 0],
                    # Center point (average of corners)
                    element.shear(
                        (xi + xj + xm + xn) / 4,
                        (yi + yj + ym + yn) / 4,
                        local,
                        load_combo.name
                    )[i, 0]
                ])

                # Track global maximum
                if Q_max is None or Q_max < Q_element:
                    Q_max = Q_element

        # Return the largest shear found, or 0.0 if nothing was found
        return 0.0 if Q_max is None else Q_max

    def min_shear(self, direction: str = 'Qx', combo_tags: str | list[str] = 'Combo 1') -> float:
        """Returns the minimum shear in the mesh.

        Checks corner and center shears in all elements contained in the mesh. The mesh
        must belong to a solved model prior to calling this method.

        :param direction: Shear component to evaluate. Use local components ``'Qx'`` or
            ``'Qy'``, or global components ``'QX'`` or ``'QY'``. Defaults to ``'Qx'``.
        :type direction: str
        :param combo_tags: Either a single load combination name (``str``) or a list of
            tags (``list[str]``). If a list is provided, any load combination with at
            least one of these tags will be considered. Defaults to ``'Combo 1'``.
        :type combo_tags: str | list[str]
        :return: The minimum shear found for the requested component (0.0 if none).
        :rtype: float
        """

        # Determine if the shear is requested in local or global axes
        if direction in ['QX', 'QY']:
            local = False
        else:
            local = True

        # Map direction string to index in element.shear() results
        if direction.upper() == 'QX':
            i = 0
        elif direction.upper() == 'QY':
            i = 1
        else:
            raise Exception(
                "Invalid direction specified for mesh shear results. "
                "Valid values are 'Qx', 'Qy', 'QX', or 'QY'."
            )

        # Initialize the minimum value
        Q_min = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine evaluation points depending on element type
            if element.type == 'Rect':
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn = 0, element.height()
            elif element.type == 'Quad':
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1
            else:
                continue  # Skip unsupported element types

            # Step through each load combination
            for load_combo in self.model.load_combos.values():

                # Decide whether this combo should be checked
                if isinstance(combo_tags, str):
                    include = (load_combo.name == combo_tags)
                elif isinstance(combo_tags, list):
                    include = any(tag in load_combo.tags for tag in combo_tags)
                else:
                    include = False

                if not include:
                    continue

                # Evaluate corner and center shears
                Q_element = min([
                    element.shear(xi, yi, local, load_combo.name)[i, 0],
                    element.shear(xj, yj, local, load_combo.name)[i, 0],
                    element.shear(xm, ym, local, load_combo.name)[i, 0],
                    element.shear(xn, yn, local, load_combo.name)[i, 0],
                    # Center point (average of corners)
                    element.shear(
                        (xi + xj + xm + xn) / 4,
                        (yi + yj + ym + yn) / 4,
                        local,
                        load_combo.name
                    )[i, 0]
                ])

                # Track global minimum
                if Q_min is None or Q_min > Q_element:
                    Q_min = Q_element

        # Return the smallest shear found, or 0.0 if nothing was found
        return 0.0 if Q_min is None else Q_min

    def max_moment(self, direction: str = 'Mx', combo_tags: str | list[str] = 'Combo 1') -> float:
        """Returns the maximum moment in the mesh.

        Checks corner and center moments in all elements contained in the mesh. The mesh
        must belong to a solved model prior to calling this method.

        :param direction: Moment component to evaluate. Use local components ``'Mx'``,
            ``'My'``, ``'Mxy'`` or global components ``'MX'``, ``'MY'``, ``'MZ'``.
            Defaults to ``'Mx'``.
        :type direction: str
        :param combo_tags: Either a single load combination name (``str``) or a list of
            tags (``list[str]``). If a list is provided, any load combination with at
            least one of these tags will be considered. Defaults to ``'Combo 1'``.
        :type combo_tags: str | list[str]
        :return: The maximum moment found for the requested component (0.0 if none).
        :rtype: float
        """

        # Determine if the moment is requested in local or global axes
        if direction in ['MX', 'MY', 'MZ']:
            local = False
        else:
            local = True

        # Map direction string to index in element.moment() results
        if direction.upper() == 'MX':
            i = 0
        elif direction.upper() == 'MY':
            i = 1
        elif direction == 'Mxy' or direction.upper() == 'MZ':
            i = 2
        else:
            raise Exception(
                "Invalid direction specified for mesh moment results. "
                "Valid values are 'Mx', 'My', 'Mxy', 'MX', 'MY', or 'MZ'."
            )

        # Initialize the maximum value
        M_max = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine evaluation points depending on element type
            if element.type == 'Rect':
                # Use the rectangle's local (x, y) coordinate system
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn = 0, element.height()
            elif element.type == 'Quad':
                # Use the quad's natural (xi, eta) coordinate system
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1
            else:
                continue  # Skip unsupported element types

            # Step through each load combination
            for load_combo in self.model.load_combos.values():

                # Decide whether this combo should be checked
                if isinstance(combo_tags, str):
                    include = (load_combo.name == combo_tags)
                elif isinstance(combo_tags, list):
                    include = any(tag in load_combo.tags for tag in combo_tags)
                else:
                    include = False

                if not include:
                    continue

                # Evaluate corner and center moments
                M_element = max([
                    element.moment(xi, yi, local, load_combo.name)[i, 0],
                    element.moment(xj, yj, local, load_combo.name)[i, 0],
                    element.moment(xm, ym, local, load_combo.name)[i, 0],
                    element.moment(xn, yn, local, load_combo.name)[i, 0],
                    # Center point (average of corners)
                    element.moment((xi + xj + xm + xn) / 4, (yi + yj + ym + yn) / 4,
                                local, load_combo.name)[i, 0]
                ])

                # Track global maximum
                if M_max is None or M_max < M_element:
                    M_max = M_element

        # Return the largest moment found, or 0.0 if nothing was found
        return 0.0 if M_max is None else M_max

    def min_moment(self, direction: str = 'Mx', combo_tags: str | list[str] = 'Combo 1') -> float:
        """Returns the minimum moment in the mesh.

        Checks corner and center moments in all elements contained in the mesh. The mesh
        must belong to a solved model prior to calling this method.

        :param direction: Moment component to evaluate. Use local components ``'Mx'``,
            ``'My'``, ``'Mxy'`` or global components ``'MX'``, ``'MY'``, ``'MZ'``.
            Defaults to ``'Mx'``.
        :type direction: str
        :param combo_tags: Either a single load combination name (``str``) or a list of
            tags (``list[str]``). If a list is provided, any load combination with at
            least one of these tags will be considered. Defaults to ``'Combo 1'``.
        :type combo_tags: str | list[str]
        :return: The minimum moment found for the requested component (0.0 if none).
        :rtype: float
        """

        # Determine if the moment is requested in local or global axes
        if direction in ['MX', 'MY', 'MZ']:
            local = False
        else:
            local = True

        # Map direction string to index in element.moment() results
        if direction.upper() == 'MX':
            i = 0
        elif direction.upper() == 'MY':
            i = 1
        elif direction == 'Mxy' or direction.upper() == 'MZ':
            i = 2
        else:
            raise Exception(
                "Invalid direction specified for mesh moment results. "
                "Valid values are 'Mx', 'My', 'Mxy', 'MX', 'MY', or 'MZ'."
            )

        # Initialize the minimum value
        M_min = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine evaluation points depending on element type
            if element.type == 'Rect':
                # Use the rectangle's local (x, y) coordinate system
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn = 0, element.height()
            elif element.type == 'Quad':
                # Use the quad's natural (xi, eta) coordinate system
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1
            else:
                continue  # Skip unsupported element types

            # Step through each load combination
            for load_combo in self.model.load_combos.values():

                # Decide whether this combo should be checked
                if isinstance(combo_tags, str):
                    include = (load_combo.name == combo_tags)
                elif isinstance(combo_tags, list):
                    include = any(tag in load_combo.tags for tag in combo_tags)
                else:
                    include = False

                if not include:
                    continue

                # Evaluate corner and center moments
                M_element = min([
                    element.moment(xi, yi, local, load_combo.name)[i, 0],
                    element.moment(xj, yj, local, load_combo.name)[i, 0],
                    element.moment(xm, ym, local, load_combo.name)[i, 0],
                    element.moment(xn, yn, local, load_combo.name)[i, 0],
                    # Center point (average of corners)
                    element.moment((xi + xj + xm + xn) / 4, (yi + yj + ym + yn) / 4,
                                local, load_combo.name)[i, 0]
                ])

                # Track global minimum
                if M_min is None or M_min > M_element:
                    M_min = M_element

        # Return the smallest moment found, or 0.0 if nothing was found
        return 0.0 if M_min is None else M_min

    def max_membrane(self, direction: str = 'Sx', combo_tags: str | list[str] = 'Combo 1') -> float:
        """Returns the maximum membrane stress in the mesh.

        Checks corner and center moments in all elements contained in the mesh. The mesh
        must belong to a solved model prior to calling this method.

        :param direction: Membrane stress component to evaluate. Use local components ``'Sx'``,
            ``'Sy'``, ``'Txy'``.
        :type direction: str
        :param combo_tags: Either a single load combination name (``str``) or a list of
            tags (``list[str]``). If a list is provided, any load combination with at
            least one of these tags will be considered. Defaults to ``'Combo 1'``.
        :type combo_tags: str | list[str]
        :return: The maximum stress found for the requested component (0.0 if none).
        :rtype: float
        """

        # Determine if local or global coordinate results have been requested
        if direction in ['SX', 'SY']:
            local = False
        else:
            local = True

        # Identify which stress component needs to be extracted
        if direction.upper() == 'SX':
            i = 0
        elif direction.upper() == 'SY':
            i = 1
        elif direction == 'Sxy':
            i = 2
        else:
            raise Exception(
                "Invalid direction specified for mesh membrane stress results. "
                "Valid values are 'Sx', 'Sy', or 'Sxy'."
            )

        # Initialize the maximum value
        S_max = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine evaluation points depending on element type
            if element.type == 'Rect':
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn = 0, element.height()
            elif element.type == 'Quad':
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1
            else:
                continue  # Skip unsupported element types

            # Step through each load combination
            for load_combo in self.model.load_combos.values():

                # Decide whether this combo should be checked
                if isinstance(combo_tags, str):
                    include = (load_combo.name == combo_tags)
                elif isinstance(combo_tags, list):
                    include = any(tag in load_combo.tags for tag in combo_tags)
                else:
                    include = False

                if not include:
                    continue

                # Evaluate corner and center membrane stresses
                S_element = max([
                    element.membrane(xi, yi, local, load_combo.name)[i, 0],
                    element.membrane(xj, yj, local, load_combo.name)[i, 0],
                    element.membrane(xm, ym, local, load_combo.name)[i, 0],
                    element.membrane(xn, yn, local, load_combo.name)[i, 0],
                    # Center point (average of corners)
                    element.membrane(
                        (xi + xj + xm + xn) / 4,
                        (yi + yj + ym + yn) / 4,
                        local,
                        load_combo.name
                    )[i, 0]
                ])

                # Track global maximum
                if S_max is None or S_max < S_element:
                    S_max = S_element

        # Return the largest membrane stress found, or 0.0 if nothing was found
        return 0.0 if S_max is None else S_max
    
    def min_membrane(self, direction: str = 'Sx', combo_tags: str | list[str] = 'Combo 1') -> float:
        """Returns the minimum membrane stress in the mesh.

        Checks corner and center moments in all elements contained in the mesh. The mesh
        must belong to a solved model prior to calling this method.

        :param direction: Membrane stress component to evaluate. Use local components ``'Sx'``,
            ``'Sy'``, ``'Txy'``.
        :type direction: str
        :param combo_tags: Either a single load combination name (``str``) or a list of
            tags (``list[str]``). If a list is provided, any load combination with at
            least one of these tags will be considered. Defaults to ``'Combo 1'``.
        :type combo_tags: str | list[str]
        :return: The minimum stress found for the requested component (0.0 if none).
        :rtype: float
        """

        # Determine if local or global coordinate results have been requested
        if direction in ['SX', 'SY']:
            local = False
        else:
            local = True

        # Identify which stress component needs to be extracted
        if direction.upper() == 'SX':
            i = 0
        elif direction.upper() == 'SY':
            i = 1
        elif direction == 'Sxy':
            i = 2
        else:
            raise Exception(
                "Invalid direction specified for mesh membrane stress results. "
                "Valid values are 'Sx', 'Sy', or 'Sxy'."
            )

        # Initialize the maximum value
        S_min = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine evaluation points depending on element type
            if element.type == 'Rect':
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn = 0, element.height()
            elif element.type == 'Quad':
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1
            else:
                continue  # Skip unsupported element types

            # Step through each load combination
            for load_combo in self.model.load_combos.values():

                # Decide whether this combo should be checked
                if isinstance(combo_tags, str):
                    include = (load_combo.name == combo_tags)
                elif isinstance(combo_tags, list):
                    include = any(tag in load_combo.tags for tag in combo_tags)
                else:
                    include = False

                if not include:
                    continue

                # Evaluate corner and center membrane stresses
                S_element = min([
                    element.membrane(xi, yi, local, load_combo.name)[i, 0],
                    element.membrane(xj, yj, local, load_combo.name)[i, 0],
                    element.membrane(xm, ym, local, load_combo.name)[i, 0],
                    element.membrane(xn, yn, local, load_combo.name)[i, 0],
                    # Center point (average of corners)
                    element.membrane(
                        (xi + xj + xm + xn) / 4,
                        (yi + yj + ym + yn) / 4,
                        local,
                        load_combo.name
                    )[i, 0]
                ])

                # Track global maximum
                if S_min is None or S_min > S_element:
                    S_min = S_element

        # Return the largest membrane stress found, or 0.0 if nothing was found
        return 0.0 if S_min is None else S_min
    
#%%
class RectangleMesh(Mesh):

    def __init__(self, mesh_size: float, width: float, height: float, thickness: float, material_name: str, model: FEModel3D, kx_mod: float = 1.0, ky_mod: float = 1.0, origin: List[float] = [0, 0, 0], plane: str = 'XY', x_control: List[float] | None = None, y_control: List[float] | None = None, start_node: str = 'N1', start_element: str = 'Q1', element_type: str = 'Quad') -> None:
        """
        A rectangular mesh of elements.

        :param mesh_size: Desired mesh size.
        :type mesh_size: float
        :param width: The overall width of the mesh measured along its local x-axis.
        :type width: float
        :param height: The overall height of the mesh measured along its local y-axis.
        :type height: float
        :param thickness: Element thickness.
        :type thickness: float
        :material_name: The name of the element material.
        :type material_name: str
        :param model: The model the mesh belongs to.
        :type model: FEModel3D
        :param kx_mod: Stiffness modification factor for in-plane stiffness in the element's local
            x-direction. Default value is 1.0 (no modification).
        :type kx_mod: float, optional
        :param ky_mod: Stiffness modification factor for in-plane stiffness in the element's local
            y-direction. Default value is 1.0 (no modification).
        :type ky_mod: float, optional
        :param origin: The origin of the rectangular mesh's local coordinate system. The default is [0.0, 0.0, 0.0].
        :type origin: list[float], optional
        :param plane: The plane the mesh will be parallel to. Options are 'XY', 'YZ', and 'XZ'. The default is 'XY'.
        :type plane: str, optional
        :param x_control: A list of control points along the mesh's local x-axis work into the mesh.
        :type x_control: list[float], optional
        :param y_control: A list of control points along the mesh's local y-axis work into the mesh.
        :type y_control: list[float], optional
        :param start_node: A unique name for the first node in the mesh. The default is 'N1'.
        :type start_node: str, optional
        :param start_element: A unique name for the first element in the mesh. The default is 'Q1' or 'R1' depending on the type of element selected.
        :type start_element: str, optional
        :param element_type: The type of element to make the mesh out of. Either 'Quad' or 'Rect'. The default is 'Quad'.
        :type element_type: str, optional
        :return: A new rectangular mesh object.
        :rtype: RectangleMesh
        """

        super().__init__(thickness, material_name, model, kx_mod, ky_mod, start_node, start_element)
        self.mesh_size = mesh_size
        self.width = width
        self.height = height
        self.origin = origin
        self.plane = plane

        if x_control is None: self.x_control = []
        else: self.x_control = x_control

        if y_control is None: self.y_control = []
        else: self.y_control = y_control

        self.element_type = element_type
        self.openings: Dict[str, RectOpening] = {}
    
    def generate(self) -> None:

        mesh_size = self.mesh_size
        width = self.width
        height = self.height
        Xo = self.origin[0]
        Yo = self.origin[1]
        Zo = self.origin[2]
        plane = self.plane
        x_control = self.x_control
        y_control = self.y_control

        element_type = self.element_type

        # Add the mesh's boundaries to the list of control points
        x_control.append(0)
        x_control.append(width)
        y_control.append(0)
        y_control.append(height)

        # Sort the control points in ascending order
        x_control = sorted(x_control)
        y_control = sorted(y_control)

        # Remove any values that are duplicates or near duplicates from `x_control`
        unique_list = []
        for i in range(len(x_control) - 1):
            # Only keep the value at `i` if it's not a duplicate or near duplicate of the next value
            if not isclose(x_control[i], x_control[i+1]):
                unique_list.append(x_control[i])
        unique_list.append(x_control[-1])
        x_control = unique_list

        # Remove any values that are duplicates or near duplicates from `y_control`
        unique_list = []
        for i in range(len(y_control) - 1):
            # Only keep the value at `i` if it's not a duplicate or near duplicate of the next value
            if not isclose(y_control[i], y_control[i+1]):
                unique_list.append(y_control[i])
        unique_list.append(y_control[-1])
        y_control = unique_list
    
        # Each node number will be increased by the offset calculated below
        node_offset = int(self.start_node[1:]) - 1

        # Each element number will be increased by the offset calculated below
        element_offset = int(self.start_element[1:]) - 1

        # Determine which prefix to assign to new elements
        if element_type == 'Quad':
            element_prefix = 'Q'
        elif element_type == 'Rect':
            element_prefix = 'R'
        else:
            raise Exception('Invalid element type specified for RectangleMesh. Select \'Quad\' or \'Rect\'.')

        # Initialize node numbering
        node_num = 1

        # Step through each y control point (except the first one which is always zero)
        num_rows = 0
        num_cols = 0
        y, h = 0, None
        for j in range(1, len(y_control), 1):
            
            # If this is not the first iteration 'y' will be too high at this point.
            if j != 1:
                y -= h

            # Determine the mesh size between this y control point and the previous one
            ny = max(1, (y_control[j] - y_control[j - 1])/mesh_size)
            h = (y_control[j] - y_control[j - 1])/ceil(ny)

            # Adjust 'y' if this is not the first iteration.
            if j != 1:
                y += h

            # Generate nodes between the y control points
            while round(y, 10) <= round(y_control[j], 10):
                
                # Count the number of rows of plates as we go
                num_rows += 1

                # Step through each x control point (except the first one which is always zero)
                x, b = 0, None
                for i in range(1, len(x_control), 1):
                    
                    # 'x' needs to be adjusted for the same reasons 'y' needed to be adjusted
                    if i != 1:
                        x -= b

                    # Determine the mesh size between this x control point and the previous one
                    nx = max(1, (x_control[i] - x_control[i - 1])/mesh_size)
                    b = (x_control[i] - x_control[i - 1])/ceil(nx)

                    if i != 1:
                        x += b

                    # Generate nodes between the x control points
                    while round(x, 10) <= round(x_control[i], 10):
                        
                        # Count the number of columns of plates as we go
                        if y == 0:
                            num_cols += 1

                        # Assign the node a name
                        node_name = 'N' + str(node_num + node_offset)

                        # Calculate the node's coordinates
                        if plane == 'XY':
                            X = Xo + x
                            Y = Yo + y
                            Z = Zo + 0
                        elif plane == 'YZ':
                            X = Xo + 0
                            Y = Yo + y
                            Z = Zo + x
                        elif plane == 'XZ':
                            X = Xo + x
                            Y = Yo + 0
                            Z = Zo + y
                        else:
                            raise Exception('Invalid plane selected for RectangleMesh.')

                        # Add the node to the mesh
                        self.nodes[node_name] = Node3D(node_name, X, Y, Z)

                        # Move to the next x coordinate
                        x += b

                        # Move to the next node number
                        node_num += 1

                # Move to the next y coordinate
                y += h
        
        # At this point `num_cols` and `num_rows` represent the number of columns and rows of
        # nodes. We'll adjust these variables to be the number of columns and rows of elements
        # instead.
        num_cols -= 1
        num_rows -= 1
        
        # Create the elements
        r = 1
        n = 1
        for i in range(1, num_cols*num_rows + 1, 1):

            # Assign the element a name
            element_name = element_prefix + str(i + element_offset)

            # Find the attached nodes
            i_node = n + (r - 1)
            j_node = i_node + 1
            m_node = j_node + (num_cols + 1)
            n_node = m_node - 1

            if i % num_cols == 0:
                r += 1
            
            n += 1
            
            if element_type == 'Quad':
                self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                                   self.nodes['N' + str(j_node + node_offset)],
                                                                   self.nodes['N' + str(m_node + node_offset)],
                                                                   self.nodes['N' + str(n_node + node_offset)],
                                                                   self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod)
            else:
                self.elements[element_name] = Plate3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                                    self.nodes['N' + str(j_node + node_offset)],
                                                                    self.nodes['N' + str(m_node + node_offset)],
                                                                    self.nodes['N' + str(n_node + node_offset)],
                                                                    self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod)

        # Initialize a list of nodes and associated elements that fall within opening boundaries
        # that will be deleted
        node_del_list = []
        element_del_list = []

        # Go back through the mesh and delete any nodes that are in the openings
        for node in self.nodes.values():
            
            # Get the node's position in the mesh's local coordinate sytem.
            x, y = self.node_local_coords(node)

            # Step through each opening in the mesh
            for opng in self.openings.values():

                # Determine if the node falls within the boundaries of the opening
                if (round(x, 10) > round(opng.x_left, 10)
                and round(x, 10) < round(opng.x_left + opng.width, 10) 
                and round(y, 10) > round(opng.y_bott, 10) 
                and round(y, 10) < round(opng.y_bott + opng.height, 10)):

                    # Mark the node for deletion if it's not already marked
                    if node.name not in node_del_list:
                        node_del_list.append(node.name)
                
        # Go back through the mesh and delete any elements that are in the openings
        for element in self.elements.values():

            # Find the top, bottom, left side and right side of the element in local coordinates
            left, top = self.node_local_coords(element.n_node)
            right, bott = self.node_local_coords(element.j_node)

            for opng in self.openings.values():

                # Determine if the element falls within the boundaries of the opening
                if ((round(opng.y_bott + opng.height, 10) >= round(top, 10))
                and (round(opng.y_bott, 10) <= round(bott, 10))
                and (round(opng.x_left, 10) <= round(left, 10))
                and (round(opng.x_left + opng.width, 10) >= round(right, 10))):

                    # Mark the element for deletion if it's not already marked
                    if element.name not in element_del_list:
                        element_del_list.append(element.name)

        # Delete the elements marked for deletion
        for element_name in element_del_list:
            del self.elements[element_name]

        # Delete the nodes marked for deletion
        for node_name in node_del_list:
            del self.nodes[node_name]
        
        # Find any remaining orphaned nodes around the perimeter of the mesh
        node_del_list = []
        for node in self.nodes.values():
            if (node not in [element.i_node for element in self.elements.values()]
            and node not in [element.j_node for element in self.elements.values()]
            and node not in [element.m_node for element in self.elements.values()]
            and node not in [element.n_node for element in self.elements.values()]):
                node_del_list.append(node.name)
        
        # Delete the orphaned nodes
        for node_name in node_del_list:
            del self.nodes[node_name]

        # Identify the last node and last element in the mesh
        self.last_node = list(self.nodes.values())[-1]
        self.last_element = list(self.elements.values())[-1]

        # At this point we have a mesh, but some of the element or node names may already be
        # being used in the model. Rename any names that are already being used.
        self._rename_duplicates()

        # Add the nodes to the model
        for key, node in self.nodes.items():
            self.model.nodes[key] = node
        
        # Add the elements to the model
        for key, element in self.elements.items():
            if element.type == 'Quad':
                self.model.quads[key] = element
            elif element.type == 'Rect':
                self.model.plates[key] = element

        # Flag the mesh as generated
        self.is_generated = True

    def node_local_coords(self, node: Node3D) -> tuple[float, float]:
        """Calculates a node's position in the mesh's local x/y coordinate system.

        :param node: The node to evaluate.
        :type node: Node3D
        :return: Local coordinates ``(x, y)`` measured from the mesh origin.
        :rtype: tuple[float, float]
        """

        if self.plane == 'XY':
            x = node.X - self.origin[0]
            y = node.Y - self.origin[1]
        elif self.plane == 'YZ':
            x = node.Z - self.origin[2]
            y = node.Y - self.origin[1]
        elif self.plane == 'XZ':
            x = node.X - self.origin[0]
            y = node.Z - self.origin[2]
        
        return x, y

    def add_rect_opening(self, name: str, x_left: float, y_bott: float, width: float, height: float) -> None:
        """Adds a rectangular opening to the mesh.

        :param name: Unique name for the opening (used as its key in ``mesh.openings``).
        :type name: str
        :param x_left: Local x-coordinate for the left edge of the opening.
        :type x_left: float
        :param y_bott: Local y-coordinate for the bottom edge of the opening.
        :type y_bott: float
        :param width: Opening width (local x-direction).
        :type width: float
        :param height: Opening height (local y-direction).
        :type height: float
        :return: None
        :rtype: NoneType
        """

        self.openings[name] = RectOpening(x_left, y_bott, width, height)
        self.x_control.append(x_left)
        self.y_control.append(y_bott)
        self.x_control.append(x_left + width)
        self.y_control.append(y_bott + height)

        # Flag the mesh as not generated yet
        self.is_generated = False

#%%
class RectOpening():
    """
    Represents a rectangular opening in a rectangular mesh.
    """

    def __init__(self, x_left: float, y_bott: float, width: float, height: float) -> None:
        """Create a rectangular opening descriptor.

        :param x_left: Local x-coordinate for the left edge of the opening.
        :type x_left: float
        :param y_bott: Local y-coordinate for the bottom edge of the opening.
        :type y_bott: float
        :param width: Opening width (local x-direction).
        :type width: float
        :param height: Opening height (local y-direction).
        :type height: float
        """

        self.x_left = x_left
        self.y_bott = y_bott
        self.width = width
        self.height = height

#%%           
class AnnulusMesh(Mesh):
    """
    A mesh of quadrilaterals forming an annulus (a donut).
    """

    def __init__(self, mesh_size: float, outer_radius: float, inner_radius: float, thickness: float, material_name: str, model: FEModel3D, kx_mod: float = 1,
        ky_mod: float = 1, origin: List[float] = [0, 0, 0], axis: str = 'Y', start_node: str = 'N1', start_element: str = 'Q1') -> None:

        """Annular (donut) mesh between inner and outer radii.

        :param mesh_size: Target element size used to seed circumferential and radial divisions.
        :type mesh_size: float
        :param outer_radius: Outer radius of the annulus.
        :type outer_radius: float
        :param inner_radius: Inner radius of the annulus.
        :type inner_radius: float
        :param thickness: Element thickness.
        :type thickness: float
        :param material_name: Name of the element material.
        :type material_name: str
        :param model: The FEModel3D this mesh will add nodes/elements to.
        :type model: FEModel3D
        :param kx_mod: In-plane stiffness modifier along the element's local x-direction. Defaults to 1.
        :type kx_mod: float, optional
        :param ky_mod: In-plane stiffness modifier along the element's local y-direction. Defaults to 1.
        :type ky_mod: float, optional
        :param origin: Local origin of the annulus in global coordinates ``[X, Y, Z]``. Defaults to ``[0, 0, 0]``.
        :type origin: list[float], optional
        :param axis: Global axis about which the mesh is generated: ``'X'``, ``'Y'``, or ``'Z'``. Defaults to ``'Y'``.
        :type axis: str, optional
        :param start_node: Name of the first node to use. Defaults to ``'N1'``.
        :type start_node: str, optional
        :param start_element: Name of the first element to use. Defaults to ``'Q1'``.
        :type start_element: str, optional
        """

        super().__init__(thickness, material_name, model, kx_mod, ky_mod, start_node, start_element)

        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.mesh_size = mesh_size
        self.origin = origin
        self.axis = axis

        self.num_quads_inner = None
        self.num_quads_outer = None
        
        # self.generate()
    
    def generate(self) -> None:
        
        mesh_size = self.mesh_size
        r_outer = self.outer_radius
        r_inner = self.inner_radius
        n = int(self.start_node[1:])
        q = int(self.start_element[1:])

        circumf = 2*pi*r_inner           # Circumference of the ring at the inner radius
        n_circ = int(circumf/mesh_size)  # Number of times `mesh_size` fits in the circumference
        self.num_quads_outer = n_circ

        # Mesh the annulus from the inside toward the outside
        while round(r_inner, 10) < round(r_outer, 10):
            
            radial = r_outer - r_inner                    # Remaining length in the radial direction to be meshed
            circumf = 2*pi*r_inner                        # Circumference of the ring at the inner radius
            b_circ = circumf/n_circ                       # Element width in the circumferential direction
            n_rad = int(radial/min(mesh_size, 3*b_circ))  # Number of times the plate width fits in the remaining unmeshed radial direction
            h_rad = radial/n_rad                          # Element height in the radial direction

            # Determine if the mesh is getting too big. If so the mesh will need to transition to a
            # finer mesh.
            if b_circ > 3*mesh_size:
                transition = True
            else:
                transition = False
        
            # Create a mesh of nodes for the ring
            if transition == True:
                ring = AnnulusTransRingMesh(r_inner + h_rad, r_inner, n_circ, self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod,
                                            self.origin, self.axis, 'N' + str(n), 'Q' + str(q))
                n += 3*n_circ
                q += 4*n_circ
                n_circ *= 3
                self.num_quads_outer = n_circ
            else:
                ring = AnnulusRingMesh(r_inner + h_rad, r_inner, n_circ, self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod, self.origin,
                                       self.axis, 'N' + str(n), 'Q' + str(q))
                n += n_circ
                q += n_circ
        
            # Add the newly generated nodes and elements to the overall mesh. Note that if duplicate
            # keys exist, the `.update()` method will overwrite them with the newly generated key value
            # pairs. This works in our favor by automatically eliminating duplicate nodes from the
            # dictionary.
            self.nodes.update(ring.nodes)
            self.elements.update(ring.elements)

            # Prepare to move to the next ring
            r_inner += h_rad

        # After calling the `.update()` method some elements are still attached to the duplicate
        # nodes that are no longer in the dictionary. Attach these plates to the nodes that are
        # still in the dictionary instead. 
        for element in self.elements.values():
            element.i_node = self.nodes[element.i_node.name]
            element.j_node = self.nodes[element.j_node.name]
            element.m_node = self.nodes[element.m_node.name]
            element.n_node = self.nodes[element.n_node.name]

        # Add the nodes to the model
        for node in self.nodes.values():
            self.model.nodes[node.name] = node
        
        # Add the elements to the model
        for element in self.elements.values():
            if element.type.upper() == 'QUAD':
                self.model.quads[element.name] = element
            elif element.type.upper() == 'RECT':
                self.model.plates[element.name] = element
        
        # Flag the mesh as generated
        self.is_generated = True

#%%
class AnnulusRingMesh(Mesh):
    """
    A mesh of quadrilaterals forming an annular ring (a donut).
    """

    def __init__(self, outer_radius: float, inner_radius: float, num_quads: int, thickness: float, material_name: str, model: FEModel3D, kx_mod: float = 1, ky_mod: float = 1,
                 origin: List[float] = [0, 0, 0], axis: str = 'Y', start_node: str = 'N1', start_element: str = 'Q1') -> None:

        """Single annular ring of quads between two radii.

        :param outer_radius: Outer radius of the ring.
        :type outer_radius: float
        :param inner_radius: Inner radius of the ring.
        :type inner_radius: float
        :param num_quads: Number of quads around the ring.
        :type num_quads: int
        :param thickness: Element thickness.
        :type thickness: float
        :param material_name: Name of the element material.
        :type material_name: str
        :param model: Owning FEModel3D.
        :type model: FEModel3D
        :param kx_mod: In-plane stiffness modifier (local x). Defaults to 1.
        :type kx_mod: float, optional
        :param ky_mod: In-plane stiffness modifier (local y). Defaults to 1.
        :type ky_mod: float, optional
        :param origin: Local origin in global coordinates ``[X, Y, Z]``. Defaults to ``[0, 0, 0]``.
        :type origin: list[float], optional
        :param axis: Global axis: ``'X'``, ``'Y'``, or ``'Z'``. Defaults to ``'Y'``.
        :type axis: str, optional
        :param start_node: First node name. Defaults to ``'N1'``.
        :type start_node: str, optional
        :param start_element: First element name. Defaults to ``'Q1'``.
        :type start_element: str, optional
        """

        super().__init__(thickness, material_name, model, kx_mod, ky_mod, start_node=start_node,
                         start_element=start_element)

        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.n = num_quads
        self.Xo = origin[0]
        self.Yo = origin[1]
        self.Zo = origin[2]

        self.axis = axis

        # Generate the nodes and elements
        self.generate()

    def generate(self) -> None:

        n = self.n  # Number of plates in the initial ring

        inner_radius = self.inner_radius  # The inner radius of the ring
        outer_radius = self.outer_radius  # The outer radius of the ring

        Xo = self.Xo  # Global X-coordinate of the center of the ring
        Yo = self.Yo  # Global Y-coordinate of the center of the ring
        Zo = self.Zo  # Global Z-coordinate of the center of the ring

        axis = self.axis

        theta = 2*pi/self.n  # Angle between nodes in the ring

        # Each node number will be increased by the offset calculated below
        node_offset = int(self.start_node[1:]) - 1

        # Each element number will be increased by the offset calculated below
        element_offset = int(self.start_element[1:]) - 1

        # Generate the nodes that make up the ring, working from the inside to the outside
        angle = 0
        for i in range(1, 2*n + 1, 1):

            # Assign the node a name
            node_name = 'N' + str(i + node_offset)

            # Generate the inner radius of nodes
            if i <= n:
                angle = theta*(i - 1)
                if axis == 'Y':
                    x = Xo + inner_radius*cos(angle)
                    y = Yo
                    z = Zo + inner_radius*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + inner_radius*sin(angle)
                    z = Zo + inner_radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + inner_radius*sin(angle)
                    y = Yo + inner_radius*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusRingMesh.')
            
            # Generate the outer radius of nodes
            else:
                angle = theta*((i - n) - 1)
                if axis == 'Y':
                    x = Xo + outer_radius*cos(angle)
                    y = Yo 
                    z = Zo + outer_radius*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + outer_radius*sin(angle)
                    z = Zo + outer_radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + outer_radius*sin(angle)
                    y = Yo + outer_radius*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusRingMesh.')
            
            self.nodes[node_name] = Node3D(node_name, x, y, z)

        # Generate the elements that make up the ring
        for i in range(1, n + 1, 1):

            # Assign the element a name
            element_name = 'Q' + str(i + element_offset)
            
            n_node = i
            i_node = i + n
            if i != n:
                m_node = i + 1
                j_node = i + 1 + n
            else:
                m_node = 1
                j_node = 1 + n

            self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                               self.nodes['N' + str(j_node + node_offset)],
                                                               self.nodes['N' + str(m_node + node_offset)],
                                                               self.nodes['N' + str(n_node + node_offset)],
                                                               self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod)

        # Add the nodes and elements to the model
        for node in self.nodes.values():
            self.model.nodes[node.name] = node
        
        # Add the elements to the model
        for element in self.elements.values():
            if element.type.upper() == 'QUAD':
                self.model.quads[element.name] = element
            if element.type.upper() == 'RECT':
                self.model.plates[element.name] = element
        
        # Flag the mesh as generated
        self.is_generated = True

#%%
class AnnulusTransRingMesh(Mesh):
    """
    A mesh of quadrilaterals forming an annular ring (a donut) with the mesh getting finer on the outer
    edge.
    """

    def __init__(self, outer_radius: float, inner_radius: float, num_inner_quads: int, thickness: float, material_name: str, model: FEModel3D,
                 kx_mod: float = 1, ky_mod: float = 1, origin: List[float] = [0, 0, 0], axis: str = 'Y', start_node: str = 'N1',
                 start_element: str = 'Q1') -> None:
        """
        Parameters
        ----------

        direction : array
            A vector indicating the direction normal to the ring.
        """

        super().__init__(thickness, material_name, model, kx_mod, ky_mod, start_node=start_node,
                         start_element=start_element)

        self.inner_radius = inner_radius
        self.outer_radius = (inner_radius + outer_radius)/2
        self.r3 = outer_radius
        self.n = num_inner_quads
        self.Xo = origin[0]
        self.Yo = origin[1]
        self.Zo = origin[2]
        self.axis = axis

        # Create the mesh
        self.generate()

    def generate(self) -> None:

        n = self.n  # Number of plates in the outside of the ring (coarse mesh)

        inner_radius = self.inner_radius  # The inner radius of the ring
        outer_radius = self.outer_radius  # The center radius of the ring
        r3 = self.r3  # The outer radius of the ring

        Xo = self.Xo  # Global X-coordinate of the center of the ring
        Yo = self.Yo  # Global Y-coordinate of the center of the ring
        Zo = self.Zo  # Global Z-coordinate of the center of the ring

        axis = self.axis

        theta1 = 2*pi/self.n      # Angle between nodes at the inner radius of the ring
        theta2 = 2*pi/(self.n*3)  # Angle between nodes at the center of the ring
        theta3 = 2*pi/(self.n*3)  # Angle between nodes at the outer radius of the ring

        # Each node number will be increased by the offset calculated below
        node_offset = int(self.start_node[1:]) - 1

        # Each element number will be increased by the offset calculated below
        element_offset = int(self.start_element[1:]) - 1

        # Generate the nodes that make up the ring, working from the inside to the outside
        angle = 0
        for i in range(1, 6*n + 1, 1):

            # Assign the node a name
            node_name = 'N' + str(i + node_offset)

            # Generate the inner radius of nodes
            if i <= n:
                angle = theta1*(i - 1)
                if axis == 'Y':
                    x = Xo + inner_radius*cos(angle)
                    y = Yo
                    z = Zo + inner_radius*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + inner_radius*sin(angle)
                    z = Zo + inner_radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + inner_radius*sin(angle)
                    y = Yo + inner_radius*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusTransRingMesh.')
            
            # Generate the center radius of nodes
            elif i <= 3*n:
                if (i - n) == 1:
                    angle = theta2
                elif (i - n) % 2 == 0:
                    angle += theta2
                else:
                    angle += 2*theta2
                if axis == 'Y':
                    x = Xo + outer_radius*cos(angle)
                    y = Yo 
                    z = Zo + outer_radius*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + outer_radius*sin(angle)
                    z = Zo + outer_radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + outer_radius*sin(angle)
                    y = Yo + outer_radius*cos(angle)
                    z = Zo
            # Generate the outer radius of nodes
            else:
                if (i - 3*n) == 1:
                    angle = 0
                else:
                    angle = theta3*((i - 3*n) - 1)
                if axis == 'Y':
                    x = Xo + r3*cos(angle)
                    y = Yo 
                    z = Zo + r3*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + r3*sin(angle)
                    z = Zo + r3*cos(angle)
                elif axis == 'Z':
                    x = Xo + r3*sin(angle)
                    y = Yo + r3*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusTransRingMesh.')
            
            self.nodes[node_name] = Node3D(node_name, x, y, z)

        # Generate the elements that make up the ring
        for i in range(1, 4*n + 1, 1):

            # Assign the element a name
            element_name = 'Q' + str(i + element_offset)

            if i <= n:
                n_node = i
                j_node = 2*i + n
                i_node = 2*i + n - 1
                if i != n:
                    m_node = i + 1
                else:
                    m_node = 1
            elif (i - n) % 3 == 1:
                n_node = 1 + (i - (n + 1))//3
                m_node = i - (i - (n + 1))//3
                j_node = i + 2*n + 1
                i_node = i + 2*n
            elif (i - n) % 3 == 2:
                n_node = i - 1 - (i - (n + 1))//3
                m_node = i - (i - (n + 1))//3
                j_node = i + 2*n + 1
                i_node = i + 2*n            
            else:
                n_node = i - 1 - (i - (n + 1))//3
                i_node = i + 2*n
                if i != 4*n:
                    m_node = 2 + (i - (n + 1))//3
                    j_node = i + 2*n + 1
                else:
                    m_node = 1
                    j_node = 1 + 3*n

            self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                               self.nodes['N' + str(j_node + node_offset)],
                                                               self.nodes['N' + str(m_node + node_offset)],
                                                               self.nodes['N' + str(n_node + node_offset)],
                                                               self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod)

        # Add the nodes and elements to the model
        for node in self.nodes.values():
            self.model.nodes[node.name] = node
        
        for element in self.elements.values():
            if element.type == 'Quad':
                self.model.quads[element.name] = element
            else:
                self.model.plates[element.name] = element
        
        # Flag the mesh as generated
        self.is_generated = True

#%%
class FrustrumMesh(AnnulusMesh):
    """
    A mesh of quadrilaterals forming a frustrum (a cone intersected by a horizontal plane).
    """

    def __init__(self, mesh_size: float, large_radius: float, small_radius: float, height: float, thickness: float, material_name: str, model: FEModel3D, kx_mod: float = 1, ky_mod: float = 1,
                 origin: List[float] = [0, 0, 0], axis: str = 'Y', start_node: str = 'N1', start_element: str = 'Q1') -> None:
        """Conical frustum mesh generated from an annulus and then tapered to height.

        :param mesh_size: Target element size for the base annulus.
        :type mesh_size: float
        :param large_radius: Large radius (base) of the frustum.
        :type large_radius: float
        :param small_radius: Small radius (top) of the frustum.
        :type small_radius: float
        :param height: Frustum height along ``axis``.
        :type height: float
        :param thickness: Element thickness.
        :type thickness: float
        :param material_name: Name of the element material.
        :type material_name: str
        :param model: Owning FEModel3D.
        :type model: FEModel3D
        :param kx_mod: In-plane stiffness modifier (local x). Defaults to 1.
        :type kx_mod: float, optional
        :param ky_mod: In-plane stiffness modifier (local y). Defaults to 1.
        :type ky_mod: float, optional
        :param origin: Local origin in global coordinates ``[X, Y, Z]``. Defaults to ``[0, 0, 0]``.
        :type origin: list[float], optional
        :param axis: Global axis: ``'X'``, ``'Y'``, or ``'Z'``. Defaults to ``'Y'``.
        :type axis: str, optional
        :param start_node: First node name. Defaults to ``'N1'``.
        :type start_node: str, optional
        :param start_element: First element name. Defaults to ``'Q1'``.
        :type start_element: str, optional
        """

        # Create an annulus mesh
        super().__init__(mesh_size, large_radius, small_radius, thickness, material_name, model, kx_mod,
                         ky_mod, origin, axis, start_node, start_element)
        
        self.height = height
    
    def generate(self) -> None:

        super().generate()
        
        Xo = self.origin[0]
        Yo = self.origin[1]
        Zo = self.origin[2]

        # Adjust the coordinates of each node to make a frustrum
        for node in self.nodes.values():

            X = node.X
            Y = node.Y
            Z = node.Z
            
            if self.axis == 'Y':
                # Calculate the radius to the node
                r = ((X - Xo)**2 + (Z - Zo)**2)**0.5
                # Adjust the height of the node
                node.Y += (r - self.outer_radius)/(self.outer_radius - self.inner_radius)*self.height
            elif self.axis == 'X':
                # Calculate the radius to the node
                r = ((Y - Yo)**2 + (Z - Zo)**2)**0.5
                # Adjust the height of the node
                node.X += (r - self.outer_radius)/(self.outer_radius - self.inner_radius)*self.height
            elif self.axis == 'Z':
                # Calculate the radius to the node
                r = ((X - Xo)**2 + (Y - Yo)**2)**0.5
                # Adjust the height of the node
                node.Z += (r - self.outer_radius)/(self.outer_radius - self.inner_radius)*self.height
            else:
                raise Exception('Invalid axis specified for frustrum mesh.')

#%%
class CylinderMesh(Mesh):

    def __init__(self, mesh_size: float, radius: float, height: float, thickness: float, material_name: str, model: FEModel3D, kx_mod: float = 1, ky_mod: float = 1,origin: List[float] = [0, 0, 0], axis: str = 'Y', start_node: str = 'N1', start_element: str = 'Q1', num_elements: int | None = None, element_type: str = 'Quad') -> None:

        """Cylindrical shell mesh.

        :param mesh_size: Target element edge size. If ``num_elements`` is provided, used only for vertical division; otherwise used to compute circumferential division too.
        :type mesh_size: float
        :param radius: Cylinder radius to element centers.
        :type radius: float
        :param height: Total height of the cylinder.
        :type height: float
        :param thickness: Element thickness.
        :type thickness: float
        :param material_name: Name of the element material.
        :type material_name: str
        :param model: Owning FEModel3D.
        :type model: FEModel3D
        :param kx_mod: In-plane stiffness modifier (local x). Defaults to 1.
        :type kx_mod: float, optional
        :param ky_mod: In-plane stiffness modifier (local y). Defaults to 1.
        :type ky_mod: float, optional
        :param origin: Local origin in global coordinates ``[X, Y, Z]``. Defaults to ``[0, 0, 0]``.
        :type origin: list[float], optional
        :param axis: Global axis about which the mesh is generated: ``'X'``, ``'Y'``, or ``'Z'``. Defaults to ``'Y'``.
        :type axis: str, optional
        :param start_node: First node name (e.g., ``'N1'``). Defaults to ``'N1'``.
        :type start_node: str, optional
        :param start_element: First element name (e.g., ``'Q1'``). Defaults to ``'Q1'``.
        :type start_element: str, optional
        :param num_elements: Optional number of quads around the circumference (overrides automatic calculation from ``mesh_size``).
        :type num_elements: int | None, optional
        :param element_type: Element family for the mesh: ``'Quad'`` or ``'Rect'``. Defaults to ``'Quad'``.
        :type element_type: str, optional
        """

        # Inherit properties and methods from the parent `Mesh` class
        super().__init__(thickness, material_name, model, kx_mod, ky_mod, start_node, start_element)

        # Define a few new additional class properties related to cylinders
        self.radius = radius
        self.h = height
        self.mesh_size = mesh_size
        self.origin = origin
        self.axis = axis

        # Check if the user has requested a specific number of elements for each course of plates. This can be useful for ensuring the mesh matches up with other meshes.
        if num_elements == None:
            # Calculate the number of elements if the user hasn't specified
            self.num_elements = int(round(2*pi*radius/mesh_size, 0))
        else:
            # Use the user specified number of elements
            self.num_elements = num_elements

        # Check which type of element the user has requested (rectangular plate or quad)
        self.element_type = element_type

        # Generate the mesh
        self.generate()
    
    def generate(self) -> None:
        
        # Get the mesh thickness and the material name
        thickness = self.thickness
        material_name = self.material_name

        mesh_size = self.mesh_size  # Desired mesh size
        num_elements = self.num_elements  # Number of quadrilaterals in each course of the ring
        n = self.num_elements  # Total number of elements in the mesh (initialized for a single ring at the moment)

        radius = self.radius
        h = self.h

        # Set the cylinder base's local y-coordinate
        if self.axis == 'Y':
            y = self.origin[1]  
        elif self.axis == 'X':
            y = self.origin[0]  
        elif self.axis == 'Z':
             y = self.origin[2]

        n = int(self.start_node[1:])
        q = int(self.start_element[1:])

        element_type = self.element_type

        # Determine the number of quads to mesh the circumference into
        if num_elements == None:
            num_elements = int(2*pi/mesh_size)

        # Mesh the cylinder from the bottom toward the top
        while round(y, 10) < round(h, 10):

            height = h - y  #Remaining height to be meshed
            # Number of times the plate height fits in the remaining unmeshed height, resulting at least one element
            n_vert = max(int(abs(height)/mesh_size),1)
            h_y = height/n_vert             # Element height in the vertical direction
            # Create a mesh of nodes for the ring
            if self.axis == 'Y':
                ring = CylinderRingMesh(radius, h_y, num_elements, thickness, material_name, self.model, 1, 1, [0, y, 0], self.axis, 'N' + str(n), 'Q' + str(q), element_type)
            elif self.axis == 'X':
                ring = CylinderRingMesh(radius, h_y, num_elements, thickness, material_name, self.model, 1, 1, [y, 0, 0], self.axis, 'N' + str(n), 'Q' + str(q), element_type)
            elif self.axis == 'Z':
                ring = CylinderRingMesh(radius, h_y, num_elements, thickness, material_name, self.model, 1, 1, [0, 0, y], self.axis, 'N' + str(n), 'Q' + str(q), element_type)

            n += num_elements
            q += num_elements
        
            # Add the newly generated nodes and elements to the overall mesh. Note that if duplicate keys exist, the `.update()` method will overwrite them with the newly generated key value pairs. This works in our favor by automatically eliminating duplicate nodes at the shared boundaries between rings.
            self.nodes.update(ring.nodes)
            self.elements.update(ring.elements)

            # Prepare to move to the next ring
            y += h_y
        
        # After calling the `.update()` method some elements are still attached to the duplicate nodes that are no longer in the dictionary. Attach these plates to the nodes that are still in the dictionary instead. 
        for element in self.elements.values():
            element.i_node = self.nodes[element.i_node.name]
            element.j_node = self.nodes[element.j_node.name]
            element.m_node = self.nodes[element.m_node.name]
            element.n_node = self.nodes[element.n_node.name]

        # Add the nodes and elements to the model
        for node in self.nodes.values():
            self.model.nodes[node.name] = node
        
        for element in self.elements.values():
            if element.type == 'Quad':
                self.model.quads[element.name] = element
            else:
                self.model.plates[element.name] = element
        
        # Flag the mesh as generated
        self.is_generated = True

#%%
class CylinderRingMesh(Mesh):
    """
    A mesh of quadrilaterals forming a cylindrical ring.

    Parameters
    ----------

    radius : number
        Radius to the center of the plates in the cylindrical ring.
    height : number
        Height of the cylindrical ring.
    num_elements : number
        Number of elements used to generate the cylindrical ring.
    thickness : number
        Element thickness.
    material_name : string
        The name of the element material.
    kx_mod : number
        Stiffness modification factor for in-plane stiffness in the element's local
        x-direction. Default value is 1.0 (no modification).
    ky_mod : number
        Stiffness modification factor for in-plane stiffness in the element's local
        y-direction. Default value is 1.0 (no modification).
    origin : list
        The location of the center of the base of the cylindrical ring. Default is [0, 0, 0].
    axis : string
        Global axis about which to revolve the ring ('X', 'Y', or 'Z'). Default is 'Y'.
    start_node : string, optional
        The name of the first node in the mesh. The name must be formatted starting with a single
        letter followed by a number (e.g. 'N12'). The mesh will begin numbering nodes from this
        number. The default is 'N1'. 
    start_element : string, optional
        The name of the first element in the mesh. The name must be formatted starting with a
        single letter followed by a number (e.g. 'Q32'). The mesh will begin numbering elements
        from this number. The default is 'Q1'.
    num_elements : number
        The number of elements to divide the circumference into.
    
    """

    def __init__(self, radius: float, height: float, num_elements: int, thickness: float, material_name: str, model: FEModel3D, kx_mod: float = 1, ky_mod: float = 1,
                 origin: List[float] = [0, 0, 0], axis: str = 'Y', start_node: str = 'N1', start_element: str = 'Q1',
                 element_type: str = 'Quad') -> None:

        super().__init__(thickness, material_name, model, kx_mod, ky_mod, start_node=start_node, start_element=start_element)

        self.radius = radius
        self.height = height
        self.num_elements = num_elements
        self.Xo = origin[0]
        self.Yo = origin[1]
        self.Zo = origin[2]
        self.axis = axis
        self.element_type = element_type

        # Generate the nodes and elements
        self.generate()

    def generate(self) -> None:
        """
        Generates the nodes and elements in the mesh.
        """

        num_elements = self.num_elements  # Number of quadrilaterals in the ring
        n = self.num_elements

        radius = self.radius  # The radius of the ring
        height = self.height  # The height of the ring

        Xo = self.Xo  # Global X-coordinate of the center of the bottom of the ring
        Yo = self.Yo  # Global Y-coordinate of the center of the bottom of the ring
        Zo = self.Zo  # Global Z-coordinate of the center of the bottom of the ring

        axis = self.axis
        
        # Calculate the angle between nodes in the circumference of the ring
        theta = 2*pi/num_elements

        # Each node number will be increased by the offset calculated below
        try:
            node_offset = int(self.start_node[1:]) - 1
        except:
            raise ValueError('Invalid node name. Enter a letter followed by a number (e.g. \'N25\')')

        # Each element number will be increased by the offset calculated below
        try:
            element_offset = int(self.start_element[1:]) - 1
        except:
            raise ValueError('Invalid element ame. Enter a letter followed by a number (e.g. \'Q83\')')

        # Generate the nodes that make up the ring
        angle = 0
        for i in range(1, 2*n + 1, 1):

            # Assign the node a name
            node_name = 'N' + str(i + node_offset)

            # Generate the bottom nodes of the ring
            if i <= n:
                angle = theta*(i - 1)
                if axis == 'Y':
                    x = Xo + radius*cos(angle)
                    y = Yo
                    z = Zo + radius*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + radius*sin(angle)
                    z = Zo + radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + radius*sin(angle)
                    y = Yo + radius*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for CylinderRingMesh.')

            # Generate the top nodes of the ring
            else:
                angle = theta*((i - n) - 1)
                if axis == 'Y':
                    x = Xo + radius*cos(angle)
                    y = Yo + height
                    z = Zo + radius*sin(angle)
                elif axis == 'X':
                    x = Xo + height
                    y = Yo + radius*sin(angle)
                    z = Zo + radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + radius*sin(angle)
                    y = Yo + radius*cos(angle)
                    z = Zo + height
                else:
                    raise Exception('Invalid axis specified for CylinderRingMesh.')
            
            self.nodes[node_name] = Node3D(node_name, x, y, z)

        # Generate the elements that make up the ring
        for i in range(1, n + 1, 1):

            # Assign the element a name
            if self.element_type == 'Quad':
                element_name = 'Q' + str(i + element_offset)
            elif self.element_type == 'Rect':
                element_name = 'R' + str(i + element_offset)
            else:
                raise Exception('Invalid element type specified for cylinder ring mesh.')
            
            # Assign nodes to the element
            n_node = i
            i_node = i + n
            if i != n:
                m_node = i + 1
                j_node = i + 1 + n
            else:
                m_node = 1
                j_node = 1 + n

            # Create the element and add it to the `elements` dictionary
            if self.element_type == 'Quad':
                self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                     self.nodes['N' + str(j_node + node_offset)],
                                                     self.nodes['N' + str(m_node + node_offset)],
                                                     self.nodes['N' + str(n_node + node_offset)],
                                                     self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod)
            elif self.element_type == 'Rect':
                self.elements[element_name] = Plate3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                      self.nodes['N' + str(j_node + node_offset)],
                                                      self.nodes['N' + str(m_node + node_offset)],
                                                      self.nodes['N' + str(n_node + node_offset)],
                                                      self.thickness, self.material_name, self.model, self.kx_mod, self.ky_mod)
        
        # Add the nodes and elements to the model
        for node in self.nodes.values():
            self.model.nodes[node.name] = node
        
        for element in self.elements.values():
            if element.type == 'Quad':
                self.model.quads[element.name] = element
            else:
                self.model.plates[element.name] = element
            
        # Flag the mesh as generated
        self.is_generated = True
        
def check_mesh_integrity(mesh: Mesh, console_log: bool = True) -> Union[str, List[str], None]:
    """Runs basic integrity checks to ensure the mesh is in sync with its model. Usually you don't
    want to run this check unless the mesh has been generated since generating the mesh is what
    syncs it to the model.

    :param mesh: A mesh of finite elements.
    :type mesh: Mesh
    :param console_log: Determines whether the results will be printed to the console or returned
                        as a list of strings. Defaults to True.
    :type console_log: bool, optional
    :return: Any errors discovered by the integrity check
    :rtype: list
    """

    # Initialize an empty list for error messages
    errors = []

    # Check that every element in the mesh is attached at all four nodes to nodes that are in the mesh
    count = 0
    for element in mesh.elements.values():

        if (element.i_node not in mesh.nodes.values() or
            element.j_node not in mesh.nodes.values() or
            element.m_node not in mesh.nodes.values() or
            element.n_node not in mesh.nodes.values()):

            count += 1
    
    # Prepare the error message
    if count != 0:
        errors.append(str(count) + ' elements are attached to nodes that are not in the mesh.')

    # Check that every element in the mesh is attached at all four nodes to nodes that are in the model
    count = 0
    for element in mesh.elements.values():

        if (element.i_node not in mesh.model.nodes.values() or
            element.j_node not in mesh.model.nodes.values() or
            element.m_node not in mesh.model.nodes.values() or
            element.n_node not in mesh.model.nodes.values()):

            count += 1
    
    # Prepare the error message
    if count != 0:
        errors.append(str(count) + ' elements are attached to nodes that are not in the model.')
    
    # Check that each element has 4 unique nodes
    count = 0
    for element in mesh.elements.values():

        if (element.i_node is element.j_node or
            element.i_node is element.m_node or
            element.i_node is element.n_node or
            element.j_node is element.m_node or
            element.j_node is element.n_node or
            element.m_node is element.n_node):

            count += 1
        
    # Prepare the error message
    if count != 0:
        errors.append(str(count) + ' elements have duplicate nodes in them.')

    # Check that each element's name in the mesh is in the model
    count = 0
    for element in mesh.elements.values():
        if element.name not in mesh.model.plates.keys() and element.name not in mesh.model.quads.keys():
            count += 1
    
    if count != 0:
        errors.append(str(count) + ' element names in the mesh were not found in the model.')
    
    # Check that each element in the mesh is in the model
    count = 0
    for element in mesh.elements.values():
        if element not in mesh.model.plates.values() and element.name not in mesh.model.quads.values():
            count += 1
    
    if count != 0:
        errors.append(str(count) + ' elements in the mesh were not found in the model.')

    # Check that each element in the mesh matches its corresponding element in the model
    count = 0
    for element in mesh.elements.values():
        
        if mesh.element_type == 'Rect':
            if element.name not in mesh.model.plates.keys() or mesh.elements[element.name] is not mesh.model.plates[element.name]:
                count +=1
        elif mesh.element_type == 'Quad':
            if element.name not in mesh.model.quads.keys() or mesh.elements[element.name] is not mesh.model.quads[element.name]:
                count += 1
    
    # Prepare the error message
    if count != 0:
        errors.append(str(count) + ' elements in the mesh do not match their respective elements in the model.')
    
    # Check that the mesh's key for each node matches the node's name
    count = 0
    for node in mesh.nodes.values():
        if node.name not in mesh.nodes.keys():
            count += 1
    
    if count != 0:
        errors.append(str(count) + ' node names in the mesh do not match the mesh\'s `nodes` dictionary\'s keys.')

    # # Check that each node in the mesh matches its corresponding node in the model
    # count = 0
    # for node in mesh.nodes.values():

    #     if node.name not in mesh.model.nodes.keys() or mesh.nodes[node.name] is not mesh.model.nodes[node.name]:
    #         count += 1
    
    # # Prepare the error message
    # if count != 0:
    #     errors.append(str(count) + ' nodes in the mesh do not match their respective nodes in the model.')

    # TODO: Add more integrity checks and error messages

    # Report if no errors were found
    if errors == []:
        errors = ['No errors detected.']

    # Return the error messages
    if console_log == True:
        for error in errors:
            print(error)
            return ''
    else:
        return errors


