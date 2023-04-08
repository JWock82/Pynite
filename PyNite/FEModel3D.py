# %%
from os import rename
import warnings
from matplotlib.pyplot import get

from numpy import array, zeros, matmul, divide, subtract, atleast_2d, nanmax
from numpy import seterr
from numpy.linalg import solve

from PyNite.Node3D import Node3D
from PyNite.Material import Material
from PyNite.PhysMember import PhysMember
from PyNite.Spring3D import Spring3D
from PyNite.Member3D import Member3D
from PyNite.Quad3D import Quad3D
from PyNite.Plate3D import Plate3D
from PyNite.LoadCombo import LoadCombo
from PyNite.Mesh import Mesh, RectangleMesh, AnnulusMesh, FrustrumMesh, CylinderMesh

# %%
class FEModel3D():
    """A 3D finite element model.
    """

    def __init__(self):
        """Creates a new 3D finite element model.
        """
        
        # Initialize the model's various dictionaries. The dictionaries will be prepopulated with
        # the data types they store, and then those types will be removed. This will give us the
        # ability to get type-based hints when using the dictionaries.

        self.Nodes = {str:Node3D}          # A dictionary of the model's nodes
        self.Nodes.pop(str)
        self.AuxNodes = {str:Node3D}       # A dictionary of the model's auxiliary nodes
        self.AuxNodes.pop(str)
        self.Materials = {str:Material}    # A dictionary of the model's materials
        self.Materials.pop(str)
        self.Springs = {str:Spring3D}      # A dictionary of the model's springs
        self.Springs.pop(str)
        self.Members = {str:PhysMember}    # A dictionary of the model's physical members
        self.Members.pop(str)
        self.Quads = {str:Quad3D}          # A dictionary of the model's quadiralterals
        self.Quads.pop(str)
        self.Plates = {str:Plate3D}        # A dictionary of the model's rectangular plates
        self.Plates.pop(str)
        self.Meshes = {str:Mesh}           # A dictionary of the model's meshes
        self.Meshes.pop(str)         
        self.LoadCombos = {str:LoadCombo}  # A dictionary of the model's load combinations
        self.LoadCombos.pop(str)
        self._D = {str:[]}                 # A dictionary of the model's nodal displacements by load combination
        self._D.pop(str)
        
        self.solution = None  # Indicates the solution type for the latest run of the model

    @property
    def LoadCases(self):
        """Returns a list of all the load cases in the model (in alphabetical order).
        """
        
        # Create an empty list of load cases
        cases = []

        # Step through each node
        for node in self.Nodes.values():
            # Step through each nodal load
            for load in node.NodeLoads:
                # Get the load case for each nodal laod
                cases.append(load[2])
        
        # Step through each member
        for member in self.Members.values():
            # Step through each member point load
            for load in member.PtLoads:
                # Get the load case for each member point load
                cases.append(load[3])
            # Step through each member distributed load
            for load in member.DistLoads:
                # Get the load case for each member distributed load
                cases.append(load[5])
        
        # Step through each plate/quad
        for plate in list(self.Plates.values()) + list(self.Quads.values()):
            # Step through each surface load
            for load in plate.pressures:
                # Get the load case for each plate/quad pressure
                cases.append(load[1])

        # Remove duplicates and return the list (sorted ascending)
        return sorted(list(dict.fromkeys(cases)))

    def add_node(self, name, X, Y, Z):
        """Adds a new node to the model.

        :param name: A unique user-defined name for the node. If set to None or "" a name will be
                     automatically assigned.
        :type name: str
        :param X: The node's global X-coordinate.
        :type X: number
        :param Y: The node's global Y-coordinate.
        :type Y: number
        :param Z: The node's global Z-coordinate.
        :type Z: number
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the node added to the model.
        :rtype: str
        """
        
        # Name the node or check it doesn't already exist
        if name:
            if name in self.Nodes:
                raise NameError(f"Node name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "N" + str(len(self.Nodes))
            count = 1
            while name in self.Nodes: 
                name = "N" + str(len(self.Nodes) + count)
                count += 1
        
        # Create a new node
        new_node = Node3D(name, X, Y, Z)
        
        # Add the new node to the list
        self.Nodes[name] = new_node
        
        # Flag the model as unsolved
        self.solution = None

        #Return the node name
        return name

    def add_auxnode(self, name, X, Y, Z):
        """Adds a new auxiliary node to the model. Together with a member's `i` and `j` nodes, an
        auxiliary node defines the plane in which the member's local z-axis lies, and the side of
        the member the z-axis points toward. If no auxiliary node is specified for a member, PyNite
        uses its own default configuration.

        :param name: A unique user-defined name for the node. If None or "", a name will be
                     automatically assigned.
        :type name: str
        :param X: The global X-coordinate of the node.
        :type X: number
        :param Y: The global Y-coordinate of the node.
        :type Y: number
        :param Z: The global Z-coordinate of the node.
        :type Z: number
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the auxiliary node that was added to the model.
        :rtype: str
        """ 
        
        # Name the node or check it doesn't already exist
        if name:
            if name in self.AuxNodes:
                raise NameError(f"Auxnode name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "AN" + str(len(self.AuxNodes))
            count = 1
            while name in self.AuxNodes: 
                name = "AN" + str(len(self.AuxNodes)+count)
                count += 1
                
        # Create a new node
        new_node = Node3D(name, X, Y, Z)
        
        # Add the new node to the list
        self.AuxNodes[name] = new_node

        # Flag the model as unsolved
        self.solution = None
        
        #Return the node name
        return name 

    def add_material(self, name, E, G, nu, rho):
        """Adds a new material to the model.

        :param name: A unique user-defined name for the material.
        :type name: str
        :param E: The modulus of elasticity of the material.
        :type E: number
        :param G: The shear modulus of elasticity of the material.
        :type G: number
        :param nu: Poisson's ratio of the material.
        :type nu: number
        :param rho: The density of the material
        :type rho: number
        :raises NameError: Occurs when the specified name already exists in the model.
        """

        # Name the material or check it doesn't already exist
        if name:
            if name in self.Materials:
                raise NameError(f"Material name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "M" + str(len(self.Materials))
            count = 1
            while name in self.Materials: 
                name = "M" + str(len(self.Materials) + count)
                count += 1
                
        # Create a new material
        new_material = Material(name, E, G, nu, rho)
        
        # Add the new material to the list
        self.Materials[name] = new_material
        
        # Flag the model as unsolved
        self.solution = None

    def add_spring(self, name, i_node, j_node, ks, tension_only=False, comp_only=False):
        """Adds a new spring to the model.

        :param name: A unique user-defined name for the member. If None or "", a name will be
                    automatically assigned
        :type name: str
        :param i_node: The name of the i-node (start node).
        :type i_node: str
        :param j_node: The name of the j-node (end node).
        :type j_node: str
        :param ks: The spring constant (force/displacement).
        :type ks: number
        :param tension_only: Indicates if the member is tension-only, defaults to False
        :type tension_only: bool, optional
        :param comp_only: Indicates if the member is compression-only, defaults to False
        :type comp_only: bool, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the spring that was added to the model.
        :rtype: str
        """ 
        
        # Name the spring or check it doesn't already exist
        if name:
            if name in self.Springs:
                raise NameError(f"Spring name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "S" + str(len(self.Springs))
            count = 1
            while name in self.Springs: 
                name = "S" + str(len(self.Springs) + count)
                count += 1
                
        # Create a new spring
        new_spring = Spring3D(name, self.Nodes[i_node], self.Nodes[j_node],
                              ks, self.LoadCombos, tension_only=tension_only,
                              comp_only=comp_only)
        
        # Add the new spring to the list
        self.Springs[name] = new_spring
        
        # Flag the model as unsolved
        self.solution = None

        # Return the spring name
        return name

    def add_member(self, name, i_node, j_node, material, Iy, Iz, J, A, auxNode=None,
                   tension_only=False, comp_only=False):
        """Adds a new physical member to the model.

        :param name: A unique user-defined name for the member. If None or "", a name will be
                    automatically assigned
        :type name: str
        :param i_node: The name of the i-node (start node).
        :type i_node: str
        :param j_node: The name of the j-node (end node).
        :type j_node: str
        :param material: The name of the material of the member.
        :type material: str
        :param Iy: The moment of inertia of the member about its local y-axis.
        :type Iy: number
        :param Iz: The moment of inertia of the member about its local z-axis.
        :type Iz: number
        :param J: The polar moment of inertia of the member.
        :type J: number
        :param A: The cross-sectional area of the member.
        :type A: number
        :param auxNode: The name of the auxiliary node used to define the local z-axis. The default
                        is None, in which case the program defines the axis instead of using an
                        auxiliary node.
        :type auxNode: str, optional
        :param tension_only: Indicates if the member is tension-only, defaults to False
        :type tension_only: bool, optional
        :param comp_only: Indicates if the member is compression-only, defaults to False
        :type comp_only: bool, optional
        :raises NameError: Occurs if the specified name already exists.
        :return: The name of the member added to the model.
        :rtype: str
        """

        # Name the member or check it doesn't already exist
        if name:
            if name in self.Members: raise NameError(f"Member name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "M" + str(len(self.Members))
            count = 1
            while name in self.Members: 
                name = "M" + str(len(self.Members)+count)
                count += 1
                
        # Create a new member
        if auxNode == None:
            new_member = PhysMember(name, self.Nodes[i_node], self.Nodes[j_node], material, Iy, Iz, J,
                                    A, model=self, tension_only=tension_only, comp_only=comp_only)
        else:
            new_member = PhysMember(name, self.Nodes[i_node], self.Nodes[j_node], material, Iy, Iz, J,
                                    A, model=self, auxnode=self.GetAuxNode(auxNode),
                                    tension_only=tension_only, comp_only=comp_only)
        
        # Add the new member to the list
        self.Members[name] = new_member
        
        # Flag the model as unsolved
        self.solution = None

        # Return the member name
        return name

    def add_plate(self, name, i_node, j_node, m_node, n_node, t, material, kx_mod=1.0, ky_mod=1.0):
        """Adds a new rectangular plate to the model. The plate formulation for in-plane (membrane)
        stiffness is based on an isoparametric formulation. For bending, it is based on a 12-term
        polynomial formulation. This element must be rectangular, and must not be used where a
        thick plate formulation is needed. For a more versatile plate element that can handle
        distortion and thick plate conditions, consider using the `add_quad` method instead.

        :param name: A unique user-defined name for the plate. If None or "", a name will be
                     automatically assigned.
        :type name: str
        :param i_node: The name of the i-node.
        :type i_node: str
        :param j_node: The name of the j-node.
        :type j_node: str
        :param m_node: The name of the m-node.
        :type m_node: str
        :param n_node: The name of the n-node.
        :type n_node: str
        :param t: The thickness of the element.
        :type t: number
        :param material: The name of the material for the element.
        :type material: str
        :param kx_mod: Stiffness modification factor for in-plane stiffness in the element's local
                       x-direction, defaults to 1 (no modification).
        :type kx_mod: number, optional
        :param ky_mod: Stiffness modification factor for in-plane stiffness in the element's local
                       y-direction, defaults to 1 (no modification).
        :type ky_mod: number, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the element added to the model.
        :rtype: str
        """
        
        # Name the plate or check it doesn't already exist
        if name:
            if name in self.Plates: raise NameError(f"Plate name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "P" + str(len(self.Plates))
            count = 1
            while name in self.Plates: 
                name = "P" + str(len(self.Plates)+count)
                count += 1
        
        # Create a new plate
        new_plate = Plate3D(name, self.Nodes[i_node], self.Nodes[j_node], self.Nodes[m_node],
                           self.Nodes[n_node], t, material, self, kx_mod, ky_mod)
        
        # Add the new plate to the list
        self.Plates[name] = new_plate

        # Flag the model as unsolved
        self.solution = None
        
        # Return the plate name
        return name

    def add_quad(self, name, i_node, j_node, m_node, n_node, t, material, kx_mod=1.0, ky_mod=1.0):
        """Adds a new quadrilateral to the model. The quad formulation for in-plane (membrane)
        stiffness is based on an isoparametric formulation. For bending, it is based on an MITC4
        formulation. This element handles distortion relatively well, and is appropriate for thick
        and thin plates. One limitation with this element is that it does a poor job of reporting
        corner stresses. Corner forces, however are very accurate. Center stresses are very
        accurate as well. For cases where corner stress results are important, consider using the
        `add_plate` method instead.

        :param name: A unique user-defined name for the quadrilateral. If None or "", a name will
                     be automatically assigned.
        :type name: str
        :param i_node: The name of the i-node.
        :type i_node: str
        :param j_node: The name of the j-node.
        :type j_node: str
        :param m_node: The name of the m-node.
        :type m_node: str
        :param n_node: The name of the n-node.
        :type n_node: str
        :param t: The thickness of the element.
        :type t: number
        :param material: The name of the material for the element.
        :type material: str
        :param kx_mod: Stiffness modification factor for in-plane stiffness in the element's local
            x-direction, defaults to 1 (no modification).
        :type kx_mod: number, optional
        :param ky_mod: Stiffness modification factor for in-plane stiffness in the element's local
            y-direction, defaults to 1 (no modification).
        :type ky_mod: number, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the element added to the model.
        :rtype: str
        """
        
        # Name the quad or check it doesn't already exist
        if name:
            if name in self.Quads: raise NameError(f"Quad name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "Q" + str(len(self.Quads))
            count = 1
            while name in self.Quads: 
                name = "Q" + str(len(self.Quads) + count)
                count += 1
        
        # Create a new member
        new_quad = Quad3D(name, self.Nodes[i_node], self.Nodes[j_node], self.Nodes[m_node],
                         self.Nodes[n_node], t, material, self, kx_mod, ky_mod)
        
        # Add the new member to the list
        self.Quads[name] = new_quad

        # Flag the model as unsolved
        self.solution = None
        
        #Return the quad name
        return name

    def add_rectangle_mesh(self, name, mesh_size, width, height, thickness, material, kx_mod=1.0, 
            ky_mod=1.0, origin=[0, 0, 0], plane='XY', x_control=None, y_control=None, start_node=None,
            start_element = None, element_type='Quad'):
        """Adds a rectangular mesh of elements to the model.

        :param name: A unique name for the mesh.
        :type name: str
        :param mesh_size: The desired mesh size.
        :type mesh_size: number
        :param width: The overall width of the rectangular mesh measured along its local x-axis.
        :type width: number
        :param height: The overall height of the rectangular mesh measured along its local y-axis.
        :type height: number
        :param thickness: The thickness of each element in the mesh.
        :type thickness: number
        :param material: The name of the material for elements in the mesh.
        :type material: str
        :param kx_mod: Stiffness modification factor for in-plane stiffness in the element's local
                       x-direction. Defaults to 1.0 (no modification).
        :type kx_mod: float, optional
        :param ky_mod: Stiffness modification factor for in-plane stiffness in the element's local
                       y-direction. Defaults to 1.0 (no modification).
        :type ky_mod: float, optional
        :param origin: The origin of the regtangular mesh's local coordinate system. Defaults to [0, 0, 0]
        :type origin: list, optional
        :param plane: The plane the mesh will be parallel to. Options are 'XY', 'YZ', and 'XZ'.
                      Defaults to 'XY'.
        :type plane: str, optional
        :param x_control: A list of control points along the mesh's local x-axis to work into the
                          mesh. Defaults to `None`.
        :type x_control: list, optional
        :param y_control: A list of control points along the mesh's local y-axis to work into the
                          mesh. Defaults to None.
        :type y_control: list, optional
        :param start_node: The name of the first node in the mesh. If set to `None` the program
                           will use the next available node name. Default is `None`.
        :type start_node: str, optional
        :param start_element: The name of the first element in the mesh. If set to `None` the
                              program will use the next available element name. Default is `None`.
        :type start_element: str, optional
        :param element_type: They type of element to make the mesh out of. Either 'Quad' or 'Rect'.
                             Defaults to 'Quad'.
        :type element_type: str, optional
        :raises NameError: Occurs when the specified name already exists in the model.
        :return: The name of the mesh added to the model.
        :rtype: str
        """
        
        # Check if a mesh name has been provided
        if name:
            # Check that the mesh name isn't already being used
            if name in self.Meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Rename the mesh if necessary
        else:
            name = self.unique_name(self.Meshes, 'MSH')
        
        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.Nodes, 'N')
        if element_type == 'Rect' and start_element is None:
            start_element = self.unique_name(self.Plates, 'R')
        elif element_type == 'Quad' and start_element is None:
            start_element = self.unique_name(self.Quads, 'Q')
        
        # Create the mesh
        new_mesh = RectangleMesh(mesh_size, width, height, thickness, material, self, kx_mod,
                                 ky_mod, origin, plane, x_control, y_control, start_node,
                                 start_element, element_type=element_type)

        # Add the new mesh to the `Meshes` dictionary
        self.Meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name
    
    def add_annulus_mesh(self, name, mesh_size, outer_radius, inner_radius, thickness, material, kx_mod=1.0, 
            ky_mod=1.0, origin=[0, 0, 0], axis='Y', start_node=None, start_element=None):
        """Adds a mesh of quadrilaterals forming an annulus (a donut).

        :param name: A unique name for the mesh.
        :type name: str
        :param mesh_size: The target mesh size.
        :type mesh_size: float
        :param outer_radius: The radius to the outside of the annulus.
        :type outer_radius: float
        :param inner_radius: The radius to the inside of the annulus.
        :type inner_radius: float
        :param thickness: Element thickness.
        :type thickness: float
        :param material: The name of the element material.
        :type material: str
        :param kx_mod: Stiffness modification factor for radial stiffness in the element's local
                       x-direction. Default is 1.0 (no modification).
        :type kx_mod: float, optional
        :param ky_mod: Stiffness modification factor for meridional stiffness in the element's
                       local y-direction. Default is 1.0 (no modification).
        :type ky_mod: float, optional
        :param origin: The origin of the mesh. The default is [0, 0, 0].
        :type origin: list, optional
        :param axis: The global axis about which the mesh will be generated. The default is 'Y'.
        :type axis: str, optional
        :param start_node: The name of the first node in the mesh. If set to `None` the program
                           will use the next available node name. Default is `None`.
        :type start_node: str, optional
        :param start_element: The name of the first element in the mesh. If set to `None` the
                              program will use the next available element name. Default is `None`.
        :type start_element: str, optional
        :raises NameError: Occurs if the specified name already exists in the model.
        :return: The name of the mesh added to the model.
        :rtype: str
        """
        
        # Check if a mesh name has been provided
        if name:
            # Check that the mesh name doesn't already exist
            if name in self.Meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Give the mesh a new name if necessary
        else:
            name = self.unique_name(self.Meshes, 'MSH')

        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.Nodes, 'N')
        if start_element is None:
            start_element = self.unique_name(self.Quads, 'Q')
        
        # Create a new mesh
        new_mesh = AnnulusMesh(mesh_size, outer_radius, inner_radius, thickness, material, self,
                               kx_mod, ky_mod, origin, axis, start_node, start_element)

        # Add the new mesh to the `Meshes` dictionary
        self.Meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name

    def add_frustrum_mesh(self, name, mesh_size, large_radius, small_radius, height, thickness,
                          material, kx_mod=1.0, ky_mod=1.0, origin=[0, 0, 0], axis='Y',
                          start_node=None, start_element=None):
        """Adds a mesh of quadrilaterals forming a frustrum (a cone intersected by a horizontal
        plane).

        :param name: A unique name for the mesh.
        :type name: str
        :param mesh_size: The target mesh size
        :type mesh_size: number
        :param large_radius: The larger of the two end radii.
        :type large_radius: number
        :param small_radius: The smaller of the two end radii.
        :type small_radius: number
        :param height: The height of the frustrum.
        :type height: number
        :param thickness: The thickness of the elements.
        :type thickness: number
        :param material: The name of the element material.
        :type material: str
        :param kx_mod: Stiffness modification factor for radial stiffness in each element's local
                       x-direction, defaults to 1 (no modification).
        :type kx_mod: number, optional
        :param ky_mod: Stiffness modification factor for meridional stiffness in each
                       element's local y-direction, defaults to 1 (no modification).
        :type ky_mod: number, optional
        :param origin: The origin of the mesh, defaults to [0, 0, 0].
        :type origin: list, optional
        :param axis: The global axis about which the mesh will be generated, defaults to 'Y'.
        :type axis: str, optional
        :param start_node: The name of the first node in the mesh. If set to None the program
                           will use the next available node name, defaults to None.
        :type start_node: str, optional
        :param start_element: The name of the first element in the mesh. If set to `None` the
                              program will use the next available element name, defaults to None
        :type start_element: str, optional
        :raises NameError: Occurs if the specified name already exists.
        :return: The name of the mesh added to the model.
        :rtype: str
        """

        # Check if a name has been provided
        if name:
            # Check that the mesh name doesn't already exist
            if name in self.Meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Give the mesh a new name if necessary
        else:
            name = self.unique_name(self.Meshes, 'MSH')

        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.Nodes, 'N')
        if start_element is None:
            start_element = self.unique_name(self.Quads, 'Q')
        
        # Create a new mesh
        new_mesh = FrustrumMesh(mesh_size, large_radius, small_radius, height, thickness, material,
                                self, kx_mod, ky_mod, origin, axis, start_node, start_element)

        # Add the new mesh to the `Meshes` dictionary
        self.Meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name
    
    def add_cylinder_mesh(self, name, mesh_size, radius, height, thickness, material, kx_mod=1,
                          ky_mod=1, origin=[0, 0, 0], axis='Y', num_elements=None, start_node=None,
                          start_element=None, element_type='Quad'):
        """
        Adds a mesh of elements forming a cylinder.

        Parameters
        ----------
        name : str
            A unique name for the mesh.
        mesh_size : number
            The target mesh size
        radius : number
            The radius of the cylinder.
        height : number
            The height of the cylinder.
        thickness : number
            Element thickness.
        material : str
            The name of the element material.
        kx_mod : number
            Stiffness modification factor for hoop stiffness in each
            element's local x-direction. Default value is 1.0 (no
            modification).
        ky_mod : number
            Stiffness modification factor for meridional stiffness in each
            element's local y-direction. Default value is 1.0 (no
            modification).
        origin : list, optional
            The origin of the mesh. The default is [0, 0, 0].
        axis : str, optional
            The global axis about which the mesh will be generated. The default is 'Y'.
        num_elements : number
            The number of elements to use to form the perimeter of each course. This is typically
            only used if you are trying to match the nodes to another mesh's nodes. If set to
            `None` the program will automatically calculate the number of elements to use based on
            the mesh size. The default is None.
        start_node : str
            The name of the first node in the mesh. If set to `None` the program will use the next
            available node name. Default is `None`
        start_element : str
            The name of the first element in the mesh. If set to `None` the program will use the
            next available element name. Default is `None`
        element_type : str, optional
            The type of element to make the mesh out of. Either 'Quad' or 'Rect'. The default is
            'Quad'.

        """
        
        # Check if a name has been provided
        if name:
            # Check that the mesh name doesn't already exist
            if name in self.Meshes: raise NameError(f"Mesh name '{name}' already exists")
        # Give the mesh a new name if necessary
        else:
            name = self.unique_name(self.Meshes, 'MSH')

        # Identify the starting node and element
        if start_node is None:
            start_node = self.unique_name(self.Nodes, 'N')
        if element_type == 'Rect' and start_element is None:
            start_element = self.unique_name(self.Plates, 'R')
        elif element_type == 'Quad' and start_element is None:
            start_element = self.unique_name(self.Quads, 'Q')
        
        # Create a new mesh
        new_mesh = CylinderMesh(mesh_size, radius, height, thickness, material, self,
                               kx_mod, ky_mod, origin, axis, start_node, start_element,
                               num_elements, element_type)

        # Add the new mesh to the `Meshes` dictionary
        self.Meshes[name] = new_mesh

        # Flag the model as unsolved
        self.solution = None
        
        #Return the mesh's name
        return name

    def merge_duplicate_nodes(self, tolerance=0.001):
        """
        Removes duplicate nodes from the model and returns a list of the removed node names.

        Parameters
        ----------
        tolerance : number
            The maximum distance between two nodes in order to consider them duplicates.

        Returns
        -------
        remove_list : list
            A list of th enames of the nodes that were removed.
        """

        # Initialize a list of nodes to be removed from the `Nodes` dictionary
        node_remove_list = []

        # Initialize a dictionary marking where each node is used
        node_lookup = {node_name: [] for node_name in self.Nodes.keys()}
        element_dicts = ('Springs', 'Members', 'Plates', 'Quads')
        node_types = ('i_node', 'j_node', 'm_node', 'n_node')

        # Step through each dictionary of elements in the model (springs, members, plates, quads)
        for element_dict in element_dicts:

            # Step through each element in the dictionary
            for element in getattr(self, element_dict).values():

                # Step through each possible node type in the element (i-node, j-node, m-node, n-node)
                for node_type in node_types:

                    # Get the current element's node having the current type
                    # Return `None` if the element doesn't have this node type
                    node = getattr(element, node_type, None)

                    # Determine if the node exists on the element
                    if node is not None:
                        # Add the element to the list of elements attached to the node
                        node_lookup[node.name].append((element, node_type))

        # Make a copy of the `Nodes` dictionary
        temp = list(self.Nodes.values())

        # Step through each node in the copy of the `Nodes` dictionary
        for i, node_1 in enumerate(temp):

            # Skip iteration if `node_1` has already been removed
            if node_lookup[node_1.name] is None:
                continue

            # There is no need to check `node_1` against itself
            for node_2 in temp[i + 1:]:

                # Skip iteration if node_2 has already been removed
                if node_lookup[node_2.name] is None:
                    continue

                # Calculate the distance between nodes
                if self.Nodes[node_1.name].distance(self.Nodes[node_2.name]) > tolerance:
                    continue

                # Replace references to `node_2` in each element with references to `node_1`
                for element, node_type in node_lookup[node_2.name]:
                    setattr(element, node_type, node_1)

                # Replace references to `node_2` in each mesh with references to `node_1`
                for mesh in self.Meshes.values():

                    if node_2.name in mesh.nodes.keys():
                        mesh.nodes[node_2.name] = node_1

                    for element in mesh.elements.values():
                        if element.i_node is node_2: element.i_node = node_1
                        if element.j_node is node_2: element.j_node = node_1
                        if element.m_node is node_2: element.m_node = node_1
                        if element.n_node is node_2: element.n_node = node_1

                # Mark `node_2` for removal from the `Nodes` dictionary
                node_remove_list.append(node_2.name)

                # Flag `node_2` as not used
                node_lookup[node_2.name] = None

                # Merge any boundary conditions
                support_cond = ('support_DX', 'support_DY', 'support_DZ', 'support_RX', 'support_RY', 'support_RZ')
                for dof in support_cond:
                    if getattr(node_2, dof) == True:
                        setattr(node_1, dof, True)
                
                # Merge any spring supports
                spring_cond = ('spring_DX', 'spring_DY', 'spring_DZ', 'spring_RX', 'spring_RY', 'spring_RZ')
                for dof in spring_cond:
                    value = getattr(node_2, dof)
                    if value != [None, None, None]:
                        setattr(node_1, dof, value)

        # Delete the duplicate nodes from the model
        for name in node_remove_list:
            self.Nodes.pop(name)
        
        # Flag the model as unsolved
        self.solution = None

        # Return a list of the names of nodes that were removed from the model
        return node_remove_list

    def delete_node(self, node_name):
        '''
        Removes a node from the model. All nodal loads associated with the
        node and elements attached to the node will also be removed.
        
        Parameters
        ----------
        node_name : str
            The name of the node to be removed.
        '''
        
        # Remove the node. Nodal loads are stored within the node, so they
        # will be deleted automatically when the node is deleted.
        self.Nodes.pop(node_name)
        
        # Find any elements attached to the node and remove them
        self.Members = {name: member for name, member in self.Members.items() if member.i_node.name != node_name and member.j_node.name != node_name}
        self.Plates = {name: plate for name, plate in self.Plates.items() if plate.i_node.name != node_name and plate.j_node.name != node_name and plate.m_node.name != node_name and plate.n_node.name != node_name}
        self.Quads = {name: quad for name, quad in self.Quads.items() if quad.i_node.name != node_name and quad.j_node.name != node_name and quad.m_node.name != node_name and quad.n_node.name != node_name}

        # Flag the model as unsolved
        self.solution = None

    def delete_auxnode(self, auxnode_name):
        '''
        Removes an auxiliary node from the model.

        Parameters
        ----------
        auxnode_name : str
            The name of the auxiliary node to be removed
        '''

        # Remove the auxiliary node
        self.AuxNodes.pop(auxnode_name)

        # Remove the auxiliary node from any members that were using it
        for member in self.Members.values():
            if member.auxNode == auxnode_name:
                member.auxNode = None
        
        # Flag the model as unsolved
        self.solution = None

    def delete_spring(self, spring_name):
        '''
        Removes a spring from the model.
        
        Parameters
        ----------
        spring_name : str
            The name of the spring to be removed.
        '''
        
        # Remove the spring
        self.Springs.pop(spring_name)

        # Flag the model as unsolved
        self.solution = None

    def delete_member(self, member_name):
        '''
        Removes a member from the model. All member loads associated with the
        member will also be removed.
        
        Parameters
        ----------
        member_name : str
            The name of the member to be removed.
        '''
        
        # Remove the member. Member loads are stored within the member, so they
        # will be deleted automatically when the member is deleted.
        self.Members.pop(member_name)

        # Flag the model as unsolved
        self.solution = None
        
    def def_support(self, node_name, support_DX=False, support_DY=False, support_DZ=False, support_RX=False, support_RY=False, support_RZ=False):
        """Defines the support conditions at a node. Nodes will default to fully unsupported
           unless specified otherwise.

        :param node_name: The name of the node where the support is being defined.
        :type node_name: str
        :param support_DX: Indicates whether the node is supported against translation in the
                           global X-direction. Defaults to False.
        :type support_DX: bool, optional
        :param support_DY: Indicates whether the node is supported against translation in the
                           global Y-direction. Defaults to False.
        :type support_DY: bool, optional
        :param support_DZ: Indicates whether the node is supported against translation in the
                           global Z-direction. Defaults to False.
        :type support_DZ: bool, optional
        :param support_RX: Indicates whether the node is supported against rotation about the
                           global X-axis. Defaults to False.
        :type support_RX: bool, optional
        :param support_RY: Indicates whether the node is supported against rotation about the
                           global Y-axis. Defaults to False.
        :type support_RY: bool, optional
        :param support_RZ: Indicates whether the node is supported against rotation about the
                           global Z-axis. Defaults to False.
        :type support_RZ: bool, optional
        """            
        
        # Get the node to be supported
        node = self.Nodes[node_name]
                
        # Set the node's support conditions
        node.support_DX = support_DX
        node.support_DY = support_DY
        node.support_DZ = support_DZ
        node.support_RX = support_RX
        node.support_RY = support_RY
        node.support_RZ = support_RZ

        # Flag the model as unsolved
        self.solution = None

    def def_support_spring(self, node_name, dof, stiffness, direction=None):
        """
        Defines a spring support at a node.

        Parameters
        ----------
        node_name : str
            The name of the node to apply the spring support to.
        dof : str ('DX', 'DY', 'DZ', 'RX', 'RY', or 'RZ')
            The degree of freedom to apply the spring support to.
        stiffness : number
            The translational or rotational stiffness of the spring support.
        direction : str or None ('+', '-', None)
            The direction in which the spring can act. '+' allows the spring
            to resist positive displacements. '-' allows the spring to resist
            negative displacements. None allows the spring to act in both
            directions. Default is None.
        """
        
        if dof in ('DX', 'DY', 'DZ', 'RX', 'RY', 'RZ'):
            if direction in ('+', '-', None):
                if dof == 'DX':
                    self.Nodes[node_name].spring_DX = [stiffness, direction, True]
                elif dof == 'DY':
                    self.Nodes[node_name].spring_DY = [stiffness, direction, True]
                elif dof == 'DZ':
                    self.Nodes[node_name].spring_DZ = [stiffness, direction, True]
                elif dof == 'RX':
                    self.Nodes[node_name].spring_RX = [stiffness, direction, True]
                elif dof == 'RY':
                    self.Nodes[node_name].spring_RY = [stiffness, direction, True]
                elif dof == 'RZ':
                    self.Nodes[node_name].spring_RZ = [stiffness, direction, True]
            else:
                raise ValueError('Invalid support spring direction. Specify \'+\', \'-\', or None.')
        else:
            raise ValueError('Invalid support spring degree of freedom. Specify \'DX\', \'DY\', \'DZ\', \'RX\', \'RY\', or \'RZ\'')
        
        # Flag the model as unsolved
        self.solution = None

    def def_node_disp(self, node_name, direction, magnitude): 
        '''
        Defines a nodal displacement at a node.

        node_name : str
            The name of the node where the nodal displacement is being applied.
        direction : str
            The global direction the nodal displacement is being applied in. Displacements are 'DX', 'DY', and 'DZ'. Rotations are 'RX', 'RY', and 'RZ'.
        magnitude : number
            The magnitude of the displacement.
        '''
        # Validate the value of Direction
        if direction not in ('DX', 'DY', 'DZ', 'RX', 'RY', 'RZ'):
            raise ValueError(f"Direction must be 'DX', 'DY', 'DZ', 'RX', 'RY', or 'RZ'. {direction} was given.")
        # Get the node
        node = self.Nodes[node_name]

        if direction == 'DX':
            node.EnforcedDX = magnitude
        if direction == 'DY':
            node.EnforcedDY = magnitude
        if direction == 'DZ':
            node.EnforcedDZ = magnitude
        if direction == 'RX':
            node.EnforcedRX = magnitude
        if direction == 'RY':
            node.EnforcedRY = magnitude
        if direction == 'RZ':
            node.EnforcedRZ = magnitude
        
        # Flag the model as unsolved
        self.solution = None

    def def_releases(self, Member, Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=False, Rzi=False, Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=False, Rzj=False):
        '''
        Defines member end releases.
        
        All member end releases will default to unreleased unless specified otherwise.
        
        Parameters
        ----------
        Member : str
            The name of the member to have its releases modified.
        Dxi : boolean
            Indicates whether the member is released axially at its start.
        Dyi : boolean
            Indicates whether the member is released for shear in the local y-axis at its start.
        Dzi : boolean
            Indicates whether the member is released for shear in the local z-axis at its start.
        Rxi : boolean
            Indicates whether the member is released for torsion at its start.
        Ryi : boolean
            Indicates whether the member is released for moment about the local y-axis at its start.
        Rzi : boolean
            Indicates whether the member is released for moment about the local z-axis at its start.
        Dxj : boolean
            Indicates whether the member is released axially at its end.
        Dyj : boolean
            Indicates whether the member is released for shear in the local y-axis at its end.
        Dzj : boolean
            Indicates whether the member is released for shear in the local z-axis.
        Rxj : boolean
            Indicates whether the member is released for torsion at its end.
        Ryj : boolean
            Indicates whether the member is released for moment about the local y-axis at its end.
        Rzj : boolean
            Indicates whether the member is released for moment about the local z-axis at its end.
        '''
        
        # Apply the end releases to the member
        self.Members[Member].Releases = [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]     

        # Flag the model as unsolved
        self.solution = None

    def add_load_combo(self, name, factors, combo_type='strength'):
        '''
        Adds a load combination to the model

        Parameters
        ----------
        name : str
            A unique name for the load combination (e.g. '1.2D+1.6L+0.5S' or 'Gravity Combo').
        factors : dictionary
            A dictionary containing load cases and their corresponding factors (e.g. {'D':1.2, 'L':1.6, 'S':0.5}).
        combo_type : str
            A description of the type of load combination (e.g. 'strength', 'service'). Currently
            this does nothing in the program, and is a placeholder for future features.
        '''

        # Create a new load combination object
        new_combo = LoadCombo(name, combo_type, factors)

        # Add the load combination to the dictionary of load combinations
        self.LoadCombos[name] = new_combo

        # Flag the model as solved
        self.solution = None

    def add_node_load(self, Node, Direction, P, case='Case 1'):
        '''
        Adds a nodal load to the model.
        
        Parameters
        ----------
        Node : str
            The name of the node where the load is being applied.
        Direction : {'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'}
            The global direction the load is being applied in. Forces are 'FX', 'FY', and 'FZ'. Moments are 'MX', 'MY', and 'MZ'.
        P : number
            The numeric value (magnitude) of the load.
        case : str
            The name of the load case the load belongs to.
        '''
        # Validate the value of Direction
        if Direction not in ('FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'):
            raise ValueError(f"Direction must be 'FX', 'FY', 'FZ', 'MX', 'MY', or 'MZ'. {Direction} was given.")
        # Add the node load to the model
        self.Nodes[Node].NodeLoads.append((Direction, P, case))

        # Flag the model as unsolved
        self.solution = None

    def add_member_pt_load(self, Member, Direction, P, x, case='Case 1'):
        '''
        Adds a member point load to the model.
        
        Parameters
        ----------
        Member : str
            The name of the member the load is being applied to.
        Direction : {'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'}
            The direction in which the force is to be applied. Note that
            typical beam sign convention is used. Transverse forces acting
            toward the beam are positive. Moments are positive if they act
            counter-clockwise relative to the beam's local coordinate system.
            Torsional point loads follow the right hand rule for sign convention.
        P : number
            The numeric value (magnitude) of the load.
        x : number
            The load's location along the member's local x-axis.
        '''

        # Validate the value of Direction
        if Direction not in ('Fx', 'Fy', 'Fz', 'FX', 'FY', 'FZ', 'Mx', 'My', 'Mz', 'MX', 'MY', 'MZ'):
            raise ValueError(f"Direction must be 'Fx', 'Fy', 'Fz', 'FX', 'FY', FZ', 'Mx', 'My', 'Mz', 'MX', 'MY', or 'MZ'. {Direction} was given.")
        
        # Add the point load to the member
        self.Members[Member].PtLoads.append((Direction, P, x, case))
        
        # Flag the model as unsolved
        self.solution = None

    def add_member_dist_load(self, Member, Direction, w1, w2, x1=None, x2=None, case='Case 1'):
        '''
        Adds a member distributed load to the model.
        
        Parameters
        ----------
        Member : str
            The name of the member the load is being appied to
        Direction : {'Fx', 'Fy', 'Fz'}
            The direction in which the load is to be applied. Note that
            typical beam sign convention is used. Forces acting toward the beam
            are positive.
        w1 : number
            The starting value (magnitude) of the load.
        w2 : number
            The ending value (magnitude) of the load.
        x1 : number
            The load's start location along the member's local x-axis. If this argument
            is not specified, the start of the member will be used.
        x2 : number
            The load's end location along the member's local x-axis. If this argument
            is not specified, the end of the member will be used.
        '''
        # Validate the value of Direction
        if Direction not in ('Fx', 'Fy', 'Fz', 'FX', 'FY', 'FZ'):
            raise ValueError(f"Direction must be 'Fx', 'Fy', 'Fz', 'FX', 'FY', or 'FZ'. {Direction} was given.")
        # Determine if a starting and ending points for the load have been specified.
        # If not, use the member start and end as defaults
        if x1 == None:
            start = 0
        else:
            start = x1
        
        if x2 == None:
            end = self.Members[Member].L()
        else:
            end = x2

        # Add the distributed load to the member
        self.Members[Member].DistLoads.append((Direction, w1, w2, start, end, case))
        
        # Flag the model as unsolved
        self.solution = None

    def add_plate_surface_pressure(self, plate_name, pressure, case='Case 1'):
        """
        Adds a surface pressure to the rectangular plate element.

        Parameters
        ----------
        plate_name : str
            The name for the rectangular plate to add the surface pressure to.
        pressure : number
            The value for the surface pressure.
        case : str, optional
            The load case to add the surface pressure to. Default is 'Case 1'.
        
        """

        # Add the surface pressure to the rectangle
        if plate_name in self.Plates.keys():
            self.Plates[plate_name].pressures.append([pressure, case])
        else:
            raise Exception('Invalid plate name specified for plate surface pressure.')
        
        # Flag the model as unsolved
        self.solution = None

    def add_quad_surface_pressure(self, quad_name, pressure, case='Case 1'):
        """
        Adds a surface pressure to the quadrilateral element.

        Parameters
        ----------
        quad_name : str
            The name for the quad to add the surface pressure to.
        pressure : number
            The value for the surface pressure.
        case : str, optional
            The load case to add the surface pressure to. Default is 'Case 1'.
        
        """

        # Add the surface pressure to the quadrilateral
        if quad_name in self.Quads.keys():
            self.Quads[quad_name].pressures.append([pressure, case])
        else:
            raise Exception('Invalid quad name specified for quad surface pressure.')
        
        # Flag the model as unsolved
        self.solution = None

    def delete_loads(self):
        '''
        Deletes all loads from the model along with any results based on the loads.
        '''

        # Delete the member loads and the calculated internal forces
        for member in self.Members.values():
            member.DistLoads = []
            member.PtLoads = []
            member.SegmentsZ = []
            member.SegmentsY = []
            member.SegmentsX = []
        
        # Delete the plate loads
        for plate in self.Plates.values():
            plate.pressures = []
        
        # Delete the quadrilateral loads
        for quad in self.Quads.values():
            quad.pressures = []
        
        # Delete the nodal loads, calculated displacements, and calculated reactions
        for node in self.Nodes.values():

            node.NodeLoads = []

            node.DX = {}
            node.DY = {}
            node.DZ = {}
            node.RX = {}
            node.RY = {}
            node.RZ = {}

            node.RxnFX = {}
            node.RxnFY = {}
            node.RxnFZ = {}
            node.RxnMX = {}
            node.RxnMY = {}
            node.RxnMZ = {}
        
        # Flag the model as unsolved
        self.solution = None

    def _aux_list(self):
        """Builds a list with known nodal displacements and with the positions in global stiffness
           matrix of known and unknown nodal displacements

        :return: A list of the global matrix indices for the unknown nodal displacements (D1_indices). A
                 list of the global matrix indices for the known nodal displacements (D2_indices). A list
                 of the known nodal displacements (D2).
        :rtype: list, list, list
        """

        D1_indices = [] # A list of the indices for the unknown nodal displacements
        D2_indices = [] # A list of the indices for the known nodal displacements
        D2 = []         # A list of the values of the known nodal displacements (D != None)

        # Create the auxiliary table
        for node in self.Nodes.values():
            
            # Unknown displacement DX
            if node.support_DX==False and node.EnforcedDX == None:
                D1_indices.append(node.ID*6 + 0)
            # Known displacement DX
            elif node.EnforcedDX != None:
                D2_indices.append(node.ID*6 + 0)
                D2.append(node.EnforcedDX)
            # Support at DX
            else:
                D2_indices.append(node.ID*6 + 0)
                D2.append(0.0)

            # Unknown displacement DY
            if node.support_DY == False and node.EnforcedDY == None:
                D1_indices.append(node.ID*6 + 1)
            # Known displacement DY
            elif node.EnforcedDY != None:
                D2_indices.append(node.ID*6 + 1)
                D2.append(node.EnforcedDY)
            # Support at DY
            else:
                D2_indices.append(node.ID*6 + 1)
                D2.append(0.0)

            # Unknown displacement DZ
            if node.support_DZ == False and node.EnforcedDZ == None:
                D1_indices.append(node.ID*6 + 2)
            # Known displacement DZ
            elif node.EnforcedDZ != None:
                D2_indices.append(node.ID*6 + 2)
                D2.append(node.EnforcedDZ)
            # Support at DZ
            else:
                D2_indices.append(node.ID*6 + 2)
                D2.append(0.0)

            # Unknown displacement RX
            if node.support_RX == False and node.EnforcedRX == None:
                D1_indices.append(node.ID*6 + 3)
            # Known displacement RX
            elif node.EnforcedRX != None:
                D2_indices.append(node.ID*6 + 3)
                D2.append(node.EnforcedRX)
            # Support at RX
            else:
                D2_indices.append(node.ID*6 + 3)
                D2.append(0.0)

            # Unknown displacement RY
            if node.support_RY == False and node.EnforcedRY == None:
                D1_indices.append(node.ID*6 + 4)
            # Known displacement RY
            elif node.EnforcedRY != None:
                D2_indices.append(node.ID*6 + 4)
                D2.append(node.EnforcedRY)
            # Support at RY
            else:
                D2_indices.append(node.ID*6 + 4)
                D2.append(0.0)

            # Unknown displacement RZ
            if node.support_RZ == False and node.EnforcedRZ == None:
                D1_indices.append(node.ID*6 + 5)
            # Known displacement RZ
            elif node.EnforcedRZ != None:
                D2_indices.append(node.ID*6 + 5)
                D2.append(node.EnforcedRZ)
            # Support at RZ
            else:
                D2_indices.append(node.ID*6 + 5)
                D2.append(0.0)

        # Return the indices and the known displacements
        return D1_indices, D2_indices, D2
               
    def K(self, combo_name='Combo 1', log=False, check_stability=True, sparse=True):
        """
        Returns the model's global stiffness matrix.

        The stiffness matrix will be returned in scipy's sparse lil format,
        which reduces memory usage and can be easily converted to other
        formats.

        Parameters
        ----------
        combo_name : str, optional
            The load combination to get the stiffness matrix for. Default is 'Combo 1'.
        log : bool, optional
            Prints updates to the console if set to True. Default is False.
        sparse : bool, optional
            Returns a sparse matrix if set to True, and a dense matrix
            otherwise. Default is True
        
        Returns
        -------
        K : ndarray or coo_matrix
            The global stiffness matrix for the structure.
        """
        
        # Determine if a sparse matrix has been requested
        if sparse == True:
            # The stiffness matrix will be stored as a scipy `coo_matrix`. Scipy's
            # documentation states that this type of matrix is ideal for efficient
            # construction of finite element matrices. When converted to another
            # format, the `coo_matrix` sums values at the same (i, j) index. We'll
            # build the matrix from three lists.
            row = []
            col = []
            data = []
        else:
            # Initialize a dense matrix of zeros
            K = zeros((len(self.Nodes)*6, len(self.Nodes)*6))
        
        # Add stiffness terms for each nodal spring in the model
        if log: print('- Adding nodal spring support stiffness terms to global stiffness matrix')
        for node in self.Nodes.values():
            
            # Determine if the node has any spring supports
            if node.spring_DX[0] != None:

                # Check for an active spring support
                if node.spring_DX[2] == True:
                    m, n = node.ID*6, node.ID*6
                    if sparse == True:
                        data.append(float(node.spring_DX[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_DX[0])

            if node.spring_DY[0] != None:

                # Check for an active spring support
                if node.spring_DY[2] == True:
                    m, n = node.ID*6 + 1, node.ID*6 + 1
                    if sparse == True:
                        data.append(float(node.spring_DY[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_DY[0])

            if node.spring_DZ[0] != None:

                # Check for an active spring support
                if node.spring_DZ[2] == True:
                    m, n = node.ID*6 + 2, node.ID*6 + 2
                    if sparse == True:
                        data.append(float(node.spring_DZ[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_DZ[0])
        
            if node.spring_RX[0] != None:

                # Check for an active spring support
                if node.spring_RX[2] == True:
                    m, n = node.ID*6 + 3, node.ID*6 + 3
                    if sparse == True:
                        data.append(float(node.spring_RX[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_RX[0])

            if node.spring_RY[0] != None:

                # Check for an active spring support
                if node.spring_RY[2] == True:
                    m, n = node.ID*6 + 4, node.ID*6 + 4
                    if sparse == True:
                        data.append(float(node.spring_RY[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_RY[0])
            
            if node.spring_RZ[0] != None:

                # Check for an active spring support
                if node.spring_RZ[2] == True:
                    m, n = node.ID*6 + 5, node.ID*6 + 5
                    if sparse == True:
                        data.append(float(node.spring_RZ[0]))
                        row.append(m)
                        col.append(n)
                    else:
                        K[m, n] += float(node.spring_RZ[0])

        # Add stiffness terms for each spring in the model
        if log: print('- Adding spring stiffness terms to global stiffness matrix')
        for spring in self.Springs.values():
            
            if spring.active[combo_name] == True:

                # Get the spring's global stiffness matrix
                # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                spring_K = spring.K()

                # Step through each term in the spring's stiffness matrix
                # 'a' & 'b' below are row/column indices in the spring's stiffness matrix
                # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
                for a in range(12):
                
                    # Determine if index 'a' is related to the i-node or j-node
                    if a < 6:
                        # Find the corresponding index 'm' in the global stiffness matrix
                        m = spring.i_node.ID*6 + a
                    else:
                        # Find the corresponding index 'm' in the global stiffness matrix
                        m = spring.j_node.ID*6 + (a-6)
                    
                    for b in range(12):
                    
                        # Determine if index 'b' is related to the i-node or j-node
                        if b < 6:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = spring.i_node.ID*6 + b
                        else:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = spring.j_node.ID*6 + (b-6)
                    
                        # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                        if sparse == True:
                            row.append(m)
                            col.append(n)
                            data.append(spring_K[a, b])
                        else:
                            K[m, n] += spring_K[a, b]

        # Add stiffness terms for each physical member in the model
        if log: print('- Adding member stiffness terms to global stiffness matrix')
        for phys_member in self.Members.values():
            
            # Check to see if the physical member is active for the given load combination
            if phys_member.active[combo_name] == True:

                # Step through each sub-member in the physical member and add terms
                for member in phys_member.sub_members.values():
                    
                    # Get the member's global stiffness matrix
                    # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                    member_K = member.K()

                    # Step through each term in the member's stiffness matrix
                    # 'a' & 'b' below are row/column indices in the member's stiffness matrix
                    # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
                    for a in range(12):
                    
                        # Determine if index 'a' is related to the i-node or j-node
                        if a < 6:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.i_node.ID*6 + a
                        else:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.j_node.ID*6 + (a-6)
                        
                        for b in range(12):
                        
                            # Determine if index 'b' is related to the i-node or j-node
                            if b < 6:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.i_node.ID*6 + b
                            else:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.j_node.ID*6 + (b-6)
                        
                            # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                            if sparse == True:
                                row.append(m)
                                col.append(n)
                                data.append(member_K[a, b])
                            else:
                                K[m, n] += member_K[a, b]
                
        # Add stiffness terms for each quadrilateral in the model
        if log: print('- Adding quadrilateral stiffness terms to global stiffness matrix')
        for quad in self.Quads.values():
            
            # Get the quadrilateral's global stiffness matrix
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            quad_K = quad.K()

            # Step through each term in the quadrilateral's stiffness matrix
            # 'a' & 'b' below are row/column indices in the quadrilateral's stiffness matrix
            # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
            for a in range(24):

                # Determine which node the index 'a' is related to
                if a < 6:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.m_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.n_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.i_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.j_node.ID*6 + (a - 18)

                for b in range(24):

                    # Determine which node the index 'b' is related to
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.m_node.ID*6 + b
                    elif b < 12:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.n_node.ID*6 + (b - 6)
                    elif b < 18:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.i_node.ID*6 + (b - 12)
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.j_node.ID*6 + (b - 18)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                    if sparse == True:
                        row.append(m)
                        col.append(n)
                        data.append(quad_K[a, b])
                    else:
                        K[m, n] += quad_K[a, b]
        
        # Add stiffness terms for each plate in the model
        if log: print('- Adding plate stiffness terms to global stiffness matrix')
        for plate in self.Plates.values():
            
            # Get the plate's global stiffness matrix
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            plate_K = plate.K()

            # Step through each term in the plate's stiffness matrix
            # 'a' & 'b' below are row/column indices in the plate's stiffness matrix
            # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
            for a in range(24):

                # Determine which node the index 'a' is related to
                if a < 6:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.i_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.j_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.m_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.n_node.ID*6 + (a - 18)

                for b in range(24):

                    # Determine which node the index 'b' is related to
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.i_node.ID*6 + b
                    elif b < 12:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.j_node.ID*6 + (b - 6)
                    elif b < 18:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.m_node.ID*6 + (b - 12)
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.n_node.ID*6 + (b - 18)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                    if sparse == True:
                        row.append(m)
                        col.append(n)
                        data.append(plate_K[a, b])
                    else:
                        K[m, n] += plate_K[a, b]

        if sparse:
            # The stiffness matrix will be stored as a scipy `coo_matrix`. Scipy's
            # documentation states that this type of matrix is ideal for efficient
            # construction of finite element matrices. When converted to another
            # format, the `coo_matrix` sums values at the same (i, j) index.
            from scipy.sparse import coo_matrix
            row = array(row)
            col = array(col)
            data = array(data)
            K = coo_matrix((data, (row, col)), shape=(len(self.Nodes)*6, len(self.Nodes)*6))

        # Check that there are no nodal instabilities
        if check_stability:
            if log: print('- Checking nodal stability')
            if sparse: self._check_stability(K.tocsr())
            else: self._check_stability(K)

        # Return the global stiffness matrix
        return K    
   
    def Kg(self, combo_name='Combo 1', log=False, sparse=True):
        """
        Returns the model's global geometric stiffness matrix.

        The model must have a static solution prior to obtaining the geometric stiffness matrix.
        Stiffness of plates is not included.

        Parameters
        ----------
        combo_name : str, optional.
            The name of the load combination to derive the matrix for (not the load combination itself).
        log : bool, optional
            Prints updates to the console if set to True. Default is False.
        sparse : bool, optional
            Returns a sparse matrix if set to True, and a dense matrix
            otherwise. Default is True.
                
        Returns
        -------
        Kg : ndarray or coo_matrix
            The global geometric stiffness matrix for the structure.
        """
        
        if sparse == True:
            # Initialize a zero matrix to hold all the stiffness terms. The matrix will be stored as a
            # scipy sparse `lil_matrix`. This matrix format has several advantages. It uses less memory
            # if the matrix is sparse, supports slicing, and can be converted to other formats (sparse
            # or dense) later on for mathematical operations.
            from scipy.sparse import lil_matrix
            Kg = lil_matrix((len(self.Nodes)*6, len(self.Nodes)*6))
        else:
            Kg = zeros(len(self.Nodes)*6, len(self.Nodes)*6)
        
        # Add stiffness terms for each physical member in the model
        if log: print('- Adding member geometric stiffness terms to global geometric stiffness matrix')
        for phys_member in self.Members.values():
            
            # Check to see if the physical member is active for the given load combination
            if phys_member.active[combo_name] == True:

                # Step through each sub-member in the physical member and add terms
                for member in phys_member.sub_members.values():

                    # Calculate the axial force in the member
                    E = member.E
                    A = member.A
                    L = member.L()
                    d = member.d(combo_name)
                    P = E*A/L*(d[6, 0] - d[0, 0])

                    # Get the member's global stiffness matrix
                    # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                    member_Kg = member.Kg(P)

                    # Step through each term in the member's stiffness matrix
                    # 'a' & 'b' below are row/column indices in the member's stiffness matrix
                    # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
                    for a in range(12):
                    
                        # Determine if index 'a' is related to the i-node or j-node
                        if a < 6:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.i_node.ID*6 + a
                        else:
                            # Find the corresponding index 'm' in the global stiffness matrix
                            m = member.j_node.ID*6 + (a-6)
                        
                        for b in range(12):
                        
                            # Determine if index 'b' is related to the i-node or j-node
                            if b < 6:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.i_node.ID*6 + b
                            else:
                                # Find the corresponding index 'n' in the global stiffness matrix
                                n = member.j_node.ID*6 + (b-6)
                        
                            # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                            Kg[m, n] += member_Kg[(a, b)]

        # Return the global geometric stiffness matrix
        return Kg
      
    def FER(self, combo_name='Combo 1'):
        '''
        Assembles and returns the global fixed end reaction vector.

        Parameters
        ----------
        combo_name : str
            The name of the load combination to get the fixed end reaction vector for (not the load combination itself).
        '''
        
        # Initialize a zero vector to hold all the terms
        FER = zeros((len(self.Nodes) * 6, 1))
        
        # Step through each physical member in the model
        for phys_member in self.Members.values():
            
            # Step through each sub-member and add terms
            for member in phys_member.sub_members.values():

                # Get the member's global fixed end reaction vector
                # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
                member_FER = member.FER(combo_name)

                # Step through each term in the member's fixed end reaction vector
                # 'a' below is the row index in the member's fixed end reaction vector
                # 'm' below is the corresponding row index in the global fixed end reaction vector
                for a in range(12):
                    
                    # Determine if index 'a' is related to the i-node or j-node
                    if a < 6:
                        # Find the corresponding index 'm' in the global fixed end reaction vector
                        m = member.i_node.ID * 6 + a
                    else:
                        # Find the corresponding index 'm' in the global fixed end reaction vector
                        m = member.j_node.ID * 6 + (a - 6)
                    
                    # Now that 'm' is known, place the term in the global fixed end reaction vector
                    FER[m, 0] += member_FER[a, 0]
        
        # Add terms for each rectangle in the model
        for plate in self.Plates.values():
            
            # Get the quadrilateral's global fixed end reaction vector
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            plate_FER = plate.FER(combo_name)

            # Step through each term in the quadrilateral's fixed end reaction vector
            # 'a' below is the row index in the quadrilateral's fixed end reaction vector
            # 'm' below is the corresponding row index in the global fixed end reaction vector
            for a in range(24):
                
                # Determine if index 'a' is related to the i-node, j-node, m-node, or n-node
                if a < 6:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.i_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.j_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.m_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.n_node.ID*6 + (a - 18)
                
                # Now that 'm' is known, place the term in the global fixed end reaction vector
                FER[m, 0] += plate_FER[a, 0]

        # Add terms for each quadrilateral in the model
        for quad in self.Quads.values():
            
            # Get the quadrilateral's global fixed end reaction vector
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            quad_FER = quad.FER(combo_name)

            # Step through each term in the quadrilateral's fixed end reaction vector
            # 'a' below is the row index in the quadrilateral's fixed end reaction vector
            # 'm' below is the corresponding row index in the global fixed end reaction vector
            for a in range(24):
                
                # Determine if index 'a' is related to the i-node, j-node, m-node, or n-node
                if a < 6:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.m_node.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.n_node.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.i_node.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.j_node.ID*6 + (a - 18)
                
                # Now that 'm' is known, place the term in the global fixed end reaction vector
                FER[m, 0] += quad_FER[a, 0]

        # Return the global fixed end reaction vector
        return FER
    
    def P(self, combo_name='Combo 1'):
        '''
        Assembles and returns the global nodal force vector.

        Parameters
        ----------
        combo_name : str
            The name of the load combination to get the force vector for (not the load combination itself).
        '''
            
        # Initialize a zero vector to hold all the terms
        P = zeros((len(self.Nodes)*6, 1))
        
        # Get the load combination for the given 'combo_name'
        combo = self.LoadCombos[combo_name]

        # Add terms for each node in the model
        for node in self.Nodes.values():
            
            # Get the node's ID
            ID = node.ID

            # Step through each load factor in the load combination
            for case, factor in combo.factors.items():

                # Add the node's loads to the global nodal load vector
                for load in node.NodeLoads:

                    if load[2] == case:

                        if load[0] == 'FX':
                            P[ID*6 + 0, 0] += factor*load[1]
                        elif load[0] == 'FY':
                            P[ID*6 + 1, 0] += factor*load[1]
                        elif load[0] == 'FZ':
                            P[ID*6 + 2, 0] += factor*load[1]
                        elif load[0] == 'MX':
                            P[ID*6 + 3, 0] += factor*load[1]
                        elif load[0] == 'MY':
                            P[ID*6 + 4, 0] += factor*load[1]
                        elif load[0] == 'MZ':
                            P[ID*6 + 5, 0] += factor*load[1]
        
        # Return the global nodal force vector
        return P

    def D(self, combo_name='Combo 1'):
        """Returns the global displacement vector for the model.

        :param combo_name: The name of the load combination to get the results for. Defaults to
                           'Combo 1'.
        :type combo_name: str, optional
        :return: The global displacement vector for the model
        :rtype: array
        """
 
        # Return the global displacement vector
        return self._D[combo_name]

    def _partition(self, unp_matrix, D1_indices, D2_indices):
        '''
        Partitions a matrix (or vector) into submatrices (or subvectors) based on degree of freedom
        boundary conditions.

        Parameters
        ----------
        unp_matrix : ndarray or lil_matrix
            The unpartitioned matrix (or vector) to be partitioned.
        D1_indices : list
            A list of the indices for degrees of freedom that have unknown displacements.
        D2_indices : list
            A list of the indices for degrees of freedom that have known displacements.
        '''

        # Determine if this is a 1D vector or a 2D matrix

        # 1D vectors
        if unp_matrix.shape[1] == 1:
            # Partition the vector into 2 subvectors
            m1 = unp_matrix[D1_indices, :]
            m2 = unp_matrix[D2_indices, :]
            return m1, m2
        # 2D matrices
        else:
            # Partition the matrix into 4 submatrices
            m11 = unp_matrix[D1_indices, :][:, D1_indices]
            m12 = unp_matrix[D1_indices, :][:, D2_indices]
            m21 = unp_matrix[D2_indices, :][:, D1_indices]
            m22 = unp_matrix[D2_indices, :][:, D2_indices]
            return m11, m12, m21, m22

    def analyze(self, log=False, check_stability=True, check_statics=False, max_iter=30, sparse=True):
        """
        Performs first-order static analysis.
        
        Iterations are performed if tension-only members or
        compression-only members are present.


        Parameters
        ----------
        log : bool, optional
            Prints the analysis log to the console if set to True. Default is False.
        check_statics : bool, optional
            When set to True, causes a statics check to be performed
        max_iter : number, optional
            The maximum number of iterations to try to get convergence
            for tension/compression-only analysis.
        sparse : bool, optional
            Indicates whether the sparse matrix solver should be used. A matrix can be considered
            sparse or dense depening on how many zero terms there are. Structural stiffness
            matrices often contain many zero terms. The sparse solver can offer faster solutions
            for such matrices. Using the sparse solver on dense matrices may lead to slower
            solution times.
        """

        if log:
            print('+-----------+')
            print('| Analyzing |')
            print('+-----------+')

        # Import `scipy` features if the sparse solver is being used
        if sparse == True:
            from scipy.sparse.linalg import spsolve

        # Ensure there is at least 1 load combination to solve if the user didn't define any
        if self.LoadCombos == {}:
            # Create and add a default load combination to the dictionary of load combinations
            self.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})
        
        # Generate all meshes
        for mesh in self.Meshes.values():
            if mesh.is_generated == False:
                mesh.generate()

        # Activate all springs and members for all load combinations
        for spring in self.Springs.values():
            for combo_name in self.LoadCombos.keys():
                spring.active[combo_name] = True
        
        # Activate all physical members for all load combinations
        for phys_member in self.Members.values():
            for combo_name in self.LoadCombos.keys():
                phys_member.active[combo_name] = True
        
        # Assign an internal ID to all nodes and elements in the model
        self._renumber()

        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = self._aux_list()

        # Convert D2 from a list to a vector
        D2 = atleast_2d(D2).T

        # Step through each load combination
        for combo in self.LoadCombos.values():

            if log:
                print('')
                print('- Analyzing load combination ' + combo.name)

            # Keep track of the number of iterations
            iter_count = 1
            convergence = False
            divergence = False

            # Iterate until convergence or divergence occurs
            while convergence == False and divergence == False:

                # Get the partitioned global stiffness matrix K11, K12, K21, K22
                if sparse == True:
                    K11, K12, K21, K22 = self._partition(self.K(combo.name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)
                else:
                    K11, K12, K21, K22 = self._partition(self.K(combo.name, log, check_stability, sparse), D1_indices, D2_indices)

                # Get the partitioned global fixed end reaction vector
                FER1, FER2 = self._partition(self.FER(combo.name), D1_indices, D2_indices)

                # Get the partitioned global nodal force vector       
                P1, P2 = self._partition(self.P(combo.name), D1_indices, D2_indices)          

                # Calculate the global displacement vector
                if log: print('- Calculating global displacement vector')
                if K11.shape == (0, 0):
                    # All displacements are known, so D1 is an empty vector
                    D1 = []
                else:
                    try:
                        # Calculate the unknown displacements D1
                        if sparse == True:
                            # The partitioned stiffness matrix is in `lil` format, which is great
                            # for memory, but slow for mathematical operations. The stiffness
                            # matrix will be converted to `csr` format for mathematical operations.
                            # The `@` operator performs matrix multiplication on sparse matrices.
                            D1 = spsolve(K11.tocsr(), subtract(subtract(P1, FER1), K12.tocsr() @ D2))
                            D1 = D1.reshape(len(D1), 1)
                        else:
                            D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))
                    except:
                        # Return out of the method if 'K' is singular and provide an error message
                        raise Exception('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')

                # Form the global displacement vector, D, from D1 and D2
                D = zeros((len(self.Nodes)*6, 1))

                for node in self.Nodes.values():

                    if D2_indices.count(node.ID*6 + 0) == 1:
                        D.itemset((node.ID*6 + 0, 0), D2[D2_indices.index(node.ID*6 + 0), 0])
                    else:
                        D.itemset((node.ID*6 + 0, 0), D1[D1_indices.index(node.ID*6 + 0), 0]) 

                    if D2_indices.count(node.ID*6 + 1) == 1:
                        D.itemset((node.ID*6 + 1, 0), D2[D2_indices.index(node.ID*6 + 1), 0])
                    else:
                        D.itemset((node.ID*6 + 1, 0), D1[D1_indices.index(node.ID*6 + 1), 0]) 

                    if D2_indices.count(node.ID*6 + 2) == 1:
                        D.itemset((node.ID*6 + 2, 0), D2[D2_indices.index(node.ID*6 + 2), 0])
                    else:
                        D.itemset((node.ID*6 + 2, 0), D1[D1_indices.index(node.ID*6 + 2), 0]) 

                    if D2_indices.count(node.ID*6 + 3) == 1:
                        D.itemset((node.ID*6 + 3, 0), D2[D2_indices.index(node.ID*6 + 3), 0])
                    else:
                        D.itemset((node.ID*6 + 3, 0), D1[D1_indices.index(node.ID*6 + 3), 0]) 

                    if D2_indices.count(node.ID*6 + 4) == 1:
                        D.itemset((node.ID*6 + 4, 0), D2[D2_indices.index(node.ID*6 + 4), 0])
                    else:
                        D.itemset((node.ID*6 + 4, 0), D1[D1_indices.index(node.ID*6 + 4), 0]) 

                    if D2_indices.count(node.ID*6 + 5) == 1:
                        D.itemset((node.ID*6 + 5, 0), D2[D2_indices.index(node.ID*6 + 5), 0])
                    else:
                        D.itemset((node.ID*6 + 5, 0), D1[D1_indices.index(node.ID*6 + 5), 0]) 

                # Save the global displacement vector
                self._D[combo.name] = D

                # Store the calculated global nodal displacements into each node
                for node in self.Nodes.values():
                    node.DX[combo.name] = D[node.ID*6 + 0, 0]
                    node.DY[combo.name] = D[node.ID*6 + 1, 0]
                    node.DZ[combo.name] = D[node.ID*6 + 2, 0]
                    node.RX[combo.name] = D[node.ID*6 + 3, 0]
                    node.RY[combo.name] = D[node.ID*6 + 4, 0]
                    node.RZ[combo.name] = D[node.ID*6 + 5, 0]
                
                # Check for divergence
                if iter_count > max_iter:
                    divergence = True
                    raise Exception('Model diverged during tension/compression-only analysis')
                
                # Assume the model has converged (to be checked below)
                convergence = True

                # Check tension/compression-only spring supports
                if log: print('- Checking for tension/compression-only support spring convergence')
                for node in self.Nodes.values():
                    
                    # Check if a tension/compression-only support spring has not yet converged
                    if node.spring_DX[2] == True and node.spring_DX[1] == '-' and node.DX[combo.name] > 0:
                        node.spring_DX[2] = False
                        convergence = False
                    if node.spring_DY[2] == True and node.spring_DY[1] == '-' and node.DY[combo.name] > 0:
                        node.spring_DY[2] = False
                        convergence = False
                    if node.spring_DZ[2] == True and node.spring_DZ[1] == '-' and node.DZ[combo.name] > 0:
                        node.spring_DZ[2] = False
                        convergence = False
                    if node.spring_RX[2] == True and node.spring_RX[1] == '-' and node.RX[combo.name] > 0:
                        node.spring_RX[2] = False
                        convergence = False
                    if node.spring_RY[2] == True and node.spring_RY[1] == '-' and node.RY[combo.name] > 0:
                        node.spring_RY[2] = False
                        convergence = False
                    if node.spring_RZ[2] == True and node.spring_RZ[1] == '-' and node.RZ[combo.name] > 0:
                        node.spring_RZ[2] = False
                        convergence = False
                    if node.spring_DX[2] == True and node.spring_DX[1] == '+' and node.DX[combo.name] < 0:
                        node.spring_DX[2] = False
                        convergence = False
                    if node.spring_DY[2] == True and node.spring_DY[1] == '+' and node.DY[combo.name] < 0:
                        node.spring_DY[2] = False
                        convergence = False
                    if node.spring_DZ[2] == True and node.spring_DZ[1] == '+' and node.DZ[combo.name] < 0:
                        node.spring_DZ[2] = False
                        convergence = False
                    if node.spring_RX[2] == True and node.spring_RX[1] == '+' and node.RX[combo.name] < 0:
                        node.spring_RX[2] = False
                        convergence = False
                    if node.spring_RY[2] == True and node.spring_RY[1] == '+' and node.RY[combo.name] < 0:
                        node.spring_RY[2] = False
                        convergence = False
                    if node.spring_RZ[2] == True and node.spring_RZ[1] == '+' and node.RZ[combo.name] < 0:
                        node.spring_RZ[2] = False
                        convergence = False

                # Check tension/compression-only springs
                if log: print('- Checking for tension/compression-only spring convergence')
                for spring in self.Springs.values():

                    if spring.active[combo.name] == True:

                        # Check if tension-only conditions exist
                        if spring.tension_only == True and spring.axial(combo.name) > 0:
                            spring.active[combo.name] = False
                            convergence = False
                        
                        # Check if compression-only conditions exist
                        elif spring.comp_only == True and spring.axial(combo.name) < 0:
                            spring.active[combo.name] = False
                            convergence = False

                # Check tension/compression only members
                if log: print('- Checking for tension/compression-only member convergence')
                for phys_member in self.Members.values():

                    # Only run the tension/compression only check if the member is still active
                    if phys_member.active[combo.name] == True:

                        # Check if tension-only conditions exist
                        if phys_member.tension_only == True and phys_member.max_axial(combo.name) > 0:
                            phys_member.active[combo.name] = False
                            convergence = False

                        # Check if compression-only conditions exist
                        elif phys_member.comp_only == True and phys_member.min_axial(combo.name) < 0:
                            phys_member.active[combo.name] = False
                            convergence = False
                
                if convergence == False:
                    if log: print('- Tension/compression-only analysis did not converge. Adjusting stiffness matrix and reanalyzing.')
                else:
                    if log: print('- Tension/compression-only analysis converged after ' + str(iter_count) + ' iteration(s)')

                # Keep track of the number of tension/compression only iterations
                iter_count += 1

        # Calculate reactions
        self._calc_reactions()

        if log:
            print('')     
            print('- Analysis complete')
            print('')

        # Check statics if requested
        if check_statics == True:
            self._check_statics()
        
        # Flag the model as solved
        self.solution = 'Linear TC'

    def analyze_linear(self, log=False, check_stability=True, check_statics=False, sparse=True):
        '''
        Performs first-order static analysis.
        
        This analysis procedure is much faster since it only assembles the global stiffness matrix
        once, rather than once for each load combination. It is not appropriate when non-linear
        behavior such as tension/compression only analysis or P-Delta analysis are required.

        Parameters
        ----------
        log : bool, optional
            Prints the analysis log to the console if set to True. Default is False.
        check_statics : bool, optional
            When set to True, causes a statics check to be performed
        sparse : bool, optional
            Indicates whether the sparse matrix solver should be used. A matrix can be considered
            sparse or dense depening on how many zero terms there are. Structural stiffness
            matrices often contain many zero terms. The sparse solver can offer faster solutions
            for such matrices. Using the sparse solver on dense matrices may lead to slower
            solution times.
        '''

        if log:
            print('+-------------------+')
            print('| Analyzing: Linear |')
            print('+-------------------+')

        # Import `scipy` features if the sparse solver is being used
        if sparse == True:
            from scipy.sparse.linalg import spsolve

        # Ensure there is at least 1 load combination to solve if the user didn't define any
        if self.LoadCombos == {}:
            # Create and add a default load combination to the dictionary of load combinations
            self.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})
                
        # Generate all meshes
        for mesh in self.Meshes.values():
            if mesh.is_generated == False:
                mesh.generate()
        
        # Activate all springs for all load combinations
        for spring in self.Springs.values():
            for combo_name in self.LoadCombos.keys():
                spring.active[combo_name] = True
        
        # Activate all physical members for all load combinations
        for phys_member in self.Members.values():
            for combo_name in self.LoadCombos.keys():
                phys_member.active[combo_name] = True
        
        # Assign an internal ID to all nodes and elements in the model
        self._renumber()

        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = self._aux_list()

        # Convert D2 from a list to a vector
        D2 = atleast_2d(D2).T

        # Get the partitioned global stiffness matrix K11, K12, K21, K22
        combo_name = list(self.LoadCombos.keys())[0]
        if sparse == True:
            K11, K12, K21, K22 = self._partition(self.K(combo_name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)
        else:
            K11, K12, K21, K22 = self._partition(self.K(combo_name, log, check_stability, sparse), D1_indices, D2_indices)

        # Step through each load combination
        for combo in self.LoadCombos.values():

            if log:
                print('')
                print('- Analyzing load combination ' + combo.name)

            # Get the partitioned global fixed end reaction vector
            FER1, FER2 = self._partition(self.FER(combo.name), D1_indices, D2_indices)

            # Get the partitioned global nodal force vector       
            P1, P2 = self._partition(self.P(combo.name), D1_indices, D2_indices)          

            # Calculate the global displacement vector
            if log: print('- Calculating global displacement vector')
            if K11.shape == (0, 0):
                # All displacements are known, so D1 is an empty vector
                D1 = []
            else:
                try:
                    # Calculate the unknown displacements D1
                    if sparse == True:
                        # The partitioned stiffness matrix is in `lil` format, which is great
                        # for memory, but slow for mathematical operations. The stiffness
                        # matrix will be converted to `csr` format for mathematical operations.
                        # The `@` operator performs matrix multiplication on sparse matrices.
                        D1 = spsolve(K11.tocsr(), subtract(subtract(P1, FER1), K12.tocsr() @ D2))
                        D1 = D1.reshape(len(D1), 1)
                    else:
                        D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))
                except:
                    # Return out of the method if 'K' is singular and provide an error message
                    raise Exception('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')

            # Form the global displacement vector, D, from D1 and D2
            D = zeros((len(self.Nodes)*6, 1))

            for node in self.Nodes.values():

                if D2_indices.count(node.ID*6 + 0) == 1:
                    D.itemset((node.ID*6 + 0, 0), D2[D2_indices.index(node.ID*6 + 0), 0])
                else:
                    D.itemset((node.ID*6 + 0, 0), D1[D1_indices.index(node.ID*6 + 0), 0]) 

                if D2_indices.count(node.ID*6 + 1) == 1:
                    D.itemset((node.ID*6 + 1, 0), D2[D2_indices.index(node.ID*6 + 1), 0])
                else:
                    D.itemset((node.ID*6 + 1, 0), D1[D1_indices.index(node.ID*6 + 1), 0]) 

                if D2_indices.count(node.ID*6 + 2) == 1:
                    D.itemset((node.ID*6 + 2, 0), D2[D2_indices.index(node.ID*6 + 2), 0])
                else:
                    D.itemset((node.ID*6 + 2, 0), D1[D1_indices.index(node.ID*6 + 2), 0]) 

                if D2_indices.count(node.ID*6 + 3) == 1:
                    D.itemset((node.ID*6 + 3, 0), D2[D2_indices.index(node.ID*6 + 3), 0])
                else:
                    D.itemset((node.ID*6 + 3, 0), D1[D1_indices.index(node.ID*6 + 3), 0]) 

                if D2_indices.count(node.ID*6 + 4) == 1:
                    D.itemset((node.ID*6 + 4, 0), D2[D2_indices.index(node.ID*6 + 4), 0])
                else:
                    D.itemset((node.ID*6 + 4, 0), D1[D1_indices.index(node.ID*6 + 4), 0]) 

                if D2_indices.count(node.ID*6 + 5) == 1:
                    D.itemset((node.ID*6 + 5, 0), D2[D2_indices.index(node.ID*6 + 5), 0])
                else:
                    D.itemset((node.ID*6 + 5, 0), D1[D1_indices.index(node.ID*6 + 5), 0]) 

            # Save the global displacement vector
            self._D[combo.name] = D

            # Store the calculated global nodal displacements into each node
            for node in self.Nodes.values():
                node.DX[combo.name] = D[node.ID*6 + 0, 0]
                node.DY[combo.name] = D[node.ID*6 + 1, 0]
                node.DZ[combo.name] = D[node.ID*6 + 2, 0]
                node.RX[combo.name] = D[node.ID*6 + 3, 0]
                node.RY[combo.name] = D[node.ID*6 + 4, 0]
                node.RZ[combo.name] = D[node.ID*6 + 5, 0]

        # Calculate reactions
        self._calc_reactions()

        if log:
            print('')     
            print('- Analysis complete')
            print('')

        # Check statics if requested
        if check_statics == True:
            self._check_statics()
        
        # Flag the model as solved
        self.solution = 'Linear'

    def analyze_PDelta(self, log=False, check_stability=True, max_iter=30, tol=0.01, sparse=True):
        """
        Performs second order (P-Delta) analysis.

        Parameters
        ----------
        log : bool, optional
            Prints updates to the console if set to True. Default is False.
        max_iter : number
            The maximum number of iterations permitted. If this value is exceeded the program will
            report divergence.
        tol : number
            The deflection tolerance (as a percentage) between iterations that will be used to
            define whether the model has converged (e.g. 0.01 = deflections must converge within 1%
            between iterations).
        sparse : bool, optional
            Indicates whether the sparse matrix solver should be used. A matrix can be considered
            sparse or dense depening on how many zero terms there are. Structural stiffness
            matrices often contain many zero terms. The sparse solver can offer faster solutions
            for such matrices. Using the sparse solver on dense matrices may lead to slower
            solution times.
        """
        
        if log:
            print('+--------------------+')
            print('| Analyzing: P-Delta |')
            print('+--------------------+')

        # Import `scipy` features if the sparse solver is being used
        if sparse == True:
            from scipy.sparse.linalg import spsolve

        # Ensure there is at least 1 load combination to solve if the user didn't define any
        if self.LoadCombos == {}:
            # Create and add a default load combination to the dictionary of load combinations
            self.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})
                
        # Generate all meshes
        for mesh in self.Meshes.values():
            if mesh.is_generated == False:
                mesh.generate()
        
        # Activate all springs for all load combinations. They can be turned inactive
        # during the course of the tension/compression-only analysis
        for spring in self.Springs.values():
            for combo_name in self.LoadCombos.keys():
                spring.active[combo_name] = True
                
        # Activate all physical members for all load combinations
        for phys_member in self.Members.values():
            for combo_name in self.LoadCombos.keys():
                phys_member.active[combo_name] = True
               
        # Assign an internal ID to all nodes and elements in the model
        self._renumber()
        
        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = self._aux_list()

        # Convert D2 from a list to a matrix
        D2 = array(D2, ndmin=2).T

        # Step through each load combination
        for combo in self.LoadCombos.values():
            
            if log:
                print('')
                print('- Analyzing load combination ' + combo.name)

            iter_count_TC = 1    # Tracks tension/compression-only iterations
            iter_count_PD = 1    # Tracks P-Delta iterations
            prev_results = None  # Used to store results from the previous iteration

            convergence_TC = False  # Tracks tension/compression-only convergence
            convergence_PD = False  # Tracks P-Delta convergence

            divergence_TC = False  # Tracks tension/compression-only divergence
            divergence_PD = False  # Tracks P-Delta divergence

            # Iterate until convergence or divergence occurs
            while ((convergence_TC == False or convergence_PD == False) 
                  and (divergence_TC == False and divergence_PD == False)):

                # Inform the user which iteration we're on
                if log:
                    print('- Beginning tension/compression-only iteration #' + str(iter_count_TC))
                    print('- Beginning P-Delta iteration #' + str(iter_count_PD))

                # On the first iteration, get all the partitioned global matrices
                if iter_count_PD == 1:

                    if sparse == True:
                        K11, K12, K21, K22 = self._partition(self.K(combo.name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)  # Initial stiffness matrix
                    else:
                        K11, K12, K21, K22 = self._partition(self.K(combo.name, log, check_stability, sparse), D1_indices, D2_indices)  # Initial stiffness matrix
                                                       
                    # Check that the structure is stable
                    if log: print('- Checking stability')
                    self._check_stability(K11)

                    # Assemble the force matrices
                    FER1, FER2 = self._partition(self.FER(combo.name), D1_indices, D2_indices)  # Fixed end reactions
                    P1, P2 = self._partition(self.P(combo.name), D1_indices, D2_indices)        # Nodal forces

                # On subsequent iterations, recalculate the stiffness matrix to account for P-Delta
                # effects
                else:

                    # Calculate the partitioned global stiffness matrices
                    if sparse == True:

                        K11, K12, K21, K22 = self._partition(self.K(combo.name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)  # Initial stiffness matrix
                        Kg11, Kg12, Kg21, Kg22 = self._partition(self.Kg(combo.name, log, sparse), D1_indices, D2_indices)                      # Geometric stiffness matrix
                        
                        # The stiffness matrices are currently `lil` format which is great for
                        # memory, but slow for mathematical operations. They will be converted to
                        # `csr` format. The `+` operator performs matrix addition on `csr`
                        # matrices.
                        K11 = K11.tocsr() + Kg11.tocsr()
                        K12 = K12.tocsr() + Kg12.tocsr()
                        K21 = K21.tocsr() + Kg21.tocsr()
                        K22 = K22.tocsr() + Kg22.tocsr()

                    else:

                        K11, K12, K21, K22 = self._partition(self.K(combo.name, log, check_stability, sparse), D1_indices, D2_indices)  # Initial stiffness matrix
                        Kg11, Kg12, Kg21, Kg22 = self._partition(self.Kg(combo.name, log, sparse), D1_indices, D2_indices)              # Geometric stiffness matrix
                        
                        K11 = K11 + Kg11
                        K12 = K12 + Kg12
                        K21 = K21 + Kg21
                        K22 = K22 + Kg22

                # Calculate the global displacement vector
                if log: print('- Calculating the global displacement vector')
                if K11.shape == (0, 0):
                    # All displacements are known, so D1 is an empty vector
                    D1 = []
                else:
                    try:
                        # Calculate the unknown displacements D1
                        if sparse == True:
                            # The partitioned stiffness matrix is already in `csr` format. The `@`
                            # operator performs matrix multiplication on sparse matrices.
                            D1 = spsolve(K11.tocsr(), subtract(subtract(P1, FER1), K12.tocsr() @ D2))
                            D1 = D1.reshape(len(D1), 1)
                        else:
                            # The partitioned stiffness matrix is in `csr` format. It will be
                            # converted to a 2D dense array for mathematical operations.
                            D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))

                    except:
                        # Return out of the method if 'K' is singular and provide an error message
                        raise ValueError('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')


                D = zeros((len(self.Nodes)*6, 1))

                for node in self.Nodes.values():
                    
                    if node.ID*6 + 0 in D2_indices:
                        D[(node.ID*6 + 0, 0)] = D2[D2_indices.index(node.ID*6 + 0), 0]
                    else:
                        D[(node.ID*6 + 0, 0)] = D1[D1_indices.index(node.ID*6 + 0), 0]

                    if node.ID*6 + 1 in D2_indices:
                        D[(node.ID*6 + 1, 0)] = D2[D2_indices.index(node.ID*6 + 1), 0]
                    else:
                        D[(node.ID*6 + 1, 0)] = D1[D1_indices.index(node.ID*6 + 1), 0]

                    if node.ID*6 + 2 in D2_indices:
                        D[(node.ID*6 + 2, 0)] = D2[D2_indices.index(node.ID*6 + 2), 0]
                    else:
                        D[(node.ID*6 + 2, 0)] = D1[D1_indices.index(node.ID*6 + 2), 0]

                    if node.ID*6 + 3 in D2_indices:
                        D[(node.ID*6 + 3, 0)] = D2[D2_indices.index(node.ID*6 + 3), 0]
                    else:
                        D[(node.ID*6 + 3, 0)] = D1[D1_indices.index(node.ID*6 + 3), 0]

                    if node.ID*6 + 4 in D2_indices:
                        D[(node.ID*6 + 4, 0)] = D2[D2_indices.index(node.ID*6 + 4), 0]
                    else:
                        D[(node.ID*6 + 4, 0)] = D1[D1_indices.index(node.ID*6 + 4), 0]

                    if node.ID*6 + 5 in D2_indices:
                        D[(node.ID*6 + 5, 0)] = D2[D2_indices.index(node.ID*6 + 5), 0]
                    else:
                        D[(node.ID*6 + 5, 0)] = D1[D1_indices.index(node.ID*6 + 5), 0]

                # Save the global displacement vector
                self._D[combo.name] = D

                # Store the calculated global nodal displacements into each node
                for node in self.Nodes.values():

                    node.DX[combo.name] = D[node.ID*6 + 0, 0]
                    node.DY[combo.name] = D[node.ID*6 + 1, 0]
                    node.DZ[combo.name] = D[node.ID*6 + 2, 0]
                    node.RX[combo.name] = D[node.ID*6 + 3, 0]
                    node.RY[combo.name] = D[node.ID*6 + 4, 0]
                    node.RZ[combo.name] = D[node.ID*6 + 5, 0]
                
                # Assume the model has converged (to be checked below)
                convergence_TC = True
                
                # Check tension/compression-only spring supports
                if log: print('- Checking for tension/compression-only support spring convergence')
                for node in self.Nodes.values():
                    
                    # Check if a tension/compression-only support spring has not yet converged
                    if node.spring_DX[2] == True and node.spring_DX[1] == '-' and node.DX[combo.name] > 0:
                        node.spring_DX[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_DY[2] == True and node.spring_DY[1] == '-' and node.DY[combo.name] > 0:
                        node.spring_DY[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_DZ[2] == True and node.spring_DZ[1] == '-' and node.DZ[combo.name] > 0:
                        node.spring_DZ[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_RX[2] == True and node.spring_RX[1] == '-' and node.RX[combo.name] > 0:
                        node.spring_RX[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_RY[2] == True and node.spring_RY[1] == '-' and node.RY[combo.name] > 0:
                        node.spring_RY[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_RZ[2] == True and node.spring_RZ[1] == '-' and node.RZ[combo.name] > 0:
                        node.spring_RZ[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_DX[2] == True and node.spring_DX[1] == '+' and node.DX[combo.name] < 0:
                        node.spring_DX[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_DY[2] == True and node.spring_DY[1] == '+' and node.DY[combo.name] < 0:
                        node.spring_DY[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_DZ[2] == True and node.spring_DZ[1] == '+' and node.DZ[combo.name] < 0:
                        node.spring_DZ[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_RX[2] == True and node.spring_RX[1] == '+' and node.RX[combo.name] < 0:
                        node.spring_RX[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_RY[2] == True and node.spring_RY[1] == '+' and node.RY[combo.name] < 0:
                        node.spring_RY[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False
                    if node.spring_RZ[2] == True and node.spring_RZ[1] == '+' and node.RZ[combo.name] < 0:
                        node.spring_RZ[2] = False
                        convergence_TC = False
                        iter_count_PD = 0
                        convergence_PD = False

                # Check for tension/compression-only springs that need to be deactivated
                if log: print('- Checking for tension/compression-only spring convergence')
                for spring in self.Springs.values():

                    # Only run the tension/compression only check if the spring is still active
                    if spring.active[combo.name] == True:

                        # Check if tension-only conditions exist
                        if spring.tension_only == True and spring.axial(combo.name) > 0:
                            
                            spring.active[combo.name] = False
                            convergence_TC = False

                            # Reset the P-Delta analysis for the new geometry
                            iter_count_PD = 0
                            convergence_PD = False

                        # Check if compression-only conditions exist
                        elif spring.comp_only == True and spring.axial(combo.name) < 0:
                            
                            spring.active[combo.name] = False
                            convergence_TC = False

                            # Reset the P-Delta analysis for the new geometry
                            iter_count_PD = 0
                            convergence_PD = False
                
                # Check for tension/compression-only members that need to be deactivated
                if log: print('- Checking for tension/compression-only member convergence')
                for phys_member in self.Members.values():

                    # Only run the tension/compression only check if the member is still active
                    if phys_member.active[combo.name] == True:

                        # Check if tension-only conditions exist
                        if phys_member.tension_only == True and phys_member.max_axial(combo.name) > 0:
                            
                            phys_member.active[combo.name] = False
                            for member in phys_member.sub_members.values():
                                member.active[combo.name] = False
                            convergence_TC = False

                            # Reset the P-Delta analysis for the new geometry
                            iter_count_PD = 0
                            convergence_PD = False

                        # Check if compression-only conditions exist
                        elif phys_member.comp_only == True and phys_member.min_axial(combo.name) < 0:
                            
                            phys_member.active[combo.name] = False
                            for member in phys_member.sub_members.values():
                                member.active[combo.name] = False
                            convergence_TC = False

                            # Reset the P-Delta analysis for the new geometry
                            iter_count_PD = 0
                            convergence_PD = False
                
                # Report on convergence of tension/compression only analysis
                if convergence_TC == False:
                    
                    if log:
                        print('- Tension/compression-only analysis did not converge on this iteration')
                        print('- Stiffness matrix will be adjusted for newly deactivated elements')
                        print('- P-Delta analysis will be restarted')
                    
                    # Increment the tension/compression-only iteration count
                    iter_count_TC += 1

                else:
                    if log: print('- Tension/compression-only analysis converged after ' + str(iter_count_TC) + ' iteration(s)')
                
                # Check for divergence in the tension/compression-only analysis
                if iter_count_TC > max_iter:
                    divergence_TC = True
                    raise Exception('- Model diverged during tension/compression-only analysis')

                # Check for P-Delta convergence
                if iter_count_PD > 1:
                
                    # Print a status update for the user
                    if log: print('- Checking for P-Delta convergence')

                    # Temporarily disable error messages for invalid values.
                    # We'll be dealing with some 'nan' values due to division by zero at supports with zero deflection.
                    seterr(invalid='ignore')

                    # Check for convergence
                    # Note: if the shape of K11 is (0, 0) then all degrees of freedom are fully
                    # restrained, and P-Delta analysis automatically converges
                    if K11.shape == (0, 0) or abs(nanmax(divide(D1, prev_results)) - 1) <= tol:
                        convergence_PD = True
                        if log: print('- P-Delta analysis converged after ' + str(iter_count_PD) + ' iteration(s)')
                    # Check for divergence
                    elif iter_count_PD > max_iter:
                        divergence_PD = True
                        if log: print('- P-Delta analysis failed to converge after ' + str(max_iter) + ' iteration(s)')

                    # Turn invalid value warnings back on
                    seterr(invalid='warn') 

                # Save the results for the next iteration
                prev_results = D1

                # Increment the P-Delta iteration count
                iter_count_PD += 1
        
        # Calculate reactions
        self._calc_reactions()

        if log:
            print('')
            print('- Analysis complete')
            print('')
        
        # Flag the model as solved
        self.solution = 'P-Delta'

    def _calc_reactions(self, log=False):
        """
        Calculates reactions internally once the model is solved.

        Parameters
        ----------
        log : bool, optional
            Prints updates to the console if set to True. Default is False.
        """

        # Print a status update to the console
        if log: print('- Calculating reactions')

        # Calculate the reactions node by node
        for node in self.Nodes.values():
            
            # Step through each load combination
            for combo in self.LoadCombos.values():
                
                # Initialize reactions for this node and load combination
                node.RxnFX[combo.name] = 0.0
                node.RxnFY[combo.name] = 0.0
                node.RxnFZ[combo.name] = 0.0
                node.RxnMX[combo.name] = 0.0
                node.RxnMY[combo.name] = 0.0
                node.RxnMZ[combo.name] = 0.0

                # Determine if the node has any supports
                if (node.support_DX or node.support_DY or node.support_DZ 
                or  node.support_RX or node.support_RY or node.support_RZ):

                    # Sum the spring end forces at the node
                    for spring in self.Springs.values():

                        if spring.i_node == node and spring.active[combo.name] == True:
                            
                            # Get the spring's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            spring_F = spring.F(combo.name)

                            node.RxnFX[combo.name] += spring_F[0, 0]
                            node.RxnFY[combo.name] += spring_F[1, 0]
                            node.RxnFZ[combo.name] += spring_F[2, 0]
                            node.RxnMX[combo.name] += spring_F[3, 0]
                            node.RxnMY[combo.name] += spring_F[4, 0]
                            node.RxnMZ[combo.name] += spring_F[5, 0]

                        elif spring.j_node == node and spring.active[combo.name] == True:
                        
                            # Get the spring's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            spring_F = spring.F(combo.name)
                        
                            node.RxnFX[combo.name] += spring_F[6, 0]
                            node.RxnFY[combo.name] += spring_F[7, 0]
                            node.RxnFZ[combo.name] += spring_F[8, 0]
                            node.RxnMX[combo.name] += spring_F[9, 0]
                            node.RxnMY[combo.name] += spring_F[10, 0]
                            node.RxnMZ[combo.name] += spring_F[11, 0]

                    # Step through each physical member in the model
                    for phys_member in self.Members.values():

                        # Sum the sub-member end forces at the node
                        for member in phys_member.sub_members.values():
                            
                            if member.i_node == node and phys_member.active[combo.name] == True:
                            
                                # Get the member's global force matrix
                                # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                                member_F = member.F(combo.name)

                                node.RxnFX[combo.name] += member_F[0, 0]
                                node.RxnFY[combo.name] += member_F[1, 0]
                                node.RxnFZ[combo.name] += member_F[2, 0]
                                node.RxnMX[combo.name] += member_F[3, 0]
                                node.RxnMY[combo.name] += member_F[4, 0]
                                node.RxnMZ[combo.name] += member_F[5, 0]

                            elif member.j_node == node and phys_member.active[combo.name] == True:
                            
                                # Get the member's global force matrix
                                # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                                member_F = member.F(combo.name)
                            
                                node.RxnFX[combo.name] += member_F[6, 0]
                                node.RxnFY[combo.name] += member_F[7, 0]
                                node.RxnFZ[combo.name] += member_F[8, 0]
                                node.RxnMX[combo.name] += member_F[9, 0]
                                node.RxnMY[combo.name] += member_F[10, 0]
                                node.RxnMZ[combo.name] += member_F[11, 0]

                    # Sum the plate forces at the node
                    for plate in self.Plates.values():

                        if plate.i_node == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[0, 0]
                            node.RxnFY[combo.name] += plate_F[1, 0]
                            node.RxnFZ[combo.name] += plate_F[2, 0]
                            node.RxnMX[combo.name] += plate_F[3, 0]
                            node.RxnMY[combo.name] += plate_F[4, 0]
                            node.RxnMZ[combo.name] += plate_F[5, 0]

                        elif plate.j_node == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[6, 0]
                            node.RxnFY[combo.name] += plate_F[7, 0]
                            node.RxnFZ[combo.name] += plate_F[8, 0]
                            node.RxnMX[combo.name] += plate_F[9, 0]
                            node.RxnMY[combo.name] += plate_F[10, 0]
                            node.RxnMZ[combo.name] += plate_F[11, 0]

                        elif plate.m_node == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[12, 0]
                            node.RxnFY[combo.name] += plate_F[13, 0]
                            node.RxnFZ[combo.name] += plate_F[14, 0]
                            node.RxnMX[combo.name] += plate_F[15, 0]
                            node.RxnMY[combo.name] += plate_F[16, 0]
                            node.RxnMZ[combo.name] += plate_F[17, 0]

                        elif plate.n_node == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[18, 0]
                            node.RxnFY[combo.name] += plate_F[19, 0]
                            node.RxnFZ[combo.name] += plate_F[20, 0]
                            node.RxnMX[combo.name] += plate_F[21, 0]
                            node.RxnMY[combo.name] += plate_F[22, 0]
                            node.RxnMZ[combo.name] += plate_F[23, 0]

                    # Sum the quad forces at the node
                    for quad in self.Quads.values():

                        if quad.m_node == node:

                            # Get the quad's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            quad_F = quad.F(combo.name)

                            node.RxnFX[combo.name] += quad_F[0, 0]
                            node.RxnFY[combo.name] += quad_F[1, 0]
                            node.RxnFZ[combo.name] += quad_F[2, 0]
                            node.RxnMX[combo.name] += quad_F[3, 0]
                            node.RxnMY[combo.name] += quad_F[4, 0]
                            node.RxnMZ[combo.name] += quad_F[5, 0]

                        elif quad.n_node == node:

                            # Get the quad's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            quad_F = quad.F(combo.name)
                    
                            node.RxnFX[combo.name] += quad_F[6, 0]
                            node.RxnFY[combo.name] += quad_F[7, 0]
                            node.RxnFZ[combo.name] += quad_F[8, 0]
                            node.RxnMX[combo.name] += quad_F[9, 0]
                            node.RxnMY[combo.name] += quad_F[10, 0]
                            node.RxnMZ[combo.name] += quad_F[11, 0]

                        elif quad.i_node == node:

                            # Get the quad's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            quad_F = quad.F(combo.name)
                    
                            node.RxnFX[combo.name] += quad_F[12, 0]
                            node.RxnFY[combo.name] += quad_F[13, 0]
                            node.RxnFZ[combo.name] += quad_F[14, 0]
                            node.RxnMX[combo.name] += quad_F[15, 0]
                            node.RxnMY[combo.name] += quad_F[16, 0]
                            node.RxnMZ[combo.name] += quad_F[17, 0]

                        elif quad.j_node == node:

                            # Get the quad's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            quad_F = quad.F(combo.name)
                    
                            node.RxnFX[combo.name] += quad_F[18, 0]
                            node.RxnFY[combo.name] += quad_F[19, 0]
                            node.RxnFZ[combo.name] += quad_F[20, 0]
                            node.RxnMX[combo.name] += quad_F[21, 0]
                            node.RxnMY[combo.name] += quad_F[22, 0]
                            node.RxnMZ[combo.name] += quad_F[23, 0]
                    
                    # Sum the joint loads applied to the node
                    for load in node.NodeLoads:

                        for case, factor in combo.factors.items():
                            
                            if load[2] == case:

                                if load[0] == 'FX':
                                    node.RxnFX[combo.name] -= load[1]*factor
                                elif load[0] == 'FY':
                                    node.RxnFY[combo.name] -= load[1]*factor
                                elif load[0] == 'FZ':
                                    node.RxnFZ[combo.name] -= load[1]*factor
                                elif load[0] == 'MX':
                                    node.RxnMX[combo.name] -= load[1]*factor
                                elif load[0] == 'MY':
                                    node.RxnMY[combo.name] -= load[1]*factor
                                elif load[0] == 'MZ':
                                    node.RxnMZ[combo.name] -= load[1]*factor
                
                # Calculate reactions due to active spring supports at the node
                elif node.spring_DX[0] != None and node.spring_DX[2] == True:
                    sign = node.spring_DX[1]
                    k = node.spring_DX[0]
                    if sign != None: k = float(sign + str(k))
                    DX = node.DX[combo.name]
                    node.RxnFX[combo.name] += k*DX
                elif node.spring_DY[0] != None and node.spring_DY[2] == True:
                    sign = node.spring_DY[1]
                    k = node.spring_DY[0]
                    if sign != None: k = float(sign + str(k))
                    DY = node.DY[combo.name]
                    node.RxnFY[combo.name] += k*DY
                elif node.spring_DZ[0] != None and node.spring_DZ[2] == True:
                    sign = node.spring_DZ[1]
                    k = node.spring_DZ[0]
                    if sign != None: k = float(sign + str(k))
                    DZ = node.DZ[combo.name]
                    node.RxnFZ[combo.name] += k*DZ
                elif node.spring_RX[0] != None and node.spring_RX[2] == True:
                    sign = node.spring_RX[1]
                    k = node.spring_RX[0]
                    if sign != None: k = float(sign + str(k))
                    RX = node.RX[combo.name]
                    node.RxnMX[combo.name] += k*RX
                elif node.spring_RY[0] != None and node.spring_RY[2] == True:
                    sign = node.spring_RY[1]
                    k = node.spring_RY[0]
                    if sign != None: k = float(sign + str(k))
                    RY = node.RY[combo.name]
                    node.RxnMY[combo.name] += k*RY
                elif node.spring_RZ[0] != None and node.spring_RZ[2] == True:
                    sign = node.spring_RZ[1]
                    k = node.spring_RZ[0]
                    if sign != None: k = float(sign + str(k))
                    RZ = node.RZ[combo.name]
                    node.RxnMZ[combo.name] += k*RZ

    def _check_stability(self, K):
        """
        Identifies nodal instabilities in the model's stiffness matrix.
        """

        # Initialize the `unstable` flag to `False`
        unstable = False

        # Step through each diagonal term in the stiffness matrix
        for i in range(K.shape[0]):
            
            # Determine which node this term belongs to
            node = [node for node in self.Nodes.values() if node.ID == int(i/6)][0]

            # Determine which degree of freedom this term belongs to
            dof = i%6

            # Check to see if this degree of freedom is supported
            if dof == 0:
                supported = node.support_DX
            elif dof == 1:
                supported = node.support_DY
            elif dof == 2:
                supported = node.support_DZ
            elif dof == 3:
                supported = node.support_RX
            elif dof == 4:
                supported = node.support_RY
            elif dof == 5:
                supported = node.support_RZ

            # Check if the degree of freedom on this diagonal is unstable
            if K[i, i] == 0 and not supported:

                # Flag the model as unstable
                unstable = True

                # Identify which direction this instability effects
                if i%6 == 0: direction = 'for translation in the global X direction.'
                if i%6 == 1: direction = 'for translation in the global Y direction.'
                if i%6 == 2: direction = 'for translation in the global Z direction.'
                if i%6 == 3: direction = 'for rotation about the global X axis.'
                if i%6 == 4: direction = 'for rotation about the global Y axis.'
                if i%6 == 5: direction = 'for rotation about the global Z axis.'

                # Print a message to the console
                print('* Nodal instability detected: node' + node.name + 'is unstable for' + direction)

        if unstable:
            raise Exception('Unstable node(s). See console output for details.')

        return
    
    def _check_statics(self):
        '''
        Checks static equilibrium and prints results to the console.

        Parameters
        ----------
        precision : number
            The number of decimal places to carry the results to.
        '''

        print('+----------------+')
        print('| Statics Check: |')
        print('+----------------+')
        print('')

        from prettytable import PrettyTable

        # Start a blank table and create a header row
        statics_table = PrettyTable()
        statics_table.field_names = ['Load Combination', 'Sum FX', 'Sum RX', 'Sum FY', 'Sum RY', 'Sum FZ', 'Sum RZ', 'Sum MX', 'Sum RMX', 'Sum MY', 'Sum RMY', 'Sum MZ', 'Sum RMZ']

        # Step through each load combination
        for combo in self.LoadCombos.values():

            # Initialize force and moment summations to zero
            SumFX, SumFY, SumFZ = 0.0, 0.0, 0.0
            SumMX, SumMY, SumMZ = 0.0, 0.0, 0.0
            SumRFX, SumRFY, SumRFZ = 0.0, 0.0, 0.0
            SumRMX, SumRMY, SumRMZ = 0.0, 0.0, 0.0

            # Get the global force vector and the global fixed end reaction vector
            P = self.P(combo.name)
            FER = self.FER(combo.name)

            # Step through each node and sum its forces
            for node in self.Nodes.values():

                # Get the node's coordinates
                X = node.X
                Y = node.Y
                Z = node.Z

                # Get the nodal forces
                FX = P[node.ID*6+0][0] - FER[node.ID*6+0][0]
                FY = P[node.ID*6+1][0] - FER[node.ID*6+1][0]
                FZ = P[node.ID*6+2][0] - FER[node.ID*6+2][0]
                MX = P[node.ID*6+3][0] - FER[node.ID*6+3][0]
                MY = P[node.ID*6+4][0] - FER[node.ID*6+4][0]
                MZ = P[node.ID*6+5][0] - FER[node.ID*6+5][0]

                # Get the nodal reactions
                RFX = node.RxnFX[combo.name]
                RFY = node.RxnFY[combo.name]
                RFZ = node.RxnFZ[combo.name]
                RMX = node.RxnMX[combo.name]
                RMY = node.RxnMY[combo.name]
                RMZ = node.RxnMZ[combo.name]

                # Sum the global forces
                SumFX += FX
                SumFY += FY
                SumFZ += FZ
                SumMX += MX - FY*Z + FZ*Y
                SumMY += MY + FX*Z - FZ*X
                SumMZ += MZ - FX*Y + FY*X

                # Sum the global reactions
                SumRFX += RFX
                SumRFY += RFY
                SumRFZ += RFZ
                SumRMX += RMX - RFY*Z + RFZ*Y
                SumRMY += RMY + RFX*Z - RFZ*X
                SumRMZ += RMZ - RFX*Y + RFY*X   

            # Add the results to the table
            statics_table.add_row([combo.name, '{:.3g}'.format(SumFX), '{:.3g}'.format(SumRFX),
                                               '{:.3g}'.format(SumFY), '{:.3g}'.format(SumRFY),
                                               '{:.3g}'.format(SumFZ), '{:.3g}'.format(SumRFZ),
                                               '{:.3g}'.format(SumMX), '{:.3g}'.format(SumRMX),
                                               '{:.3g}'.format(SumMY), '{:.3g}'.format(SumRMY),
                                               '{:.3g}'.format(SumMZ), '{:.3g}'.format(SumRMZ)])

        # Print the static check table
        print(statics_table)
        print('')
    
    def _renumber(self):
        """
        Assigns node and element ID numbers to be used internally by the program. Numbers are
        assigned according to the order in which they occur in each dictionary.
        """
        
        # Number each node in the model
        for id, node in enumerate(self.Nodes.values()):
            node.ID = id
        
        # Number each spring in the model
        for id, spring in enumerate(self.Springs.values()):
            spring.ID = id

        # Descritize all the physical members and number each member in the model
        id = 0
        for phys_member in self.Members.values():
            phys_member.descritize()
            for member in phys_member.sub_members.values():
                member.ID = id
                id += 1
        
        # Number each plate in the model
        for id, plate in enumerate(self.Plates.values()):
            plate.ID = id
        
        # Number each quadrilateral in the model
        for id, quad in enumerate(self.Quads.values()):
            quad.ID = id
    
    def unique_name(self, dictionary, prefix):
        """Returns the next available unique name for a dictionary of objects.

        :param dictionary: The dictionary to get a unique name for.
        :type dictionary: dict
        :param prefix: The prefix to use for the unique name.
        :type prefix: str
        :return: A unique name for the dictionary.
        :rtype: str
        """

        # Select a trial value for the next available name
        name = prefix + str(len(dictionary) + 1)
        i = 1
        while name in dictionary.keys():
            name = prefix + str(len(dictionary) + i)
            i += 1
        
        # Return the next available name
        return name
    
    def rename(self):
        """
        Renames all the nodes and elements in the model.
        """

        # Rename each node in the model
        temp = self.Nodes.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'N' + str(id)
            self.Nodes[new_key] = self.Nodes.pop(old_key)
            self.Nodes[new_key].name = new_key
            id += 1
        
        # Rename each spring in the model
        temp = self.Springs.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'S' + str(id)
            self.Springs[new_key] = self.Springs.pop(old_key)
            self.Springs[new_key].name = new_key
            id += 1

        # Rename each member in the model
        temp = self.Members.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'M' + str(id)
            self.Members[new_key] = self.Members.pop(old_key)
            self.Members[new_key].name = new_key
            id += 1
        
        # Rename each plate in the model
        temp = self.Plates.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'P' + str(id)
            self.Plates[new_key] = self.Plates.pop(old_key)
            self.Plates[new_key].name = new_key
            id += 1
        
        # Rename each quad in the model
        temp = self.Quads.copy()
        id = 1
        for old_key in temp.keys():
            new_key = 'Q' + str(id)
            self.Quads[new_key] = self.Quads.pop(old_key)
            self.Quads[new_key].name = new_key
            id += 1

    def orphaned_nodes(self):
        """
        Returns a list of the names of nodes that are not attached to any elements.
        """

        # Initialize a list of orphaned nodes
        orphans = []

        # Step through each node in the model
        for node in self.Nodes.values():

            orphaned = False

            # Check to see if the node is attached to any elements
            quads = [quad.name for quad in self.Quads.values() if quad.i_node == node or quad.j_node == node or quad.m_node == node or quad.n_node == node]
            plates = [plate.name for plate in self.Plates.values() if plate.i_node == node or plate.j_node == node or plate.m_node == node or plate.n_node == node]
            members = [member.name for member in self.Members.values() if member.i_node == node or member.j_node == node]
            springs = [spring.name for spring in self.Springs.values() if spring.i_node == node or spring.j_node == node]

            # Determine if the node is orphaned
            if quads == [] and plates == [] and members == [] and springs == []:
                orphaned = True
            
            # Add the orphaned nodes to the list of orphaned nodes
            if orphaned == True:
                orphans.append(node.name)
        
        return orphans
      