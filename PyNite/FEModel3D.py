# %%
from numpy import array, matrix, zeros, matmul, divide, add, subtract, atleast_2d
from numpy import nanmax, seterr
from numpy.linalg import solve

from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix, csr_matrix
from PyNite.Node3D import Node3D
from PyNite.Spring3D import Spring3D
from PyNite.Member3D import Member3D
from PyNite.Quad3D import Quad3D
from PyNite.Plate3D import Plate3D
from PyNite.LoadCombo import LoadCombo

import warnings

# %%
class FEModel3D():
    '''
    A class representing a 3D finite element model.
    '''

#%% 
    def __init__(self):
        '''
        Initializes a new 3D finite element model.
        '''
        
        self.Nodes = {}      # A dictionary of the structure's nodes
        self.AuxNodes = {}   # A dictionary of the structure's auxiliary nodes
        self.Springs = {}    # A dictionary of the structure's springs
        self.Members = {}    # A dictionary of the structure's members
        self.Quads = {}      # A dictionary of the structure's quadiralterals
        self.Plates = {}     # A dictionary of the structure's rectangular plates
        self.__D = {}        # A dictionary of the structure's nodal displacements by load combination
        self.LoadCombos = {} # A dictionary of the structure's load combinations

#%%
    def AddNode(self, name, X, Y, Z):
        warnings.warn('`AddNode` will be replaced with `add_node` in a future version of PyNite.', FutureWarning)
        self.add_node(name, X, Y, Z)

    def add_node(self, name, X, Y, Z):
        '''
        Adds a new node to the model. The node name will be returned
        
        Parameters
        ----------
        name : string
            A unique user-defined name for the node. If None or "", a name will be automatically assigned
        X : number
            The global X-coordinate of the node.
        Y : number
            The global Y-coordinate of the node.
        Z : number
            The global Z-coordinate of the node.
        '''
        
        # Name the node or check it doesn't already exist
        if name:
            if name in self.Nodes: raise NameError(f"Node name '{name}' already exists")
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
        
        #Return the node name
        return name

#%%
    def AddAuxNode(self, name, X, Y, Z):
        warnings.warn('`AddAuxNode` will be replaced with `add_auxnode` in a future version of PyNite.', FutureWarning)
        self.add_auxnode(name, X, Y, Z)

    def add_auxnode(self, name, X, Y, Z):
        '''
        Adds a new auxiliary node to the model. The node name will be returned
        
        Together with a member's `i` and `j` nodes, an auxiliary node defines the plane in which
        the member's local z-axis lies, and the side of the member the z-axis point towards. If no
        auxiliary node is specified for a member, PyNite uses its own defaut configuration.

        Parameters
        ----------
        name : string
            A unique user-defined name for the node. If None or "", a name will be automatically assigned
        X : number
            The global X-coordinate of the node.
        Y : number
            The global Y-coordinate of the node.
        Z : number
            The global Z-coordinate of the node.
        '''
        
        # Name the node or check it doesn't already exist
        if name:
            if name in self.AuxNodes: raise NameError(f"Auxnode name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "AN" + str(len(self.AuxNodes))
            count = 1
            while name in self.AuxNodes: 
                name = "AN" + str(len(self.AuxNodes)+count)
                count += 1
                
        # Create a new node
        newNode = Node3D(name, X, Y, Z)
        
        # Add the new node to the list
        self.AuxNodes[name] = newNode
        
        #Return the node name
        return name
  
#%%
    def AddSpring(self, name, iNode, jNode, ks, tension_only=False, comp_only=False):
        warnings.warn('`AddSpring` will be replaced with `add_spring` in a future version of PyNite.', FutureWarning)
        self.add_spring(name, iNode, jNode, ks, tension_only, comp_only)

    def add_spring(self, name, iNode, jNode, ks, tension_only=False, comp_only=False):
        '''
        Adds a new spring to the model. The spring name will be returned.
        
        Parameters
        ----------
        name : string
            A unique user-defined name for the member. If None or "", a name will be automatically assigned
        iNode : string
            The name of the i-node (start node).
        jNode : string
            The name of the j-node (end node).
        ks : number
            The spring constant (force/displacement).
        tension_only : bool, optional
            Indicates if the member is tension-only. Default is False.
        comp_only : bool, optional
            Indicates if the member is compression-only. Default is False.
        '''

        # Name the spring or check it doesn't already exist
        if name:
            if name in self.Springs: raise NameError(f"Spring name '{name}' already exists")
        else:
            # As a guess, start with the length of the dictionary
            name = "S" + str(len(self.Springs))
            count = 1
            while name in self.Springs: 
                name = "S" + str(len(self.Springs)+count)
                count += 1
                
        # Create a new spring
        newSpring = Spring3D(name, self.Nodes[iNode], self.Nodes[jNode], ks,
                             self.LoadCombos, tension_only=tension_only, comp_only=comp_only)
        
        # Add the new member to the list
        self.Springs[name] = newSpring
        
        #Return the spring name
        return name

#%%
    def AddMember(self, name, iNode, jNode, E, G, Iy, Iz, J, A, auxNode=None,
                  tension_only=False, comp_only=False):
        warnings.warn('`AddMember` will be replaced with `add_member` in a future version of PyNite.', FutureWarning)
        self.add_member(name, iNode, jNode, E, G, Iy, Iz, J, A, auxNode, tension_only, comp_only)

    def add_member(self, name, iNode, jNode, E, G, Iy, Iz, J, A, auxNode=None,
                   tension_only=False, comp_only=False):
        '''
        Adds a new member to the model. The member name will be returned.
        
        Parameters
        ----------
        name : string
            A unique user-defined name for the member. If None or "", a name will be automatically assigned
        iNode : string
            The name of the i-node (start node).
        jNode : string
            The name of the j-node (end node).
        E : number
            The modulus of elasticity of the member.
        G : number
            The shear modulus of the member.
        Iy : number
            The moment of inertia of the member about its local y-axis.
        Iz : number
            The moment of inertia of the member about its local z-axis.
        J : number
            The polar moment of inertia of the member.
        A : number
            The cross-sectional area of the member.
        auxNode : string, optional
            The name of the auxialary node used to define the local z-axis.
            The default is for the program to define the axis instead of
            using an auxiliary node.
        tension_only : bool, optional
            Indicates if the member is tension-only. Default is False.
        comp_only : bool, optional
            Indicates if the member is compression-only. Default is False.
        '''

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
            newMember = Member3D(name, self.Nodes[iNode],
            self.Nodes[jNode], E, G, Iy, Iz, J, A,
            LoadCombos=self.LoadCombos, tension_only=tension_only, comp_only=comp_only)
        else:
            newMember = Member3D(name, self.Nodes[iNode],
            self.Nodes[jNode], E, G, Iy, Iz, J, A, self.GetAuxNode(auxNode),
            self.LoadCombos, tension_only=tension_only, comp_only=comp_only)
        
        # Add the new member to the list
        self.Members[name] = newMember
        
        #Return the member name
        return name

#%%
    def AddPlate(self, name, iNode, jNode, mNode, nNode, t, E, nu):
        warnings.warn('`AddPlate` will be replaced with `add_plate` in a future version of PyNite.', FutureWarning)
        self.add_plate(name, iNode, jNode, mNode, nNode, t, E, nu)

    def add_plate(self, name, iNode, jNode, mNode, nNode, t, E, nu):
        '''
        Adds a new plate to the model. The plate name will be returned
        
        Parameters
        ----------
        name : string
            A unique user-defined name for the plate. If None or "", a name will be automatically assigned
        iNode : string
            The name of the i-node (1st node definded in clockwise order).
        jNode : string
            The name of the j-node (2nd node defined in clockwise order).
        mNode : string
            The name of the m-node (3rd node defined in clockwise order).
        nNode : string
            The name of the n-node (4th node defined in clockwise order).
        t : number
            The thickness of the plate.
        E : number
            The modulus of elasticity of the plate.
        mew : number
            Posson's ratio for the plate.
        '''
        
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
        
        # Create a new member
        newPlate = Plate3D(name, self.Nodes[iNode], self.Nodes[jNode], self.Nodes[mNode], self.Nodes[nNode], t, E, nu, self.LoadCombos)
        
        # Add the new member to the list
        self.Plates[name] = newPlate
        
        #Return the plate name
        return name

#%%
    def AddQuad(self, name, iNode, jNode, mNode, nNode, t, E, nu):
        warnings.warn('`AddQuad` will be replaced with `add_quad` in a future version of PyNite.', FutureWarning)
        self.add_quad(name, iNode, jNode, mNode, nNode, t, E, nu)

    def add_quad(self, name, iNode, jNode, mNode, nNode, t, E, nu):
        '''
        Adds a new quadrilateral to the model. The quad name will be returned

        Quadrilaterals are similar to plates, except they do not have to be
        rectangular.
        
        Parameters
        ----------
        name : string
            A unique user-defined name for the quadrilateral. If None or "", a name will be automatically assigned
        iNode : string
            The name of the i-node (1st node definded in counter-clockwise order).
        jNode : string
            The name of the j-node (2nd node defined in counter-clockwise order).
        mNode : string
            The name of the m-node (3rd node defined in counter-clockwise order).
        nNode : string
            The name of the n-node (4th node defined in counter-clockwise order).
        t : number
            The thickness of the quadrilateral.
        E : number
            The modulus of elasticity of the quadrilateral.
        nu : number
            Posson's ratio for the quadrilateral.
        x_axis : list, Default is None
            A vector in the direction of the plate's local x-axis, for example [1, 0, 0]. If None
            the program defines the x_axis as runing from the i-node to the j-node.
        '''
        
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
        newQuad = Quad3D(name, self.Nodes[iNode], self.Nodes[jNode], self.Nodes[mNode], self.Nodes[nNode], t, E, nu, self.LoadCombos)
        
        # Add the new member to the list
        self.Quads[name] = newQuad
        
        #Return the quad name
        return name

#%%
    def AddMesh(self, mesh):
        warnings.warn('`AddMesh` will be replaced with `add_mesh` in a future version of PyNite.', FutureWarning)
        self.add_mesh(mesh)

    def add_mesh(self, mesh):
        """
        Adds a predefined mesh to the model.

        Parameters
        ----------
        mesh : Mesh
            A mesh object.
        """

        # Add the mesh's nodes to the finite element model
        self.Nodes.update(mesh.nodes)

        # Add the mesh's elements to the finite element model
        if list(mesh.elements.values())[0].type == 'Quad':
            self.Quads.update(mesh.elements)
        elif list(mesh.elements.values())[0].type == 'Rect':
            self.Plates.update(mesh.elements)

        # If the mesh contained duplicate keys that were already in the `Nodes` dictionary, the old
        # nodes will have been abandoned in favor of the new nodes. Reattach any elements that are
        # still referencing the old nodes to the new nodes.
        for element in self.Quads.values():
            element.iNode = self.Nodes[element.iNode.Name]
            element.jNode = self.Nodes[element.jNode.Name]
            element.mNode = self.Nodes[element.mNode.Name]
            element.nNode = self.Nodes[element.nNode.Name]
        
        for element in self.Plates.values():
            element.iNode = self.Nodes[element.iNode.Name]
            element.jNode = self.Nodes[element.jNode.Name]
            element.mNode = self.Nodes[element.mNode.Name]
            element.nNode = self.Nodes[element.nNode.Name]

        for element in self.Members.values():
            element.iNode = self.Nodes[element.iNode.Name]
            element.jNode = self.Nodes[element.jNode.Name]

        # Add the mesh's elements to the finite element model
        for element in mesh.elements.values():
            element.LoadCombos = self.LoadCombos

#%%
    def MergeDuplicateNodes(self, tolerance=0.001):
        warnings.warn('`MergeDuplicateNodes` will be replaced with `merge_duplicate_nodes` in a future version of PyNite.', FutureWarning)
        self.merge_duplicate_nodes(tolerance)

    def merge_duplicate_nodes(self, tolerance=0.001):

        # Initialize a list of nodes to be removed from the `Nodes` dictionary
        remove_list = []

        # Make a copy of the `Nodes` dictionary
        temp = self.Nodes.copy()

        # Step through each node in the `Nodes` dictionary
        for name_1, node_1 in self.Nodes.items():

            # Make a temporary copy of the dictionary that's missing all the previously checked nodes
            temp = temp.copy()
            temp.pop(name_1)

            # Compare all the other nodes to `node_1`
            for name_2, node_2 in temp.items():

                # Calculate the distance between nodes
                dist = ((node_2.X - node_1.X)**2 + (node_2.Y - node_1.Y)**2 + (node_2.Z - node_1.Z)**2)**0.5

                # Determine if the nodes need to be merged
                if dist < tolerance:

                    # Attach any elements that were using `node_2` to `node_1` instead
                    for element in self.Quads.values():
                        if element.iNode.Name == name_2:
                            element.iNode = self.Nodes[name_1]
                        if element.jNode.Name == name_2:
                            element.jNode = self.Nodes[name_1]
                        if element.mNode.Name == name_2:
                            element.mNode = self.Nodes[name_1]
                        if element.nNode.Name == name_2:
                            element.nNode = self.Nodes[name_1]
        
                    for element in self.Plates.values():
                        if element.iNode.Name == name_2:
                            element.iNode = self.Nodes[name_1]
                        if element.jNode.Name == name_2:
                            element.jNode = self.Nodes[name_1]
                        if element.mNode.Name == name_2:
                            element.mNode = self.Nodes[name_1]
                        if element.nNode.Name == name_2:
                            element.nNode = self.Nodes[name_1]

                    for element in self.Members.values():
                        if element.iNode.Name == name_2:
                            element.iNode = self.Nodes[name_1]
                        if element.jNode.Name == name_2:
                            element.jNode = self.Nodes[name_1]
                    
                    # Add `name_2` to the list of nodes to be removed
                    remove_list.append(name_2)

        for name in remove_list:
            self.Nodes.pop(name)
            
#%%
    def RemoveNode(self, node_name):
        warnings.warn('`RemoveNode` will be replaced with `delete_node` in a future version of PyNite.', FutureWarning)
        self.delete_node(node_name)

    def delete_node(self, node_name):
        '''
        Removes a node from the model. All nodal loads associated with the
        node and elements attached to the node will also be removed.
        
        Parameters
        ----------
        node_name : string
            The name of the node to be removed.
        '''
        
        # Remove the node. Nodal loads are stored within the node, so they
        # will be deleted automatically when the node is deleted.
        self.Nodes.pop(node_name)
        
        # Find any elements attached to the node and remove them
        self.Members = {name: member for name, member in self.Members.items() if member.iNode.Name != node_name and member.jNode.Name != node_name}
        self.Plates = {name: plate for name, plate in self.Plates.items() if plate.iNode.Name != node_name and plate.jNode.Name != node_name and plate.mNode.Name != node_name and plate.nNode.Name != node_name}
        self.Quads = {name: quad for name, quad in self.Quads.items() if quad.iNode.Name != node_name and quad.jNode.Name != node_name and quad.mNode.Name != node_name and quad.nNode.Name != node_name}

#%%
    def RemoveAuxNode(self, auxnode_name):
        warnings.warn('`RemoveAuxNode` will be replaced with `delete_auxnode` in a future version of PyNite.', FutureWarning)
        self.delete_auxnode(auxnode_name)

    def delete_auxnode(self, auxnode_name):
        '''
        Removes an auxiliary node from the model.

        Parameters
        ----------
        auxnode_name : string
            The name of the auxiliary node to be removed
        '''

        # Remove the auxiliary node
        self.AuxNodes.pop(auxnode_name)

        # Remove the auxiliary node from any members that were using it
        for member in self.Members.values():
            if member.auxNode == auxnode_name:
                member.auxNode = None

#%%
    def RemoveSpring(self, spring_name):
        warnings.warn('`RemoveSpring` will be replaced with `delete_spring` in a future version of PyNite.', FutureWarning)
        self.delete_spring(spring_name)

    def delete_spring(self, spring_name):
        '''
        Removes a spring from the model.
        
        Parameters
        ----------
        spring_name : string
            The name of the spring to be removed.
        '''
        
        # Remove the spring
        self.Springs.pop(spring_name)

#%%
    def RemoveMember(self, member_name):
        warnings.warn('`RemoveMember` will be replaced with `delete_member` in a future version of PyNite.', FutureWarning)
        self.delete_member(member_name)

    def delete_member(self, member_name):
        '''
        Removes a member from the model. All member loads associated with the
        member will also be removed.
        
        Parameters
        ----------
        member_name : string
            The name of the member to be removed.
        '''
        
        # Remove the member. Member loads are stored within the member, so they
        # will be deleted automatically when the member is deleted.
        self.Members.pop(member_name)
        
#%%
    def DefineSupport(self, node_name, SupportDX=False, SupportDY=False, SupportDZ=False, SupportRX=False, SupportRY=False, SupportRZ=False):
        warnings.warn('`DefineSupport` will be replaced with `def_support` in a future version of PyNite.', FutureWarning)
        self.def_support(node_name, SupportDX, SupportDY, SupportDZ, SupportRX, SupportRY, SupportRZ)

    def def_support(self, node_name, SupportDX=False, SupportDY=False, SupportDZ=False, SupportRX=False, SupportRY=False, SupportRZ=False):
        '''
        Defines the support conditions at a node.
        
        Nodes will default to fully unsupported unless specified otherwise.
        
        Parameters
        ----------
        node_name : string
            The name of the node where the support is being defined
        SupportDX : bool
            Indicates whether the node is supported against translation in the global X-direction.
        SupportDY : bool
            Indicates whether the node is supported against translation in the global Y-direction.
        SupportDZ : bool
            Indicates whether the node is supported against translation in the global Z-direction.
        SupportRX : bool
            Indicates whether the node is supported against rotation about the global X-axis.
        SupportRY : bool
            Indicates whether the node is supported against rotation about the global Y-axis.
        SupportRZ : bool
            Indicates whether the node is supported against rotation about the global Z-axis.
        '''
        
        # Get the node to be supported
        node = self.Nodes[node_name]
                
        # Set the node's support conditions
        node.SupportDX = SupportDX
        node.SupportDY = SupportDY
        node.SupportDZ = SupportDZ
        node.SupportRX = SupportRX
        node.SupportRY = SupportRY
        node.SupportRZ = SupportRZ

#%%    
    def AddNodeDisplacement(self, Node, Direction, Magnitude):
        warnings.warn('`AddNodeDisplacement` will be replaced with `def_node_disp` in a future version of PyNite.', FutureWarning)
        self.def_node_disp(Node, Direction, Magnitude)
        
    def def_node_disp(self, Node, Direction, Magnitude): 
        '''
        Defines a nodal displacement at a node.

        Node : string
            The name of the node where the nodal displacement is being applied.
        Direction : {'DX', 'DY', 'DZ', 'RX', 'RY', 'RZ'}
            The global direction the nodal displacement is being applied in. Displacements are 'DX', 'DY', and 'DZ'. Rotations are 'RX', 'RY', and 'RZ'.
            Sign convention follows the model's global coordinate system.
        Magnitude : number
            The magnitude of the displacement.
        '''
        # Validate the value of Direction
        if Direction not in ('DX', 'DY', 'DZ', 'RX', 'RY', 'RZ'):
            raise ValueError(f"Direction must be 'DX', 'DY', 'DZ', 'RX', 'RY', or 'RZ'. {Direction} was given.")
        # Get the node
        node = self.Nodes[Node]

        if Direction == 'DX':
            node.EnforcedDX = Magnitude
        if Direction == 'DY':
            node.EnforcedDY = Magnitude
        if Direction == 'DZ':
            node.EnforcedDZ = Magnitude
        if Direction == 'RX':
            node.EnforcedRX = Magnitude
        if Direction == 'RY':
            node.EnforcedRY = Magnitude
        if Direction == 'RZ':
            node.EnforcedRZ = Magnitude

#%%
    def DefineReleases(self, Member, Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=False, Rzi=False, Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=False, Rzj=False):
        warnings.warn('`DefineReleases` will be replaced with `def_releases` in a future version of PyNite.', FutureWarning)
        self.def_releases(Member, Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj)

    def def_releases(self, Member, Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=False, Rzi=False, Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=False, Rzj=False):
        '''
        Defines member end releases.
        
        All member end releases will default to unreleased unless specified otherwise.
        
        Parameters
        ----------
        Member : string
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

#%%
    def AddLoadCombo(self, name, factors, combo_type='strength'):
        warnings.warn('`AddLoadCombo` will be replaced with `add_load_combo` in a future version of PyNite.', FutureWarning)
        self.add_load_combo(name, factors, combo_type)

    def add_load_combo(self, name, factors, combo_type='strength'):
        '''
        Adds a load combination to the model

        Parameters
        ----------
        name : string
            A unique name for the load combination (e.g. '1.2D+1.6L+0.5S' or 'Gravity Combo').
        factors : dictionary
            A dictionary containing load cases and their corresponding factors (e.g. {'D':1.2, 'L':1.6, 'S':0.5}).
        combo_type : string
            A description of the type of load combination (e.g. 'strength', 'service'). Currently
            this does nothing in the program, and is a placeholder for future features.
        '''

        # Create a new load combination object
        new_combo = LoadCombo(name, combo_type, factors)

        # Add the load combination to the dictionary of load combinations
        self.LoadCombos[name] = new_combo

#%%
    def add_node_load(self, Node, Direction, P, case='Case 1'):
        warnings.warn('`AddNodeLoad` will be replaced with `add_node_load` in a future version of PyNite.', FutureWarning)
        self.add_node_load(Node, Direction, P, case)

    def add_node_load(self, Node, Direction, P, case='Case 1'):
        '''
        Adds a nodal load to the model.
        
        Parameters
        ----------
        Node : string
            The name of the node where the load is being applied.
        Direction : {'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'}
            The global direction the load is being applied in. Forces are 'FX', 'FY', and 'FZ'. Moments are 'MX', 'MY', and 'MZ'.
        P : number
            The numeric value (magnitude) of the load.
        case : string
            The name of the load case the load belongs to.
        '''
        # Validate the value of Direction
        if Direction not in ('FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'):
            raise ValueError(f"Direction must be 'FX', 'FY', 'FZ', 'MX', 'MY', or 'MZ'. {Direction} was given.")
        # Add the node load to the model
        self.Nodes[Node].NodeLoads.append((Direction, P, case))

#%%
    def AddMemberPtLoad(self, Member, Direction, P, x, case='Case 1'):
        warnings.warn('`AddMemberPtLoad` will be replaced with `add_member_pt_load` in a future version of PyNite.', FutureWarning)
        self.add_member_pt_load(Member, Direction, P, x, case)

    def add_member_pt_load(self, Member, Direction, P, x, case='Case 1'):
        '''
        Adds a member point load to the model.
        
        Parameters
        ----------
        Member : string
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
        if Direction not in ('Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'):
            raise ValueError(f"Direction must be 'Fx', 'Fy', 'Fz', 'Mx', 'My', or 'Mz'. {Direction} was given.")
        # Add the point load to the member
        self.Members[Member].PtLoads.append((Direction, P, x, case))

#%%
    def AddMemberDistLoad(self, Member, Direction, w1, w2, x1=None, x2=None, case='Case 1'):
        warnings.warn('`AddMemberDistLoad` will be replaced with `add_member_dist_load` in a future version of PyNite.', FutureWarning)
        self.add_member_dist_load(Member, Direction, w1, w2, x1, x2, case)

    def add_member_dist_load(self, Member, Direction, w1, w2, x1=None, x2=None, case='Case 1'):
        '''
        Adds a member distributed load to the model.
        
        Parameters
        ----------
        Member : string
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
        if Direction not in ('Fx', 'Fy', 'Fz'):
            raise ValueError(f"Direction must be 'Fx', 'Fy', 'Fz'. {Direction} was given.")
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
        
#%%
    def add_plate_surface_pressure(self, plate_ID, pressure, case='Case 1'):
        warnings.warn('`AddPlateSurfacePressure` will be replaced with `add_plate_surface_pressure` in a future version of PyNite.', FutureWarning)
        self.add_plate_surface_pressure(plate_ID, pressure, case)

    def add_plate_surface_pressure(self, plate_ID, pressure, case='Case 1'):
        '''
        Adds a surface pressure to the rectangular plate element.
        '''

        # Add the surface pressure to the rectangle
        self.Plates[plate_ID].pressures.append([pressure, case])

#%%
    def AddQuadSurfacePressure(self, quad_ID, pressure, case='Case 1'):
        warnings.warn('`AddQuadSurfacePressure` will be replaced with `add_quad_surface_pressure` in a future version of PyNite.', FutureWarning)
        self.add_quad_surface_pressure(quad_ID, pressure, case)

    def add_quad_surface_pressure(self, quad_ID, pressure, case='Case 1'):
        '''
        Adds a surface pressure to the quadrilateral element.
        '''

        # Add the surface pressure to the quadrilateral
        self.Quads[quad_ID].pressures.append([pressure, case])

#%%
    def ClearLoads(self):
        warnings.warn('`ClearLoads` will be replaced with `delete_loads` in a future version of PyNite.', FutureWarning)
        self.delete_loads()

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

#%%
    def GetNode(self, name):
        '''
        Returns the node with the given name.
        
        Parameters
        ----------
        name : string
            The name of the node to be returned.
        '''
        
        warnings.warn('`GetNode` will be deleted in a future version of PyNite. Use `Nodes[node_name]` instead.', FutureWarning)

        try:
            return self.Nodes[name]
        except:
            # If the node name is not found
            raise ValueError(f"Node '{name}' was not found in the model.")

#%%         
    def GetAuxNode(self, name):
        '''
        Returns the auxiliary node with the given name.
        
        Parameters
        ----------
        name : string
            The name of the auxiliary node to be returned.
        '''
        
        warnings.warn('`GetAuxNode` will be deleted in a future version of PyNite. Use `AuxNodes[auxnode_name]` instead.', FutureWarning)

        try:
            return self.AuxNodes[name]
        except:
            # If the auxiliary node name is not found
            raise ValueError(f"Auxiliary node '{name}' was not found in the model.")

#%%
    def GetSpring(self, name):
        '''
        Returns the spring with the given name.
        
        Parameters
        ----------
        name : string
            The name of the spring to be returned.
        '''
        
        warnings.warn('`GetSpring` will be deleted in a future version of PyNite. Use `Springs[spring_name]` instead.', FutureWarning)
        
        try:
            return self.Springs[name]
        except:
            # If the auxiliary node name is not found
            raise ValueError(f"Spring '{name}' was not found in the model.")

#%%
    def GetMember(self, name):
        '''
        Returns the member with the given name.
        
        Parameters
        ----------
        name : string
            The name of the member to be returned.
        '''

        warnings.warn('`GetMember` will be deleted in a future version of PyNite. Use `Members[member_name]` instead.', FutureWarning)

        try:
            return self.Members[name]
        except:
            # If the member name is not found
            raise ValueError(f"Member '{name}' was not found in the model.")

#%%
    def GetPlate(self, name):
        '''
        Returns the plate with the given name.
        
        Parameters
        ----------
        name : string
            The name of the plate to be returned.
        '''
        
        warnings.warn('`GetPlate` will be deleted in a future version of PyNite. Use `Plates[plate_name]` instead.', FutureWarning)

        try:
            return self.Plates[name]
        except:
            # If the plate name is not found
            raise ValueError(f"Plate '{name}' was not found in the model")

#%%
    def GetQuad(self, name):
        '''
        Returns the quadrilateral with the given name.
        
        Parameters
        ----------
        name : string
            The name of the quadrilateral to be returned.
        '''
        
        warnings.warn('`GetQuad` will be deleted in a future version of PyNite. Use `Quads[quad_name]` instead.', FutureWarning)
        
        try:
            return self.Quads[name]
        except:
            # If the quad name is not found
            raise ValueError(f"Quad '{name}' was not found in the model")

#%%
    def _renumber(self):
        '''
        Assigns node, spring, member, and plate member ID numbers to be used internally by the
        program. Numbers are assigned according to the order nodes, springs, members, and plates
        were added to the model.
        '''
        
        # Number each node in the model
        for id, node in enumerate(self.Nodes.values()):
            node.ID = id
        
        # Number each spring in the model
        for id, spring in enumerate(self.Springs.values()):
            spring.ID = id

        # Number each member in the model
        for id, member in enumerate(self.Members.values()):
            member.ID = id
        
        # Number each plate in the model
        for id, plate in enumerate(self.Plates.values()):
            plate.ID = id
        
        # Number each quadrilateral in the model
        for id, quad in enumerate(self.Quads.values()):
            quad.ID = id

#%%
    def _aux_list(self):
        '''
        Builds a list with known nodal displacements and with the positions in global stiffness matrix of known 
        and unknown nodal displacements

        Returns
        -------
        D1_indices : number
            A list of the global matrix indices for the unknown nodal displacements
        D2_indices : number
            A list of the global matrix indices for the known nodal displacements
        D2 : number
            A list of the known nodal displacements
        '''

        D1_indices = [] # A list of the indices for the unknown nodal displacements
        D2_indices = [] # A list of the indices for the known nodal displacements
        D2 = []         # A list of the values of the known nodal displacements (D != None)

        # Create the auxiliary table
        for node in self.Nodes.values():
            
            # Unknown displacement DX
            if ((node.SupportDX == False and node.EnforcedDX == None)
               or type(node.SupportDX) == int or type(node.SupportDX) == float):
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
            if ((node.SupportDY == False and node.EnforcedDY == None) 
               or type(node.SupportDY) == int or type(node.SupportDY) == float):
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
            if ((node.SupportDZ == False and node.EnforcedDZ == None)
               or type(node.SupportDZ) == int or type(node.SupportDZ) == float):
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
            if ((node.SupportRX == False and node.EnforcedRX == None)
               or type(node.SupportRX) == int or type(node.SupportRX) == float):
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
            if ((node.SupportRY == False and node.EnforcedRY == None)
               or type(node.SupportRY) == int or type(node.SupportRY) == float):
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
            if ((node.SupportRZ == False and node.EnforcedRZ == None)
               or type(node.SupportRZ) == int or type(node.SupportRZ) == float):
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
            
#%%    
    def K(self, combo_name='Combo 1', log=False):
        '''
        Returns the model's global stiffness matrix.

        Parameters
        ----------
        combo_name : string, optional
            The load combination to get the stiffness matrix for. Default is 'Combo 1'.
        log : bool, optional
            Prints updates to the console if set to True. Default is False.
        '''
        
        # Initialize a zero matrix to hold all the stiffness terms. The matrix will be stored as a
        # scipy sparse `lil_matrix`. This matrix format has several advantages. It uses less memory
        # if the matrix is sparse, supports slicing, and can be converted to other formats (sparse
        # or dense) later on for mathematical operations.
        K = lil_matrix((len(self.Nodes)*6, len(self.Nodes)*6))
        
        # Add stiffness terms for each nodal spring in the model
        if log: print('- Adding nodal spring support stiffness terms to global stiffness matrix')
        for node in self.Nodes.values():
            
            # Determine if the node has any spring supports
            if type(node.SupportDX) == float or type(node.SupportDX) == int:
                m, n = node.ID*6, node.ID*6
                K[(m, n)] += node.SupportDX
            if type(node.SupportDY) == float or type(node.SupportDY) == int:
                m, n = node.ID*6 + 1, node.ID*6 + 1
                K[(m, n)] += node.SupportDY
            if type(node.SupportDZ) == float or type(node.SupportDZ) == int:
                m, n = node.ID*6 + 2, node.ID*6 + 2
                K[(m, n)] += node.SupportDZ
            if type(node.SupportRX) == float or type(node.SupportRX) == int:
                m, n = node.ID*6 + 3, node.ID*6 + 3
                K[(m, n)] += node.SupportRX
            if type(node.SupportRY) == float or type(node.SupportRY) == int:
                m, n = node.ID*6 + 4, node.ID*6 + 4
                K[(m, n)] += node.SupportRY
            if type(node.SupportRZ) == float or type(node.SupportRZ) == int:
                m, n = node.ID*6 + 5, node.ID*6 + 5
                K[(m, n)] += node.SupportRZ

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
                        m = spring.iNode.ID*6 + a
                    else:
                        # Find the corresponding index 'm' in the global stiffness matrix
                        m = spring.jNode.ID*6 + (a-6)
                    
                    for b in range(12):
                    
                        # Determine if index 'b' is related to the i-node or j-node
                        if b < 6:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = spring.iNode.ID*6 + b
                        else:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = spring.jNode.ID*6 + (b-6)
                    
                        # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                        K[(m, n)] += spring_K[(a, b)]

        # Add stiffness terms for each member in the model
        if log: print('- Adding member stiffness terms to global stiffness matrix')
        for member in self.Members.values():
            
            if member.active[combo_name] == True:

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
                        m = member.iNode.ID*6 + a
                    else:
                        # Find the corresponding index 'm' in the global stiffness matrix
                        m = member.jNode.ID*6 + (a-6)
                    
                    for b in range(12):
                    
                        # Determine if index 'b' is related to the i-node or j-node
                        if b < 6:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = member.iNode.ID*6 + b
                        else:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = member.jNode.ID*6 + (b-6)
                    
                        # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                        K[(m, n)] += member_K[(a, b)]
                
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
                    m = quad.mNode.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.nNode.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.iNode.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = quad.jNode.ID*6 + (a - 18)

                for b in range(24):

                    # Determine which node the index 'b' is related to
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.mNode.ID*6 + b
                    elif b < 12:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.nNode.ID*6 + (b - 6)
                    elif b < 18:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.iNode.ID*6 + (b - 12)
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = quad.jNode.ID*6 + (b - 18)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
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
                    m = plate.iNode.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.nNode.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.mNode.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.jNode.ID*6 + (a - 18)

                for b in range(24):

                    # Determine which node the index 'b' is related to
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.iNode.ID*6 + b
                    elif b < 12:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.nNode.ID*6 + (b - 6)
                    elif b < 18:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.mNode.ID*6 + (b - 12)
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.jNode.ID*6 + (b - 18)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                    K[m, n] += plate_K[a, b]

        # Return the global stiffness matrix
        return K      

#%%    
    def Kg(self, combo_name='Combo 1', log=False):
        '''
        Returns the model's global geometric stiffness matrix.

        The model must have a static solution prior to obtaining the geometric stiffness matrix.
        Stiffness of plates is not included.

        Parameters
        ----------
        combo_name : string, optional.
            The name of the load combination to derive the matrix for (not the load combination itself).
        log : bool, optional
            Prints updates to the console if set to True. Default is False.
        '''
        
        # Initialize a zero matrix to hold all the stiffness terms. The matrix will be stored as a
        # scipy sparse `lil_matrix`. This matrix format has several advantages. It uses less memory
        # if the matrix is sparse, supports slicing, and can be converted to other formats (sparse
        # or dense) later on for mathematical operations.
        Kg = lil_matrix((len(self.Nodes)*6, len(self.Nodes)*6))
        
        # Add stiffness terms for each member in the model
        if log: print('- Adding member geometric stiffness terms to global geometric stiffness matrix')
        for member in self.Members.values():
            
            if member.active[combo_name] == True:

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
                        m = member.iNode.ID*6 + a
                    else:
                        # Find the corresponding index 'm' in the global stiffness matrix
                        m = member.jNode.ID*6 + (a-6)
                    
                    for b in range(12):
                    
                        # Determine if index 'b' is related to the i-node or j-node
                        if b < 6:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = member.iNode.ID*6 + b
                        else:
                            # Find the corresponding index 'n' in the global stiffness matrix
                            n = member.jNode.ID*6 + (b-6)
                    
                        # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                        Kg[m, n] += member_Kg[(a, b)]

        # Return the global geometric stiffness matrix
        return Kg
     
#%%    
    def FER(self, combo_name='Combo 1'):
        '''
        Assembles and returns the global fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the fixed end reaction vector for (not the load combination itself).
        '''
        
        # Initialize a zero vector to hold all the terms
        FER = zeros((len(self.Nodes) * 6, 1))
        
        # Add terms for each member in the model
        for member in self.Members.values():
            
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
                    m = member.iNode.ID * 6 + a
                else:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = member.jNode.ID * 6 + (a - 6)
                
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
                    m = plate.iNode.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.nNode.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.mNode.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = plate.jNode.ID*6 + (a - 18)
                
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
                    m = quad.mNode.ID*6 + a
                elif a < 12:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.nNode.ID*6 + (a - 6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.iNode.ID*6 + (a - 12)
                else:
                    # Find the corresponding index 'm' in the global fixed end reaction vector
                    m = quad.jNode.ID*6 + (a - 18)
                
                # Now that 'm' is known, place the term in the global fixed end reaction vector
                FER[m, 0] += quad_FER[a, 0]

        # Return the global fixed end reaction vector
        return FER
    
#%%
    def P(self, combo_name='Combo 1'):
        '''
        Assembles and returns the global nodal force vector.

        Parameters
        ----------
        combo_name : string
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

#%%
    def D(self, combo_name='Combo 1'):
        '''
        Returns the global displacement vector for the model.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the displacements for (not the load combination itself).
        '''
 
        # Return the global displacement vector
        return self.__D[combo_name]

#%%
    def _partition(self, unp_matrix, D1_indices, D2_indices):
        '''
        Partitions a matrix (or vector) into submatrices (or subvectors) based on degree of freedom
        boundary conditions.

        Parameters
        ----------
        unp_matrix : array
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

#%%
    def Analyze(self, log=False, check_statics=False, max_iter=30, sparse=True):
        warnings.warn('`Analyze` will be replaced with `analyze` in a future version of PyNite.', FutureWarning)
        self.analyze(log, check_statics, max_iter, sparse)

    def analyze(self, log=False, check_statics=False, max_iter=30, sparse=True):
        '''
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
        '''

        if log:
            print('+-----------+')
            print('| Analyzing |')
            print('+-----------+')

        # Assign an ID to all nodes and elements in the model
        self._renumber()

        # Ensure there is at least 1 load combination to solve if the user didn't define any
        if self.LoadCombos == {}:
            # Create and add a default load combination to the dictionary of load combinations
            self.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})
        
        # Activate all springs and members for all load combinations
        for spring in self.Springs.values():
            for combo_name in self.LoadCombos.keys():
                spring.active[combo_name] = True

        for member in self.Members.values():
            for combo_name in self.LoadCombos.keys():
                member.active[combo_name] = True 

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
                K11, K12, K21, K22 = self._partition(self.K(combo.name, log), D1_indices, D2_indices)

                # Get the partitioned global fixed end reaction vector
                FER1, FER2 = self._partition(self.FER(combo.name), D1_indices, D2_indices)

                # Get the partitioned global nodal force vector       
                P1, P2 = self._partition(self.P(combo.name), D1_indices, D2_indices)          

                # Calculate the global displacement vector
                if log: print('- Calculating global displacement vector for load combination', combo.name)
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
                            # The partitioned stiffness matrix is in `lil` format. It will be
                            # converted to a 2D dense array for mathematical operations.
                            D1 = solve(K11.toarray(), subtract(subtract(P1, FER1), matmul(K12.toarray(), D2)))
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
                self.__D[combo.name] = D

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
                    raise Exception('- Model diverged during tension/compression-only analysis')
                
                # Assume the model has converged (to be checked below)
                convergence = True

                # Check tension-only and compression-only springs
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

                # Check tension-only and compression-only members
                if log: print('- Checking for tension/compression-only member convergence')
                for member in self.Members.values():

                    # Only run the tension/compression only check if the member is still active
                    if member.active[combo.name] == True:

                        # Check if tension-only conditions exist
                        if member.tension_only == True and member.max_axial(combo.name) > 0:
                            member.active[combo.name] = False
                            convergence = False

                        # Check if compression-only conditions exist
                        elif member.comp_only == True and member.min_axial(combo.name) < 0:
                            member.active[combo.name] = False
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

#%%
    def Analyze_PDelta(self, log=False, max_iter=30, tol=0.01, sparse=True):
        warnings.warn('`Analyze_PDelta` will be replaced with `analyze_P_Delta` in a future version of PyNite.', FutureWarning)
        self.analyze_PDelta(log, max_iter, tol, sparse)

    def analyze_PDelta(self, log=False, max_iter=30, tol=0.01, sparse=True):
        '''
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
        '''
        
        if log:
            print('+--------------------+')
            print('| Analyzing: P-Delta |')
            print('+--------------------+')

        # Assign an ID to all nodes and elements in the model
        self._renumber()

        # Ensure there is at least 1 load combination to solve if the user didn't define any
        if self.LoadCombos == {}:

            # Create and add a default load combination to the dictionary of load combinations
            self.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})
    
        # Activate all springs and members for all load combinations. They can be turned inactive
        # during the course of the tension/compression-only analysis
        for spring in self.Springs.values():
            for combo_name in self.LoadCombos.keys():
                spring.active[combo_name] = True

        for member in self.Members.values():
            for combo_name in self.LoadCombos.keys():
                member.active[combo_name] = True 
        
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

                # Get the partitioned global matrices
                if iter_count_PD == 1:
                    
                    K11, K12, K21, K22 = self._partition(self.K(combo.name, log), D1_indices, D2_indices) # Initial stiffness matrix
                    FER1, FER2 = self._partition(self.FER(combo.name), D1_indices, D2_indices)       # Fixed end reactions
                    P1, P2 = self._partition(self.P(combo.name), D1_indices, D2_indices)             # Nodal forces

                else:

                    # Calculate the global stiffness matrices (partitioned)
                    K11, K12, K21, K22 = self._partition(self.K(combo.name, log), D1_indices, D2_indices)      # Initial stiffness matrix
                    Kg11, Kg12, Kg21, Kg22 = self._partition(self.Kg(combo.name, log), D1_indices, D2_indices) # Geometric stiffness matrix

                    # Combine the partitioned stiffness matrices. They are currently `lil` format
                    # which is great for memory, but slow for mathematical operations. They will be
                    # converted to `csr` format. The `+` operator performs matrix addition on `csr`
                    # matrices.
                    K11 = K11.tocsr() + Kg11.tocsr()
                    K12 = K12.tocsr() + Kg12.tocsr()
                    K21 = K21.tocsr() + Kg21.tocsr()
                    K22 = K22.tocsr() + Kg22.tocsr()                

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
                            D1 = solve(K11.toarray(), subtract(subtract(P1, FER1), matmul(K12.toarray(), D2)))

                    except:
                        # Return out of the method if 'K' is singular and provide an error message
                        raise ValueError('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')
            
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
                self.__D[combo.name] = D

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
                for member in self.Members.values():

                    # Only run the tension/compression only check if the member is still active
                    if member.active[combo.name] == True:

                        # Check if tension-only conditions exist
                        if member.tension_only == True and member.max_axial(combo.name) > 0:
                            
                            member.active[combo.name] = False
                            convergence_TC = False

                            # Reset the P-Delta analysis for the new geometry
                            iter_count_PD = 0
                            convergence_PD = False

                        # Check if compression-only conditions exist
                        elif member.comp_only == True and member.min_axial(combo.name) < 0:
                            
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
                    if log: print('- Checking for convergence')

                    # Temporarily disable error messages for invalid values.
                    # We'll be dealing with some 'nan' values due to division by zero at supports with zero deflection.
                    seterr(invalid='ignore')

                    # Check for convergence
                    if abs(1 - nanmax(divide(prev_results, D1))) <= tol:
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

#%%
    def _calc_reactions(self, log=False):
        '''
        Calculates reactions once the model is solved.

        Parameters
        ----------
        log : bool, optional
            Prints updates to the console if set to True. Default is False.
        '''

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
                if (node.SupportDX == True) \
                or (node.SupportDY == True) \
                or (node.SupportDZ == True) \
                or (node.SupportRX == True) \
                or (node.SupportRY == True) \
                or (node.SupportRZ == True):

                    # Sum the spring end forces at the node
                    for spring in self.Springs.values():

                        if spring.iNode == node and spring.active[combo.name] == True:
                            
                            # Get the spring's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            spring_F = spring.F(combo.name)

                            node.RxnFX[combo.name] += spring_F[0, 0]
                            node.RxnFY[combo.name] += spring_F[1, 0]
                            node.RxnFZ[combo.name] += spring_F[2, 0]
                            node.RxnMX[combo.name] += spring_F[3, 0]
                            node.RxnMY[combo.name] += spring_F[4, 0]
                            node.RxnMZ[combo.name] += spring_F[5, 0]

                        elif spring.jNode == node and spring.active[combo.name] == True:
                        
                            # Get the spring's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            spring_F = spring.F(combo.name)
                        
                            node.RxnFX[combo.name] += spring_F[6, 0]
                            node.RxnFY[combo.name] += spring_F[7, 0]
                            node.RxnFZ[combo.name] += spring_F[8, 0]
                            node.RxnMX[combo.name] += spring_F[9, 0]
                            node.RxnMY[combo.name] += spring_F[10, 0]
                            node.RxnMZ[combo.name] += spring_F[11, 0]

                    # Sum the member end forces at the node
                    for member in self.Members.values():
                    
                        if member.iNode == node and member.active[combo.name] == True:
                        
                            # Get the member's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            member_F = member.F(combo.name)

                            node.RxnFX[combo.name] += member_F[0, 0]
                            node.RxnFY[combo.name] += member_F[1, 0]
                            node.RxnFZ[combo.name] += member_F[2, 0]
                            node.RxnMX[combo.name] += member_F[3, 0]
                            node.RxnMY[combo.name] += member_F[4, 0]
                            node.RxnMZ[combo.name] += member_F[5, 0]

                        elif member.jNode == node and member.active[combo.name] == True:
                        
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

                        if plate.iNode == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[0, 0]
                            node.RxnFY[combo.name] += plate_F[1, 0]
                            node.RxnFZ[combo.name] += plate_F[2, 0]
                            node.RxnMX[combo.name] += plate_F[3, 0]
                            node.RxnMY[combo.name] += plate_F[4, 0]
                            node.RxnMZ[combo.name] += plate_F[5, 0]

                        elif plate.nNode == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[6, 0]
                            node.RxnFY[combo.name] += plate_F[7, 0]
                            node.RxnFZ[combo.name] += plate_F[8, 0]
                            node.RxnMX[combo.name] += plate_F[9, 0]
                            node.RxnMY[combo.name] += plate_F[10, 0]
                            node.RxnMZ[combo.name] += plate_F[11, 0]

                        elif plate.mNode == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[12, 0]
                            node.RxnFY[combo.name] += plate_F[13, 0]
                            node.RxnFZ[combo.name] += plate_F[14, 0]
                            node.RxnMX[combo.name] += plate_F[15, 0]
                            node.RxnMY[combo.name] += plate_F[16, 0]
                            node.RxnMZ[combo.name] += plate_F[17, 0]

                        elif plate.jNode == node:

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

                        if quad.mNode == node:

                            # Get the quad's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            quad_F = quad.F(combo.name)
                    
                            node.RxnFX[combo.name] += quad_F[0, 0]
                            node.RxnFY[combo.name] += quad_F[1, 0]
                            node.RxnFZ[combo.name] += quad_F[2, 0]
                            node.RxnMX[combo.name] += quad_F[3, 0]
                            node.RxnMY[combo.name] += quad_F[4, 0]
                            node.RxnMZ[combo.name] += quad_F[5, 0]

                        elif quad.nNode == node:

                            # Get the quad's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            quad_F = quad.F(combo.name)
                    
                            node.RxnFX[combo.name] += quad_F[6, 0]
                            node.RxnFY[combo.name] += quad_F[7, 0]
                            node.RxnFZ[combo.name] += quad_F[8, 0]
                            node.RxnMX[combo.name] += quad_F[9, 0]
                            node.RxnMY[combo.name] += quad_F[10, 0]
                            node.RxnMZ[combo.name] += quad_F[11, 0]

                        elif quad.iNode == node:

                            # Get the quad's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            quad_F = quad.F(combo.name)
                    
                            node.RxnFX[combo.name] += quad_F[12, 0]
                            node.RxnFY[combo.name] += quad_F[13, 0]
                            node.RxnFZ[combo.name] += quad_F[14, 0]
                            node.RxnMX[combo.name] += quad_F[15, 0]
                            node.RxnMY[combo.name] += quad_F[16, 0]
                            node.RxnMZ[combo.name] += quad_F[17, 0]

                        elif quad.jNode == node:

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

#%%
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
