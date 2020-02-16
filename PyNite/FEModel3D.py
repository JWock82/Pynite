# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 21:11:20 2017

@author: D. Craig Brinck, SE
"""
# %%
from numpy import matrix, zeros, delete, insert, matmul, divide, add, subtract, nanmax, seterr, shape
from numpy.linalg import inv, matrix_rank
from math import isclose
from PyNite.Node3D import Node3D
from PyNite.Member3D import Member3D
from PyNite.Plate3D import Plate3D

# %%
class FEModel3D():
    """
    A class representing a 3D finite element model.
    """

#%% 
    def __init__(self):
        """
        Initializes a new 3D finite element model.
        """
        
        self.Nodes = []    # A list of the structure's nodes
        self.auxNodes = [] # A list of the structure's auxiliary nodes
        self.Members = []  # A list of the structure's members
        self.Plates = []   # A list of the structure's plates
        self.__D = []      # A list of the structure's nodal displacements

#%%
    def AddNode(self, Name, X, Y, Z):
        """
        Adds a new node to the model.
        
        Parameters
        ----------
        Name : string
            A unique user-defined name for the node.
        X : number
            The global X-coordinate of the node.
        Y : number
            The global Y-coordinate of the node.
        Z : number
            The global Z-coordinate of the node.
        """
        
        # Create a new node
        newNode = Node3D(Name, X, Y, Z)
        
        # Add the new node to the list
        self.Nodes.append(newNode)
        

    def AddAuxNode(self, Name, X, Y, Z):
        """
        Adds a new auxiliary node to the model.
        
        Parameters
        ----------
        Name : string
            A unique user-defined name for the node.
        X : number
            The global X-coordinate of the node.
        Y : number
            The global Y-coordinate of the node.
        Z : number
            The global Z-coordinate of the node.
        """
        
        # Create a new node
        newNode = Node3D(Name, X, Y, Z)
        
        # Add the new node to the list
        self.auxNodes.append(newNode)
        
        
#%%
    def AddMember(self, Name, iNode, jNode, E, G, Iy, Iz, J, A, auxNode=None):
        """
        Adds a new member to the model.
        
        Parameters
        ----------
        Name : string
            A unique user-defined name for the member.
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
        """
        
        # Create a new member
        if auxNode == None:
            newMember = Member3D(Name, self.GetNode(iNode), self.GetNode(jNode), E, G, Iy, Iz, J, A)
        else:
            newMember = Member3D(Name, self.GetNode(iNode), self.GetNode(jNode), E, G, Iy, Iz, J, A, self.GetAuxNode(auxNode))
        
        # Add the new member to the list
        self.Members.append(newMember)

#%%
    def AddPlate(self, Name, iNode, jNode, mNode, nNode, t, E, nu):
        """
        Adds a new plate to the model.
        
        Parameters
        ----------
        Name : string
            A unique user-defined name for the plate.
        iNode : string
            The name of the i-node (1st node definded in counter-clockwise order).
        jNode : string
            The name of the j-node (2nd node defined in counter-clockwise order).
        mNode : string
            The name of the m-node (3rd node defined in counter-clockwise order).
        nNode : string
            The name of the n-node (4th node defined in counter-clockwise order).
        t : number
            The thickness of the plate.
        E : number
            The modulus of elasticity of the plate.
        mew : number
            Posson's ratio for the plate.
        """
        
        # Create a new member
        newPlate = Plate3D(Name, self.GetNode(iNode), self.GetNode(jNode), self.GetNode(mNode), self.GetNode(nNode), t, E, nu)
        
        # Add the new member to the list
        self.Plates.append(newPlate)

#%%
    def RemoveNode(self, Node):
        """
        Removes a node from the model. All nodal loads associated with the
        node and members attached to the node will also be removed.
        
        Parameters
        ----------
        Node : string
            The name of the node to be removed.
        """
        
        # Remove the node. Nodal loads are stored within the node, so they
        # will be deleted automatically when the node is deleted.
        self.Nodes.remove(self.GetNode(Node))
        
        # Find any members attached to the node and remove them
        self.Members = [member for member in self.Members if member.iNode.Name != Node and member.jNode.Name != Node]
        
#%%
    def RemoveMember(self, Member):
        """
        Removes a member from the model. All member loads associated with the
        member will also be removed.
        
        Parameters
        ----------
        Member : string
            The name of the member to be removed.
        """
        
        # Remove the member. Member loads are stored within the member, so they
        # will be deleted automatically when the member is deleted.
        self.Members.remove(self.GetMember(Member))
        
#%%
    def DefineSupport(self, Node, SupportDX=False, SupportDY=False, SupportDZ=False, SupportRX=False, SupportRY=False, SupportRZ=False):
        """
        Defines the support conditions at a node.
        
        Nodes will default to fully unsupported unless specified otherwise.
        
        Parameters
        ----------
        Node : string
            The name of the node where the support is being defined
        SupportDX : number
            Indicates whether the node is supported against translation in the global X-direction or if there is a support's settlement.
        SupportDY : number
            Indicates whether the node is supported against translation in the global Y-direction or if there is a support's settlement.
        SupportDZ : number
            Indicates whether the node is supported against translation in the global Z-direction or if there is a support's settlement.
        SupportRX : number
            Indicates whether the node is supported against rotation about the global X-axis or if there is a support's settlement.
        SupportRY : number
            Indicates whether the node is supported against rotation about the global Y-axis or if there is a support's settlement.
        SupportRZ : number
            Indicates whether the node is supported against rotation about the global Z-axis or if there is a support's settlement.
        """
        
        # Get the node to be supported
        node = self.GetNode(Node)
                
        # Set the node's support conditions
        if SupportDX == True:
            node.DX = 0.0
        elif SupportDX != False:
            node.DX = SupportDX
        
        if SupportDY == True:
            node.DY = 0.0
        elif SupportDY != False:
            node.DY = SupportDY

        if SupportDZ == True:
            node.DZ = 0.0
        elif SupportDZ != False:
            node.DZ = SupportDZ

        if SupportRX == True:
            node.RX = 0.0
        elif SupportRX != False:
            node.RX = SupportRX

        if SupportRY == True:
            node.RY = 0.0
        elif SupportRY != False:
            node.RY = SupportRY

        if SupportRZ == True:
            node.RZ = 0.0
        elif SupportRZ != False:
            node.RZ = SupportRZ

#%%            
    def AddNodalDisplacement (self, Node, Direction, Magnitude): 
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

        # Get the node
        node = self.GetNode(Node)

        if Direction == 'DX':
            node.DX = value
        if Direction == 'DY':
            node.DY = value
        if Direction == 'DZ':
            node.DZ = value
        if Direction == 'RX':
            node.RX = value
        if Direction == 'RY':
            node.RY = value
        if Direction == 'RZ':
            node.RZ = value

#%%
    def DefineReleases(self, Member, Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=False, Rzi=False, Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=False, Rzj=False):
        """
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
        """
        
        # Apply the end releases to the member
        self.GetMember(Member).Releases = [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]     
            
#%%
    def AddNodeLoad(self, Node, Direction, P):
        """
        Adds a nodal load to the model.
        
        Parameters
        ----------
        Node : string
            The name of the node where the load is being applied.
        Direction : {'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'}
            The global direction the load is being applied in. Forces are 'FX', 'FY', and 'FZ'. Moments are 'MX', 'MY', and 'MZ'.
        P : number
            The numeric value (magnitude) of the load.
        """
        
        # Add the node load to the model
        self.GetNode(Node).NodeLoads.append((Direction, P))

#%%      
    def AddMemberPtLoad(self, Member, Direction, P, x):
        """
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
        """
        
        # Add the point load to the member
        self.GetMember(Member).PtLoads.append((Direction, P, x))

#%%
    def AddMemberDistLoad(self, Member, Direction, w1, w2, x1=None, x2=None):
        """
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
        """
        
        # Determine if a starting and ending points for the load have been specified.
        # If not, use the member start and end as defaults
        if x1 == None:
            start = 0
        else:
            start = x1
        
        if x2 == None:
            end = self.GetMember(Member).L()
        else:
            end = x2

        # Add the distributed load to the member
        self.GetMember(Member).DistLoads.append((Direction, w1, w2, start, end))

#%%
    def ClearLoads(self):
        """
        Clears all loads from the model along with any results based on the loads.
        """

        # Clear out the member loads and the calculated internal forces
        for member in self.Members:
            member.DistLoads = []
            member.PtLoads = []
            member.SegmentsZ = []
            member.SegmentsY = []
            member.SegmentsX = []
        
        # Clear out the nodal loads and the nodal displacements
        for node in self.Nodes:
            node.NodeLoads = []
            node.SupportDX = None
            node.SupportDY = None
            node.SupportDZ = None
            node.SupportRX = None
            node.SupportRY = None
            node.SupportRZ = None       

#%%
    def GetNode(self, Name):
        """
        Returns the node with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the node to be returned.
        """
        
        # Step through each node in the 'Nodes' list
        for node in self.Nodes:
            
            # Check the name of the node
            if node.Name == Name:
                
                # Return the node of interest
                return node

            
    def GetAuxNode(self, Name):
        """
        Returns the auxiliary node with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the auxiliary node to be returned.
        """
        
        # Step through each node in the 'Nodes' list
        for node in self.auxNodes:
            
            # Check the name of the node
            if node.Name == Name:
                
                # Return the node of interest
                return node            
            
#%%
    def GetMember(self, Name):
        """
        Returns the member with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the member to be returned.
        """
        
        # Step through each member in the 'Members' list
        for member in self.Members:
            
            # Check the name of the member
            if member.Name == Name:
                
                # Return the member of interest
                return member

#%%
    def GetPlate(self, Name):
        """
        Returns the plate with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the plate to be returned.
        """
        
        # Step through each plate in the 'Plates' list
        for plate in self.Plates:
            
            # Check the name of the plate
            if plate.Name == Name:
                
                # Return the plate of interest
                return plate

#%%
    def __Renumber(self):
        """
        Assigns node, plate, and member ID numbers to be used internally by the
        program. Numbers are assigned according to the order nodes, members, and plates
        were added to the model.
        
        """
        
        # Number each node in the model
        i = 0
        for node in self.Nodes:
            node.ID = i
            i += 1
        
        # Number each member in the model
        i = 0
        for member in self.Members:
            member.ID = i
            i += 1
        
        # Number each plate in the model
        i = 0
        for plate in self.Plates:
            plate.ID = i
            i += 1

#%%
    def __AuxTable(self):
        '''
        Builds a table with known nodal displacements and with the positions in global stiffness matrix of known 
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
        for node in self.Nodes:
            
            # Known displacement
            if node.DX != None:
                D2_indices.append((node.ID*6) + 0)
                D2.append(node.DX)
            # Unknown displacement
            else:
                D1_indices.append((node.ID*6) + 0)

            # Known displacement
            if node.DY != None:
                D2_indices.append((node.ID*6) + 1)
                D2.append(node.DY)
            # Unknown displacement
            else:
                D1_indices.append((node.ID*6) + 1)

            # Known displacement
            if node.DZ != None:
                D2_indices.append((node.ID*6) + 2)
                D2.append(node.DZ)
            # Unknown displacement
            else:
                D1_indices.append((node.ID*6) + 2)

            # Known displacement
            if node.RX != None:
                D2_indices.append((node.ID*6) + 3)
                D2.append(node.RX)
            # Unknown displacement
            else:
                D1_indices.append((node.ID*6) + 3)

            # Known displacement
            if node.RY != None:
                D2_indices.append((node.ID*6) + 4)
                D2.append(node.RY)
            # Unknown displacement
            else:
                D1_indices.append((node.ID*6) + 4)

            # Known displacement
            if node.RZ != None:
                D2_indices.append((node.ID*6) + 5)
                D2.append(node.RZ)
            # Unknown displacement
            else:
                D1_indices.append((node.ID*6) + 5)

        # Return the indices and the known displacements
        return D1_indices, D2_indices, D2
            
#%%    
    def K(self):
        '''
        Assembles and returns the global stiffness matrix.
        '''
        
        # Initialize a zero matrix to hold all the stiffness terms
        K = zeros((len(self.Nodes)*6, len(self.Nodes)*6))
        
        # Add stiffness terms for each member in the model
        print('...Adding member stiffness terms to global stiffness matrix')
        for member in self.Members:
            
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
                    K.itemset((m, n), K.item((m, n)) + member_K.item((a, b)))
        
        # Add stiffness terms for each plate in the model
        print('...Adding plate stiffness terms to global stiffness matrix')
        for plate in self.Plates:
            
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
                    m = plate.jNode.ID*6 + (a-6)
                elif a < 18:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.mNode.ID*6 + (a-12)
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = plate.nNode.ID*6 + (a-18)

                for b in range(24):

                    # Determine which node the index 'b' is related to
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.iNode.ID*6 + b
                    elif b < 12:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.jNode.ID*6 + (b-6)
                    elif b < 18:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.mNode.ID*6 + (b-12)
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = plate.nNode.ID*6 + (b-18)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                    K.itemset((m, n), K.item((m, n)) + plate_K.item((a, b)))

        # Return the global stiffness matrix
        return K

#%%
    def __K_Partition(self, K):  
        '''
        Partitions global stiffness matrix in preparation for analysis
        '''
        
        # Get the auxiliary table to help with partitioning the matrix by known & unknown DOFs
        D1_indices, D2_indices, D2 = self.__AuxTable()

        # Initialize each partitioned matrix     
        K11 = zeros((K.shape[0] - len(D2), K.shape[0] - len(D2)))
        K12 = zeros((K.shape[0] - len(D2), len(D2)))
        K21 = zeros((len(D2), K.shape[0] - len(D2)))
        K22 = zeros((len(D2), len(D2)))

        # Initialize variables used to track rows and columns as the matrix is partitioned

        # m = Row
        # n = Column
        # 1 = Unknown DOF
        # 2 = Known DOF

        m11, n11 = 0, 0
        m12, n12 = 0, 0
        m21, n21 = 0, 0
        m22, n22 = 0, 0

        for m in range(K.shape[0]):

            for n in range(K.shape[0]):

                # K11 ---> DOF in row and column are both unknown
                if D2_indices.count(m) == 0 and D2_indices.count(n) == 0:
                    
                    K11.itemset((m11, n11), K[m, n])
                    n11 += 1
                    if n11 == K.shape[0] - len(D2):
                        n11 = 0
                        m11 += 1

                # K12 ---> DOF in row is unknown, DOF in column is known
                elif D2_indices.count(m) == 0 and D2_indices.count(n) == 1:
                    
                    K12.itemset((m12, n12), K[m, n])
                    n12 += 1
                    if n12 == len(D2):
                        n12 = 0
                        m12 += 1

                # K21 ---> DOF in row is known, DOF in column is unknown
                elif D2_indices.count(m) == 1 and D2_indices.count(n) == 0:
                    
                    K21.itemset((m21, n21), K[m, n])
                    n21 += 1
                    if n21 == K.shape[0] - len(D2):
                        n21 = 0
                        m21 += 1

                # K22 --> DOF in row and column are both known 
                elif D2_indices.count(m) == 1 and D2_indices.count(n) == 1:
                    
                    K22.itemset((m22, n22), K[m, n])
                    n22 += 1
                    if n22 == len(D2):
                        n22 = 0
                        m22 += 1

        # Return the matrix partitioned into 4 submatrices
        return K11, K12, K21, K22       

#%%    
    def Kg(self):
        '''
        Assembles and returns the global geometric stiffness matrix.

        The model must have a static solution prior to obtaining the geometric stiffness matrix.
        Geometric stiffness of plates is not included.
        '''
        
        # Initialize a zero matrix to hold all the stiffness terms
        Kg = zeros((len(self.Nodes) * 6, len(self.Nodes) * 6))
        
        # Add stiffness terms for each member in the model
        print('...Adding member geometric stiffness terms to global geometric stiffness matrix')
        for member in self.Members:
            
            # Calculate the axial force in the member
            E = member.E
            A = member.A
            L = member.L()
            d = member.d()
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
                    Kg.itemset((m, n), Kg.item((m, n)) + member_Kg.item((a, b)))

        # Return the global geometric stiffness matrix
        return Kg
     
#%%    
    def FER(self):
        '''
        Assembles and returns the global fixed end reaction vector.
        '''
        
        # Initialize a zero vector to hold all the terms
        FER = zeros((len(self.Nodes) * 6, 1))
        
        # Add terms for each member in the model
        for member in self.Members:
            
            # Get the member's global fixed end reaction vector
            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed
            member_FER = member.FER()

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
                FER.itemset((m, 0), FER[m, 0] + member_FER[a, 0])
        
        # Return the global fixed end reaction vector
        return FER

#%%    
    def __FER_Partition(self):
        '''
        Partitions global fixed end reaction vector prior to analysis
        '''
    
        D1_indices, D2_indices, D2 = self.__AuxTable()        
        FER = self.FER()

        # Initialize each partitioned matrix          
        FER1 = zeros((FER.shape[0] - len(D2), 1))
        FER2 = zeros((len(D2), 1))

        # Variables used to track indexing during partitioning

        # m = Row
        # n = Column
        # 1 = Unknown DOF
        # 2 = Known DOF

        m1 = 0
        m2 = 0

        for m in range(FER.shape[0]):
            
            #FER1 --> unknown dof
            if D2_indices.count(m) == 0:
                FER1.itemset((m1, 0), FER[m, 0])
                m1 += 1
            #FER2 --> known dof 
            else:
                FER2.itemset((m2, 0), FER[m, 0])
                m2 += 1

        # Return the vector partitioned as 2 subvectors
        return FER1, FER2
    
#%%
    def P(self):
        '''
        Assembles and returns the global nodal force vector.
        '''
            
        # Initialize a zero vector to hold all the terms
        Pvector = zeros((len(self.Nodes)*6, 1))
        
        # Add terms for each node in the model
        for node in self.Nodes:
            
            # Get the node's ID
            ID = node.ID
            
            # Add the node's loads to the global nodal load vector
            for load in node.NodeLoads:
                
                if load[0] == 'FX':
                    Pvector.itemset((ID*6 + 0, 0), Pvector[ID*6 + 0, 0] + load[1])
                elif load[0] == 'FY':
                    Pvector.itemset((ID*6 + 1, 0), Pvector[ID*6 + 1, 0] + load[1])
                elif load[0] == 'FZ':
                    Pvector.itemset((ID*6 + 2, 0), Pvector[ID*6 + 2, 0] + load[1])
                elif load[0] == 'MX':
                    Pvector.itemset((ID*6 + 3, 0), Pvector[ID*6 + 3, 0] + load[1])
                elif load[0] == 'MY':
                    Pvector.itemset((ID*6 + 4, 0), Pvector[ID*6 + 4, 0] + load[1])
                elif load[0] == 'MZ':
                    Pvector.itemset((ID*6 + 5, 0), Pvector[ID*6 + 5, 0] + load[1])
        
        # Return the global nodal force vector
        return Pvector
    
#%%
    def __P_Partition(self):
        '''
        Partitions the global nodal force vector prior to analysis 
        '''
        
        D1_indices, D2_indices, D2= self.__AuxTable()             
        P = self.P()

        # Initialize each partitioned matrix         
        P1 = zeros((P.shape[0] - len(D2), 1))
        P2 = zeros((len(D2), 1))

        # Variables used to track indexing during partitioning

        # m = row
        # n = column
        
        # 1 = unknown dof
        # 2 = known dof

        m1 = 0
        m2 = 0

        for m in range(P.shape[0]):

                #P1 --> unknown dof 
                if D2_indices.count(m) == 0:

                    P1.itemset((m1, 0), P[m, 0])
                    m1 += 1

                #P2 --> known dof 
                else:

                    P2.itemset((m2, 0), P[m, 0])
                    m2 += 1

        # Return the vector partitioned as 2 subvectors
        return P1, P2
    
#%%
    def D(self):
        """
        Returns the global displacement vector for the model.
        """
        
        # Return the global displacement vector
        return self.__D
        
#%%  
    def Analyze(self, check_statics=True):
        '''
        Analyzes the model.
        '''
        
        print('**Analyzing**')

        # Assign an ID to all nodes and elements in the model
        self.__Renumber()

        # Get the partitioned global stiffness matrix K11, K12, K21, K22
        K11, K12, K21, K22 = self.__K_Partition(self.K())
        
        # Get the partitioned global fixed end reaction vector
        FER1, FER2 = self.__FER_Partition()
        
        # Get the partitioned global nodal force vector
        P1, P2 = self.__P_Partition()
        
        # The global displacement vector is composed of two vectors, D1 and D2:

        # D1: Unknown displacements
        # D2: Known displacements
        
        # Get vector D2 wich contains known nodal displacements (D != None)
        D1_indices, D2_indices, D2 = self.__AuxTable()             

        # Convert D2 from a list to a matrix
        D2 = matrix(D2).T

        # Check for global stability by determining if 'K11' is singular
        print('...Checking global stability')
        if matrix_rank(K11) < min(K11.shape):
            # Return out of the method if 'K' is singular and provide an error message
            print('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')
            return
        else:
            # Calculate the unknown displacements D1
            print('...Calculating global displacement vector')
            D1 = matmul(inv(K11), subtract(subtract(P1, FER1), matmul(K12, D2)))

        # Form the global displacement vector, D, from D1 and D2
        D = zeros((len(self.Nodes)*6, 1))

        for node in self.Nodes:

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

        # Save the displacements as a local variable for easier reference below
        self.__D = D

        # Store the calculated global nodal displacements into each node
        for node in self.Nodes:
            node.DX = D.item((node.ID * 6 + 0, 0))
            node.DY = D.item((node.ID * 6 + 1, 0))
            node.DZ = D.item((node.ID * 6 + 2, 0))
            node.RX = D.item((node.ID * 6 + 3, 0))
            node.RY = D.item((node.ID * 6 + 4, 0))
            node.RZ = D.item((node.ID * 6 + 5, 0))
        
        # Segment all members in the model to make member results available
        print('...Calculating member internal forces')
        for member in self.Members:
            member.SegmentMember()
        
        # Calculate reactions
        self.__CalcReactions()
    
        # Check statics if requested
        if check_statics == True:
            self.__CheckStatics()

#%%
    def Analyze_PDelta(self, max_iter=30, tol=0.01):
        """
        Runs a second order (P-Delta) analysis on the structure.
        """
        print("**Running P-Delta analysis**")

        # Keep track of the number of iterations
        iter_count = 1
        convergence = False
        divergence = False

        # Iterate until convergence or divergence occurs
        while convergence == False and divergence == False:
            
            # Inform the user which iteration we're on
            print('...Beginning P-Delta iteration #' + str(iter_count))

            # Get the global stiffness matrix and renumber the nodes & members
            # in the process of creating it
            if iter_count == 1:
                K = self.K(True)      # Initial stiffness matrix
                FER = self.FER(False) # Fixed end reactions
                P = self.P(False)     # Nodal forces
            else:
                # Adjust the stiffness matrix
                K = add(self.K(False), self.Kg())        

            # Eliminate supported degrees of freedom from each of the matrices/vectors
            # Work backwards through the node list so that the relationship between
            # the DOF's and node ID's is unnafected by the matrices/vectors
            # shrinking
            for node in reversed(self.Nodes):
                
                if node.SupportRZ == True:
                    K = delete(K, node.ID * 6 + 5, axis = 0)
                    K = delete(K, node.ID * 6 + 5, axis = 1)
                    if iter_count == 1:
                        FER = delete(FER, node.ID * 6 + 5, axis = 0)
                        P = delete(P, node.ID * 6 + 5, axis = 0)
                
                if node.SupportRY == True:
                    K = delete(K, node.ID * 6 + 4, axis = 0)
                    K = delete(K, node.ID * 6 + 4, axis = 1)
                    if iter_count == 1:
                        FER = delete(FER, node.ID * 6 + 4, axis = 0)
                        P = delete(P, node.ID * 6 + 4, axis = 0)

                if node.SupportRX == True:
                    K = delete(K, node.ID * 6 + 3, axis = 0)
                    K = delete(K, node.ID * 6 + 3, axis = 1)
                    if iter_count == 1:
                        FER = delete(FER, node.ID * 6 + 3, axis = 0)
                        P = delete(P, node.ID * 6 + 3, axis = 0)

                if node.SupportDZ == True:
                    K = delete(K, node.ID * 6 + 2, axis = 0)
                    K = delete(K, node.ID * 6 + 2, axis = 1)
                    if iter_count == 1:
                        FER = delete(FER, node.ID * 6 + 2, axis = 0)
                        P = delete(P, node.ID * 6 + 2, axis = 0)

                if node.SupportDY == True:
                    K = delete(K, node.ID * 6 + 1, axis = 0)
                    K = delete(K, node.ID * 6 + 1, axis = 1)
                    if iter_count == 1:
                        FER = delete(FER, node.ID * 6 + 1, axis = 0)
                        P = delete(P, node.ID * 6 + 1, axis = 0)

                if node.SupportDX == True:
                    K = delete(K, node.ID * 6 + 0, axis = 0)
                    K = delete(K, node.ID * 6 + 0, axis = 1)
                    if iter_count == 1:
                        FER = delete(FER, node.ID * 6 + 0, axis = 0)
                        P = delete(P, node.ID * 6 + 0, axis = 0)
        
            # Determine if 'K' is singular
            print('...Checking global stability')
            if matrix_rank(K) < min(K.shape):
                # Return out of the method if 'K' is singular and provide an error message
                print('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')
                return
            else:
                # Calculate the global displacement vector
                print('...Calculating global displacement vector')
                D = matmul(inv(K), subtract(P, FER))
        
            # Expand the global displacement vector to include supported degrees of freedom
            # Work forwards through the node list so that the relationship between
            # the DOF's and node ID's is unnafected by the vector expanding
            for node in self.Nodes:
                if node.SupportDX == True:
                    D = insert(D, node.ID * 6 + 0, 0, axis = 0)
                if node.SupportDY == True:
                    D = insert(D, node.ID * 6 + 1, 0, axis = 0)
                if node.SupportDZ == True:
                    D = insert(D, node.ID * 6 + 2, 0, axis = 0)
                if node.SupportRX == True:
                    D = insert(D, node.ID * 6 + 3, 0, axis = 0)
                if node.SupportRY == True:
                    D = insert(D, node.ID * 6 + 4, 0, axis = 0)
                if node.SupportRZ == True:
                    D = insert(D, node.ID * 6 + 5, 0, axis = 0)

            # Save the expanded global displacement vector
            self.__D = D

            # Store the calculated global nodal displacements into each node
            for node in self.Nodes:
                node.DX = D.item((node.ID * 6 + 0, 0))
                node.DY = D.item((node.ID * 6 + 1, 0))
                node.DZ = D.item((node.ID * 6 + 2, 0))
                node.RX = D.item((node.ID * 6 + 3, 0))
                node.RY = D.item((node.ID * 6 + 4, 0))
                node.RZ = D.item((node.ID * 6 + 5, 0))
            
            if iter_count != 1:
                
                # Print a status update for the user
                print("...Checking for convergence.")

                # Temporarily disable error messages for invalid values.
                # We'll be dealing with some 'nan' values due to division by zero at supports with zero deflection.
                seterr(invalid='ignore')

                # Check for convergence
                if abs(1 - nanmax(divide(prev_results, D))) <= tol:
                    convergence = True
                    print("...P-Delta analysis converged after "+str(iter_count)+" iterations.")
                # Check for divergence
                elif iter_count > max_iter:
                    divergence = True
                    print("...P-Delta analysis failed to converge after 30 iterations.")

                # Turn invalid value warnings back on
                seterr(invalid='warn') 

            # Save the results for the next iteration
            prev_results = D

            # Increment the iteration count
            iter_count += 1
        
        # Calculate reactions
        self.__CalcReactions()
                
        # Segment all members in the model to make member results available
        print('...Calculating member internal forces')
        for member in self.Members:
            member.SegmentMember()

#%%
    def __CalcReactions(self):
        """
        Calculates reactions once the model is solved.
        """

        # Print a status update to the console
        print('...Calculating reactions')

        # Calculate the reactions, node by node
        for node in self.Nodes:
            
            # Determine if the node has any supports
            if isclose(node.DX, 0.0) or isclose(node.DY, 0.0) or isclose(node.DZ, 0.0) or isclose(node.RX, 0.0) or isclose(node.RY, 0) or isclose(node.RZ, 0):

                # Sum the member end forces at the node
                for member in self.Members:
                    
                    if member.iNode == node:
                        
                        # Get the member's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        member_F = member.F()

                        node.RxnFX += member_F[0, 0]
                        node.RxnFY += member_F[1, 0]
                        node.RxnFZ += member_F[2, 0]
                        node.RxnMX += member_F[3, 0]
                        node.RxnMY += member_F[4, 0]
                        node.RxnMZ += member_F[5, 0]

                    elif member.jNode == node:
                        
                        # Get the member's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        member_F = member.F()
                        
                        node.RxnFX += member_F[6, 0]
                        node.RxnFY += member_F[7, 0]
                        node.RxnFZ += member_F[8, 0]
                        node.RxnMX += member_F[9, 0]
                        node.RxnMY += member_F[10, 0]
                        node.RxnMZ += member_F[11, 0]

                # Sum the plate forces at the node
                for plate in self.Plates:

                    if plate.iNode == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F()
                    
                        node.RxnFX += plate_F[0, 0]
                        node.RxnFY += plate_F[1, 0]
                        node.RxnFZ += plate_F[2, 0]
                        node.RxnMX += plate_F[3, 0]
                        node.RxnMY += plate_F[4, 0]
                        node.RxnMZ += plate_F[5, 0]

                    elif plate.jNode == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F()
                    
                        node.RxnFX += plate_F[6, 0]
                        node.RxnFY += plate_F[7, 0]
                        node.RxnFZ += plate_F[8, 0]
                        node.RxnMX += plate_F[9, 0]
                        node.RxnMY += plate_F[10, 0]
                        node.RxnMZ += plate_F[11, 0]

                    elif plate.mNode == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F()
                    
                        node.RxnFX += plate_F[12, 0]
                        node.RxnFY += plate_F[13, 0]
                        node.RxnFZ += plate_F[14, 0]
                        node.RxnMX += plate_F[15, 0]
                        node.RxnMY += plate_F[16, 0]
                        node.RxnMZ += plate_F[17, 0]

                    elif plate.nNode == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F()
                    
                        node.RxnFX += plate_F[18, 0]
                        node.RxnFY += plate_F[19, 0]
                        node.RxnFZ += plate_F[20, 0]
                        node.RxnMX += plate_F[21, 0]
                        node.RxnMY += plate_F[22, 0]
                        node.RxnMZ += plate_F[23, 0]

                # Sum the joint forces at the node
                for load in node.NodeLoads:
                    
                    if load[0] == "FX":
                        node.RxnFX -= load[1]
                    elif load[0] == "FY":
                        node.RxnFY -= load[1]
                    elif load[0] == "FZ":
                        node.RxnFZ -= load[1]
                    elif load[0] == "MX":
                        node.RxnMX -= load[1]
                    elif load[0] == "MY":
                        node.RxnMY -= load[1]
                    elif load[0] == "MZ":
                        node.RxnMZ -= load[1]

#%%
    def __CheckStatics(self):
        """
        Sums global forces and reactions and prints them to the console.
        """
        
        # Print a status update to the console
        print('...Checking statics')

        # Initialize force summations to zero
        SumFX = 0
        SumFY = 0
        SumFZ = 0
        SumMX = 0
        SumMY = 0
        SumMZ = 0
        SumRFX = 0
        SumRFY = 0
        SumRFZ = 0
        SumRMX = 0
        SumRMY = 0
        SumRMZ = 0

        # Get the global force vector and the global fixed end reaction vector
        P = self.P()
        FER = self.FER()

        # Step through each node and sum its forces
        for node in self.Nodes:

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
            RFX = node.RxnFX
            RFY = node.RxnFY
            RFZ = node.RxnFZ
            RMX = node.RxnMX
            RMY = node.RxnMY
            RMZ = node.RxnMZ

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
        
        # Print the load summation
        print('**Applied Loads**')
        print('Sum Forces X: ', SumFX, ', Sum Forces Y: ', SumFY, ', Sum Forces Z: ', SumFZ)
        print('Sum Moments MX: ', SumMX, ', Sum Moments MY: ', SumMY, ', Sum Moments MZ: ', SumMZ)

        # Print the reaction summation
        print('**Reactions**')
        print('Sum Forces X: ', SumRFX, ', Sum Forces Y: ', SumRFY, ', Sum Forces Z: ', SumRFZ)
        print('Sum Moments MX: ', SumRMX, ', Sum Moments MY: ', SumRMY, ', Sum Moments MZ: ', SumRMZ)

        return SumFX, SumFY, SumFZ, SumMX, SumMY, SumMZ
