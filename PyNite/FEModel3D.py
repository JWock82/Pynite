# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 21:11:20 2017

@author: D. Craig Brinck, SE
"""
# %%
from numpy import zeros, delete, insert, matmul, subtract
from numpy.linalg import inv, matrix_rank
from PyNite.Node3D import Node3D
from PyNite.Member3D import Member3D

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
        
        self.Nodes = []           # A list of the structure's nodes
        self.Members = []         # A list of the structure's members
        self.__D = []               # A list of the structure's nodal displacements

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

#%%
    def AddMember(self, Name, iNode, jNode, E, G, Iy, Iz, J, A):
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
        newMember = Member3D(Name, self.GetNode(iNode), self.GetNode(jNode), E, G, Iy, Iz, J, A)
        
        # Add the new member to the list
        self.Members.append(newMember)

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
    def DefineSupport(self, Node, SupportDX = False, SupportDY = False, SupportDZ = False, SupportRX = False, SupportRY = False, SupportRZ = False):
        """
        Defines the support conditions at a node.
        
        Nodes will default to fully unsupported unless specified otherwise.
        
        Parameters
        ----------
        Node : string
            The name of the node where the support is being defined
        SupportDX : boolean
            Indicates whether the node is supported against translation in the global X-direction.
        SupportDY : boolean
            Indicates whether the node is supported against translation in the global Y-direction.
        SupportDZ : boolean
            Indicates whether the node is supported against translation in the global Z-direction.
        SupportRX : boolean
            Indicates whether the node is supported against rotation about the global X-axis.
        SupportRY : boolean
            Indicates whether the node is supported against rotation about the global Y-axis.
        SupportRZ : boolean
            Indicates whether the node is supported against rotation about the global Z-axis.
        """
        
        # Get the node to be supported
        node = self.GetNode(Node)
                
        # Set the node's supports
        node.SupportDX = SupportDX
        node.SupportDY = SupportDY
        node.SupportDZ = SupportDZ
        node.SupportRX = SupportRX
        node.SupportRY = SupportRY
        node.SupportRZ = SupportRZ

#%%
    def DefineReleases(self, Member, Dxi = False, Dyi = False, Dzi = False, Rxi = False, Ryi = False, Rzi = False, Dxj = False, Dyj = False, Dzj = False, Rxj = False, Ryj = False, Rzj = False):
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
        Direction : {'Fx', 'Fy', 'Fz', 'My', 'Mz'}
            The direction in which the force is to be applied. Note that
            typical beam sign convention is used. Transverse forces acting
            toward the beam are positive. Moments are positive if they act
            counter-clockwise relative to the beam's local coordinate system.
            Torsional point loads are not supported at this time.
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
            end = self.GetMember(Member).L
        else:
            end = x2

        # Add the distributed load to the member
        self.GetMember(Member).DistLoads.append((Direction, w1, w2, start, end))

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
    def __Renumber(self):
        """
        Assigns node and member ID numbers to be used internally by the
        program. Numbers are assigned according to the order nodes and members
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
            
#%%    
    def K(self, Renumber=False):
        """
        Assembles and returns the global stiffness matrix.
        
        Parameters
        ----------
        Renumber : boolean
            Indicates whether nodes and members should be renumbered prior to
            calculating the stiffness matrix. This may be necessary if a model
            is being solved for the first time, or if it has been changed since
            the last run, potentially creating a gap in the numbering.
        """
        
        # Renumber the nodes and members in the model if requested
        if Renumber == True:
            self.__Renumber()
        
        # Initialize a zero matrix to hold all the stiffness terms
        K = zeros((len(self.Nodes) * 6, len(self.Nodes) * 6))
        
        # Add stiffness terms for each member in the model
        for member in self.Members:
            
            # Step through each term in the member's stiffness matrix
            # 'a' & 'b' below are row/column indices in the member's stiffness matrix
            # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
            for a in range(12):
                
                # Determine if index 'a' is related to the i-node or j-node
                if a < 6:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = member.iNode.ID * 6 + a
                else:
                    # Find the corresponding index 'm' in the global stiffness matrix
                    m = member.jNode.ID * 6 + (a - 6)
                    
                for b in range(12):
                    
                    # Determine if index 'b' is related to the i-node or j-node
                    if b < 6:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = member.iNode.ID * 6 + b
                    else:
                        # Find the corresponding index 'n' in the global stiffness matrix
                        n = member.jNode.ID * 6 + (b - 6)
                    
                    # Now that 'm' and 'n' are known, place the term in the global stiffness matrix
                    K.itemset((m, n), K.item((m, n)) + member.K().item((a, b)))
        
        # Return the global stiffness matrix
        return K
    
#%%    
    def FER(self, Renumber=False):
        """
        Assembles and returns the global fixed end reaction vector.
        
        Parameters
        ----------
        Renumber : boolean
            Indicates whether nodes and members should be renumbered prior to
            calculating the fixed end reaction vector. This may be necessary if
            a model is being solved for the first time, or if it has been
            changed since the last run, potentially creating a gap in the
            numbering.
        """
        
        # Renumber the nodes and members in the model if requested
        if Renumber == True:
            self.__Renumber()
        
        # Initialize a zero vector to hold all the terms
        FER = zeros((len(self.Nodes) * 6, 1))
        
        # Add terms for each member in the model
        for member in self.Members:
            
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
                FER.itemset((m, 0), FER[m, 0] + member.FER()[a, 0])
        
        # Return the global fixed end reaction vector
        return FER
    
#%%
    def P(self, Renumber=False):
        """
        Assembles and returns the global nodal force vector.
        
        Parameters
        ----------
        Renumber : boolean
            Indicates whether nodes and members should be renumbered prior to
            calculating the fixed end reaction vector. This may be necessary if
            a model is being solved for the first time, or if it has been
            changed since the last run, potentially creating a gap in the
            numbering.
        """
        
        # Renumber the nodes and members in the model if requested
        if Renumber == True:
            self.__Renumber()
            
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
    def D(self):
        """
        Returns the global displacement vector for the model.
        """
        
        # Return the global displacement vector
        return self.__D
        
#%%  
    def Analyze(self, check_statics=True):
        """
        Analyzes the model.
        """
            
        # Get the global stiffness matrix and renumber the nodes & members
        # in the process of creating it
        K = self.K(True)
        
        # Get the global fixed end reaction vector
        FER = self.FER(False)
        
        # Get the global nodal force vector
        P = self.P(False)
        
        # Eliminate supported degrees of freedom from each of the matrices/vectors
        # Work backwards through the node list so that the relationship between
        # the DOF's and node ID's is unnafected by the matrices/vectors
        # shrinking
        for node in reversed(self.Nodes):
            
            if node.SupportRZ == True:
                K = delete(K, node.ID * 6 + 5, axis = 0)
                K = delete(K, node.ID * 6 + 5, axis = 1)
                FER = delete(FER, node.ID * 6 + 5, axis = 0)
                P = delete(P, node.ID * 6 + 5, axis = 0)
                
            if node.SupportRY == True:
                K = delete(K, node.ID * 6 + 4, axis = 0)
                K = delete(K, node.ID * 6 + 4, axis = 1)
                FER = delete(FER, node.ID * 6 + 4, axis = 0)
                P = delete(P, node.ID * 6 + 4, axis = 0)
                
            if node.SupportRX == True:
                K = delete(K, node.ID * 6 + 3, axis = 0)
                K = delete(K, node.ID * 6 + 3, axis = 1)
                FER = delete(FER, node.ID * 6 + 3, axis = 0)
                P = delete(P, node.ID * 6 + 3, axis = 0)
                
            if node.SupportDZ == True:
                K = delete(K, node.ID * 6 + 2, axis = 0)
                K = delete(K, node.ID * 6 + 2, axis = 1)
                FER = delete(FER, node.ID * 6 + 2, axis = 0)
                P = delete(P, node.ID * 6 + 2, axis = 0)
                
            if node.SupportDY == True:
                K = delete(K, node.ID * 6 + 1, axis = 0)
                K = delete(K, node.ID * 6 + 1, axis = 1)
                FER = delete(FER, node.ID * 6 + 1, axis = 0)
                P = delete(P, node.ID * 6 + 1, axis = 0)
                
            if node.SupportDX == True:
                K = delete(K, node.ID * 6 + 0, axis = 0)
                K = delete(K, node.ID * 6 + 0, axis = 1)
                FER = delete(FER, node.ID * 6 + 0, axis = 0)
                P = delete(P, node.ID * 6 + 0, axis = 0)
                        
        # Determine if 'K' is singular
        if matrix_rank(K) < min(K.shape):
            # Return out of the method if 'K' is singular and provide an error message
            print('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')
            return
        else:
            # Calculate the global displacement vector
            self.__D = matmul(inv(K), subtract(P, FER))
        
        # Save the displacements as a local variable for easier reference below
        D = self.__D
        
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

        # Store the calculated global nodal displacements into each node
        for node in self.Nodes:
            node.DX = D.item((node.ID * 6 + 0, 0))
            node.DY = D.item((node.ID * 6 + 1, 0))
            node.DZ = D.item((node.ID * 6 + 2, 0))
            node.RX = D.item((node.ID * 6 + 3, 0))
            node.RY = D.item((node.ID * 6 + 4, 0))
            node.RZ = D.item((node.ID * 6 + 5, 0))
        
        # Calculate and store the reactions at each node
        for node in self.Nodes:
            
            # Sum the member end forces at the node
            for member in self.Members:
                
                if member.iNode == node:
                    
                    node.RxnFX += member.F()[0, 0]
                    node.RxnFY += member.F()[1, 0]
                    node.RxnFZ += member.F()[2, 0]
                    node.RxnMX += member.F()[3, 0]
                    node.RxnMY += member.F()[4, 0]
                    node.RxnMZ += member.F()[5, 0]
                
                elif member.jNode == node:
                    
                    node.RxnFX += member.F()[6, 0]
                    node.RxnFY += member.F()[7, 0]
                    node.RxnFZ += member.F()[8, 0]
                    node.RxnMX += member.F()[9, 0]
                    node.RxnMY += member.F()[10, 0]
                    node.RxnMZ += member.F()[11, 0]
                
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
                
        # Segment all members in the model to make member results available
        for member in self.Members:
            member.SegmentMember()
        
        # Check statics if requested
        if check_statics == True:
            self.__CheckStatics()

#%%
    def __CheckStatics(self):

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
        P = self.P(False)
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
        print('**Reactions**')
        print('Sum Forces X: ', SumRFX, ', Sum Forces Y: ', SumRFY, ', Sum Forces Z: ', SumRFZ)
        print('Sum Moments MX: ', SumRMX, ', Sum Moments MY: ', SumRMY, ', Sum Moments MZ: ', SumRMZ)

        return SumFX, SumFY, SumFZ, SumMX, SumMY, SumMZ