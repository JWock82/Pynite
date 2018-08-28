# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 21:11:20 2017

@author: D. Craig Brinck, SE
"""
# %%
from numpy import zeros, delete, insert, matmul, subtract
from numpy.linalg import inv
from PyNite.Node3D import Node3D
from PyNite.Member3D import Member3D

# %%
class FEModel3D():
    """
    A 3D finite element model class
    """  

#%% 
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        
        self.Nodes = []     # A list of nodes in the structure
        self.Members = []   # A list of members in the structure
        self.NodeLoads = [] # A list of nodal loads in the structure
        self.D = []
        self.Solved = False

#%%
    # Adds a node to the finite element mdoel
    def AddNode(self, Name, X, Y, Z):
        """
        Adds a new node to the model
        
        Parameters
        ----------
        Name : string
            A unique name for the node
        X : number
            The global X-coordinate of the node
        Y : number
            The global Y-coordinate of the node
        Z : number
            The global Z-coordinate of the node
        """
        
        # Create the new node
        newNode = Node3D(Name, len(self.Nodes), X, Y, Z)
        
        # Add the new node to the 'Nodes' list
        self.Nodes.append(newNode)
        
        # Flag the model as unsolved due to the new node
        self.Solved = False

#%%
    # Adds a new member element to the finite element model
    def AddMember(self, Name, iNode, jNode, E, G, Iy, Iz, J, A):
        """
        Adds a new member to the finite element model
        
        Parameters
        ----------
        Name : string
            A unique name for the member
        iNode : string
            The name of the i-node (start node)
        jNode : string
            The name of the j-node (end node)
        E : number
            The modulus of elasticity of the member
        G : number
            The shear modulus of the member
        Iy : number
            The moment of inertia of the member about its local y-axis (usually the weak axis)
        Iz : number
            The moment of inertia of the member about its local z-axis (usually the strong axis)
        J : number
            The polar moment of inertia of the member
        A : number
            The area of the member
        """
        
        # Find the node object associated with the user's given node name
        for node in self.Nodes:
            if node.Name == iNode:
                iNode = node
            if node.Name == jNode:
                jNode = node
        
        # Create the new member
        newMember = Member3D(Name, len(self.Members), iNode, jNode, E, G, Iy, Iz, J, A)
        
        # Add the new member to the 'Members' list
        self.Members.append(newMember)
        
        # Flag the model as unsolved due to the new member
        self.Solved = False

#%%
    # Defines the support conditions at a node
    def DefineSupport(self, Node, SupportX=False, SupportY=False, SupportZ=False, SupportRX=False, SupportRY=False, SupportRZ=False):
        """
        Defines the support conditions at a node. Nodes default to unsupported unless specified otherwise.
        
        Parameters
        ----------
        Node : string
            The name of the node where the support is being defined
        SupportX : boolean
            Flag indicating whether the node is supported against translation in the global X direction
        SupportY : boolean
            Flag indicating whether the node is supported against translation in the global Y direction
        SupportZ : boolean
            Flag indicating whether the node is supported against translation in the global Z direction
        SupportRX : boolean
            Flag indicating whether the node is supported against rotation about the global X axis
        SupportRY : boolean
            Flag indicating whether the node is supported against rotation about the global Y axis
        SupportRZ : boolean
            Flag indicating whether the node is supported against rotation about the global Z axis
        """
        
        # Find the node object associated with the user's given node name
        for node in self.Nodes:
            
            if node.Name == Node:
                
                # Define the supports
                node.SupportDX = SupportX
                node.SupportDY = SupportY
                node.SupportDZ = SupportZ
                node.SupportRX = SupportRZ
                node.SupportRY = SupportRY
                node.SupportRZ = SupportRZ
                break
        
        # Flag the model as unsolved due to support changes
        self.Solved = False
        
#%%
    # Adds a nodal load to the model
    def AddNodeLoad(self, Node, Direction, Load):
        """
        Adds a nodal load to the model
        
        Parameters
        ----------
        Node : string
            The name of the node where the load is being applied
        Direction : string
            The direction the load is being applied in. Must be one of the following values:
            'FX' = Force in the global X-direction
            'FY' = Force in the global Y-direction
            'FZ' = Force in the global Z-direction
            'MX' = Moment about the global X-axis
            'MY' = Moment about the global Y-axis
            'MZ' = Moment about the global Z-axis
        Load : number
            The magnitude (value) of the load
        """
        
        # Find the node object associated with the user's given node name
        for node in self.Nodes:
            
            if node.Name == Node:
                self.NodeLoads.append([node, Direction, Load])
        
        # Flag the model as unsolved due to the new nodal load
        self.Solved = False

#%%      
    # Adds a member point load to the model
    def AddMemberPtLoad(self, Member, Direction, P, x):
        """
        Adds a member point load to the model
        
        Parameters
        ----------
        Member : string
            The name of the member the load is being applied to
        Direction : string
            Note that beam sign convention is used:
            'Fx' = Axial force in the member's local x-direction (positive is +x)
            'Fy' = Transverse force in the member's local y-direction (positive is -y)
            'Fz' = Transverse force in the member's local z-direction (positive is -z)
            'Mx' = Moment about the member's local x-axis
            'My' = Moment about the member's local y-axis
            'Mz' = Moment about the member's local z-axis    
        P : number
            The load magnitude (value)
        x : number
            The load's location along the member's local x-axis
        """
        
        # Find the member
        for member in self.Members:  
            if member.Name == Member:
                # Add the distributed load to the member
                member.PtLoads.append((Direction, P, x))

#%%
    # Adds a member distributed load to the model
    def AddMemberDistLoad(self, Member, Direction, w1, w2, x1, x2):
        """
        Adds a member distributed load to the model
        
        Parameters
        ----------
        Member : string
            The name of the member the load is being appied to
        Direction : string
            Note that beam sign convention is used:
            'Fx' = Axial force in the member's local x-direction (positive is +x)
            'Fy' = Transverse force in the member's local y-direction (positive is -y)
            'Fz' = Transverse force in the member's local z-direction (positive is -z)
        w1 : number
            The starting magnitude (value) of the load
        w2 : number
            The ending magnitude (value) of the load
        x1 : number
            The load's start location along the member's local x-axis
        x2 : number
            The load's end location along the member's local x-axis
        """
        
        # Find the member
        for member in self.Members:  
            if member.Name == Member:
                # Add the distributed load to the member
                member.DistLoads.append((Direction, w1, w2, x1, x2))

#%%    
    # Returns the global stiffness matrix for the structure
    def K(self):
        """
        Assembles and returns the global stiffness matrix
        """
        
        # Initialize a zero matrix to house all the stiffness terms
        stiffMatrix = zeros((len(self.Nodes)*6, len(self.Nodes)*6))
        
        # Add stiffness terms for each member in the model
        for member in self.Members:
            
            # Step through each term in the member's stiffness matrix
            # 'a' & 'b' below are row/column indices in the member's stiffness matrix
            # 'm' & 'n' are corresponding row/column indices in the global stiffness matrix
            
            for a in range(12):
                                    
                if a < 6:
                    m = member.iNode.ID*6 + a
                else:
                    m = member.jNode.ID*6 + (a - 6)
                    
                for b in range(12):
                    
                    if b < 6:
                        n = member.iNode.ID*6 + b
                    else:
                        n = member.jNode.ID*6 + (b - 6)

                    stiffMatrix.itemset((m,n), stiffMatrix.item((m,n)) + member.K().item((a,b)))
        
        # Return the global stiffness matrix
        return stiffMatrix
    
#%%    
    # Returns the global fixed end reaction vecotr
    def FER(self):
        """
        Assembles and returns the global fixed end reaction vector
        """
        
        # Initialize the fixed end reaction vector
        ferVector = zeros((len(self.Nodes)*6, 1))
        
        # Step through each member in the model
        for member in self.Members:
            
            # Step through each term in the member's fixed end reaction vector
            for a in range(12):
                
                if a < 6:
                    m = member.iNode.ID*6 + a
                else:
                    m = member.jNode.ID*6 + (a - 6)

                ferVector.itemset((m,0), ferVector.item((m,0)) + member.fer().item((a,0)))
        
        # Return the fixed end reaction vector
        return ferVector
    
#%%
    # Returns the global nodal force vector
    def P(self):
        """
        Assembles and returns the global nodal force vector for the structure
        """
        # Initialize the global nodal force vector
        pVector = zeros((len(self.Nodes)*6, 1))
        
        # Populate the vector
        # Step through each nodal load in 'NodeLoads'
        for nodeLoad in self.NodeLoads:
            
            # Find the node ID corresponding to this node load
            for node in self.Nodes:
                if node == nodeLoad[0]:
                    ID = node.ID
                    break
            
            # Add the load to the global nodal load vector
            if nodeLoad[1] == 'FX':
                pVector.itemset((ID*6+0, 0), pVector.item((ID*6+0, 0))+nodeLoad[2])
            elif nodeLoad[1] == 'FY':
                pVector.itemset((ID*6+1, 0), pVector.item((ID*6+1, 0))+nodeLoad[2])
            elif nodeLoad[1] == 'FZ':
                pVector.itemset((ID*6+2, 0), pVector.item((ID*6+2, 0))+nodeLoad[2])
            elif nodeLoad[1] == 'MX':
                pVector.itemset((ID*6+3, 0), pVector.item((ID*6+3, 0))+nodeLoad[2])
            elif nodeLoad[1] == 'MY':
                pVector.itemset((ID*6+4, 0), pVector.item((ID*6+4, 0))+nodeLoad[2])
            elif nodeLoad[1] == 'MZ':
                pVector.itemset((ID*6+5, 0), pVector.item((ID*6+5, 0))+nodeLoad[2])
        
        # Return the global nodal force vector
        return pVector
#%%
    def D(self):
        
        if self.Solved == True:
            return self.D
        else:
            self.Analyze()
            return self.D
#%%  
    # Analyzes the model
    def Analyze(self):
        """
        Analyzes the model
        """
        
        # Get the global stiffness matrix
        K = self.K()
        
        # Get the global fixed end reaction vector
        FER = self.FER()
        
        # Get the global nodal force vector
        P = self.P()
        
        # Eliminate supported degrees of freedom from each of the matrices/vectors
        # Work backwards through the node list so that the relationship between the DOF's and node ID's is unnafected by the matrices/vectors shrinking
        for node in reversed(self.Nodes):
            
            if node.SupportRZ == True:
                K = delete(K, node.ID*6+5, axis = 0)
                K = delete(K, node.ID*6+5, axis = 1)
                FER = delete(FER, node.ID*6+5, axis = 0)
                P = delete(P, node.ID*6+5, axis = 0)
                
            if node.SupportRY == True:
                K = delete(K, node.ID*6+4, axis = 0)
                K = delete(K, node.ID*6+4, axis = 1)
                FER = delete(FER, node.ID*6+4, axis = 0)
                P = delete(P, node.ID*6+4, axis = 0)
                
            if node.SupportRX == True:
                K = delete(K, node.ID*6+3, axis = 0)
                K = delete(K, node.ID*6+3, axis = 1)
                FER = delete(FER, node.ID*6+3, axis = 0)
                P = delete(P, node.ID*6+3, axis = 0)
                
            if node.SupportDZ == True:
                K = delete(K, node.ID*6+2, axis = 0)
                K = delete(K, node.ID*6+2, axis = 1)
                FER = delete(FER, node.ID*6+2, axis = 0)
                P = delete(P, node.ID*6+2, axis = 0)
                
            if node.SupportDY == True:
                K = delete(K, node.ID*6+1, axis = 0)
                K = delete(K, node.ID*6+1, axis = 1)
                FER = delete(FER, node.ID*6+1, axis = 0)
                P = delete(P, node.ID*6+1, axis = 0)
                
            if node.SupportDX == True:
                K = delete(K, node.ID*6+0, axis = 0)
                K = delete(K, node.ID*6+0, axis = 1)
                FER = delete(FER, node.ID*6+0, axis = 0)
                P = delete(P, node.ID*6+0, axis = 0)
                        
        # Calculate the global displacement vector
        self.D = matmul(inv(K), subtract(P, FER))
        
        # Expand the global displacement vector to include supported degrees of freedom
        # Work forwards through the node list so that the relationship between the DOF's and node ID's is unnafected by the vector expanding
        for node in self.Nodes:
            if node.SupportDX == True:
                self.D = insert(self.D, node.ID*6+0, 0, axis = 0)
            if node.SupportDY == True:
                self.D = insert(self.D, node.ID*6+1, 0, axis = 0)
            if node.SupportDZ == True:
                self.D = insert(self.D, node.ID*6+2, 0, axis = 0)
            if node.SupportRX == True:
                self.D = insert(self.D, node.ID*6+3, 0, axis = 0)
            if node.SupportRY == True:
                self.D = insert(self.D, node.ID*6+4, 0, axis = 0)
            if node.SupportRZ == True:
                self.D = insert(self.D, node.ID*6+5, 0, axis = 0)

        # Store the calculated global nodal displacements into each node
        for node in self.Nodes:
            node.DX = self.D.item((node.ID*6+0, 0))
            node.DY = self.D.item((node.ID*6+1, 0))
            node.DZ = self.D.item((node.ID*6+2, 0))
            node.RX = self.D.item((node.ID*6+3, 0))
            node.RY = self.D.item((node.ID*6+4, 0))
            node.RZ = self.D.item((node.ID*6+5, 0))

        # Flag the model as solved
        self.Solved = True
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        