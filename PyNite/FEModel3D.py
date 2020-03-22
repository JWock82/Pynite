# %%
from numpy import matrix, zeros, delete, insert, matmul, divide, add, subtract, nanmax, seterr, shape
from numpy.linalg import solve, matrix_rank
from math import isclose
from PyNite.Node3D import Node3D
from PyNite.Member3D import Member3D
from PyNite.Plate3D import Plate3D
from PyNite.LoadCombo import LoadCombo

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
        
        self.Nodes = []      # A list of the structure's nodes
        self.auxNodes = []   # A list of the structure's auxiliary nodes
        self.Members = []    # A list of the structure's members
        self.Plates = []     # A list of the structure's plates
        self.__D = {}        # A dictionary of the structure's nodal displacements by load combination
        self.LoadCombos = {} # A dictionary of the structure's load combinations

#%%
    def AddNode(self, Name, X, Y, Z):
        '''
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
        '''
        
        # Create a new node
        newNode = Node3D(Name, X, Y, Z)
        
        # Add the new node to the list
        self.Nodes.append(newNode)

#%%
    def AddAuxNode(self, Name, X, Y, Z):
        '''
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
        '''
        
        # Create a new node
        newNode = Node3D(Name, X, Y, Z)
        
        # Add the new node to the list
        self.auxNodes.append(newNode)
  
#%%
    def AddMember(self, Name, iNode, jNode, E, G, Iy, Iz, J, A, auxNode=None):
        '''
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
        '''
        
        # Create a new member
        if auxNode == None:
            newMember = Member3D(Name, self.GetNode(iNode), self.GetNode(jNode), E, G, Iy, Iz, J, A, LoadCombos=self.LoadCombos)
        else:
            newMember = Member3D(Name, self.GetNode(iNode), self.GetNode(jNode), E, G, Iy, Iz, J, A, self.GetAuxNode(auxNode), self.LoadCombos)
        
        # Add the new member to the list
        self.Members.append(newMember)

#%%
    def AddPlate(self, Name, iNode, jNode, mNode, nNode, t, E, nu):
        '''
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
        '''
        
        # Create a new member
        newPlate = Plate3D(Name, self.GetNode(iNode), self.GetNode(jNode), self.GetNode(mNode), self.GetNode(nNode), t, E, nu)
        
        # Add the new member to the list
        self.Plates.append(newPlate)

#%%
    def RemoveNode(self, Node):
        '''
        Removes a node from the model. All nodal loads associated with the
        node and members attached to the node will also be removed.
        
        Parameters
        ----------
        Node : string
            The name of the node to be removed.
        '''
        
        # Remove the node. Nodal loads are stored within the node, so they
        # will be deleted automatically when the node is deleted.
        self.Nodes.remove(self.GetNode(Node))
        
        # Find any members attached to the node and remove them
        self.Members = [member for member in self.Members if member.iNode.Name != Node and member.jNode.Name != Node]
        
#%%
    def RemoveMember(self, Member):
        '''
        Removes a member from the model. All member loads associated with the
        member will also be removed.
        
        Parameters
        ----------
        Member : string
            The name of the member to be removed.
        '''
        
        # Remove the member. Member loads are stored within the member, so they
        # will be deleted automatically when the member is deleted.
        self.Members.remove(self.GetMember(Member))
        
#%%
    def DefineSupport(self, Node, SupportDX=False, SupportDY=False, SupportDZ=False, SupportRX=False, SupportRY=False, SupportRZ=False):
        '''
        Defines the support conditions at a node.
        
        Nodes will default to fully unsupported unless specified otherwise.
        
        Parameters
        ----------
        Node : string
            The name of the node where the support is being defined
        SupportDX : number
            Indicates whether the node is supported against translation in the global X-direction.
        SupportDY : number
            Indicates whether the node is supported against translation in the global Y-direction.
        SupportDZ : number
            Indicates whether the node is supported against translation in the global Z-direction.
        SupportRX : number
            Indicates whether the node is supported against rotation about the global X-axis.
        SupportRY : number
            Indicates whether the node is supported against rotation about the global Y-axis.
        SupportRZ : number
            Indicates whether the node is supported against rotation about the global Z-axis.
        '''
        
        # Get the node to be supported
        node = self.GetNode(Node)
                
        # Set the node's support conditions
        node.SupportDX = SupportDX
        node.SupportDY = SupportDY
        node.SupportDZ = SupportDZ
        node.SupportRX = SupportRX
        node.SupportRY = SupportRY
        node.SupportRZ = SupportRZ

#%%            
    def AddNodeDisplacement (self, Node, Direction, Magnitude): 
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
            node.EnforceDX = Magnitude
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
        self.GetMember(Member).Releases = [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]     

#%%
    def AddLoadCombo(self, name, factors, combo_type='strength'):
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
    def AddNodeLoad(self, Node, Direction, P, case='Case 1'):
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
        
        # Add the node load to the model
        self.GetNode(Node).NodeLoads.append((Direction, P, case))

#%%      
    def AddMemberPtLoad(self, Member, Direction, P, x, case='Case 1'):
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
        
        # Add the point load to the member
        self.GetMember(Member).PtLoads.append((Direction, P, x, case))

#%%
    def AddMemberDistLoad(self, Member, Direction, w1, w2, x1=None, x2=None, case='Case 1'):
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
        self.GetMember(Member).DistLoads.append((Direction, w1, w2, start, end, case))

#%%
    def ClearLoads(self):
        '''
        Clears all loads from the model along with any results based on the loads.
        '''

        # Clear out the member loads and the calculated internal forces
        for member in self.Members:
            member.DistLoads = []
            member.PtLoads = []
            member.SegmentsZ = []
            member.SegmentsY = []
            member.SegmentsX = []
        
        # Clear out the nodal loads, calculated displacements, and calculated reactions
        for node in self.Nodes:

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
    def GetNode(self, Name):
        '''
        Returns the node with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the node to be returned.
        '''
        
        # Step through each node in the 'Nodes' list
        for node in self.Nodes:
            
            # Check the name of the node
            if node.Name == Name:
                
                # Return the node of interest
                return node

            
    def GetAuxNode(self, Name):
        '''
        Returns the auxiliary node with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the auxiliary node to be returned.
        '''
        
        # Step through each node in the 'Nodes' list
        for node in self.auxNodes:
            
            # Check the name of the node
            if node.Name == Name:
                
                # Return the node of interest
                return node            
            
#%%
    def GetMember(self, Name):
        '''
        Returns the member with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the member to be returned.
        '''
        
        # Step through each member in the 'Members' list
        for member in self.Members:
            
            # Check the name of the member
            if member.Name == Name:
                
                # Return the member of interest
                return member

#%%
    def GetPlate(self, Name):
        '''
        Returns the plate with the given name.
        
        Parameters
        ----------
        Name : string
            The name of the plate to be returned.
        '''
        
        # Step through each plate in the 'Plates' list
        for plate in self.Plates:
            
            # Check the name of the plate
            if plate.Name == Name:
                
                # Return the plate of interest
                return plate

#%%
    def __Renumber(self):
        '''
        Assigns node, plate, and member ID numbers to be used internally by the
        program. Numbers are assigned according to the order nodes, members, and plates
        were added to the model.
        
        '''
        
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
    def __AuxList(self):
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
        for node in self.Nodes:
            
            # Unknown displacement DX
            if node.SupportDX == False and node.EnforcedDX == None:
                D1_indices.append((node.ID*6) + 0)
            # Known displacement DX
            elif node.EnforcedDX != None:
                D2_indices.append((node.ID*6) + 0)
                D2.append(node.EnforcedDX)
            # Support at DX
            else:
                D2_indices.append((node.ID*6) + 0)
                D2.append(0.0)

            # Unknown displacement DY
            if node.SupportDY == False and node.EnforcedDY == None:
                D1_indices.append((node.ID*6) + 1)
            # Known displacement DY
            elif node.EnforcedDY != None:
                D2_indices.append((node.ID*6) + 1)
                D2.append(node.EnforcedDY)
            # Support at DY
            else:
                D2_indices.append((node.ID*6) + 1)
                D2.append(0.0)

            # Unknown displacement DZ
            if node.SupportDZ == False and node.EnforcedDZ == None:
                D1_indices.append((node.ID*6) + 2)
            # Known displacement DZ
            elif node.EnforcedDZ != None:
                D2_indices.append((node.ID*6) + 2)
                D2.append(node.EnforcedDZ)
            # Support at DZ
            else:
                D2_indices.append((node.ID*6) + 2)
                D2.append(0.0)

            # Unknown displacement RX
            if node.SupportRX == False and node.EnforcedRX == None:
                D1_indices.append((node.ID*6) + 3)
            # Known displacement RX
            elif node.EnforcedRX != None:
                D2_indices.append((node.ID*6) + 3)
                D2.append(node.EnforcedRX)
            # Support at RX
            else:
                D2_indices.append((node.ID*6) + 3)
                D2.append(0.0)

            # Unknown displacement RY
            if node.SupportRY == False and node.EnforcedRY == None:
                D1_indices.append((node.ID*6) + 4)
            # Known displacement RY
            elif node.EnforcedRY != None:
                D2_indices.append((node.ID*6) + 4)
                D2.append(node.EnforcedRY)
            # Support at RY
            else:
                D2_indices.append((node.ID*6) + 4)
                D2.append(0.0)

            # Unknown displacement RZ
            if node.SupportRZ == False and node.EnforcedRZ == None:
                D1_indices.append((node.ID*6) + 5)
            # Known displacement RZ
            elif node.EnforcedRZ != None:
                D2_indices.append((node.ID*6) + 5)
                D2.append(node.EnforcedRZ)
            # Support at RZ
            else:
                D2_indices.append((node.ID*6) + 5)
                D2.append(0.0)

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
    def Kg(self, combo_name='Combo 1'):
        '''
        Assembles and returns the global geometric stiffness matrix.

        The model must have a static solution prior to obtaining the geometric stiffness matrix.
        Geometric stiffness of plates is not included.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to derive the matrix for (not the load combination itself).
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
                    Kg.itemset((m, n), Kg.item((m, n)) + member_Kg.item((a, b)))

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
        for member in self.Members:
            
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
                FER.itemset((m, 0), FER[m, 0] + member_FER[a, 0])
        
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
        
        # Add terms for each node in the model
        for node in self.Nodes:
            
            # Get the node's ID
            ID = node.ID
            
            # Get the load combination for the given 'combo_name'
            combo = self.LoadCombos[combo_name]

            # Step through each load factor in the load combination
            for case, factor in combo.factors.items():

                # Add the node's loads to the global nodal load vector
                for load in node.NodeLoads:

                    if load[2] == case:

                        if load[0] == 'FX':
                            P.itemset((ID*6 + 0, 0), P[ID*6 + 0, 0] + factor*load[1])
                        elif load[0] == 'FY':
                            P.itemset((ID*6 + 1, 0), P[ID*6 + 1, 0] + factor*load[1])
                        elif load[0] == 'FZ':
                            P.itemset((ID*6 + 2, 0), P[ID*6 + 2, 0] + factor*load[1])
                        elif load[0] == 'MX':
                            P.itemset((ID*6 + 3, 0), P[ID*6 + 3, 0] + factor*load[1])
                        elif load[0] == 'MY':
                            P.itemset((ID*6 + 4, 0), P[ID*6 + 4, 0] + factor*load[1])
                        elif load[0] == 'MZ':
                            P.itemset((ID*6 + 5, 0), P[ID*6 + 5, 0] + factor*load[1])
        
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
    def __Partition(self, unp_matrix, D1_indices, D2_indices):
        '''
        Partitions a matrix into submatrices based on degree of freedom boundary conditions

        Parameters
        ----------
        unp_matrix : matrix
            The unpartitioned matrix to be partitioned.
        '''

        if unp_matrix.shape[1] == 1:
            m1 = unp_matrix[D1_indices, :]
            m2 = unp_matrix[D2_indices, :]
            return m1, m2
        else:
            m11 = unp_matrix[D1_indices, :][:, D1_indices]
            m12 = unp_matrix[D1_indices, :][:, D2_indices]
            m21 = unp_matrix[D2_indices, :][:, D1_indices]
            m22 = unp_matrix[D2_indices, :][:, D2_indices]
            return m11, m12, m21, m22

#%%  
    def Analyze(self, check_statics=True):
        '''
        Analyzes the model.
        '''
        
        print('**Analyzing**')

        # Assign an ID to all nodes and elements in the model
        self.__Renumber()

        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = self.__AuxList()

        # Convert D2 from a list to a matrix
        D2 = matrix(D2).T

        # Get the partitioned global stiffness matrix K11, K12, K21, K22
        K11, K12, K21, K22 = self.__Partition(self.K(), D1_indices, D2_indices)

        # Ensure there is at least 1 load combination to solve if the user didn't define any
        if self.LoadCombos == {}:
            # Create and add a default load combination to the dictionary of load combinations
            self.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})

        # Step through each load combination
        for combo in self.LoadCombos.values():

            # Get the partitioned global fixed end reaction vector
            FER1, FER2 = self.__Partition(self.FER(combo.name), D1_indices, D2_indices)

            # Get the partitioned global nodal force vector       
            P1, P2 = self.__Partition(self.P(combo.name), D1_indices, D2_indices)          

            # Check for global stability by determining if 'K11' is singular
            print('...Checking global stability')
            if K11.shape == (0, 0):
                # All displacements are known, so D1 is an empty vector
                D1 = []
            elif matrix_rank(K11) < min(K11.shape):
                # Return out of the method if 'K' is singular and provide an error message
                print('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')
                return
            else:
                # Calculate the unknown displacements D1
                print('...Calculating global displacement vector')
                D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))

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

            # Save the global displacement vector
            self.__D[combo.name] = D

            # Store the calculated global nodal displacements into each node
            for node in self.Nodes:
                node.DX[combo.name] = D[node.ID*6 + 0, 0]
                node.DY[combo.name] = D[node.ID*6 + 1, 0]
                node.DZ[combo.name] = D[node.ID*6 + 2, 0]
                node.RX[combo.name] = D[node.ID*6 + 3, 0]
                node.RY[combo.name] = D[node.ID*6 + 4, 0]
                node.RZ[combo.name] = D[node.ID*6 + 5, 0]
        
        # Calculate reactions
        self.__CalcReactions()
    
        # Check statics if requested
        if check_statics == True:
            self.__CheckStatics()

#%%
    def Analyze_PDelta(self, max_iter=30, tol=0.01):
        '''
        Runs a second order (P-Delta) analysis on the structure.

        Parameters
        ----------
        max_iter : number
            The maximum number of iterations permitted. If this value is exceeded the program will
            report divergence.
        tol : number
            The deflection tolerance (as a percentage) between iterations that will be used to define whether the model
            has converged (e.g. 0.01 = deflections must converge within 1% between iterations).
        '''
        
        print('**Running P-Delta analysis**')

        # Assign an ID to all nodes and elements in the model
        self.__Renumber()

        # Get the auxiliary list used to determine how the matrices will be partitioned
        D1_indices, D2_indices, D2 = self.__AuxList()

        # Convert D2 from a list to a matrix
        D2 = matrix(D2).T    

        # Ensure there is at least 1 load combination to solve if the user didn't define any
        if self.LoadCombos == {}:
            # Create and add a default load combination to the dictionary of load combinations
            self.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})

        # Step through each load combination
        for combo in self.LoadCombos.values():

            # Keep track of the number of iterations
            iter_count = 1
            convergence = False
            divergence = False

            # Iterate until convergence or divergence occurs
            while convergence == False and divergence == False:
            
                # Inform the user which iteration we're on
                print('...Beginning P-Delta iteration #' + str(iter_count))

                # Get the partitioned global matrices
                if iter_count == 1:
                    
                    K11, K12, K21, K22 = self.__Partition(self.K(), D1_indices, D2_indices) # Initial stiffness matrix
                    FER1, FER2 = self.__Partition(self.FER(combo.name), D1_indices, D2_indices)  # Fixed end reactions
                    P1, P2 = self.__Partition(self.P(combo.name), D1_indices, D2_indices)        # Nodal forces

                else:

                    # Calculate the global stiffness matrices (partitioned)
                    K11, K12, K21, K22 = self.__Partition(self.K(), D1_indices, D2_indices)           # Initial stiffness matrix
                    Kg11, Kg12, Kg21, Kg22 = self.__Partition(self.Kg(combo.name), D1_indices, D2_indices) # Geometric stiffness matrix

                    # Combine the stiffness matrices
                    K11 = add(K11, Kg11)
                    K12 = add(K12, Kg12)
                    K21 = add(K21, Kg21)
                    K22 = add(K22, Kg22)                     

                # Determine if 'K' is singular
                print('...Checking global stability')
                if K11.shape == (0, 0):
                    # All displacements are known, so D1 is an empty vector
                    D1 = []
                elif matrix_rank(K11) < min(K11.shape):
                    # Return out of the method if 'K' is singular and provide an error message
                    print('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')
                    return
                else:
                    # Calculate the global displacement vector
                    print('...Calculating global displacement vector')
                    D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))
            
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

                # Save the global displacement vector
                self.__D[combo.name] = D

                # Store the calculated global nodal displacements into each node
                for node in self.Nodes:

                    node.DX[combo.name] = D[node.ID*6 + 0, 0]
                    node.DY[combo.name] = D[node.ID*6 + 1, 0]
                    node.DZ[combo.name] = D[node.ID*6 + 2, 0]
                    node.RX[combo.name] = D[node.ID*6 + 3, 0]
                    node.RY[combo.name] = D[node.ID*6 + 4, 0]
                    node.RZ[combo.name] = D[node.ID*6 + 5, 0]

                if iter_count != 1:
                
                    # Print a status update for the user
                    print('...Checking for convergence.')

                    # Temporarily disable error messages for invalid values.
                    # We'll be dealing with some 'nan' values due to division by zero at supports with zero deflection.
                    seterr(invalid='ignore')

                    # Check for convergence
                    if abs(1 - nanmax(divide(prev_results, D1))) <= tol:
                        convergence = True
                        print('...P-Delta analysis converged after ' + str(iter_count) + ' iterations.')
                    # Check for divergence
                    elif iter_count > max_iter:
                        divergence = True
                        print('...P-Delta analysis failed to converge after 30 iterations.')

                    # Turn invalid value warnings back on
                    seterr(invalid='warn') 

                # Save the results for the next iteration
                prev_results = D1

                # Increment the iteration count
                iter_count += 1
        
        # Calculate reactions
        self.__CalcReactions()

#%%
    def __CalcReactions(self):
        '''
        Calculates reactions once the model is solved.
        '''

        # Print a status update to the console
        print('...Calculating reactions')

        # Calculate the reactions, node by node
        for node in self.Nodes:
            
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
                if (node.SupportDX != None and isclose(node.SupportDX, 0.0)) \
                or (node.SupportDY != None and isclose(node.SupportDY, 0.0)) \
                or (node.SupportDZ != None and isclose(node.SupportDZ, 0.0)) \
                or (node.SupportRX != None and isclose(node.SupportRX, 0.0)) \
                or (node.SupportRY != None and isclose(node.SupportRY, 0.0)) \
                or (node.SupportRZ != None and isclose(node.SupportRZ, 0.0)):

                    # Sum the member end forces at the node
                    for member in self.Members:
                    
                        if member.iNode == node:
                        
                            # Get the member's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            member_F = member.F(combo.name)

                            node.RxnFX[combo.name] += member_F[0, 0]
                            node.RxnFY[combo.name] += member_F[1, 0]
                            node.RxnFZ[combo.name] += member_F[2, 0]
                            node.RxnMX[combo.name] += member_F[3, 0]
                            node.RxnMY[combo.name] += member_F[4, 0]
                            node.RxnMZ[combo.name] += member_F[5, 0]

                        elif member.jNode == node:
                        
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
                    for plate in self.Plates:

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

                        elif plate.jNode == node:

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

                        elif plate.nNode == node:

                            # Get the plate's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            plate_F = plate.F(combo.name)
                    
                            node.RxnFX[combo.name] += plate_F[18, 0]
                            node.RxnFY[combo.name] += plate_F[19, 0]
                            node.RxnFZ[combo.name] += plate_F[20, 0]
                            node.RxnMX[combo.name] += plate_F[21, 0]
                            node.RxnMY[combo.name] += plate_F[22, 0]
                            node.RxnMZ[combo.name] += plate_F[23, 0]

                    # Sum the joint forces at the node
                    for load in node.NodeLoads:
                    
                        if load[0] == 'FX':
                            node.RxnFX[combo.name] -= load[1]
                        elif load[0] == 'FY':
                            node.RxnFY[combo.name] -= load[1]
                        elif load[0] == 'FZ':
                            node.RxnFZ[combo.name] -= load[1]
                        elif load[0] == 'MX':
                            node.RxnMX[combo.name] -= load[1]
                        elif load[0] == 'MY':
                            node.RxnMY[combo.name] -= load[1]
                        elif load[0] == 'MZ':
                            node.RxnMZ[combo.name] -= load[1]

#%%
    def __CheckStatics(self):
        '''
        Sums global forces and reactions and prints them to the console.


        '''
        
        # Print a status update to the console
        print('...Checking statics')

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
        
            # Print the load summation
            print('**Load Combination:', combo.name + '**')
            print('**Applied Loads**')
            print('Sum Forces X:', SumFX, ', Sum Forces Y:', SumFY, ', Sum Forces Z:', SumFZ, ', Sum Moments MX:', SumMX, ', Sum Moments MY:', SumMY, ', Sum Moments MZ:', SumMZ)

            # Print the reaction summation
            print('**Reactions**')
            print('Sum Forces X:', SumRFX, ', Sum Forces Y:', SumRFY, ', Sum Forces Z:', SumRFZ, ', Sum Moments MX:', SumRMX, ', Sum Moments MY:', SumRMY, ', Sum Moments MZ:', SumRMZ)
