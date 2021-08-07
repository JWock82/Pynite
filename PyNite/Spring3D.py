# %%
from numpy import zeros, array, transpose, add, subtract, matmul, insert, cross, divide
from numpy.linalg import inv
from math import isclose
import PyNite.FixedEndReactions
from PyNite.LoadCombo import LoadCombo

# %%
class Spring3D():
    '''
    A class representing a 3D spring element in a finite element model.
    '''

    # '__plt' is used to store the 'pyplot' from matplotlib once it gets imported. Setting it to 'None' for now allows
    # us to defer importing it until it's actually needed.
    __plt = None

#%%
    def __init__(self, Name, iNode, jNode, ks, LoadCombos={'Combo 1':LoadCombo('Combo 1', factors={'Case 1':1.0})},
                 tension_only=False, comp_only=False):
        '''
        Initializes a new spring.
        '''
        self.Name = Name    # A unique name for the spring given by the user
        self.ID = None      # Unique index number for the spring assigned by the program
        self.iNode = iNode  # The spring's i-node
        self.jNode = jNode  # The spring's j-node
        self.ks = ks        # The spring constant (force/displacement)
        self.LoadCombos = LoadCombos # The dictionary of load combinations in the model this spring belongs to
        self.tension_only = tension_only # Indicates whether the spring is tension-only
        self.comp_only = comp_only # Indicates whether the spring is compression-only

        # Springs need to track whether they are active or not for any given load combination.
        # They may become inactive for a load combination during a tension/compression-only
        # analysis. This dictionary will be used when the model is solved.
        self.active = {} # Key = load combo name, Value = True or False

#%%
    def L(self):
        '''
        Returns the length of the spring.
        '''

        # Get the i-node and the j-node for the spring
        iNode = self.iNode
        jNode = self.jNode

        # Return the distance between the two nodes
        return ((jNode.X-iNode.X)**2+(jNode.Y-iNode.Y)**2+(jNode.Z-iNode.Z)**2)**0.5

#%%
    def k(self):
        '''
        Returns the local stiffness matrix for the spring.
        '''

        # Get the spring constant
        ks = self.ks
               
        # Calculate the local stiffness matrix
        k = array([
            [ks,  0, 0, 0, 0, 0, -ks, 0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [-ks, 0, 0, 0, 0, 0, ks,  0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0],
            [0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0]])

        # Return the local stiffness matrix
        return k

#%%   
    def f(self, combo_name='Combo 1'):
        '''
        Returns the spring's local end force vector for the given load combination.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the local end force vector for (not the load combination itself).
        '''
        
        # Calculate and return the spring's local end force vector
        return matmul(self.k(), self.d(combo_name))

#%%
    def d(self, combo_name='Combo 1'):
        '''
        Returns the spring's local displacement vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to construct the displacement vector for (not the load combination itself).
        '''
        
        # Calculate and return the local displacement vector
        return matmul(self.T(), self.D(combo_name))
        
#%%
    def T(self):
        '''
        Returns the transformation matrix for the spring.
        '''

        x1 = self.iNode.X
        x2 = self.jNode.X
        y1 = self.iNode.Y
        y2 = self.jNode.Y
        z1 = self.iNode.Z
        z2 = self.jNode.Z
        L = self.L()
        
        # Calculate the direction cosines for the local x-axis
        x = [(x2-x1)/L, (y2-y1)/L, (z2-z1)/L]
            
        # Calculate the remaining direction cosines. The local z-axis will be kept parallel to the global XZ plane in all cases
        # Vertical springs
        if isclose(x1, x2) and isclose(z1, z2):
                
            # For vertical springs, keep the local y-axis in the XY plane to make 2D problems easier to solve in the XY plane
            if y2 > y1:
                y = [-1, 0, 0]
                z = [0, 0, 1]
            else:
                y = [1, 0, 0]
                z = [0, 0, 1]

        # Horizontal springs
        elif isclose(y1, y2):
            
            # Find a vector in the direction of the local z-axis by taking the cross-product
            # of the local x-axis and the local y-axis. This vector will be perpendicular to
            # both the local x-axis and the local y-axis.
            y = [0, 1, 0]
            z = cross(x, y)

            # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
            z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)

        # Members neither vertical or horizontal
        else:

            # Find the projection of x on the global XZ plane
            proj = [x2-x1, 0, z2-z1]

            # Find a vector in the direction of the local z-axis by taking the cross-product
            # of the local x-axis and its projection on a plane parallel to the XZ plane. This
            # produces a vector perpendicular to both the local x-axis and its projection. This
            # vector will always be horizontal since it's parallel to the XZ plane. The order
            # in which the vectors are 'crossed' has been selected to ensure the y-axis always
            # has an upward component (i.e. the top of the beam is always on top).
            if y2 > y1:
                z = cross(proj, x)
            else:
                z = cross(x, proj)

            # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
            z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)
            
            # Find the direction cosines for the local y-axis
            y = cross(z, x)
            y = divide(y, (y[0]**2 + y[1]**2 + y[2]**2)**0.5)

        # Create the direction cosines matrix
        dirCos = array([x, y, z])
      
        # Build the transformation matrix
        transMatrix = zeros((12, 12))
        transMatrix[0:3, 0:3] = dirCos
        transMatrix[3:6, 3:6] = dirCos
        transMatrix[6:9, 6:9] = dirCos
        transMatrix[9:12, 9:12] = dirCos
        
        return transMatrix

#%%
    def K(self):
        '''
        Spring global stiffness matrix
        '''
        
        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(transpose(self.T()), self.k()), self.T())

#%%
    def F(self, combo_name='Combo 1'):
        '''
        Returns the spring's global end force vector for the given load combination.
        '''
        
        # Calculate and return the global force vector
        return matmul(inv(self.T()), self.f(combo_name))

#%%
    def D(self, combo_name='Combo 1'):
        '''
        Returns the spring's global displacement vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to construct the global
            displacement vector for (not the load combination itelf).
        '''
        
        # Initialize the displacement vector
        D = zeros((12, 1))
        
        # Read in the global displacements from the nodes
        # Apply axial displacements only if the spring is active
        if self.active[combo_name] == True:
            D[0, 0] = self.iNode.DX[combo_name]
            D[6, 0] = self.jNode.DX[combo_name]

        # Apply the remaining displacements
        D[1, 0] = self.iNode.DY[combo_name]
        D[2, 0] = self.iNode.DZ[combo_name]
        D[3, 0] = self.iNode.RX[combo_name]
        D[4, 0] = self.iNode.RY[combo_name]
        D[5, 0] = self.iNode.RZ[combo_name]
        D[7, 0] = self.jNode.DY[combo_name]
        D[8, 0] = self.jNode.DZ[combo_name]
        D[9, 0] = self.jNode.RX[combo_name]
        D[10, 0] = self.jNode.RY[combo_name]
        D[11, 0] = self.jNode.RZ[combo_name]      

        # Return the global displacement vector
        return D

#%%
    def axial(self, combo_name='Combo 1'):
        '''
        Returns the axial force in the spring.
        
        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
            
        # Calculate the axial force
        return self.f(combo_name)[0, 0]