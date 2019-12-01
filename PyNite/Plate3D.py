from numpy import zeros, matrix, matmul, transpose, insert, cross, divide, add
from numpy.linalg import inv

# A 3D isotropic plate bending element
# Membrane stresses are not yet supported
class Plate3D():

    def __init__(self, Name, iNode, jNode, mNode, nNode, t, E, mew):

        self.Name = Name
        self.ID = None

        self.iNode = iNode
        self.jNode = jNode
        self.mNode = mNode
        self.nNode = nNode

        self.t = t
        self.E = E
        self.mew = mew
    
    # Creates the local stiffness matrix
    def k(self):

        x2 = ((self.jNode.X-self.iNode.X)**2 + (self.jNode.Y-self.iNode.Y)**2 + (self.jNode.Z-self.iNode.Z)**2)**0.5
        y3 = ((self.mNode.X-self.jNode.X)**2 + (self.mNode.Y-self.jNode.Y)**2 + (self.mNode.Z-self.jNode.Z)**2)**0.5

        E = self.E
        t = self.t
        m = self.mew

        # Stiffness matrix for plate bending and out-of-plane displacement based on a 12-term polynomial
        b = x2/y3
        g = y3/x2

        k = matrix([[120*(b**2+g**2)-24*m+84,   0,                       0,                      0,                         0,                      0,                      0,                        0,                      0,                      0,                       0,                      0],
                    [(10*b**2+(1+4*m))*6*y3,    40*x2**2+8*(1-m)*y3**2,  0,                      0,                         0,                      0,                      0,                        0,                      0,                      0,                       0,                      0],
                    [-(10*g**2+(1+4*m))*6*x2,   -30*m*x2*y3,             40*y3**2+8*(1-m)*x2**2, 0,                         0,                      0,                      0,                        0,                      0,                      0,                       0,                      0],
                    [60*(g**2-2*b**2)+24*m-84,  -(10*b**2+(1-m))*6*y3,   (-5*g**2+(1+4*m))*6*x2, 120*(b**2+g**2)-24*m+84,   0,                      0,                      0,                        0,                      0,                      0,                       0,                      0],
                    [(10*b**2+(1-m))*6*y3,      20*x2**2-2*(1-m)*y3**2,  0,                      -(10*b**2+(1+4*m))*6*y3,   40*x2**2+8*(1-m)*y3**2, 0,                      0,                        0,                      0,                      0,                       0,                      0],
                    [(-5*g**2+(1+4*m))*6*x2,    0,                       20*y3**2-8*(1-m)*x2**2, -(10*g**2+(1+4*m))*6*x2,   30*m*x2*y3,             40*y3**2+8*(1-m)*x2**2, 0,                        0,                      0,                      0,                       0,                      0],
                    [-60*(g**2+b**2)-24*m+84,   (-5*b**2+(1-m))*6*y3,    (5*g**2-(1-m))*6*x2,    -60*(2*g**2-b**2)+24*m-84, (-5*b**2+(1+4*m))*6*y3, (10*g**2+(1-m))*6*x2,   120*(b**2+g**2)-24*m+84,  0,                      0,                      0,                       0,                      0],
                    [(5*b**2-(1-m))*6*y3,       10*x2**2+2*(1-m)*y3**2,  0,                      (-5*b**2+(1+4*m))*6*y3,    20*x2**2-8*(1-m)*y3**2, 0,                      -(10*b**2+(1+4*m))*6*y3,  40*x2**2+8*(1-m)*y3**2, 0,                      0,                       0,                      0],
                    [(-5*g**2+(1-m))*6*x2,      0,                       10*y3**2+2*(1-m)*x2**2, -(10*g**2+(1-m))*6*x2,     0,                      20*y3**2-2*(1-m)*x2**2, (10*g**2+(1+4*m))*6*x2,   -30*m*x2*y3,            40*y3**2+8*(1-m)*x2**2, 0,                       0,                      0],
                    [-60*(2*g**2-b**2)+24*m-84, (-5*b**2+(1+4*m))*6*y3,  (10*g**2+(1-m))*6*x2,   -60*(b**2+g**2)-24*m+84,   (5*b**2-(1-m))*6*y3,    (5*g**2-(1-m))*6*x2,    60*(g**2-2*b**2)+24*m-84, (10*b**2+(1-m))*6*y3,   (5*g**2-(1+4*m))*6*x2,  120*(b**2+g**2)-24*m+84, 0,                      0],
                    [(5*b**2-(1+4*m))*6*y3,     20*x2**2-8*(1-m)*y3**2,  0,                      (-5*b**2+(1-m))*6*y3,      10*x2**2+2*(1-m)*y3**2, 0,                      -(10*b**2+(1-m))*6*y3,    20*x2**2-2*(1-m)*y3**2, 0,                      (10*b**2+(1+4*m))*6*y3,  40*x2**2+8*(1-m)*y3**2, 0],
                    [-(10*g**2+(1-m))*6*x2,     0,                       20*y3**2-2*(1-m)*x2**2, (-5*g**2+(1-m))*6*x2,      0,                      10*y3**2+2*(1-m)*x2**2, (5*g**2-(1+4*m))*6*x2,    0,                      20*y3**2-8*(1-m)*x2**2, (10*g**2+(1+4*m))*6*x2,  30*m*x2*y3,             40*y3**2+8*(1-m)*x2**2]])

        k = k*(E*t**3/(360*(1-m**2)*x2*y3)) 

        # Apply symmetry to the matrix
        for i in range(12):
            for j in range(i, 12):
                k[i, j] = k[j, i]
        
        # Insert rows of zeros for degrees of freedom not included in the matrix above
        k = insert(k, 12, 0, axis=0)
        k = insert(k, 12, 0, axis=1)

        k = insert(k, 9, 0, axis=0)
        k = insert(k, 9, 0, axis=1)
        k = insert(k, 9, 0, axis=0)
        k = insert(k, 9, 0, axis=1)

        k = insert(k, 9, 0, axis=0)
        k = insert(k, 9, 0, axis=1)

        k = insert(k, 6, 0, axis=0)
        k = insert(k, 6, 0, axis=1)
        k = insert(k, 6, 0, axis=0)
        k = insert(k, 6, 0, axis=1)

        k = insert(k, 6, 0, axis=0)
        k = insert(k, 6, 0, axis=1)

        k = insert(k, 3, 0, axis=0)
        k = insert(k, 3, 0, axis=1)
        k = insert(k, 3, 0, axis=0)
        k = insert(k, 3, 0, axis=1)

        k = insert(k, 3, 0, axis=0)
        k = insert(k, 3, 0, axis=1)

        k = insert(k, 0, 0, axis=0)
        k = insert(k, 0, 0, axis=1)
        k = insert(k, 0, 0, axis=0)
        k = insert(k, 0, 0, axis=1)

        # Return the local stiffness matrix
        return k

#%%   
    def f(self):
        """
        Returns the plate's local end force vector
        """
        
        # Calculate and return the plate's local end force vector
        return add(matmul(self.k(), self.d()), self.fer())

#%%
    def fer(self):
        """
        Returns the plate's local fixed end reaction vector (zero's for now until surface loads get added)
        """

        return zeros((24, 1))

#%%
    def d(self):
       """
       Returns the plate's local displacement vector
       """

       # Calculate and return the local displacement vector
       return matmul(self.T(), self.D())

#%%
    def F(self):
        
        # Calculate and return the global force vector
        return matmul(inv(self.T()), self.f())

#%%
    def D(self):
        """
        Returns the plate's global displacement vector
        """
        
        # Initialize the displacement vector
        Dvector = zeros((24, 1))
        
        # Read in the global displacements from the nodes
        Dvector.itemset((0, 0), self.iNode.DX)
        Dvector.itemset((1, 0), self.iNode.DY)
        Dvector.itemset((2, 0), self.iNode.DZ)
        Dvector.itemset((3, 0), self.iNode.RX)
        Dvector.itemset((4, 0), self.iNode.RY)
        Dvector.itemset((5, 0), self.iNode.RZ)

        Dvector.itemset((6, 0), self.jNode.DX)
        Dvector.itemset((7, 0), self.jNode.DY)
        Dvector.itemset((8, 0), self.jNode.DZ)
        Dvector.itemset((9, 0), self.jNode.RX)
        Dvector.itemset((10, 0), self.jNode.RY)
        Dvector.itemset((11, 0), self.jNode.RZ)

        Dvector.itemset((12, 0), self.mNode.DX)
        Dvector.itemset((13, 0), self.mNode.DY)
        Dvector.itemset((14, 0), self.mNode.DZ)
        Dvector.itemset((15, 0), self.mNode.RX)
        Dvector.itemset((16, 0), self.mNode.RY)
        Dvector.itemset((17, 0), self.mNode.RZ)

        Dvector.itemset((18, 0), self.nNode.DX)
        Dvector.itemset((19, 0), self.nNode.DY)
        Dvector.itemset((20, 0), self.nNode.DZ)
        Dvector.itemset((21, 0), self.nNode.RX)
        Dvector.itemset((22, 0), self.nNode.RY)
        Dvector.itemset((23, 0), self.nNode.RZ)
        
        # Return the global displacement vector
        return Dvector

#%%  
    # Transformation matrix
    def T(self):

        # Calculate the direction cosines for the local x-axis
        # The local x-axis will run from the i-node to the j-node
        xi = self.iNode.X
        xj = self.jNode.X
        yi = self.iNode.Y
        yj = self.jNode.Y
        zi = self.iNode.Z
        zj = self.jNode.Z
        Lx = ((xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)**0.5
        x = [(xj-xi)/Lx, (yj-yi)/Lx, (zj-zi)/Lx]
        
        # The local y-axis will be in the plane of the plate
        # Find a vector in the plate's local xy plane
        xm = self.mNode.X
        ym = self.mNode.Y
        zm = self.mNode.Z
        xy = [xm-xj, ym-yj, zm-zj]

        # Find a vector perpendicular to the plate surface to get the orientation of the local z-axis
        z = cross(x, xy)
        
        # Divide the vector by its magnitude to produce a unit z-vector of direction cosines
        z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)

        # Calculate the local y-axis as a vector perpendicular to the local z and x-axes
        y = cross(z, x)
        
        # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
        y = divide(y, (y[0]**2 + y[1]**2 + y[2]**2)**0.5)

        # Create the direction cosines matrix
        dirCos = matrix([x, y, z])
        
        # Build the transformation matrix
        transMatrix = zeros((24, 24))
        transMatrix[0:3, 0:3] = dirCos
        transMatrix[3:6, 3:6] = dirCos
        transMatrix[6:9, 6:9] = dirCos
        transMatrix[9:12, 9:12] = dirCos
        transMatrix[12:15, 12:15] = dirCos
        transMatrix[15:18, 15:18] = dirCos
        transMatrix[18:21, 18:21] = dirCos
        transMatrix[21:24, 21:24] = dirCos
        
        return transMatrix

#%%
    # Plate global stiffness matrix
    def K(self):
        
        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(transpose(self.T()), self.k()), self.T())

#%%
    # CODE BELOW NOT USED - INCOMPLETE ALTERNATE STIFFNESS MATRIX DERIVATION
    # **********************************************************************

    # # Calculates and returns the [C] coefficient matrix
    # def __C(self):

    #     # Find the x and y coordinates at each node
    #     xi = iNode.X
    #     yi = iNode.Y
    #     xj = jNode.X
    #     yj = jNode.Y
    #     xm = mNode.X
    #     ym = mNode.Y
    #     xn = nNode.X
    #     yn = nNode.Y

    #     # Calculate the [C] coefficient matrix
    #     C = matrix([[1, xi, yi, xi**2, xi*yi, yi**2, xi**3, xi**2*yi, xi*yi**2, yi**3, xi**3*yi, xi*yi**3],
    #                 [0, 0, 1, 0, xi, 2*yi, 0, xi**2, 2*xi*yi, 3*yi**2, xi**3, 3*xi*yi**2],
    #                 [0, -1, 0, -2*xi, -yi, 0, -3*xi**2, -2*xi*yi, -yi**2, 0, -3*xi**2*yi, -yi**3],

    #                 [1, xj, yj, xj**2, xj*yj, yj**2, xj**3, xj**2*yj, xj*yj**2, yj**3, xj**3*yj, xj*yj**3],
    #                 [0, 0, 1, 0, xj, 2*yj, 0, xj**2, 2*xj*yj, 3*yj**2, xj**3, 3*xj*yj**2],
    #                 [0, -1, 0, -2*xj, -yj, 0, -3*xj**2, -2*xj*yj, -yj**2, 0, -3*xj**2*yj, -yj**3],

    #                 [1, xm, ym, xm**2, xm*ym, ym**2, xm**3, xm**2*ym, xm*ym**2, ym**3, xm**3*ym, xm*ym**3],
    #                 [0, 0, 1, 0, xm, 2*ym, 0, xm**2, 2*xm*ym, 3*ym**2, xm**3, 3*xm*ym**2],
    #                 [0, -1, 0, -2*xm, -ym, 0, -3*xm**2, -2*xm*ym, -ym**2, 0, -3*xm**2*ym, -ym**3],
                    
    #                 [1, xn, yn, xn**2, xn*yn, yn**2, xn**3, xn**2*yn, xn*yn**2, yn**3, xn**3*yn, xn*yn**3],
    #                 [0, 0, 1, 0, xn, 2*yn, 0, xn**2, 2*xn*yn, 3*yn**2, xn**3, 3*xn*yn**2],
    #                 [0, -1, 0, -2*xn, -yn, 0, -3*xn**2, -2*xn*yn, -yn**2, 0, -3*xn**2*yn, -yn**3]])

    #     # Return the coefficient matrix
    #     return C

    # # Calculates and returns the [Q] coefficient matrix
    # def __Q(self):

    #     # Find the x and y coordinates at each node
    #     xi = iNode.X
    #     yi = iNode.Y
    #     xj = jNode.X
    #     yj = jNode.Y
    #     xm = mNode.X
    #     ym = mNode.Y
    #     xn = nNode.X
    #     yn = nNode.Y

    #     # Calculate the [Q] coefficient matrix
    #     Q = matrix([[0, 0, 0, -2, 0, 0, -6*xi, -2*yi, 0, 0, -6*xi*yi, 0],
    #                 [0, 0, 0, 0, 0, -2, 0, 0, -2*xi, -6*yi, 0, -6*xi*yi],
    #                 [0, 0, 0, 0, -2, 0, 0, -4*xi, -4*yi, 0, -6*xi**2, -6*yi**2],

    #                 [0, 0, 0, -2, 0, 0, -6*xj, -2*yj, 0, 0, -6*xj*yj, 0],
    #                 [0, 0, 0, 0, 0, -2, 0, 0, -2*xj, -6*yj, 0, -6*xj*yj],
    #                 [0, 0, 0, 0, -2, 0, 0, -4*xj, -4*yj, 0, -6*xj**2, -6*yj**2],

    #                 [0, 0, 0, -2, 0, 0, -6*xm, -2*ym, 0, 0, -6*xm*ym, 0],
    #                 [0, 0, 0, 0, 0, -2, 0, 0, -2*xm, -6*ym, 0, -6*xm*ym],
    #                 [0, 0, 0, 0, -2, 0, 0, -4*xm, -4*ym, 0, -6*xm**2, -6*ym**2],
                    
    #                 [0, 0, 0, -2, 0, 0, -6*xn, -2*yn, 0, 0, -6*xn*yn, 0],
    #                 [0, 0, 0, 0, 0, -2, 0, 0, -2*xn, -6*yn, 0, -6*xn*yn],
    #                 [0, 0, 0, 0, -2, 0, 0, -4*xn, -4*yn, 0, -6*xn**2, -6*yn**2]])
        
    #     # Return the [Q] coefficient matrix
    #     return Q
    
    # # Calculates and returns the gradient matrix [B]
    # def __B(self):

    #     # Calculate the gradient matrix [B]
    #     B = matmul(self.__Q(), inv(self.__C()))

    #     # Return the gradient matrix [B]
    #     return B
    
    # # Calculates and returns the constitutive matrix for isotropic materials [D]
    # def __D(self):

    #     # Calculate the coefficient for the constitutive matrix [D]
    #     C = self.E*self.t**3/(12*(1-self.nu**2))

    #     # Calculate the constitutive matrix [D]
    #     D = matrix([[C, C*self.nu, 0],
    #                 [C*self.nu, C, 0],
    #                 [0, 0, C*(1-self.nu)/2]])

    #     # Return the constitutive matrix [D]
    #     return D