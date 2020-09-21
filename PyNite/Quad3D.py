# This is an isoparametric general quad element based on the MITC4
# formulation. This element performs well for most basic structural and
# mechanical engineering problems, even when distorted, and eliminates
# the "shear locking" problems that can occur with some other plate bending
# elements.

# References used to derive this element:
# 1. "Finite Element Procedures, 2nd Edition", Klaus-Jurgen Bathe
# 2. "A First Course in the Finite Element Method, 4th Edition",
#    Daryl L. Logan

from numpy import array, arccos, dot, cross, matmul, add, zeros, insert
from numpy.linalg import inv, det, norm
from math import atan, sin, cos

class Quad3D():

#%%
    def __init__(self, Name, iNode, jNode, mNode, nNode, t, E, nu):

        self.Name = Name
        self.ID = None

        self.iNode = iNode
        self.jNode = jNode
        self.mNode = mNode
        self.nNode = nNode

        self.t = t
        self.E = E
        self.nu = nu

#%%
    def J(self, r, s):
        '''
        Returns the Jacobian matrix for the element
        '''
        x1, y1 = self.mNode.X, self.mNode.Y
        x2, y2 = self.nNode.X, self.nNode.Y
        x3, y3 = self.iNode.X, self.iNode.Y
        x4, y4 = self.jNode.X, self.jNode.Y

        # Return the Jacobian matrix
        return 1/4*array([[x1*(s + 1) - x2*(s + 1) + x3*(s - 1) - x4*(s - 1), y1*(s + 1) - y2*(s + 1) + y3*(s - 1) - y4*(s - 1)],
                          [x1*(r + 1) - x2*(r - 1) + x3*(r - 1) - x4*(r + 1), y1*(r + 1) - y2*(r - 1) + y3*(r - 1) - y4*(r + 1)]])

#%%
    def B_kappa(self, r, s):

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to r and s to derivatives with respect to x and y
        dH = matmul(inv(self.J(r, s)), 1/4*array([[s + 1, -s - 1, s - 1, -s + 1],                 
                                                  [r + 1, -r + 1, r - 1, -r - 1]]))
        
        # Row 1 = d(beta_x)/dx divided by the local displacement vector 'u'
        # Row 2 = d(beta_y)/dy divided by the local displacement vector 'u'
        # Row 3 = d(beta_x)/dy + d(beta_y)/dx divided by the local displacement vector 'u'
        # Note that beta_x is a function of -theta_y and beta_y is a function of +theta_x (Equations 5.99, p. 423)
        B_kappa = array([[0,    0,     -dH[0, 0], 0,    0,     -dH[0, 1], 0,    0,     -dH[0, 2], 0,    0,     -dH[0, 3]],
                         [0, dH[1, 0],     0,     0, dH[1, 1],     0,     0, dH[1, 2],     0,     0, dH[1, 3],     0   ],
                         [0, dH[0, 0], -dH[1, 0], 0, dH[0, 1], -dH[1, 1], 0, dH[0, 2], -dH[1, 2], 0, dH[0, 3], -dH[1, 3]]])

        return B_kappa
    
#%%
    def B_gamma_MITC4(self, r, s):

        # The derivation of this matrix is adapted from "Finite Element Procedures"
        X1, Y1, Z1 = self.mNode.X, self.mNode.Y, self.mNode.Z
        X2, Y2, Z2 = self.nNode.X, self.nNode.Y, self.nNode.Z
        X3, Y3, Z3 = self.iNode.X, self.iNode.Y, self.iNode.Z
        X4, Y4, Z4 = self.jNode.X, self.jNode.Y, self.jNode.Z

        # Following Reference 1, Figure 5.26, node 3 will be used as the
        # origin of the plate's local (x, y) coordinate system. Find the
        # vector from the origin to each node.
        vector_32 = array([X2 - X3, Y2 - Y3, Z2 - Z3]).T
        vector_31 = array([X1 - X3, Y1 - Y3, Z1 - Z3]).T
        vector_34 = array([X4 - X3, Y4 - Y3, Z4 - Z3]).T

        # Define the plate's local x, y, and z axes
        x_axis = vector_34
        z_axis = cross(x_axis, vector_32)
        y_axis = cross(z_axis, x_axis)

        # Convert the axes into unit vectors
        x_axis = x_axis/norm(x_axis)
        y_axis = y_axis/norm(y_axis)
        z_axis = z_axis/norm(z_axis)

        # Calculate the local (x, y) coordinates for each node
        x1 = dot(vector_31, x_axis)
        x2 = dot(vector_32, x_axis)
        x3 = 0
        x4 = dot(vector_34, x_axis)
        y1 = dot(vector_31, y_axis)
        y2 = dot(vector_32, y_axis)
        y3 = 0
        y4 = dot(vector_34, y_axis)

        # Reference 1, Equations 5.105
        Ax = x1 - x2 - x3 + x4
        Bx = x1 - x2 + x3 - x4
        Cx = x1 + x2 - x3 - x4
        Ay = y1 - y2 - y3 + y4
        By = y1 - y2 + y3 - y4
        Cy = y1 + y2 - y3 - y4

        # Find the angles between the axes of the natural coordinate system and
        # the local x-axis.
        r_axis = array([(x1 + x4)/2 - (x2 + x3)/2, (y1 + y4)/2 - (y2 + y3)/2, 0]).T
        s_axis = array([(x2 + x1)/2 - (x3 + x4)/2, (y2 + y1)/2 - (y3 + y4)/2, 0]).T

        r_axis = r_axis/norm(r_axis)
        s_axis = s_axis/norm(s_axis)

        alpha = arccos(dot(r_axis, x_axis))
        beta = arccos(dot(s_axis, x_axis))
        
        # Reference 1, Equations 5.103 and 5.104 (p. 426)
        J = self.J(r, s)

        gr = ((Cx + r*Bx)**2 + (Cy + r*By)**2)**0.5/(8*det(J)**2)
        gs = ((Ax + s*Bx)**2 + (Ay + s*By)**2)**0.5/(8*det(J)**2)

        # See Jupyter Notebook derivation for this next part
        gamma_rz = gr/4*array([[2*(s + 1), -s*y1 + s*y2 - y1 + y2, s*x1 - s*x2 + x1 - x2, 2*(-s - 1), -s*y1 + s*y2 - y1 + y2,  s*x1 - s*x2 + x1 - x2, 2*(s - 1), -s*y3 + s*y4 + y3 - y4,  s*x3 - s*x4 - x3 + x4, 2*(1 - s),  -s*y3 + s*y4 + y3 - y4, s*x3 - s*x4 - x3 + x4]])
        gamma_sz = gs/4*array([[2*(r + 1), -r*y1 + r*y4 - y1 + y4, r*x1 - r*x4 + x1 - x4, 2*(1 - r),   r*y2 - r*y3 - y2 + y3, -r*x2 + r*x3 + x2 - x3, 2*(r - 1),  r*y2 - r*y3 - y2 + y3, -r*x2 + r*x3 + x2 - x3, 2*(-r - 1), -r*y1 + r*y4 - y1 + y4, r*x1 - r*x4 + x1 - x4]])

        # Reference 1, Equations 5.102
        B_gamma_MITC4 = zeros((2, 12))
        B_gamma_MITC4[0, :] = gamma_rz*sin(beta) - gamma_sz*sin(alpha)
        B_gamma_MITC4[1, :] = -gamma_rz*cos(beta) + gamma_sz*cos(alpha)
        
        # Return the [B] matrix for shear
        return B_gamma_MITC4

#%%
    def B_m(self, r, s):

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to r and s to derivatives with respect to x and y
        dH = matmul(inv(self.J(r, s)), 1/4*array([[s + 1, -s - 1, s - 1, -s + 1],                 
                                                  [r + 1, -r + 1, r - 1, -r - 1]]))

        # Reference 1, Example 5.5 (page 353)
        B_m = 1/4*array([[dH[0, 0],    0,     dH[0, 1],    0,     dH[0, 2],    0,     dH[0, 3],    0    ],
                         [   0,     dH[1, 0],    0,     dH[1, 1],    0,     dH[1, 2],    0,     dH[1, 3]],
                         [dH[1, 0], dH[0, 0], dH[1, 1], dH[0, 1], dH[1, 2], dH[0, 2], dH[1, 3], dH[0, 3]]])

        return B_m

#%%
    def Cb(self):
        '''
        Returns the stress-strain matrix for plate bending.
        '''

        # Referemce 1, Table 4.3, page 194
        nu = self.nu
        E = self.E
        h = self.t

        Cb = E*h**3/(12*(1 - nu**2))*array([[1,  nu,      0    ],
                                            [nu, 1,       0    ],
                                            [0,  0,  (1 - nu)/2]])
        
        return Cb

#%%
    def Cs(self):
        '''
        Returns the stress-strain matrix for shear.
        '''

        k = 5/6
        E = self.E
        h = self.t
        nu = self.nu

        Cs = E*h*k/(2*(1 + nu))*array([[1, 0],
                                       [0, 1]])

        return Cs

#%%
    def C(self):
        '''
        Returns the stress-strain matrix for a plane stress.
        '''
        
        E = self.E
        nu = self.nu

        # Reference 1, Table 4.3, page 194
        C = E/(1 - nu**2)*array([[1,  nu,     0    ],
                                 [nu, 1,      0    ],
                                 [0,  0, (1 - nu)/2]])
        
        return C

#%%
    def k_b(self):
        '''
        Returns the local stiffness matrix for bending stresses
        '''

        Cb = self.Cb()
        Cs = self.Cs()

        # See Reference 1, Equation 5.94
        k = (matmul(self.B_kappa(-0.5773, -0.5773).T, matmul(Cb, self.B_kappa(-0.5773, -0.5773)))*det(self.J(-0.5773, -0.5773)) +
             matmul(self.B_kappa(-0.5773, 0.5773).T, matmul(Cb, self.B_kappa(-0.5773, 0.5773)))*det(self.J(-0.5773, 0.5773)) +
             matmul(self.B_kappa(0.5773, 0.5773).T, matmul(Cb, self.B_kappa(0.5773, 0.5773)))*det(self.J(0.5773, 0.5773)) +
             matmul(self.B_kappa(0.5773, -0.5773).T, matmul(Cb, self.B_kappa(0.5773, -0.5773)))*det(self.J(0.5773, -0.5773)))
        
        k += (matmul(self.B_gamma_MITC4(-0.5773, -0.5773).T, matmul(Cs, self.B_gamma_MITC4(-0.5773, -0.5773)))*det(self.J(-0.5773, -0.5773)) +
              matmul(self.B_gamma_MITC4(-0.5773, 0.5773).T, matmul(Cs, self.B_gamma_MITC4(-0.5773, 0.5773)))*det(self.J(-0.5773, 0.5773)) +
              matmul(self.B_gamma_MITC4(0.5773, 0.5773).T, matmul(Cs, self.B_gamma_MITC4(0.5773, 0.5773)))*det(self.J(0.5773, 0.5773)) +
              matmul(self.B_gamma_MITC4(0.5773, -0.5773).T, matmul(Cs, self.B_gamma_MITC4(0.5773, -0.5773)))*det(self.J(0.5773, -0.5773)))
        
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

        return k

#%%
    def k_m(self):
        '''
        Returns the local stiffness matrix for membrane (in-plane) stresses.

        Plane stress is assumed
        '''

        t = self.t
        C = self.C()

        # See reference 1 at the bottom of page 353, and reference 2 page 466
        k = t*(matmul(self.B_m(-0.5773, -0.5773).T, matmul(C, self.B_m(-0.5773, -0.5773)))*det(self.J(-0.5773, -0.5773)) +
               matmul(self.B_m(-0.5773, 0.5773).T, matmul(C, self.B_m(-0.5773, 0.5773)))*det(self.J(-0.5773, 0.5773)) +
               matmul(self.B_m(0.5773, 0.5773).T, matmul(C, self.B_m(0.5773, 0.5773)))*det(self.J(0.5773, 0.5773)) +
               matmul(self.B_m(0.5773, -0.5773).T, matmul(C, self.B_m(0.5773, -0.5773)))*det(self.J(0.5773, -0.5773)))
        
        # Insert rows of zeros for degrees of freedom not included in the matrix above
        k = insert(k, 8, 0, axis=0)
        k = insert(k, 8, 0, axis=1)
        k = insert(k, 8, 0, axis=0)
        k = insert(k, 8, 0, axis=1)
        k = insert(k, 8, 0, axis=0)
        k = insert(k, 8, 0, axis=1)
        k = insert(k, 8, 0, axis=0)
        k = insert(k, 8, 0, axis=1)

        k = insert(k, 6, 0, axis=0)
        k = insert(k, 6, 0, axis=1)
        k = insert(k, 6, 0, axis=0)
        k = insert(k, 6, 0, axis=1)
        k = insert(k, 6, 0, axis=0)
        k = insert(k, 6, 0, axis=1)
        k = insert(k, 6, 0, axis=0)
        k = insert(k, 6, 0, axis=1)

        k = insert(k, 4, 0, axis=0)
        k = insert(k, 4, 0, axis=1)
        k = insert(k, 4, 0, axis=0)
        k = insert(k, 4, 0, axis=1)
        k = insert(k, 4, 0, axis=0)
        k = insert(k, 4, 0, axis=1)
        k = insert(k, 4, 0, axis=0)
        k = insert(k, 4, 0, axis=1)
        
        k = insert(k, 2, 0, axis=0)
        k = insert(k, 2, 0, axis=1)
        k = insert(k, 2, 0, axis=0)
        k = insert(k, 2, 0, axis=1)
        k = insert(k, 2, 0, axis=0)
        k = insert(k, 2, 0, axis=1)
        k = insert(k, 2, 0, axis=0)
        k = insert(k, 2, 0, axis=1)
        
        return k

#%%
    def k(self):
        '''
        Returns the quad element's local stiffness matrix.
        '''

        return self.k_b() + self.k_m()

#%%   
    def f(self, combo_name='Combo 1'):
        '''
        Returns the quad element's local end force vector
        '''
        
        # Calculate and return the plate's local end force vector
        return add(matmul(self.k(), self.d(combo_name)), self.fer(combo_name))

#%%
    def fer(self, combo_name='Combo 1'):
        '''
        Returns the quad element's local fixed end reaction vector (zero's for
        now until surface loads get added).
        '''

        return zeros((24, 1))

#%%
    def d(self, combo_name='Combo 1'):
       '''
       Returns the quad element's local displacement vector
       '''

       # Calculate and return the local displacement vector
       return matmul(self.T(), self.D(combo_name))

#%%
    def F(self, combo_name='Combo 1'):
        '''
        Returns the quad element's global force vector
        '''
        
        # Calculate and return the global force vector
        return matmul(inv(self.T()), self.f(combo_name))

#%%
    def D(self, combo_name='Combo 1'):
        '''
        Returns the quad element's global displacement vector for the given
        load combination.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the displacement vector
            for (not the load combination itself).
        '''
        
        # Initialize the displacement vector
        D = zeros((24, 1))
        
        # Read in the global displacements from the nodes
        D.itemset((0, 0), self.mNode.DX[combo_name])
        D.itemset((1, 0), self.mNode.DY[combo_name])
        D.itemset((2, 0), self.mNode.DZ[combo_name])
        D.itemset((3, 0), self.mNode.RX[combo_name])
        D.itemset((4, 0), self.mNode.RY[combo_name])
        D.itemset((5, 0), self.mNode.RZ[combo_name])

        D.itemset((6, 0), self.nNode.DX[combo_name])
        D.itemset((7, 0), self.nNode.DY[combo_name])
        D.itemset((8, 0), self.nNode.DZ[combo_name])
        D.itemset((9, 0), self.nNode.RX[combo_name])
        D.itemset((10, 0), self.nNode.RY[combo_name])
        D.itemset((11, 0), self.nNode.RZ[combo_name])

        D.itemset((12, 0), self.iNode.DX[combo_name])
        D.itemset((13, 0), self.iNode.DY[combo_name])
        D.itemset((14, 0), self.iNode.DZ[combo_name])
        D.itemset((15, 0), self.iNode.RX[combo_name])
        D.itemset((16, 0), self.iNode.RY[combo_name])
        D.itemset((17, 0), self.iNode.RZ[combo_name])

        D.itemset((18, 0), self.jNode.DX[combo_name])
        D.itemset((19, 0), self.jNode.DY[combo_name])
        D.itemset((20, 0), self.jNode.DZ[combo_name])
        D.itemset((21, 0), self.jNode.RX[combo_name])
        D.itemset((22, 0), self.jNode.RY[combo_name])
        D.itemset((23, 0), self.jNode.RZ[combo_name])
        
        # Return the global displacement vector
        return D

#%%
    def K(self):
        '''
        Returns the quad element's global stiffness matrix
        '''
        
        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(transpose(self.T()), self.k()), self.T())

#%%  
    def T(self):
        '''
        Returns the coordinate transformation matrix for the quad element.
        '''

        # Calculate the direction cosines for the local x-axis.The local
        # x-axis will run from the i-node to the j-node.
        xi = self.iNode.X
        xj = self.jNode.X
        yi = self.iNode.Y
        yj = self.jNode.Y
        zi = self.iNode.Z
        zj = self.jNode.Z
        x = [xj - xi, yj - yi, zj - zi]

        # Divide the vector by its magnitude to produce a unit x-vector of
        # direction cosines
        x = x/norm(x)
        
        # The local y-axis will be in the plane of the plate. Find a vector in
        # the plate's local xy plane.
        xn = self.nNode.X
        yn = self.nNode.Y
        zn = self.nNode.Z
        xy = [xn-xi, yn-yi, zn-zi]

        # Find a vector perpendicular to the plate surface to get the
        # orientation of the local z-axis.
        z = cross(x, xy)
        
        # Divide the vector by its magnitude to produce a unit z-vector of
        # direction cosines.
        z = z/norm(z)

        # Calculate the local y-axis as a vector perpendicular to the local z
        # and x-axes.
        y = cross(z, x)
        
        # Divide the z-vector by its magnitude to produce a unit vector of
        # direction cosines.
        y = y/norm(y)

        # Create the direction cosines matrix.
        dirCos = array([x,
                        y,
                        z])
        
        # Build the transformation matrix.
        T = zeros((24, 24))
        T[0:3, 0:3] = dirCos
        T[3:6, 3:6] = dirCos
        T[6:9, 6:9] = dirCos
        T[9:12, 9:12] = dirCos
        T[12:15, 12:15] = dirCos
        T[15:18, 15:18] = dirCos
        T[18:21, 18:21] = dirCos
        T[21:24, 21:24] = dirCos
        
        # Return the transformation matrix.
        return T

#%%
    def shear(self, r=0, s=0, combo_name='Combo 1'):
        '''
        Returns the interal shears at any point in the quad element.

        Internal shears are reported as a 2D array [[Qx], [Qy]] at the
        specified location in the (r, s) natural coordinate system.

        Parameters
        ----------
        r : number
            The r-coordinate. Default is 0.
        s : number
            The s-coordinate. Default is 0.
        
        Returns
        -------
        Internal shear force per unit length of the quad element.
        '''

        # Get the plate's local displacement vector
        # Slice out terms not related to plate bending
        d = self.d(combo_name)[[2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22], :]

        # Calculate and return internal shears
        return matmul(self.Cs(), matmul(self.B_gamma_MITC4(r, s), d))

#%%   
    def moment(self, r=0, s=0, combo_name='Combo 1'):
        '''
        Returns the interal moments at any point in the quad element.

        Internal moments are reported as a 2D array [[Mx], [My], [Mxy]] at the
        specified location in the (r, s) natural coordinate system.

        Parameters
        ----------
        r : number
            The r-coordinate. Default is 0.
        s : number
            The s-coordinate. Default is 0.
        
        Returns
        -------
        Internal moment per unit length of the quad element.
        '''

        # Get the plate's local displacement vector
        # Slice out terms not related to plate bending
        d = self.d(combo_name)[[2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22], :]

        # Calculate and return internal moments
        return matmul(self.Cb(), matmul(self.B_kappa(r, s), d))

#%%
    def membrane(self, r=0, s=0, combo_name='Combo 1'):
        
        # Get the plate's local displacement vector
        # Slice out terms not related to membrane forces
        d = self.d(combo_name)[[0, 1, 6, 7, 12, 13, 18, 19], :]

from Node3D import Node3D
i = Node3D('i', 0, 0, 0)
j = Node3D('j', 2, 2, 0)
m = Node3D('m', 3, 5, 0)
n = Node3D('n', 1, 3, 0)
pl = Quad3D('pl', i, j, m, n, 0.25, 29000, 0.3)
print(pl.k())