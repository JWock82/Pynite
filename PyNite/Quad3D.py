# This is an isoparametric general quad element. The bending portion is based on the MITC4
# formulation. This element performs well for most basic structural and mechanical engineering
# problems, even when distorted, and eliminates the "shear locking" problems that can occur with
# some other plate bending elements.

# References used to derive this element:
# 1. "Finite Element Procedures, 2nd Edition", Klaus-Jurgen Bathe
# 2. "A First Course in the Finite Element Method, 4th Edition", Daryl L. Logan

from numpy import array, arccos, dot, cross, matmul, add, zeros, transpose
from numpy.linalg import inv, det, norm
from math import sin, cos
from PyNite.LoadCombo import LoadCombo

class Quad3D():

#%%
    def __init__(self, Name, i_node, j_node, m_node, n_node, t, E, nu,
                 LoadCombos={'Combo 1':LoadCombo('Combo 1', factors={'Case 1':1.0})}):

        self.Name = Name
        self.ID = None
        self.type = 'Quad'

        self.i_node = i_node
        self.j_node = j_node
        self.m_node = m_node
        self.n_node = n_node

        self.t = t
        self.E = E
        self.nu = nu

        self.pressures = []  # A list of surface pressures [pressure, case='Case 1']
        self.LoadCombos = LoadCombos

#%%
    def _local_coords(self):
        '''
        Calculates or recalculates and stores the local (x, y) coordinates for each node of the
        quadrilateral.
        '''

        # Get the global coordinates for each node
        X1, Y1, Z1 = self.m_node.X, self.m_node.Y, self.m_node.Z
        X2, Y2, Z2 = self.n_node.X, self.n_node.Y, self.n_node.Z
        X3, Y3, Z3 = self.i_node.X, self.i_node.Y, self.i_node.Z
        X4, Y4, Z4 = self.j_node.X, self.j_node.Y, self.j_node.Z

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

        # Convert the x and y axes into unit vectors
        x_axis = x_axis/norm(x_axis)
        y_axis = y_axis/norm(y_axis)

        # Calculate the local (x, y) coordinates for each node
        self.x1 = dot(vector_31, x_axis)
        self.x2 = dot(vector_32, x_axis)
        self.x3 = 0
        self.x4 = dot(vector_34, x_axis)
        self.y1 = dot(vector_31, y_axis)
        self.y2 = dot(vector_32, y_axis)
        self.y3 = 0
        self.y4 = dot(vector_34, y_axis)

#%%
    def J(self, r, s):
        '''
        Returns the Jacobian matrix for the element
        '''
        
        # Get the local coordinates for the element
        x1, y1, x2, y2, x3, y3, x4, y4 = self.x1, self.y1, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4

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
        dH = matmul(inv(self.J(r, s)), 1/4*array([[1 + s, -1 - s, -1 + s,  1 - s],                 
                                                  [1 + r,  1 - r, -1 + r, -1 - r]]))
        
        # Row 1 = d(beta_x)/dx divided by the local displacement vector 'u'
        # Row 2 = d(beta_y)/dy divided by the local displacement vector 'u'
        # Row 3 = d(beta_x)/dy + d(beta_y)/dx divided by the local displacement vector 'u'
        # Note that beta_x is a function of -theta_y and beta_y is a function of +theta_x (Equations 5.99, p. 423)
        B_kappa = array([[0,    0,     -dH[0, 0], 0,    0,     -dH[0, 1], 0,    0,     -dH[0, 2], 0,    0,     -dH[0, 3]],
                         [0, dH[1, 0],     0,     0, dH[1, 1],     0,     0, dH[1, 2],     0,     0, dH[1, 3],     0    ],
                         [0, dH[0, 0], -dH[1, 0], 0, dH[0, 1], -dH[1, 1], 0, dH[0, 2], -dH[1, 2], 0, dH[0, 3], -dH[1, 3]]])
        
        # Below is the matrix derived from the 1984 version of the MITC4 element. It appears to be
        # the same, but with a different sign convention for the section curvatures.
        # B_kappa = array([[0,     0,     dH[0, 0],  0,     0,     dH[0, 1],  0,     0,     dH[0, 2],  0,     0,     dH[0, 3]],
        #                  [0,  dH[1, 0],     0,     0,  dH[1, 1],     0,     0,  dH[1, 2],     0,     0,  dH[1, 3],     0   ],
        #                  [0, -dH[0, 0], dH[1, 0],  0, -dH[0, 1], dH[1, 1],  0, -dH[0, 2], dH[1, 2],  0, -dH[0, 3], dH[1, 3]]])

        return B_kappa

#%%
    def B_gamma(self, r, s):
        '''
        Returns the [B] matrix for shear.

        This is provided for reference only and is not actually used by
        PyNite. This is the theoretical solution, but it is known to
        produce spurious shear forces. It is prone to a phenomenon called
        shear locking. Instead of this matrix, the MITC4 [B] matrix is used,
        which eliminates shear-locking and can be used for thick and thin
        plates.
        '''

        H = 1/4*array([(1 + r)*(1 + s), (1 - r)*(1 + s), (1 - r)*(1 - s), (1 + r)*(1 - s)])

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with respect to r and s
        # to derivatives with respect to x and y
        dH = matmul(inv(self.J(r, s)), 1/4*array([[1 + s, -1 - s, -1 + s,  1 - s],                 
                                                  [1 + r,  1 - r, -1 + r, -1 - r]]))

        # Row 1 = d(beta_x)/dx divided by the local displacement vector 'u'
        # Row 2 = d(beta_y)/dy divided by the local displacement vector 'u'
        # Row 3 = d(beta_x)/dy + d(beta_y)/dx divided by the local displacement vector 'u'
        # Note that beta_x is a function of -theta_y and beta_y is a function of +theta_x (Equations 5.99, p. 423)
        B_gamma = array([[dH[0, 0],   0,   H[0], dH[0, 1],   0,   H[1], dH[0, 2],   0,   H[2], dH[0, 3],   0,   H[3]],
                         [dH[1, 0], -H[0],  0,   dH[1, 1], -H[1],  0,   dH[1, 2], -H[2],  0,   dH[1, 3], -H[3],  0  ]])
        
        return B_gamma

#%%
    def B_gamma_MITC4(self, r, s):
        '''
        Returns the [B] matrix for shear.

        MITC stands for mixed interpolation tensoral components. MITC elements
        are used in many programs and are known to perform well for thick and
        thin plates, and for distorted plate geometries.
        '''

        # Get the local coordinates for the element
        x1, y1, x2, y2, x3, y3, x4, y4 = self.x1, self.y1, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4
        x_axis = array([1, 0, 0]).T

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
        s_axis = array([(x1 + x2)/2 - (x3 + x4)/2, (y1 + y2)/2 - (y3 + y4)/2, 0]).T

        r_axis = r_axis/norm(r_axis)
        s_axis = s_axis/norm(s_axis)

        alpha = arccos(dot(r_axis, x_axis))
        beta = arccos(dot(s_axis, x_axis))
        # alpha = atan(Ay/Ax)
        # beta = pi/2 - atan(Cx/Cy)
        
        # Reference 1, Equations 5.103 and 5.104 (p. 426)
        det_J = det(self.J(r, s))

        gr = ((Cx + r*Bx)**2 + (Cy + r*By)**2)**0.5/(8*det_J)
        gs = ((Ax + s*Bx)**2 + (Ay + s*By)**2)**0.5/(8*det_J)

        # d        =           [    w1           theta_x1             theta_y1             w2            theta_x2              theta_y2            w3             theta_x3             theta_y3         w4             theta_x4             theta_y4       ]
        B_gamma_rz = gr*array([[(1 + s)/2, -(y1 - y2)/4*(1 + s), (x1 - x2)/4*(1 + s), -(1 + s)/2,  -(y1 - y2)/4*(1 + s), (x1 - x2)/4*(1 + s), -(1 - s)/2, -(y4 - y3)/4*(1 - s), (x4 - x3)/4*(1 - s), (1 - s)/2,  -(y4 - y3)/4*(1 - s), (x4 - x3)/4*(1 - s)]])
        B_gamma_sz = gs*array([[(1 + r)/2, -(y1 - y4)/4*(1 + r), (x1 - x4)/4*(1 + r),  (1 - r)/2,  -(y2 - y3)/4*(1 - r), (x2 - x3)/4*(1 - r), -(1 - r)/2, -(y2 - y3)/4*(1 - r), (x2 - x3)/4*(1 - r), -(1 + r)/2, -(y1 - y4)/4*(1 + r), (x1 - x4)/4*(1 + r)]])
        
        # Reference 1, Equations 5.102
        B_gamma_MITC4 = zeros((2, 12))
        B_gamma_MITC4[0, :] = B_gamma_rz*sin(beta) - B_gamma_sz*sin(alpha)
        B_gamma_MITC4[1, :] = -B_gamma_rz*cos(beta) + B_gamma_sz*cos(alpha)
        
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
        # Reference 1, Equations (5.97), page 422
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

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the determinant of the Jacobian matrix for each gauss pointing 
        # Doing this now will save us from doing it twice below
        J1 = det(self.J(gp, gp))
        J2 = det(self.J(-gp, gp))
        J3 = det(self.J(-gp, -gp))
        J4 = det(self.J(gp, -gp))

        # Get the bending B matrices for each gauss point
        B1 = self.B_kappa(gp, gp)
        B2 = self.B_kappa(-gp, gp)
        B3 = self.B_kappa(-gp, -gp)
        B4 = self.B_kappa(gp, -gp)

        # Create the stiffness matrix with bending stiffness terms
        # See Reference 1, Equation 5.94
        k = (matmul(B1.T, matmul(Cb, B1))*J1 +
             matmul(B2.T, matmul(Cb, B2))*J2 +
             matmul(B3.T, matmul(Cb, B3))*J3 +
             matmul(B4.T, matmul(Cb, B4))*J4)

        # Get the MITC4 shear B matrices for each gauss point
        B1 = self.B_gamma_MITC4(gp, gp)
        B2 = self.B_gamma_MITC4(-gp, gp)
        B3 = self.B_gamma_MITC4(-gp, -gp)
        B4 = self.B_gamma_MITC4(gp, -gp)
        
        # Alternatively the shear B matrix below could be used. However, this matrix is prone to
        # shear locking and will overestimate the stiffness.
        # B1 = self.B_gamma(gp, gp)
        # B2 = self.B_gamma(-gp, gp)
        # B3 = self.B_gamma(-gp, -gp)
        # B4 = self.B_gamma(gp, -gp)

        # Add shear stiffness terms to the stiffness matrix
        k += (matmul(B1.T, matmul(Cs, B1))*J1 +
              matmul(B2.T, matmul(Cs, B2))*J2 +
              matmul(B3.T, matmul(Cs, B3))*J3 +
              matmul(B4.T, matmul(Cs, B4))*J4)
        
        # Calculate the stiffness of a weak spring for the drilling degree of freedom (rotation about local z)
        k_rz = min(abs(k[1, 1]), abs(k[2, 2]), abs(k[4, 4]), abs(k[5, 5]),
                   abs(k[7, 7]), abs(k[8, 8]), abs(k[10, 10]), abs(k[11, 11])
                   )/1000
        
        # Initialize the expanded stiffness matrix to all zeros
        k_exp = zeros((24, 24))

        # Step through each term in the unexpanded stiffness matrix
        # i = Unexpanded matrix row
        for i in range(12):

            # j = Unexpanded matrix column
            for j in range(12):
                
                # Find the corresponding term in the expanded stiffness
                # matrix

                # m = Expanded matrix row
                if i in [0, 3, 6, 9]:  # indices associated with deflection in z
                    m = 2*i + 2
                if i in [1, 4, 7, 10]:  # indices associated with rotation about x
                    m = 2*i + 1
                if i in [2, 5, 8, 11]:  # indices associated with rotation about y
                    m = 2*i

                # n = Expanded matrix column
                if j in [0, 3, 6, 9]:  # indices associated with deflection in z
                    n = 2*j + 2
                if j in [1, 4, 7, 10]:  # indices associated with rotation about x
                    n = 2*j + 1
                if j in [2, 5, 8, 11]:  # indices associated with rotation about y
                    n = 2*j
                
                # Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                # Add the term from the unexpanded matrix into the expanded
                # matrix
                k_exp[m, n] = k[i, j]

        # Add the drilling degree of freedom's weak spring
        k_exp[5, 5] = k_rz
        k_exp[11, 11] = k_rz
        k_exp[17, 17] = k_rz
        k_exp[23, 23] = k_rz

        return k_exp

#%%
    def k_m(self):
        '''
        Returns the local stiffness matrix for membrane (in-plane) stresses.

        Plane stress is assumed
        '''

        t = self.t
        C = self.C()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the membrane B matrices for each gauss point
        # Doing this now will save us from doing it twice below
        B1 = self.B_m(gp, gp)
        B2 = self.B_m(-gp, gp)
        B3 = self.B_m(-gp, -gp)
        B4 = self.B_m(gp, -gp)

        # See reference 1 at the bottom of page 353, and reference 2 page 466
        k = t*(matmul(B1.T, matmul(C, B1))*det(self.J(gp, gp)) +
               matmul(B2.T, matmul(C, B2))*det(self.J(-gp, gp)) +
               matmul(B3.T, matmul(C, B3))*det(self.J(-gp, -gp)) +
               matmul(B4.T, matmul(C, B4))*det(self.J(gp, -gp)))
        
        k_exp = zeros((24, 24))

        # Step through each term in the unexpanded stiffness matrix
        # i = Unexpanded matrix row
        for i in range(8):

            # j = Unexpanded matrix column
            for j in range(8):
                
                # Find the corresponding term in the expanded stiffness
                # matrix

                # m = Expanded matrix row
                if i in [0, 2, 4, 6]:  # indices associated with displacement in x
                    m = i*3
                if i in [1, 3, 5, 7]:  # indices associated with displacement in y
                    m = i*3 - 2

                # n = Expanded matrix column
                if j in [0, 2, 4, 6]:  # indices associated with displacement in x
                    n = j*3
                if j in [1, 3, 5, 7]:  # indices associated with displacement in y
                    n = j*3 - 2
                
                # Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                # Add the term from the unexpanded matrix into the expanded matrix
                k_exp[m, n] = k[i, j]
        
        return k_exp

#%%
    def k(self):
        '''
        Returns the quad element's local stiffness matrix.
        '''

        # Recalculate the local coordinate system
        self._local_coords()

        # Sum the bending and membrane stiffness matrices
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
        Returns the quadrilateral's local fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the consistent load vector for.
        '''
        
        Hw = lambda r, s : 1/4*array([[(1 + r)*(1 + s), 0, 0, (1 - r)*(1 + s), 0, 0, (1 - r)*(1 - s), 0, 0, (1 + r)*(1 - s), 0, 0]])

        # Initialize the fixed end reaction vector
        fer = zeros((12,1))

        # Get the requested load combination
        combo = self.LoadCombos[combo_name]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Initialize the element's surface pressure to zero
        p = 0

        # Loop through each load case and factor in the load combination
        for case, factor in combo.factors.items():

            # Sum the pressures
            for pressure in self.pressures:

                # Check if the current pressure corresponds to the current load case
                if pressure[1] == case:

                    # Sum the pressures
                    p -= factor*pressure[0]
        
        fer = (Hw(-gp, -gp).T*p*det(self.J(-gp, -gp))
               + Hw(-gp, gp).T*p*det(self.J(-gp, gp))
               + Hw(gp, gp).T*p*det(self.J(gp, gp))
               + Hw(gp, -gp).T*p*det(self.J(gp, -gp)))

        # Initialize the expanded vector to all zeros
        fer_exp = zeros((24, 1))

        # Step through each term in the unexpanded vector
        # i = Unexpanded vector row
        for i in range(12):
                
            # Find the corresponding term in the expanded vector

            # m = Expanded vector row
            if i in [0, 3, 6, 9]:   # indices associated with deflection in z
                m = 2*i + 2
            if i in [1, 4, 7, 10]:  # indices associated with rotation about x
                m = 2*i + 1
            if i in [2, 5, 8, 11]:  # indices associated with rotation about y
                m = 2*i
                
            # Ensure the index is an integer rather than a float
            m = round(m)

            # Add the term from the unexpanded vector into the expanded vector
            fer_exp[m, 0] = fer[i, 0]

        return fer_exp

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
        D.itemset((0, 0), self.m_node.DX[combo_name])
        D.itemset((1, 0), self.m_node.DY[combo_name])
        D.itemset((2, 0), self.m_node.DZ[combo_name])
        D.itemset((3, 0), self.m_node.RX[combo_name])
        D.itemset((4, 0), self.m_node.RY[combo_name])
        D.itemset((5, 0), self.m_node.RZ[combo_name])

        D.itemset((6, 0), self.n_node.DX[combo_name])
        D.itemset((7, 0), self.n_node.DY[combo_name])
        D.itemset((8, 0), self.n_node.DZ[combo_name])
        D.itemset((9, 0), self.n_node.RX[combo_name])
        D.itemset((10, 0), self.n_node.RY[combo_name])
        D.itemset((11, 0), self.n_node.RZ[combo_name])

        D.itemset((12, 0), self.i_node.DX[combo_name])
        D.itemset((13, 0), self.i_node.DY[combo_name])
        D.itemset((14, 0), self.i_node.DZ[combo_name])
        D.itemset((15, 0), self.i_node.RX[combo_name])
        D.itemset((16, 0), self.i_node.RY[combo_name])
        D.itemset((17, 0), self.i_node.RZ[combo_name])

        D.itemset((18, 0), self.j_node.DX[combo_name])
        D.itemset((19, 0), self.j_node.DY[combo_name])
        D.itemset((20, 0), self.j_node.DZ[combo_name])
        D.itemset((21, 0), self.j_node.RX[combo_name])
        D.itemset((22, 0), self.j_node.RY[combo_name])
        D.itemset((23, 0), self.j_node.RZ[combo_name])
        
        # Return the global displacement vector
        return D

#%%
    def K(self):
        '''
        Returns the quad element's global stiffness matrix
        '''

        # Get the transpose matrix
        T = self.T()

        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(transpose(T), self.k()), T)

#%% 
    # Global fixed end reaction vector
    def FER(self, combo_name='Combo 1'):
        '''
        Returns the global fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the fixed end
            reaction vector for (not the load combination itself).
        '''
        
        # Calculate and return the fixed end reaction vector
        return matmul(inv(self.T()), self.fer(combo_name))

#%%  
    def T(self):
        '''
        Returns the coordinate transformation matrix for the quad element.
        '''

        xi = self.i_node.X
        xj = self.j_node.X
        yi = self.i_node.Y
        yj = self.j_node.Y
        zi = self.i_node.Z
        zj = self.j_node.Z

        # Calculate the direction cosines for the local x-axis.The local x-axis will run from
        # the i-node to the j-node
        x = [xj - xi, yj - yi, zj - zi]

        # Divide the vector by its magnitude to produce a unit x-vector of
        # direction cosines
        mag = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
        x = [x[0]/mag, x[1]/mag, x[2]/mag]
        
        # The local y-axis will be in the plane of the plate. Find a vector in
        # the plate's local xy plane.
        xn = self.n_node.X
        yn = self.n_node.Y
        zn = self.n_node.Z
        xy = [xn - xi, yn - yi, zn - zi]

        # Find a vector perpendicular to the plate surface to get the
        # orientation of the local z-axis.
        z = cross(x, xy)
        
        # Divide the z-vector by its magnitude to produce a unit z-vector of
        # direction cosines.
        mag = (z[0]**2 + z[1]**2 + z[2]**2)**0.5
        z = [z[0]/mag, z[1]/mag, z[2]/mag]

        # Calculate the local y-axis as a vector perpendicular to the local z
        # and x-axes.
        y = cross(z, x)
        
        # Divide the y-vector by its magnitude to produce a unit vector of
        # direction cosines.
        mag = (y[0]**2 + y[1]**2 + y[2]**2)**0.5
        y = [y[0]/mag, y[1]/mag, y[2]/mag]

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

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Define extrapolated r and s points
        r_ex = r/gp
        s_ex = s/gp

        # Define the interpolation functions
        H = 1/4*array([(1 + r_ex)*(1 + s_ex), (1 - r_ex)*(1 + s_ex), (1 - r_ex)*(1 - s_ex), (1 + r_ex)*(1 - s_ex)])

        # Get the stress-strain matrix
        Cs = self.Cs()

        # Calculate the internal shears [Qx, Qy] at each gauss point
        q1 = matmul(Cs, matmul(self.B_gamma_MITC4(gp, gp), d))
        q2 = matmul(Cs, matmul(self.B_gamma_MITC4(-gp, gp), d))
        q3 = matmul(Cs, matmul(self.B_gamma_MITC4(-gp, -gp), d))
        q4 = matmul(Cs, matmul(self.B_gamma_MITC4(gp, -gp), d))

        # Extrapolate to get the value at the requested location
        Qx = H[0]*q1[0] + H[1]*q2[0] + H[2]*q3[0] + H[3]*q4[0]
        Qy = H[0]*q1[1] + H[1]*q2[1] + H[2]*q3[1] + H[3]*q4[1]

        return array([Qx,
                      Qy])

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

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # # Define extrapolated r and s points
        r_ex = r/gp
        s_ex = s/gp

        # Define the interpolation functions
        H = 1/4*array([(1 + r_ex)*(1 + s_ex), (1 - r_ex)*(1 + s_ex), (1 - r_ex)*(1 - s_ex), (1 + r_ex)*(1 - s_ex)])

        # Get the stress-strain matrix
        Cb = self.Cb()

        # Calculate the internal moments [Mx, My, Mxy] at each gauss point
        m1 = matmul(Cb, matmul(self.B_kappa(gp, gp), d))
        m2 = matmul(Cb, matmul(self.B_kappa(-gp, gp), d))
        m3 = matmul(Cb, matmul(self.B_kappa(-gp, -gp), d))
        m4 = matmul(Cb, matmul(self.B_kappa(gp, -gp), d))

        # Extrapolate to get the value at the requested location
        Mx = H[0]*m1[0] + H[1]*m2[0] + H[2]*m3[0] + H[3]*m4[0]
        My = H[0]*m1[1] + H[1]*m2[1] + H[2]*m3[1] + H[3]*m4[1]
        Mxy = H[0]*m1[2] + H[1]*m2[2] + H[2]*m3[2] + H[3]*m4[2]
        
        return array([Mx,
                      My,
                      Mxy])

#%%
    def membrane(self, r=0, s=0, combo_name='Combo 1'):
        
        # Get the plate's local displacement vector
        # Slice out terms not related to membrane forces
        d = self.d(combo_name)[[0, 1, 6, 7, 12, 13, 18, 19], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Define extrapolated r and s points
        r_ex = r/gp
        s_ex = s/gp

        # Define the interpolation functions
        H = 1/4*array([(1 + r_ex)*(1 + s_ex), (1 - r_ex)*(1 + s_ex), (1 - r_ex)*(1 - s_ex), (1 + r_ex)*(1 - s_ex)])

        # Get the stress-strain matrix
        C = self.C()
        
        # Calculate the internal stresses [Sx, Sy, Txy] at each gauss point
        s1 = matmul(C, matmul(self.B_m(gp, gp), d))
        s2 = matmul(C, matmul(self.B_m(-gp, gp), d))
        s3 = matmul(C, matmul(self.B_m(-gp, -gp), d))
        s4 = matmul(C, matmul(self.B_m(gp, -gp), d))

        # Extrapolate to get the value at the requested location
        Sx = H[0]*s1[0] + H[1]*s2[0] + H[2]*s3[0] + H[3]*s4[0]
        Sy = H[0]*s1[1] + H[1]*s2[1] + H[2]*s3[1] + H[3]*s4[1]
        Txy = H[0]*s1[2] + H[1]*s2[2] + H[2]*s3[2] + H[3]*s4[2]

        return array([Sx,
                      Sy,
                      Txy])
