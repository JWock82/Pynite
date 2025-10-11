# References used to derive this element:
# 1. "A Comparative Formulation of DKMQ, DSQ and MITC4 Quadrilateral Plate Elements with New Numerical Results Based on s-norm Tests", Irwan Katili, 
# 2. "Finite Element Procedures, 2nd Edition", Klaus-Jurgen Bathe
# 3. "A First Course in the Finite Element Method, 4th Edition", Daryl L. Logan
# 4. "Finite Element Analysis Fundamentals", Richard H. Gallagher

from __future__ import annotations  # Allows more recent type hints features
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy import add
from numpy.linalg import inv, det, norm
import warnings

if TYPE_CHECKING:
    from typing import List, Tuple, Optional
    from numpy import float64
    from numpy.typing import NDArray
    from Pynite.FEModel3D import FEModel3D
    from Pynite.Node3D import Node3D


class Quad3D():
    """
    An isoparametric general quadrilateral element, formulated by superimposing an isoparametric DKMQ bending element with an isoparametric plane stress element. Drilling stability is provided by adding a weak rotational spring stiffness at each node. Isotropic behavior is the default, but orthotropic in-plane behavior can be modeled by specifying stiffness modification factors for the element's local x and y axes.

    This element performs well for thick and thin plates, and for skewed plates. Minor errors are introduced into the solution due to the drilling approximation. Orthotropic behavior is limited to acting along the plate's local axes.
    """

    def __init__(self, name: str, i_node: Node3D, j_node: Node3D, m_node: Node3D, n_node: Node3D, 
                 t: float, material_name: str, model: FEModel3D, kx_mod: float = 1.0,
                 ky_mod: float = 1.0):

        self.name: str = name
        self.ID: Optional[int] = None
        self.type: str = 'Quad'

        self.i_node: Node3D = i_node
        self.j_node: Node3D = j_node
        self.m_node: Node3D = m_node
        self.n_node: Node3D = n_node

        self.t: float = t
        self.kx_mod: float = kx_mod
        self.ky_mod: float = ky_mod

        self.pressures: List[Tuple[float, str]] = []  # A list of surface pressures [pressure, case='Case 1']

        # Quads need a link to the model they belong to
        self.model: FEModel3D = model

        # Get material properties for the plate from the model
        try:
            self.E: float = self.model.materials[material_name].E
            self.nu: float = self.model.materials[material_name].nu
        except:
            raise KeyError('Please define the material ' + str(material_name) + ' before assigning it to plates.')

    # def _local_coords(self):
    #     """
    #     Calculates or recalculates and stores the local (x, y) coordinates for each node of the
    #     quadrilateral.
    #     """

    #     # Get the global coordinates for each node
    #     X1, Y1, Z1 = self.i_node.X, self.i_node.Y, self.i_node.Z
    #     X2, Y2, Z2 = self.j_node.X, self.j_node.Y, self.j_node.Z
    #     X3, Y3, Z3 = self.m_node.X, self.m_node.Y, self.m_node.Z
    #     X4, Y4, Z4 = self.n_node.X, self.n_node.Y, self.n_node.Z

    #     # Node 1 will be used as the origin of the plate's local (x, y) coordinate system. Find the
    #     # vector from the origin to each node.
    #     Xi1, Yi1, Zi1 = (X1 + X4)/2, (Y1 + Y4)/2, (Z1 + Z4)/2
    #     Xi2, Yi2, Zi2 = (X2 + X3)/2, (Y2 + Y3)/2, (Z2 + Z3)/2
    #     Xo, Yo, Zo = (Xi1 + Xi2)/2, (Yi1 + Yi2)/2, (Zi1 + Zi2)/2

    #     x_axis = np.array([Xi2 - Xi1, Yi2 - Yi1, Zi2 - Zi1]).T

    #     vector_01 = np.array([X1 - Xo, Y1 - Yo, Z1 - Zo]).T
    #     vector_02 = np.array([X2 - Xo, Y2 - Yo, Z2 - Zo]).T
    #     vector_03 = np.array([X3 - Xo, Y3 - Yo, Z3 - Zo]).T
    #     vector_04 = np.array([X4 - Xo, Y4 - Yo, Z4 - Zo]).T

    #     # Define the plate's local y, and z axes
    #     vector_x3 = np.array([X3 - Xi1, Y3 - Yi1, Z3 - Zi1]).T
    #     z_axis = np.cross(x_axis, vector_x3)
    #     y_axis = np.cross(z_axis, x_axis)

    #     # Convert the x, y and z axes into unit vectors
    #     x_axis = x_axis/norm(x_axis)
    #     y_axis = y_axis/norm(y_axis)
    #     z_axis = z_axis/norm(z_axis)

    #     # Calculate the local (x, y) coordinates for each node
    #     self.x1 = np.dot(vector_01, x_axis)
    #     self.x2 = np.dot(vector_02, x_axis)
    #     self.x3 = np.dot(vector_03, x_axis)
    #     self.x4 = np.dot(vector_04, x_axis)
    #     self.y1 = np.dot(vector_01, y_axis)
    #     self.y2 = np.dot(vector_02, y_axis)
    #     self.y3 = np.dot(vector_03, y_axis)
    #     self.y4 = np.dot(vector_04, y_axis)

    def _local_coords(self):
        """
        Calculates or recalculates and stores the local (x, y) coordinates for each node of the quadrilateral.
        """

        # Get the global coordinates for each node
        X1, Y1, Z1 = self.i_node.X, self.i_node.Y, self.i_node.Z
        X2, Y2, Z2 = self.j_node.X, self.j_node.Y, self.j_node.Z
        X3, Y3, Z3 = self.m_node.X, self.m_node.Y, self.m_node.Z
        X4, Y4, Z4 = self.n_node.X, self.n_node.Y, self.n_node.Z

        # Node 1 will be used as the origin of the plate's local (x, y) coordinate system. Find the
        # vector from the origin to each node.
        vector_12 = np.array([X2 - X1, Y2 - Y1, Z2 - Z1]).T
        vector_13 = np.array([X3 - X1, Y3 - Y1, Z3 - Z1]).T
        vector_14 = np.array([X4 - X1, Y4 - Y1, Z4 - Z1]).T

        # Define the plate's local x, y, and z axes
        x_axis = vector_12
        z_axis = np.cross(x_axis, vector_13)
        y_axis = np.cross(z_axis, x_axis)

        # Convert the x and y axes into unit vectors
        x_axis = x_axis/norm(x_axis)
        y_axis = y_axis/norm(y_axis)

        # Calculate the local (x, y) coordinates for each node
        self.x1: float = 0.0
        self.x2: float = np.dot(vector_12, x_axis)
        self.x3: float = np.dot(vector_13, x_axis)
        self.x4: float = np.dot(vector_14, x_axis)
        self.y1: float = 0.0
        self.y2: float = np.dot(vector_12, y_axis)
        self.y3: float = np.dot(vector_13, y_axis)
        self.y4: float = np.dot(vector_14, y_axis)

    def L_k(self, k: Literal[5, 6, 7, 8]) -> float:

        # Figures 3 and 5
        if k == 5:
            return ((self.x2 - self.x1)**2 + (self.y2 - self.y1)**2)**0.5
        elif k == 6:
            return ((self.x3 - self.x2)**2 + (self.y3 - self.y2)**2)**0.5
        elif k == 7:
            return ((self.x4 - self.x3)**2 + (self.y4 - self.y3)**2)**0.5
        elif k == 8:
            return ((self.x1 - self.x4)**2 + (self.y1 - self.y4)**2)**0.5
        else:
            raise Exception('Invalid value for k. k must be 5, 6, 7, or 8.')

    def dir_cos(self, k: Literal[5, 6, 7, 8]) -> Tuple[float, float]:

        L_k = self.L_k(k)

        # Figures 3 and 5
        if k == 5:
            C = (self.x2 - self.x1)/L_k
            S = (self.y2 - self.y1)/L_k
        elif k == 6:
            C = (self.x3 - self.x2)/L_k
            S = (self.y3 - self.y2)/L_k
        elif k == 7:
            C = (self.x4 - self.x3)/L_k
            S = (self.y4 - self.y3)/L_k
        elif k == 8:
            C = (self.x1 - self.x4)/L_k
            S = (self.y1 - self.y4)/L_k
        else:
            raise Exception('Invalid value for k. k must be 5, 6, 7, or 8.')

        return C, S

    def phi_k(self, k: Literal[5, 6, 7, 8]) -> float:

        kappa = 5/6

        # Equation 74
        return 2/(kappa*(1-self.nu))*(self.t/self.L_k(k))**2

    def N_i(self, i: Literal[1, 2, 3, 4], xi: float, eta: float) -> float:
        """
        Returns the interpolation function for any given coordinate in the natural (xi, eta) coordinate system
        """

        if i == 1:
            return 1/4*(1 - xi)*(1 - eta)
        elif i == 2:
            return 1/4*(1 + xi)*(1 - eta)
        elif i == 3:
            return 1/4*(1 + xi)*(1 + eta)
        elif i == 4:
            return 1/4*(1 - xi)*(1 + eta)
        else:
            raise Exception('Unable to calculate interpolation function. Invalid value specifed for i.')

    def P_k(self, k: Literal[5, 6, 7, 8], xi: float, eta: float) -> float:

        if k == 5:
            return 1/2*(1 - xi**2)*(1 - eta)
        elif k == 6:
            return 1/2*(1 + xi)*(1 - eta**2)
        elif k == 7:
            return 1/2*(1 - xi**2)*(1 + eta)
        elif k == 8:
            return 1/2*(1 - xi)*(1 - eta**2)
        else:
            raise Exception('Unable to calculate shape function. Invalid value specified for k.')

    def Co(self, xi: float, eta: float) -> NDArray[float64]:
        """
        This alternate calculation of the Jacobian matrix follows "The development of DKMQ plate bending element for thick to thin shell analysis based on the Naghdi/Reissner/Mindlin shell theory" by Katili, Batoz, Maknun and Hamdouni (2015). In the reference "C^o" is used instead of "J" to refer to the Jacobian. This method does not seem to produce incorrect results, but will be kept for future reference. It is helpful for understanding this plate element and may prove a useful simplification to the code base if implemented correctly someday.
        """

        # Nodal global coordinates (i=1, j=2, m=3, n=4)
        x_1 = np.array([[self.i_node.X, self.i_node.Y, self.i_node.Z]])
        x_2 = np.array([[self.j_node.X, self.j_node.Y, self.j_node.Z]])
        x_3 = np.array([[self.m_node.X, self.m_node.Y, self.m_node.Z]])
        x_4 = np.array([[self.n_node.X, self.n_node.Y, self.n_node.Z]])

        # Derivatives of the bilinear interpolation functions
        N1_xi = 0.25*(eta - 1)
        N2_xi = -0.25*(eta - 1)
        N3_xi = 0.25*(eta + 1)
        N4_xi = -0.25*(eta + 1)
        N1_eta = 0.25*(xi - 1)
        N2_eta = -0.25*(xi + 1)
        N3_eta = 0.25*(xi + 1)
        N4_eta = -0.25*(xi - 1)

        # Equation 4 - Katili 2015
        a_1 = N1_xi*x_1 + N2_xi*x_2 + N3_xi*x_3 + N4_xi*x_4
        a_2 = N1_eta*x_1 + N2_eta*x_2 + N3_eta*x_3 + N4_eta*x_4

        # Normal vector
        n = np.cross(a_1, a_2)/np.linalg.norm(np.cross(a_1, a_2))

        # Global unit vectors
        i = np.array([[1.0, 0.0, 0.0]])
        k = np.array([[0.0, 0.0, 1.0]])

        # Equation (13)
        if np.array_equal(n, k) or np.array_equal(n, -k):
            t_1 = i
        else:
            t_1 = np.cross(n, k)

        t_2 = np.cross(n, t_1)

        # Equation (7)
        a_11 = np.dot(a_1, a_1.T)[0, 0]
        a_12 = np.dot(a_1, a_2.T)[0, 0]
        a_21 = np.dot(a_2, a_1.T)[0, 0]
        a_22 = np.dot(a_2, a_2.T)[0, 0]

        # Matrix tensor of the middle surface
        a = np.array([[a_11, a_12],
                      [a_21, a_22]])

        a_det = np.linalg.det(a)

        # Contravariant vectors
        a1 = 1/a_det*(a_22*a_1 - a_12*a_2)
        a2 = 1/a_det*(-a_21*a_1 + a_11*a_2)

        Co = np.array([[np.dot(a1, t_1.T)[0, 0], np.dot(a1, t_2.T)[0, 0]],
                       [np.dot(a2, t_1.T)[0, 0], np.dot(a2, t_2.T)[0, 0]]])

        return Co

    def J(self, xi: float, eta: float) -> NDArray[float64]:
        """
        Returns the Jacobian matrix for the element
        """

        # Get the local coordinates for the element
        x1, y1, x2, y2, x3, y3, x4, y4 = self.x1, self.y1, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4

        # Return the Jacobian matrix
        return 1/4*np.array([[x1*(eta - 1) - x2*(eta - 1) + x3*(eta + 1) - x4*(eta + 1), y1*(eta - 1) - y2*(eta - 1) + y3*(eta + 1) - y4*(eta + 1)],
                             [x1*(xi - 1)  - x2*(xi + 1)  + x3*(xi + 1)  - x4*(xi - 1),  y1*(xi - 1)  - y2*(xi + 1)  + y3*(xi + 1)  - y4*(xi - 1)]])

    def N_gamma(self, xi: float, eta: float) -> NDArray[float64]:

        # Equation 44
        return np.array([[1/2*(1 - eta),       0,      1/2*(1 + eta),       0     ],
                         [     0,        1/2*(1 + xi),       0,       1/2*(1 - xi)]])

    def A_gamma(self) -> NDArray[float64]:

        L5 = self.L_k(5)
        L6 = self.L_k(6)
        L7 = self.L_k(7)
        L8 = self.L_k(8)

        # Equation 46
        return np.array([[L5/2,   0,     0,     0 ],
                         [  0,  L6/2,    0,     0 ],
                         [  0,    0,  -L7/2,    0 ],
                         [  0,    0,     0,  -L8/2]])

    def A_u(self) -> NDArray[float64]:

        # Calculate the length of each side of the quad
        L5 = self.L_k(5)
        L6 = self.L_k(6)
        L7 = self.L_k(7)
        L8 = self.L_k(8)

        # Get the direction cosines for each side of the quad
        C5, S5 = self.dir_cos(5)
        C6, S6 = self.dir_cos(6)
        C7, S7 = self.dir_cos(7)
        C8, S8 = self.dir_cos(8)

        # Return the [A_u] matrix
        return 1/2*np.array([[-2/L5, C5, S5,  2/L5, C5, S5,   0,    0,  0,   0,    0,  0],
                             [  0,    0,  0, -2/L6, C6, S6, 2/L6,  C6, S6,   0,    0,  0],
                             [  0,    0,  0,   0,    0,  0, -2/L7, C7, S7,  2/L7, C7, S7],
                             [ 2/L8, C8, S8,   0,    0,  0,   0,    0,  0, -2/L8, C8, S8]])

    def A_Delta_inv_DKMQ(self) -> NDArray[float64]:

        phi5 = self.phi_k(5)
        phi6 = self.phi_k(6)
        phi7 = self.phi_k(7)
        phi8 = self.phi_k(8)

        return -3/2*np.array([[1/(1+phi5),     0,          0,          0     ],
                              [   0,       1/(1+phi6),     0,          0     ],
                              [   0,           0,      1/(1+phi7),     0     ],
                              [   0,           0,          0,      1/(1+phi8)]])

    def A_phi_Delta(self) -> NDArray[float64]:

        phi5 = self.phi_k(5)
        phi6 = self.phi_k(6)
        phi7 = self.phi_k(7)
        phi8 = self.phi_k(8)

        return np.array([[phi5/(1+phi5),       0,             0,             0      ],
                         [      0,       phi6/(1+phi6),       0,             0      ],
                         [      0,             0,       phi7/(1+phi7),       0      ],
                         [      0,             0,             0,       phi8/(1+phi8)]])
    
    def B_b_beta(self, xi: float, eta: float) -> NDArray[float64]:

        # Get the inverse of the Jacobian matrix
        J_inv = inv(self.J(xi, eta))

        # Get the individual terms for the Jacobian inverse
        j11 = J_inv[0, 0]
        j12 = J_inv[0, 1]
        j21 = J_inv[1, 0]
        j22 = J_inv[1, 1]

        # Derivatives of the bilinear interpolation functions
        N1_xi = 0.25*(eta - 1)
        N2_xi = -0.25*(eta - 1)
        N3_xi = 0.25*(eta + 1)
        N4_xi = -0.25*(eta + 1)
        N1_eta = 0.25*(xi - 1)
        N2_eta = -0.25*(xi + 1)
        N3_eta = 0.25*(xi + 1)
        N4_eta = -0.25*(xi - 1)

        N1x = j11*N1_xi + j12*N1_eta
        N1y = j21*N1_xi + j22*N1_eta
        N2x = j11*N2_xi + j12*N2_eta
        N2y = j21*N2_xi + j22*N2_eta
        N3x = j11*N3_xi + j12*N3_eta
        N3y = j21*N3_xi + j22*N3_eta
        N4x = j11*N4_xi + j12*N4_eta
        N4y = j21*N4_xi + j22*N4_eta

        return np.array([[0, N1x,  0,  0, N2x,  0,  0, N3x,  0,  0, N4x,  0 ],
                         [0,  0,  N1y, 0,  0,  N2y, 0,  0,  N3y, 0,  0,  N4y],
                         [0, N1y, N1x, 0, N2y, N2x, 0, N3y, N3x, 0, N4y, N4x]])
    
    def B_b_Delta_beta(self, xi: float, eta: float) -> NDArray[float64]:

        # Get the inverse of the Jacobian matrix
        J_inv = inv(self.J(xi, eta))

        # Get the individual terms for the Jacobian inverse
        j11 = J_inv[0, 0]
        j12 = J_inv[0, 1]
        j21 = J_inv[1, 0]
        j22 = J_inv[1, 1]

        # Derivatives of the quadratic interpolation functions
        P5_xi = xi*(eta - 1)
        P6_xi = -0.5*(eta - 1)*(eta + 1)
        P7_xi = -xi*(eta + 1)
        P8_xi = 0.5*(eta - 1)*(eta + 1)
        P5_eta = 0.5*(xi - 1)*(xi + 1)
        P6_eta = -eta*(xi + 1)
        P7_eta = -0.5*(xi - 1)*(xi + 1)
        P8_eta = eta*(xi - 1)

        P5x = j11*P5_xi + j12*P5_eta
        P5y = j21*P5_xi + j22*P5_eta
        P6x = j11*P6_xi + j12*P6_eta
        P6y = j21*P6_xi + j22*P6_eta
        P7x = j11*P7_xi + j12*P7_eta
        P7y = j21*P7_xi + j22*P7_eta
        P8x = j11*P8_xi + j12*P8_eta
        P8y = j21*P8_xi + j22*P8_eta

        C5, S5 = self.dir_cos(5)
        C6, S6 = self.dir_cos(6)
        C7, S7 = self.dir_cos(7)
        C8, S8 = self.dir_cos(8)

        return np.array([[    P5x*C5,          P6x*C6,          P7x*C7,          P8x*C8     ],
                         [    P5y*S5,          P6y*S6,          P7y*S7,          P8y*S8,    ],
                         [P5y*C5 + P5x*S5, P6y*C6 + P6x*S6, P7y*C7 + P7x*S7, P8y*C8 + P8x*S8]])
    
    def B_b(self, xi: float, eta: float) -> NDArray[float64]:
        """
        Returns the [B_b] matrix for bending
        """

        # Return the [B] matrix for bending
        return add(self.B_b_beta(xi, eta), self.B_b_Delta_beta(xi, eta) @ self.A_Delta_inv_DKMQ() @ self.A_u())

    def B_s(self, xi: float, eta: float) -> NDArray[float64]:
        """
        Returns the [B_s] matrix for shear
        """
        
        # Return the [B] matrix for shear
        return inv(self.J(xi, eta)) @ self.N_gamma(xi, eta) @ self.A_gamma() @ self.A_phi_Delta() @ self.A_u()

    def B_s_gamma(self, xi:float , eta: float) -> None:
        """Returns the [B_s_gamma] matrix for shear (Equation 39 in Reference 1)

        :param xi: Isoparametric coordinate xi (-1 to 1)
        :type xi: _type_
        :param eta: Isoparametric coordinate eta (-1 to 1)
        :type eta: _type_
        """
        raise NotImplementedError('This function is not implemented yet. It is not needed for the current implementation of the Quad3D element.')

    def B_m(self, xi: float, eta: float) -> NDArray[float64]:

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to xi and eta to derivatives with respect to x and y
        dH = np.matmul(inv(self.J(xi, eta)), 1/4*np.array([[eta - 1, -eta + 1, eta + 1, -eta - 1],                 
                                                           [xi - 1,  -xi - 1,  xi + 1,  -xi + 1 ]]))

        # Reference 2, Example 5.5 (page 353)
        B_m = np.array([[dH[0, 0],    0,     dH[0, 1],    0,     dH[0, 2],    0,     dH[0, 3],    0    ],
                        [   0,     dH[1, 0],    0,     dH[1, 1],    0,     dH[1, 2],    0,     dH[1, 3]],
                        [dH[1, 0], dH[0, 0], dH[1, 1], dH[0, 1], dH[1, 2], dH[0, 2], dH[1, 3], dH[0, 3]]])

        return B_m

    def Hb(self) -> NDArray[float64]:
        '''
        Returns the stress-strain matrix for plate bending.
        '''

        # Referemce 1, Table 4.3, page 194
        nu = self.nu
        E = self.E
        h = self.t

        Hb = E*h**3/(12*(1 - nu**2))*np.array([[1,  nu,      0    ],
                                               [nu, 1,       0    ],
                                               [0,  0,  (1 - nu)/2]])
        
        return Hb

    def Hs(self) -> NDArray[float64]:
        '''
        Returns the stress-strain matrix for shear.
        '''
        # Reference 2, Equations (5.97), page 422
        k = 5/6
        E = self.E
        h = self.t
        nu = self.nu

        Hs = E*h*k/(2*(1 + nu))*np.array([[1, 0],
                                          [0, 1]])

        return Hs

    def Cm(self) -> NDArray[float64]:
        """
        Returns the stress-strain matrix for an isotropic or orthotropic plane stress element
        """
        
        # Apply the stiffness modification factors for each direction to obtain orthotropic
        # behavior. Stiffness modification factors of 1.0 in each direction (the default) will
        # model isotropic behavior. Orthotropic behavior is limited to the element's local
        # coordinate system.
        Ex = self.E*self.kx_mod
        Ey = self.E*self.ky_mod
        nu_xy = self.nu
        nu_yx = self.nu

        # The shear modulus will be unafected by orthotropic behavior
        # Logan, Appendix C.3, page 750
        G = self.E/(2*(1 + self.nu))

        # Gallagher, Equation 9.3, page 251
        Cm = 1/(1 - nu_xy*nu_yx)*np.array([[   Ex,    nu_yx*Ex,           0         ],
                                           [nu_xy*Ey,    Ey,              0         ],
                                           [    0,        0,     (1 - nu_xy*nu_yx)*G]])

        return Cm

    def k_b(self) -> NDArray[float64]:
        '''
        Returns the local stiffness matrix for bending and shear stresses
        '''

        Hb = self.Hb()
        Hs = self.Hs()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the determinant of the Jacobian matrix for each gauss point. Doing this now will save us from doing it twice below.
        det_J1 = det(self.J(-gp, -gp))
        det_J2 = det(self.J( gp, -gp))
        det_J3 = det(self.J( gp,  gp))
        det_J4 = det(self.J(-gp,  gp))

        # Get the bending [B_b] matrices for each gauss point
        B1 = self.B_b(-gp, -gp)
        B2 = self.B_b( gp, -gp)
        B3 = self.B_b( gp,  gp)
        B4 = self.B_b(-gp,  gp)

        # Create the stiffness matrix with bending stiffness terms
        # See 2, Equation 5.94
        k = ((B1.T @ Hb @ B1)*det_J1 +
             (B2.T @ Hb @ B2)*det_J2 +
             (B3.T @ Hb @ B3)*det_J3 +
             (B4.T @ Hb @ B4)*det_J4)

        # Get the shear [B_s] matrices for each gauss point
        B1 = self.B_s(-gp, -gp)
        B2 = self.B_s( gp, -gp)
        B3 = self.B_s( gp,  gp)
        B4 = self.B_s(-gp,  gp)

        # Add shear stiffness terms to the stiffness matrix
        k += ((B1.T @ Hs @ B1)*det_J1 +
              (B2.T @ Hs @ B2)*det_J2 +
              (B3.T @ Hs @ B3)*det_J3 +
              (B4.T @ Hs @ B4)*det_J4)
        
        # Following Bathe's recommendation for the drilling degree of freedom
        # from Example 4.19 in "Finite Element Procedures, 2nd Ed.", calculate
        # the drilling stiffness as 1/1000 of the smallest diagonal term in
        # the element's stiffness matrix. This is not theoretically correct,
        # but it allows the model to solve without singularities, and should
        # have a minimal effect on the final solution. Bathe recommends 1/1000
        # as a value that is weak enough but not so small that it affect the
        # results. Bathe recommends looking at all the diagonals in the
        # combined bending plus membrane stiffness matrix. Some of those terms
        # relate to translational stiffness. It seems more rational to only
        # look at the terms relating to rotational stiffness. That will be
        # Pynite's approach.
        k_rz = min(abs(k[1, 1]), abs(k[2, 2]), abs(k[4, 4]), abs(k[5, 5]),
                   abs(k[7, 7]), abs(k[8, 8]), abs(k[10, 10]), abs(k[11, 11])
                   )/1000
        
        # Initialize the expanded stiffness matrix to all zeros
        k_exp = np.zeros((24, 24))

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

        # Invert the local +y bending sign convention to match Pynite's
        k_exp[[4, 10, 16, 22], :] *= -1
        k_exp[:, [4, 10, 16, 22]] *= -1

        # The way the DKMQ element was derived, the positions relating to x
        # and y in the element's stiffness matrix are swapped from Pynite's.
        # Swap them to match Pynite.
        k_exp[[3, 4, 9, 10, 15, 16, 21, 22], :] = k_exp[[4, 3, 10, 9, 16, 15, 22, 21], :]
        k_exp[:, [3, 4, 9, 10, 15, 16, 21, 22]] = k_exp[:, [4, 3, 10, 9, 16, 15, 22, 21]]

        return k_exp

    def k_m(self) -> NDArray[float64]:
        '''
        Returns the local stiffness matrix for membrane (in-plane) stresses.

        Plane stress is assumed
        '''

        t = self.t
        Cm = self.Cm()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the membrane B matrices for each gauss point
        # Doing this now will save us from doing it twice below
        B1 = self.B_m(-gp, -gp)
        B2 = self.B_m( gp, -gp)
        B3 = self.B_m( gp,  gp)
        B4 = self.B_m(-gp,  gp)

        # Calculate the determinant of the Jacobian matrix at each gauss point
        det_J1 = det(self.J(-gp, -gp))
        det_J2 = det(self.J(gp, -gp))
        det_J3 = det(self.J(gp, gp))
        det_J4 = det(self.J(-gp, gp))

        if det_J1 <= 0 or det_J2 <= 0 or det_J3 <= 0 or det_J4 <= 0:
            warnings.warn(f'The Jacobian matrix for quad element {self.name} is less than or equal to zero, indicating the element is invalid or badly distorted.')

        # See reference 2 at the bottom of page 353, and reference 2 page 466
        k = t*((B1.T @ Cm @ B1)*det_J1 +
               (B2.T @ Cm @ B2)*det_J2 +
               (B3.T @ Cm @ B3)*det_J3 +
               (B4.T @ Cm @ B4)*det_J4)
        
        k_exp = np.zeros((24, 24))

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

    def k(self) -> NDArray[float64]:
        '''
        Returns the quad element's local stiffness matrix.
        '''

        # Recalculate the local coordinate system
        self._local_coords()

        # Sum the bending and membrane stiffness matrices
        return np.add(self.k_b(), self.k_m())

    def f(self, combo_name: str='Combo 1') -> NDArray[float64]:
        """
        Returns the quad element's local end force vector
        """

        # Calculate and return the plate's local end force vector
        return np.add(self.k() @ self.d(combo_name), self.fer(combo_name))

    def fer(self, combo_name: str='Combo 1') -> NDArray[float64]:
        """
        Returns the quadrilateral's local fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the consistent load vector for.
        """

        # Update the local coordinate system
        self._local_coords()

        Hw = lambda xi, eta : 1/4*np.array([[(1 - xi)*(1 - eta), 0, 0, (1 + xi)*(1 - eta), 0, 0, (1 + xi)*(1 + eta), 0, 0, (1 - xi)*(1 + eta), 0, 0]])

        # Initialize the fixed end reaction vector
        fer = np.zeros((12, 1))

        # Get the requested load combination
        combo = self.model.load_combos[combo_name]

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
             + Hw( gp, -gp).T*p*det(self.J( gp, -gp))
             + Hw( gp,  gp).T*p*det(self.J( gp,  gp))
             + Hw(-gp,  gp).T*p*det(self.J(-gp,  gp)))

        # Initialize the expanded vector to all zeros
        fer_exp = np.zeros((24, 1))

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

    def d(self, combo_name='Combo 1') -> NDArray[float64]:
       """
       Returns the quad element's local displacement vector
       """

       # Calculate and return the local displacement vector
       return self.T() @ self.D(combo_name)

    def F(self, combo_name: str='Combo 1') -> NDArray[float64]:
        """
        Returns the quad element's global force vector

        Parameters
        ----------
        combo_name : string
            The load combination to get results for.
        """

        # Calculate and return the global force vector
        return inv(self.T()) @ self.f(combo_name)

    def D(self, combo_name:str='Combo 1') -> NDArray[float64]:
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
        D = np.zeros((24, 1))
        
        # Read in the global displacements from the nodes
        D[0, 0] = self.i_node.DX[combo_name]
        D[1, 0] = self.i_node.DY[combo_name]
        D[2, 0] = self.i_node.DZ[combo_name]
        D[3, 0] = self.i_node.RX[combo_name]
        D[4, 0] = self.i_node.RY[combo_name]
        D[5, 0] = self.i_node.RZ[combo_name]

        D[6, 0] = self.j_node.DX[combo_name]
        D[7, 0] = self.j_node.DY[combo_name]
        D[8, 0] = self.j_node.DZ[combo_name]
        D[9, 0] = self.j_node.RX[combo_name]
        D[10, 0] = self.j_node.RY[combo_name]
        D[11, 0] = self.j_node.RZ[combo_name]

        D[12, 0] = self.m_node.DX[combo_name]
        D[13, 0] = self.m_node.DY[combo_name]
        D[14, 0] = self.m_node.DZ[combo_name]
        D[15, 0] = self.m_node.RX[combo_name]
        D[16, 0] = self.m_node.RY[combo_name]
        D[17, 0] = self.m_node.RZ[combo_name]

        D[18, 0] = self.n_node.DX[combo_name]
        D[19, 0] = self.n_node.DY[combo_name]
        D[20, 0] = self.n_node.DZ[combo_name]
        D[21, 0] = self.n_node.RX[combo_name]
        D[22, 0] = self.n_node.RY[combo_name]
        D[23, 0] = self.n_node.RZ[combo_name]

        # Return the global displacement vector
        return D

    def K(self) -> NDArray[float64]:
        '''
        Returns the quad element's global stiffness matrix
        '''

        # Get the transformation matrix
        T = self.T()

        # Calculate and return the stiffness matrix in global coordinates
        return inv(T) @ self.k() @ T

    # Global fixed end reaction vector
    def FER(self, combo_name:str='Combo 1') -> NDArray[float64]:
        '''
        Returns the global fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the fixed end
            reaction vector for (not the load combination itself).
        '''
        
        # Calculate and return the fixed end reaction vector
        return inv(self.T()) @ self.fer(combo_name)
  
    def T(self) -> NDArray[float64]:
        """
        Returns the coordinate transformation matrix for the quad element.
        """

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
        # mag = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
        # x = [x[0]/mag, x[1]/mag, x[2]/mag]
        x = x/norm(x)
        
        # The local y-axis will be in the plane of the plate. Find a vector in
        # the plate's local xy plane.
        xn = self.n_node.X
        yn = self.n_node.Y
        zn = self.n_node.Z
        xy = [xn - xi, yn - yi, zn - zi]

        # Find a vector perpendicular to the plate surface to get the
        # orientation of the local z-axis.
        z = np.cross(x, xy)
        
        # Divide the z-vector by its magnitude to produce a unit z-vector of
        # direction cosines.
        # mag = (z[0]**2 + z[1]**2 + z[2]**2)**0.5
        # z = [z[0]/mag, z[1]/mag, z[2]/mag]
        z = z/norm(z)

        # Calculate the local y-axis as a vector perpendicular to the local z
        # and x-axes.
        y = np.cross(z, x)
        
        # Divide the y-vector by its magnitude to produce a unit vector of
        # direction cosines.
        # mag = (y[0]**2 + y[1]**2 + y[2]**2)**0.5
        # y = [y[0]/mag, y[1]/mag, y[2]/mag]
        y = y/norm(y)

        # Create the direction cosines matrix.
        dir_cos = np.array([x,
                            y,
                            z])
        
        # Build the transformation matrix.
        T = np.zeros((24, 24))
        T[0:3, 0:3] = dir_cos
        T[3:6, 3:6] = dir_cos
        T[6:9, 6:9] = dir_cos
        T[9:12, 9:12] = dir_cos
        T[12:15, 12:15] = dir_cos
        T[15:18, 15:18] = dir_cos
        T[18:21, 18:21] = dir_cos
        T[21:24, 21:24] = dir_cos
        
        # Return the transformation matrix.
        return T
    
    # def T(self):
    #     """
    #     Returns the coordinate transformation matrix for the quad element.
    #     """

    #     xi = self.i_node.X
    #     xj = self.j_node.X
    #     yi = self.i_node.Y
    #     yj = self.j_node.Y
    #     zi = self.i_node.Z
    #     zj = self.j_node.Z

    #     # Calculate the direction cosines for the local x-axis.The local x-axis will run from
    #     # the i-node to the j-node
    #     x = [xj - xi, yj - yi, zj - zi]

    #     # Divide the vector by its magnitude to produce a unit x-vector of
    #     # direction cosines
    #     mag = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
    #     x = [x[0]/mag, x[1]/mag, x[2]/mag]
        
    #     # The local y-axis will be in the plane of the plate. Find a vector in
    #     # the plate's local xy plane.
    #     xn = self.n_node.X
    #     yn = self.n_node.Y
    #     zn = self.n_node.Z
    #     xy = [xn - xi, yn - yi, zn - zi]

    #     # Find a vector perpendicular to the plate surface to get the
    #     # orientation of the local z-axis.
    #     z = np.cross(x, xy)
        
    #     # Divide the z-vector by its magnitude to produce a unit z-vector of
    #     # direction cosines.
    #     mag = (z[0]**2 + z[1]**2 + z[2]**2)**0.5
    #     z = [z[0]/mag, z[1]/mag, z[2]/mag]

    #     # Calculate the local y-axis as a vector perpendicular to the local z
    #     # and x-axes.
    #     y = np.cross(z, x)
        
    #     # Divide the y-vector by its magnitude to produce a unit vector of
    #     # direction cosines.
    #     mag = (y[0]**2 + y[1]**2 + y[2]**2)**0.5
    #     y = [y[0]/mag, y[1]/mag, y[2]/mag]

    #     # Create the direction cosines matrix.
    #     dir_cos = np.array([x,
    #                         y,
    #                         z])
        
    #     # Build the transformation matrix.
    #     T = np.zeros((24, 24))
    #     T[0:3, 0:3] = dir_cos
    #     T[3:6, 3:6] = dir_cos
    #     T[6:9, 6:9] = dir_cos
    #     T[9:12, 9:12] = dir_cos
    #     T[12:15, 12:15] = dir_cos
    #     T[15:18, 15:18] = dir_cos
    #     T[18:21, 18:21] = dir_cos
    #     T[21:24, 21:24] = dir_cos
        
    #     # Return the transformation matrix.
    #     return T
    
    def shear(self, xi:float=0.0, eta:float=0.0, local:bool=True, combo_name:str='Combo 1') -> NDArray[float64]:
        """
        Returns the interal shears at any point in the quad element.

        Internal shears are reported as a 2D array [[Qx], [Qy]] at the
        specified location in the (xi, eta) natural coordinate system.

        Parameters
        ----------
        xi : float
            The xi-coordinate. Default is 0.
        eta : float
            The eta-coordinate. Default is 0.
        
        Returns
        -------
        Internal shear force per unit length of the quad element: [[Qx], [Qy]]
        """

        # Get the plate's local displacement vector
        d = self.d(combo_name)
                
        # Correct the sign convention for x-axis rotation - note that +x bending and +x rotation are opposite in the DKMQ derivation. Hence when correcting d we correct the x terms, but when correcting k we correct the y terms
        d[[3, 9, 15, 21], :] *= -1

        # Slice out terms not related to plate bending, and swap the local x and y to match the DKMQ derivation
        d = d[[2, 4, 3, 8, 10, 9, 14, 16, 15, 20, 22, 21], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Define extrapolated r and s points
        xi_ex = xi/gp
        eta_ex = eta/gp

        # Define the interpolation functions
        H = 1/4*np.array([(1 - xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 + eta_ex), (1 - xi_ex)*(1 + eta_ex)])

        # Get the stress-strain matrix
        Hs = self.Hs()

        # Calculate the internal shears [Qx, Qy] at each gauss point
        q1 = np.matmul(Hs, np.matmul(self.B_s(-gp, -gp), d))
        q2 = np.matmul(Hs, np.matmul(self.B_s(gp, -gp), d))
        q3 = np.matmul(Hs, np.matmul(self.B_s(gp, gp), d))
        q4 = np.matmul(Hs, np.matmul(self.B_s(-gp, gp), d))

        # Extrapolate to get the value at the requested location
        Qx = H[0]*q1[0] + H[1]*q2[0] + H[2]*q3[0] + H[3]*q4[0]
        Qy = H[0]*q1[1] + H[1]*q2[1] + H[2]*q3[1] + H[3]*q4[1]

        if local:
            
            return np.array([Qx,
                             Qy])
        
        else:
            
            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            # Transform the results to a global vector
            Qx = float(Qx)
            Qy = float(Qy)
            Q_global = np.matmul(dir_cos.T, np.array([[Qx, 0,  0],
                                                      [0,  Qy, 0],
                                                      [0,  0,  0]]))

            # Extract results acting along each global axis
            Qx_global = Q_global[0, 0]
            Qy_global = Q_global[1, 1]
            Qz_global = Q_global[2, 1]

            # Return the results as a 2D matrix
            return np.array([[Qx_global],
                             [Qy_global],
                             [Qz_global]])

    def moment(self, xi:float=0.0, eta:float=0.0, local:bool=True, combo_name:str='Combo 1') -> NDArray[float64]:
        """
        Returns the interal moments at any point in the quad element.

        Internal moments are reported as a 2D array [[Mx], [My], [Mxy]] at the
        specified location in the (xi, eta) natural coordinate system.

        Parameters
        ----------
        xi : float
            The xi-coordinate. Default is 0.
        eta : float
            The eta-coordinate. Default is 0.

        Returns
        -------
        Internal moment per unit length of the quad element: [[Mx], [My], [Mxy]]
        """

        # Get the plate's local displacement vector
        d = self.d(combo_name)

        # Correct the sign convention for x-axis rotation - note that +x bending and +x rotation are opposite in the DKMQ derivation. Hence when correcting d we correct the x terms, but when correcting k we correct the y terms
        d[[3, 9, 15, 21], :] *= -1

        # Slice out terms not related to plate bending, and swap the local x and y to match the DKMQ derivation
        d = d[[2, 4, 3, 8, 10, 9, 14, 16, 15, 20, 22, 21], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # # Define extrapolated `xi` and `eta` points
        xi_ex = xi/gp
        eta_ex = eta/gp

        # Define the interpolation functions
        H = 1/4*np.array([(1 - xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 + eta_ex), (1 - xi_ex)*(1 + eta_ex)])

        # Get the stress-strain matrix
        Hb = self.Hb()

        # Calculate the internal moments [-My, Mx, Mxy] at each gauss point
        m1 = np.matmul(Hb, np.matmul(self.B_b(-gp, -gp), d))
        m2 = np.matmul(Hb, np.matmul(self.B_b( gp, -gp), d))
        m3 = np.matmul(Hb, np.matmul(self.B_b( gp,  gp), d))
        m4 = np.matmul(Hb, np.matmul(self.B_b(-gp,  gp), d))

        # Extrapolate to get the value at the requested location
        Mx = H[0]*m1[0] + H[1]*m2[0] + H[2]*m3[0] + H[3]*m4[0]
        My = H[0]*m1[1] + H[1]*m2[1] + H[2]*m3[1] + H[3]*m4[1]
        Mxy = H[0]*m1[2] + H[1]*m2[2] + H[2]*m3[2] + H[3]*m4[2]

        if local:

            return np.array([Mx,
                             My,
                             Mxy])

        else:

            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            # Convert the local results to global results
            Mx = float(Mx)
            My = float(-My)
            Mxy = float(Mxy)
            M_local = np.array([
                [Mx, Mxy, 0.0],
                [Mxy, My, 0.0],
                [0.0, 0.0, 0.0]
            ])

            M_global_tensor = dir_cos @ M_local @ dir_cos.T

            # Extract the results for each direction
            Mx_global = M_global_tensor[0, 0]
            My_global = M_global_tensor[1, 1]
            Mxy_global = M_global_tensor[0, 1]

            return np.array([[Mx_global],
                             [My_global],
                             [Mxy_global]])


    def membrane(self, xi:float=0, eta: float=0, local:bool=True, combo_name:str='Combo 1') -> NDArray[float64]:

        # Get the plate's local displacement vector. Slice out terms not related to membrane stresses.
        d = self.d(combo_name)[[0, 1, 6, 7, 12, 13, 18, 19], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Define extrapolated r and s points
        xi_ex = xi/gp
        eta_ex = eta/gp

        # Define the interpolation functions
        H = 1/4*np.array([(1 - xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 + eta_ex), (1 - xi_ex)*(1 + eta_ex)])

        # Get the stress-strain matrix
        Cm = self.Cm()

        # Calculate the internal stresses [Sx, Sy, Txy] at each gauss point
        s1 = np.matmul(Cm, np.matmul(self.B_m(-gp, -gp), d))
        s2 = np.matmul(Cm, np.matmul(self.B_m(gp, -gp), d))
        s3 = np.matmul(Cm, np.matmul(self.B_m(gp, gp), d))
        s4 = np.matmul(Cm, np.matmul(self.B_m(-gp, gp), d))

        # Extrapolate to get the value at the requested location
        Sx = H[0]*s1[0] + H[1]*s2[0] + H[2]*s3[0] + H[3]*s4[0]
        Sy = H[0]*s1[1] + H[1]*s2[1] + H[2]*s3[1] + H[3]*s4[1]
        Txy = H[0]*s1[2] + H[1]*s2[2] + H[2]*s3[2] + H[3]*s4[2]

        if local:

            return np.array([Sx,
                             Sy,
                             Txy])

        else:

            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            # Convert the local results to global results
            Sx = float(Sx)
            Sy = float(Sy)
            Txy = float(Txy)
            S_local = np.array([
                [Sx, Txy, 0.0],
                [Txy, Sy, 0.0],
                [0.0, 0.0, 0.0]
            ])

            S_global_tensor = dir_cos @ S_local @ dir_cos.T

            # Extract the results for each direction
            Sx_global = S_global_tensor[0, 0]
            Sy_global = S_global_tensor[1, 1]
            Sxy_global = S_global_tensor[0, 1]

            return np.array([[Sx_global],
                             [Sy_global],
                             [Sxy_global]])
