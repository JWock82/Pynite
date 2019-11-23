from numpy import matrix, matmul, transpose
from numpy.linalg import inv

# A 2D isotropic plate bending element
# Note: 2D plates can be used in a 3D model, but they must lie in a plane parallel to the global X-Y plane.
#       Any node Z-coordinates will be ignored by this element.
class Plate2D():

    def __init__(self, iNode, jNode, mNode, nNode, t, E, nu):

        self.iNode = iNode
        self.jNode = jNode
        self.mNode = mNode
        self.nNode = nNode

        self.t = t

        self.E = E
        self.nu = nu
    
    # Calculates and returns the [C] coefficient matrix
    def __C(self):

        # Find the x and y coordinates at each node
        xi = iNode.X
        yi = iNode.Y
        xj = jNode.X
        yj = jNode.Y
        xm = mNode.X
        ym = mNode.Y
        xn = nNode.X
        yn = nNode.Y

        # Calculate the [C] coefficient matrix
        C = matrix([[1, xi, yi, xi**2, xi*yi, yi**2, xi**3, xi**2*yi, xi*yi**2, yi**3, xi**3*yi, xi*yi**3],
                    [0, 0, 1, 0, xi, 2*yi, 0, xi**2, 2*xi*yi, 3*yi**2, xi**3, 3*xi*yi**2],
                    [0, -1, 0, -2*xi, -yi, 0, -3*xi**2, -2*xi*yi, -yi**2, 0, -3*xi**2*yi, -yi**3],

                    [1, xj, yj, xj**2, xj*yj, yj**2, xj**3, xj**2*yj, xj*yj**2, yj**3, xj**3*yj, xj*yj**3],
                    [0, 0, 1, 0, xj, 2*yj, 0, xj**2, 2*xj*yj, 3*yj**2, xj**3, 3*xj*yj**2],
                    [0, -1, 0, -2*xj, -yj, 0, -3*xj**2, -2*xj*yj, -yj**2, 0, -3*xj**2*yj, -yj**3],

                    [1, xm, ym, xm**2, xm*ym, ym**2, xm**3, xm**2*ym, xm*ym**2, ym**3, xm**3*ym, xm*ym**3],
                    [0, 0, 1, 0, xm, 2*ym, 0, xm**2, 2*xm*ym, 3*ym**2, xm**3, 3*xm*ym**2],
                    [0, -1, 0, -2*xm, -ym, 0, -3*xm**2, -2*xm*ym, -ym**2, 0, -3*xm**2*ym, -ym**3],
                    
                    [1, xn, yn, xn**2, xn*yn, yn**2, xn**3, xn**2*yn, xn*yn**2, yn**3, xn**3*yn, xn*yn**3],
                    [0, 0, 1, 0, xn, 2*yn, 0, xn**2, 2*xn*yn, 3*yn**2, xn**3, 3*xn*yn**2],
                    [0, -1, 0, -2*xn, -yn, 0, -3*xn**2, -2*xn*yn, -yn**2, 0, -3*xn**2*yn, -yn**3]])

        # Return the coefficient matrix
        return C

    # Calculates and returns the [Q] coefficient matrix
    def __Q(self):

        # Find the x and y coordinates at each node
        xi = iNode.X
        yi = iNode.Y
        xj = jNode.X
        yj = jNode.Y
        xm = mNode.X
        ym = mNode.Y
        xn = nNode.X
        yn = nNode.Y

        # Calculate the [Q] coefficient matrix
        Q = matrix([[0, 0, 0, -2, 0, 0, -6*xi, -2*yi, 0, 0, -6*xi*yi, 0],
                    [0, 0, 0, 0, 0, -2, 0, 0, -2*xi, -6*yi, 0, -6*xi*yi],
                    [0, 0, 0, 0, -2, 0, 0, -4*xi, -4*yi, 0, -6*xi**2, -6*yi**2],

                    [0, 0, 0, -2, 0, 0, -6*xj, -2*yj, 0, 0, -6*xj*yj, 0],
                    [0, 0, 0, 0, 0, -2, 0, 0, -2*xj, -6*yj, 0, -6*xj*yj],
                    [0, 0, 0, 0, -2, 0, 0, -4*xj, -4*yj, 0, -6*xj**2, -6*yj**2],

                    [0, 0, 0, -2, 0, 0, -6*xm, -2*ym, 0, 0, -6*xm*ym, 0],
                    [0, 0, 0, 0, 0, -2, 0, 0, -2*xm, -6*ym, 0, -6*xm*ym],
                    [0, 0, 0, 0, -2, 0, 0, -4*xm, -4*ym, 0, -6*xm**2, -6*ym**2],
                    
                    [0, 0, 0, -2, 0, 0, -6*xn, -2*yn, 0, 0, -6*xn*yn, 0],
                    [0, 0, 0, 0, 0, -2, 0, 0, -2*xn, -6*yn, 0, -6*xn*yn],
                    [0, 0, 0, 0, -2, 0, 0, -4*xn, -4*yn, 0, -6*xn**2, -6*yn**2]])
        
        # Return the [Q] coefficient matrix
        return Q
    
    # Calculates and returns the gradient matrix [B]
    def __B(self):

        # Calculate the gradient matrix [B]
        B = matmul(self.__Q(), inv(self.__C()))

        # Return the gradient matrix [B]
        return B
    
    # Calculates and returns the constitutive matrix for isotropic materials [D]
    def __D(self):

        # Calculate the coefficient for the constitutive matrix [D]
        C = self.E*self.t**3/(12*(1-self.nu**2))

        # Calculate the constitutive matrix [D]
        D = matrix([[C, C*self.nu, 0],
                    [C*self.nu, C, 0],
                    [0, 0, C*(1-self.nu)/2]])

        # Return the constitutive matrix [D]
        return D