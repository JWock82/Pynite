# %%
# 3D frame element
class Frame3D():

    # Constructor
    def __init__(self, iNode, jNode, E, G, Iy, Iz, J, A):
        
        self.iNode = iNode
        self.jNode = jNode
        self.E = E
        self.G = G
        self.Iy = Iy
        self.Iz = Iz
        self.J = J
        self.A = A
        self.L = ((jNode.X-iNode.X)**2+(jNode.Y-iNode.Y)**2+(jNode.Z-iNode.Z)**2)**0.5

    # Local stiffness matrix
    def k(self):
        
        E = self.E
        G = self.G
        Iy = self.Iy
        Iz = self.Iz
        J = self.J
        A = self.A
        L = self.L
        
        return np.matrix([[A*E/L, 0, 0 , 0, 0, 0, -A*E/L, 0, 0, 0, 0, 0],
                          [0, 12*E*Iz/L**3, 0, 0, 0, 6*E*Iz/L**2, 0, -12*E*Iz/L**3, 0, 0, 0, 6*E*Iz/L**2],
                          [0, 0, 12*E*Iy/L**3, 0, -6*E*Iy/L**2, 0, 0, 0, -12*E*Iy/L**3, 0, -6*E*Iy/L**2, 0],
                          [0, 0, 0, G*J/L, 0, 0, 0, 0, 0, -G*J/L, 0, 0],
                          [0, 0, -6*E*Iy/L**2, 0, 4*E*Iy/L, 0, 0, 0, 6*E*Iy/L**2, 0, 2*E*Iy/L, 0],
                          [0, 6*E*Iz/L**2, 0, 0, 0, 4*E*Iz/L, 0, -6*E*Iz/L**2, 0, 0, 0, 2*E*Iz/L],
                          [-A*E/L, 0, 0, 0, 0, 0, A*E/L, 0, 0, 0, 0, 0],
                          [0, -12*E*Iz/L**3, 0, 0, 0, -6*E*Iz/L**2, 0, 12*E*Iz/L**3, 0, 0, 0, -6*E*Iz/L**2],
                          [0, 0, -12*E*Iy/L**3, 0, 6*E*Iy/L**2, 0, 0, 0, 12*E*Iy/L**3, 0, 6*E*Iy/L**2, 0],
                          [0, 0, 0, -G*J/L, 0, 0, 0, 0, 0, G*J/L, 0, 0],
                          [0, 0, -6*E*Iy/L**2, 0, 2*E*Iy/L, 0, 0, 0, 6*E*Iy/L**2, 0, 4*E*Iy/L, 0],
                          [0, 6*E*Iz/L**2, 0, 0, 0, 2*E*Iz/L, 0, -6*E*Iz/L**2, 0, 0, 0, 4*E*Iz/L]])
    
    # Transformation matrix
    def T(self):
    
        x1 = self.iNode.X
        x2 = self.jNode.X
        y1 = self.iNode.Y
        y2 = self.jNode.Y
        z1 = self.iNode.Z
        z2 = self.jNode.Z
        L = self.L
        
        # Calculate direction cosines for the transformation matrix
        l = (x2-x1)/L
        m = (y2-y1)/L
        n = (z2-z1)/L
        D = (l**2+m**2)**0.5
        
        if l == 0 and m == 0 and n > 0:
            dirCos = np.matrix([[0, 0, 1],
                                [0, 1, 0],
                                [-1, 0, 0]])
        elif l ==0 and m == 0 and n < 0:
            dirCos = np.matrix([[0, 0, -1],
                                [0, 1, 0],
                                [1, 0, 0]])
        else:
            dirCos = np.matrix([[l, m, n],
                                [-m/D, l/D, 0],
                                [-l*n/D, -m*n/D, D]])
        
        # Build the transformation matrix
        transMatrix = np.zeros((12, 12))
        transMatrix[0:3, 0:3] = dirCos
        transMatrix[3:6, 3:6] = dirCos
        transMatrix[6:9, 6:9] = dirCos
        transMatrix[9:12, 9:12] = dirCos
        
        return transMatrix
        
