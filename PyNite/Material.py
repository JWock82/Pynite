class Material():

    def __init__(self, name, E, G, nu, rho, fy=None):

        self.name = name
        self.E = E
        self.G = G
        self.nu = nu
        self.rho = rho
        self.fy = fy