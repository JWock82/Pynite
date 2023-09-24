
import numpy as np

class Section():

    def __init__(self, model, name, A, Iy, Iz, J, material):
        
        self.model = model
        self.name = name
        self.A = A
        self.Iy = Iy
        self.Iz = Iz
        self.J = J
        self.material = material
        self.fy = model.Materials[material].fy
    
    def Phi(self):
        pass

    def G(self):
        pass

class SteelSection(Section):

    def __init__(self, model, name, A, Iy, Iz, J, Zy, Zz, material):

        # Basic section properties
        super().__init__(model, name, A, Iy, Iz, J, material)

        # Additional section properties for steel
        self.ry = (Iy/A)**0.5
        self.rz = (Iz/A)**0.5
        self.Zy = Zy
        self.Zz = Zz
    
    def Phi(self, fx, my, mz):
        """A method used to determine whether the cross section is elastic or plastic. Values less than 1 indicate the section is elastic.

        :param fx: Axial force divided by axial strength.
        :type fx: float
        :param my: Weak axis moment divided by weak axis strength.
        :type my: float
        :param mz: Strong axis moment divided by strong axis strength.
        :type mz: float
        :return: The total stress ratio for the cross section.
        :rtype: float
        """
                
        # Plastic strengths for material nonlinearity
        Py = self.fy*self.A
        Mpy = self.fy*self.Zy
        Mpz = self.fy*self.Zz

        # Values for p, my, and mz based on actual loads
        p = fx/Py
        m_y = my/Mpy
        m_z = mz/Mpz

        # "Matrix Structural Analysis, 2nd Edition", Equation 10.18
        return p**2 + m_z**2 + m_y**4 + 3.5*p**2*m_z**2 + 3*p**6*m_y**2 + 4.5*m_z**4*m_y**2

    def G(self, fx, my, mz):
        """Returns the gradient to the yield surface at a given point
        """
        
        # Small increment for numerical differentiation
        epsilon = 1e-6

        # Calculate the central differences for each parameter
        dPhi_dfx = (self.Phi(fx + epsilon, my, mz) - self.Phi(fx - epsilon, my, mz)) / (2 * epsilon)
        dPhi_dmy = (self.Phi(fx, my + epsilon, mz) - self.Phi(fx, my - epsilon, mz)) / (2 * epsilon)
        dPhi_dmz = (self.Phi(fx, my, mz + epsilon) - self.Phi(fx, my, mz - epsilon)) / (2 * epsilon)

        # Return the gradient
        return np.array([dPhi_dfx, 0, 0, 0, dPhi_dmy, dPhi_dmz])

# test_section = SteelSection('W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 50)
# print(test_section.G(15, 30, 50))



