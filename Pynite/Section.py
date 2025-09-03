from __future__ import annotations # Allows more recent type hints features

import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from numpy import float64
    from numpy.typing import NDArray
    from Pynite.FEModel3D import FEModel3D

class Section():
    """
    A class representing a section assigned to a Member3D element in a finite element model.

    This class stores all properties related to the geometry of the member
    """
    def __init__(self, model: 'FEModel3D', name: str, A: float, Iy: float, Iz: float, J: float) -> None:
        """
        :param model: The finite element model to which this section belongs
        :type model: FEModel3D
        :param name: Name of the section
        :type name: str
        :param A: Cross-sectional area of the section
        :type A: float
        :param Iy: The second moment of area the section about the Y (minor) axis
        :type Iy: float
        :param Iz: The second moment of area the section about the Z (major) axis
        :type Iz: float
        :param J: The torsion constant of the section
        :type J: float
        """        
        self.model: 'FEModel3D' = model
        self.name: str = name
        self.A: float = A
        self.Iy: float = Iy
        self.Iz: float = Iz
        self.J: float = J
    
    def Phi(self, fx: float = 0, my: float = 0, mz: float = 0):
        """
        Method to be overridden by subclasses for determining whether the cross section is
        elastic or plastic.
        
        :param fx: Axial force
        :type fx: float
        :param my: y-axis (weak) moment
        :type my: float
        :param mz: z-axis (strong) moment
        :type mz: float
        :return: The stress ratio
        :rtype: float
        """
        raise NotImplementedError("Phi method must be implemented in subclasses.")

    def G(self, fx: float, my: float, mz: float) -> NDArray:
        """
        Returns the gradient to the yield surface at a given point using numerical differentiation. This is a default solution. For a better solution, overwrite this method with a more precise one in the material/shape specific child class that inherits from this class.

        :param fx: Axial force at the cross-section
        :type fx: float
        :param my: y-axis (weak) moment at the cross-section
        :type my: float
        :param mz: z-axis (strong) moment at the cross-section
        :type mz: float
        :return: The gradient to the yield surface at the cross-section
        :rtype: NDArray
        """

        # Small increment for numerical differentiation
        epsilon = 1e-6

        # Calculate the central differences for each parameter
        dPhi_dfx = (self.Phi(fx + epsilon, my, mz) - self.Phi(fx - epsilon, my, mz)) / (2 * epsilon)
        dPhi_dmy = (self.Phi(fx, my + epsilon, mz) - self.Phi(fx, my - epsilon, mz)) / (2 * epsilon)
        dPhi_dmz = (self.Phi(fx, my, mz + epsilon) - self.Phi(fx, my, mz - epsilon)) / (2 * epsilon)

        # Return the gradient
        return np.array([[dPhi_dfx],
                         [dPhi_dmy],
                         [dPhi_dmz]])

class SteelSection(Section):

    def __init__(self, model: 'FEModel3D', name: str, A: float, Iy: float, Iz: float, J: float, 
                 Zy: float, Zz: float, material_name: str) -> None:
        """
        Initialize a steel section

        :param model: The finite element model to which this section belongs
        :type model: FEModel3D
        :param name: Name of the section
        :type name: str
        :param A: Cross-sectional area of the section
        :type A: float
        :param Iy: The second moment of area the section about the Y (minor) axis
        :type Iy: float
        :param Iz: The second moment of area the section about the Z (major) axis
        :type Iz: float
        :param J: The torsion constant of the section
        :type J: float
        :param Zy: Plastic section modulus about the Y (minor) axis
        :type Zy: float
        :param Zz: Plastic section modulus about the Z (major) axis
        :type Zz: float
        :param material_name: Name of the material used for this section
        :type material_name: str
        """

        # Basic section properties
        super().__init__(model, name, A, Iy, Iz, J)

        # Additional section properties for steel
        self.ry: float = (Iy/A)**0.5
        self.rz: float = (Iz/A)**0.5
        self.Zy: float = Zy
        self.Zz: float = Zz

        self.material = model.materials[material_name]

    def Phi(self, fx: float = 0, my: float = 0, mz: float = 0) -> float:
        """
        A method used to determine whether the cross section is elastic or plastic.
        Values less than 1 indicate the section is elastic.

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
        Py = self.material.fy*self.A
        Mpy = self.material.fy*self.Zy
        Mpz = self.material.fy*self.Zz

        # Values for p, my, and mz based on actual loads
        p = fx/Py
        m_y = my/Mpy
        m_z = mz/Mpz

        # "Matrix Structural Analysis, 2nd Edition", Equation 10.18
        return p**2 + m_z**2 + m_y**4 + 3.5*p**2*m_z**2 + 3*p**6*m_y**2 + 4.5*m_z**4*m_y**2

    def G(self, fx: float, my: float, mz: float) -> NDArray[float64]:
        """Returns the gradient to the material's yield surface for the given load. Used to construct the plastic reduction matrix for nonlinear behavior.

        :param fx: Axial force at the cross-section.
        :type fx: float
        :param my: y-axis (weak) moment at the cross-section.
        :type my: float
        :param mz: z-axis (strong) moment at the cross-section.
        :type mz: float
        :return: The gradient to the material's yield surface at the cross-section.
        :rtype: NDArray
        """

        # Calculate `Phi` which is essentially a stress check indicating how close to yield we are
        Phi = self.Phi(fx, my, mz)

        # If Phi is less than 1.0 the member is still elastic and there is no gradient to the yield surface
        if Phi < 1.0:

            # G = zero vector
            return np.array([[0],
                             [0],
                             [0],
                             [0],
                             [0],
                             [0]])

        # Plastic action is occuring
        else:

            # Plastic strengths for material nonlinearity
            Py = self.material.fy*self.A
            Mpy = self.material.fy*self.Zy
            Mpz = self.material.fy*self.Zz

            # Partial derivatives of Phi
            dPhi_dfx = 18*fx**5*my**2/(Mpy**2*Py**6) + 2*fx/Py**2 + 7.0*fx*mz**2/(Mpz**2*Py**2)
            dPhi_dmy = 6*fx**6*my/(Mpy**2*Py**6) + 2*my/Mpy**2 + 9.0*my*mz**4/(Mpy**2*Mpz**4)
            dPhi_dmz = 7.0*fx**2*mz/(Mpz**2*Py**2) + 2*mz/Mpz**2 + 18.0*my**2*mz**3/(Mpy**2*Mpz**4)

            # Return the gradient
            return np.array([[dPhi_dfx],
                             [0],
                             [0],
                             [0],
                             [dPhi_dmy],
                             [dPhi_dmz]])

# test_section = SteelSection('W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 50)
# print(test_section.G(15, 30, 50))
