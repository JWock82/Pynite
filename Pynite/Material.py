from __future__ import annotations # Allows more recent type hints features
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from Pynite.FEModel3D import FEModel3D

class Material():
    """
    A class representing a material assigned to a Member3D, Plate or Quad in a finite element model.

    This class stores all properties related to the physical material of the element
    """
    def __init__(self, model: FEModel3D, name: str, E: float, G: float, nu: float, rho: float, fy: float | None = None) -> None:
        """Initialize a material object.

        :param model: The finite element model this material belongs to
        :type model: FEModel3D
        :param name: A unique name for the material
        :type name: str
        :param E: Modulus of elasticity
        :type E: float
        :param G: Shear modulus
        :type G: float
        :param nu: Poisson's ratio
        :type nu: float
        :param rho: Density of the material
        :type rho: float
        :param fy: Yield strength of the material, defaults to None
        :type fy: float, optional
        """

        self.model: FEModel3D = model
        self.name: str = name
        self.E: float = E
        self.G: float = G
        self.nu: float = nu
        self.rho: float = rho
        self.fy: float | None = fy