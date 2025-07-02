# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:58:03 2017

@author: D. Craig Brinck, SE
"""
# %%
from __future__ import annotations # Allows more recent type hints features
from typing import TYPE_CHECKING

from numpy import zeros

if TYPE_CHECKING:
    from typing import Literal
    from numpy import float64
    from numpy.typing import NDArray


# %%
def FER_PtLoad(P: float, x: float, L: float, Direction: Literal["Fy", "Fz"]) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a point load

    Parameters:
    -----------
    P : float
        The magnitude of the point load
    x : float
        The location of the point load relative to the start of the member
    L : float
        The length of the member
    Direction : Literal["Fy", "Fz"]
        The direction of the point load. Must be one of the following:
            "Fy" = Force on the member's local y-axis
            "Fz" = Force on the member's local z-axis
    """
    # Define variables
    b = L - x

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    if Direction == "Fy":
        FER[1, 0] = -P*b**2*(L+2*x)/L**3
        FER[5, 0] = -P*x*b**2/L**2
        FER[7, 0] = -P*x**2*(L+2*b)/L**3
        FER[11, 0] = P*x**2*b/L**2
    elif Direction == "Fz":
        FER[2, 0] = -P*b**2*(L+2*x)/L**3
        FER[4, 0] = P*x*b**2/L**2
        FER[8, 0] = -P*x**2*(L+2*b)/L**3
        FER[10, 0] = -P*x**2*b/L**2

    return FER


def FER_Moment(M: float, x: float, L: float, Direction: Literal["My", "Mz"]) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a concentrated moment

    Parameters
    ----------
    M : float
        The magnitude of the moment
    x : float
        The location of the moment relative to the start of the member
    L : float
        The length of the member
    Direction : Literal["My", "Mz"]
        The direction of the moment. Must be one of the following:
            "My" = Moment applied about the local y-axis
            "Mz" = Moment applied about the local z-axis
    """

    # Define variables
    b = L - x

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    if Direction == "Mz":
        FER[1, 0] = 6*M*x*b/L**3
        FER[5, 0] = M*b*(2*x-b)/L**2
        FER[7, 0] = -6*M*x*b/L**3
        FER[11, 0] = M*x*(2*b-x)/L**2
    elif Direction == "My":
        FER[2, 0] = -6*M*x*b/L**3
        FER[4, 0] = M*b*(2*x-b)/L**2
        FER[8, 0] = 6*M*x*b/L**3
        FER[10, 0] = M*x*(2*b-x)/L**2
    return FER


# Returns the fixed end reaction vector for a linear distributed load
def FER_LinLoad(w1: float, w2: float, x1: float, x2: float, L: float, Direction: Literal["Fy", "Fz"]) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a linear distributed load

    Parameters
    ----------
    w1 : float
        The load magnitude at the start location
    w2 : float
        The load magnitude at the end location
    x1 : float
        The start location of the distributed load
    x2 : float
        The end location of the distributed load
    L : float
        The length of the member
    Direction : Literal["Fy", "Fz"]
        The direction of the distributed load. Must be one of the following:
            "Fy" = Force on the member's local y-axis
            "Fz" = Force on the member's local z-axis
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    if Direction == 'Fy':
        FER[1, 0] = (x1 - x2)*(10*L**3*w1 + 10*L**3*w2 - 15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3)
        FER[5, 0] = (x1 - x2)*(20*L**2*w1*x1 + 10*L**2*w1*x2 + 10*L**2*w2*x1 + 20*L**2*w2*x2 - 30*L*w1*x1**2 - 20*L*w1*x1*x2 - 10*L*w1*x2**2 - 10*L*w2*x1**2 - 20*L*w2*x1*x2 - 30*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2)
        FER[7, 0] = -(x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3)
        FER[11, 0] = (x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2)
    elif Direction == 'Fz':
        FER[2, 0] = (x1 - x2)*(10*L**3*w1 + 10*L**3*w2 - 15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3)
        FER[4, 0] = -(x1 - x2)*(20*L**2*w1*x1 + 10*L**2*w1*x2 + 10*L**2*w2*x1 + 20*L**2*w2*x2 - 30*L*w1*x1**2 - 20*L*w1*x1*x2 - 10*L*w1*x2**2 - 10*L*w2*x1**2 - 20*L*w2*x1*x2 - 30*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2)
        FER[8, 0] = -(x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3)
        FER[10, 0] = -(x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2)

    return FER


# Returns the fixed end reaction vector for an axial point load
def FER_AxialPtLoad(P: float, x: float, L: float) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for an axial point load

    Parameters
    ----------
    P : float
        The magnitude of the axial point load
    x : float
        The location of the axial point load relative to the start of the member
    L : float
        The length of the member
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    FER[0, 0] = -P*(L-x)/L
    FER[6, 0] = -P*x/L

    return FER


# Returns the fixed end reaction vector for a distributed axial load
def FER_AxialLinLoad(p1: float, p2: float, x1: float, x2: float, L: float) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a distributed axial load

    Parameters
    ----------
    p1 : float
        The axial load magnitude at the start location
    p2 : float
        The axial load magnitude at the end location
    x1 : float
        The start location of the distributed axial load
    x2 : float
        The end location of the distributed axial load
    L : float
        The length of the member
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    FER[0, 0] = 1/(6*L)*(x1-x2)*(3*L*p1+3*L*p2-2*p1*x1-p1*x2-p2*x1-2*p2*x2)
    FER[6, 0] = 1/(6*L)*(x1-x2)*(2*p1*x1+p1*x2+p2*x1+2*p2*x2)

    return FER


def FER_Torque(T: float, x: float, L: float) -> NDArray[float64]:
    """
    Returns the fixed end reaction vector for a concentrated torque

    Parameters
    ----------
    T : float
        The magnitude of the torque
    x : float
        The location of the torque relative to the start of the member
    L : float
        The length of the member
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))

    # Populate the fixed end reaction vector
    FER[3, 0] = -T*(L - x)/L
    FER[9, 0] = -T*x/L

    return FER
