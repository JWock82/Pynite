# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:04:56 2017

@author: D. Craig Brinck, SE
"""
from __future__ import annotations # Allows more recent type hints features
# %%      
from typing import List, Tuple, Dict, Optional

class Node3D():
    """
    A class representing a node in a 3D finite element model.
    """
    
    def __init__(self, name: str, X: float, Y: float, Z: float):
        
        self.name = name                 # A unique name for the node assigned by the user
        self.ID: Optional[int] = None    # A unique index number for the node assigned by the program
        
        self.X = X          # Global X coordinate
        self.Y = Y          # Global Y coordinate
        self.Z = Z          # Global Z coordinate
        
        self.NodeLoads: List[Tuple[str, float, str]] = [] # A list of loads applied to the node (Direction, P, case) or (Direction, M, case)
        
        # Initialize the dictionaries of calculated node displacements
        self.DX: Dict[str, float] = {}
        self.DY: Dict[str, float] = {}
        self.DZ: Dict[str, float] = {}
        self.RX: Dict[str, float] = {}
        self.RY: Dict[str, float] = {}
        self.RZ: Dict[str, float] = {}
        
        # Initialize the dictionaries of calculated node reactions
        self.RxnFX: Dict[str, float] = {}
        self.RxnFY: Dict[str, float] = {}
        self.RxnFZ: Dict[str, float] = {}
        self.RxnMX: Dict[str, float] = {}
        self.RxnMY: Dict[str, float] = {}
        self.RxnMZ: Dict[str, float] = {}

        # Initialize all support conditions to `False`
        self.support_DX: bool = False
        self.support_DY: bool = False
        self.support_DZ: bool = False
        self.support_RX: bool = False
        self.support_RY: bool = False
        self.support_RZ: bool = False

        # Inititialize all support springs
        self.spring_DX: List[float | str | bool | None] = [None, None, None]  # [stiffness, direction, active]
        self.spring_DY: List[float | str | bool | None] = [None, None, None]
        self.spring_DZ: List[float | str | bool | None] = [None, None, None]
        self.spring_RX: List[float | str | bool | None] = [None, None, None]
        self.spring_RY: List[float | str | bool | None] = [None, None, None]
        self.spring_RZ: List[float | str | bool | None] = [None, None, None]

        # Initialize all enforced displacements to `None`
        self.EnforcedDX: float | None = None
        self.EnforcedDY: float | None = None
        self.EnforcedDZ: float | None = None
        self.EnforcedRX: float | None = None
        self.EnforcedRY: float | None = None
        self.EnforcedRZ: float | None = None

        # Initialize the color contour value for the node. This will be used for contour smoothing.
        self.contour: List[float] = []

    def distance(self, other: 'Node3D') -> float:
        """
        Returns the distance to another node.

        Parameters
        ----------
        other : Node3D
            A node object to compare coordinates with.
        """
        return ((self.X - other.X)**2 + (self.Y - other.Y)**2 + (self.Z - other.Z)**2)**0.5
