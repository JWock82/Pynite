# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:04:56 2017

@author: D. Craig Brinck, SE
"""
# %%      
class Node3D():
    """
    A class representing a node in a 3D finite element model.
    """
    
    def __init__(self, Name, X, Y, Z):
        
        self.Name = Name    # A unique name for the node assigned by the user
        self.ID = None      # A unique index number for the node assigned by the program
        
        self.X = X          # Global X coordinate
        self.Y = Y          # Global Y coordinate
        self.Z = Z          # Global Z coordinate
        
        self.NodeLoads = [] # A list of loads applied to the node (Direction, P, case) or (Direction, M, case)
        
        # Initialize the dictionaries of calculated node displacements
        self.DX = {}
        self.DY = {}
        self.DZ = {}
        self.RX = {}
        self.RY = {}
        self.RZ = {}
        
        # Initialize the dictionaries of calculated node reactions
        self.RxnFX = {}
        self.RxnFY = {}
        self.RxnFZ = {}
        self.RxnMX = {}
        self.RxnMY = {}
        self.RxnMZ = {}

        # Initialize all support conditions to `False`
        self.support_DX = False
        self.support_DY = False
        self.support_DZ = False
        self.support_RX = False
        self.support_RY = False
        self.support_RZ = False

        # Inititialize all support springs
        self.spring_DX = [None, None, None]  # [stiffness, direction, active]
        self.spring_DY = [None, None, None]
        self.spring_DZ = [None, None, None]
        self.spring_RX = [None, None, None]
        self.spring_RY = [None, None, None]
        self.spring_RZ = [None, None, None]

        # Initialize all enforced displacements to `None`
        self.EnforcedDX = None
        self.EnforcedDY = None
        self.EnforcedDZ = None
        self.EnforcedRX = None
        self.EnforcedRY = None
        self.EnforcedRZ = None

        # Initialize the color contour value for the node. This will be used for contour smoothing.
        self.contour = []
