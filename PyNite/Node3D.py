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
        """
        Initializes a new node.
        """
        
        self.Name = Name    # A unique name for the node assigned by the user
        self.ID = None      # A unique index number for the node assigned by the program
        
        self.X = X          # Global X coordinate
        self.Y = Y          # Global Y coordinate
        self.Z = Z          # Global Z coordinate
        
        self.NodeLoads = []     # A list of loads applied to the node (Direction, P) or (Direction, M)
        
        # Initialize nodal displacements to 'None'
        self.DX = None
        self.DY = None
        self.DZ = None
        self.RX = None
        self.RY = None
        self.RZ = None
        
        # Initialize Reactions to zero
        self.RxnFX = 0
        self.RxnFY = 0
        self.RxnFZ = 0
        self.RxnMX = 0
        self.RxnMY = 0
        self.RxnMZ = 0
