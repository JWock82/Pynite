# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:04:56 2017

@author: D. Craig Brinck, SE
"""
  # %%      
class Node3D():
    
    # Constructor
    def __init__(self, Name, ID, X, Y, Z):
        
        self.Name = Name # A unique name for the node assigned by the user
        self.ID = ID # A unique index number assigned by the program
        self.X = X # Global X coordinate
        self.Y = Y # Global Y coordinate
        self.Z = Z # Global Z coordinate
        self.SupportX = False # Flag indicating support in the X direction
        self.SupportY = False # Flag indicating support in the Y direction
        self.SupportZ = False # Flag indicating support in the Z direction
        self.SupportRX = False # Flag indicating rotational support about the X axis
        self.SupportRY = False # Flag indicating rotational support about the Y axis
        self.SupportRZ = False # Flag indicating rotational support about the Z axis
        
        # Initialize nodal displacements to zero
        self.DX = 0
        self.DY = 0
        self.DZ = 0
        self.RX = 0
        self.RY = 0
        self.RZ = 0
