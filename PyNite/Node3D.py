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
        
        self.NodeLoads = []     # A list of loads applied to the node (Direction, P, case) or (Direction, M, case)
        
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

        # Initialize all support conditions to 'None'
        self.SupportDX = None
        self.SupportDY = None
        self.SupportDZ = None
        self.SupportRX = None
        self.SupportRY = None
        self.SupportRZ = None
