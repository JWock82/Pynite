# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import math
import sys
from io import StringIO

class Test_2D_Frame(unittest.TestCase):
    ''' Tests of analyzing 2D frames. '''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_XY_gravity_load(self):
        # Create a new model
        frame = FEModel3D()

        # Define the nodes
        frame.AddNode('N1', 0, 0, 0)
        frame.AddNode('N2', 0, 30*12, 0)
        frame.AddNode('N3', 15*12, 40*12, 0)
        frame.AddNode('N4', 35*12, 40*12, 0)
        frame.AddNode('N5', 50*12, 30*12, 0)
        frame.AddNode('N6', 50*12, 0, 0)

        # Define the supports
        frame.DefineSupport('N1', True, True, True, True, True, True)
        frame.DefineSupport('N6', True, True, True, True, True, True)

        # Create members (all members will have the same properties in this example)
        J = 250
        Iy = 250
        Iz = 200
        E = 30000
        G = 250
        A = 12

        frame.AddMember('M1', 'N1', 'N2', E, G, Iy, Iz, J, A)
        frame.AddMember('M2', 'N2', 'N3', E, G, Iy, Iz, J, A)
        frame.AddMember('M3', 'N3', 'N4', E, G, Iy, Iz, J, A)
        frame.AddMember('M4', 'N4', 'N5', E, G, Iy, Iz, J, A)
        frame.AddMember('M5', 'N5', 'N6', E, G, Iy, Iz, J, A)

        # Add nodal loads
        frame.AddNodeLoad('N3', 'FY', -30)
        frame.AddNodeLoad('N4', 'FY', -30)

        # Analyze the model
        frame.Analyze()

        node1 = frame.GetNode('N1')
        node6 = frame.GetNode('N6')
        self.assertAlmostEqual(node1.RxnFX['Combo 1'], 11.6877, 4)
        self.assertAlmostEqual(node1.RxnFY['Combo 1'], 30, 4)
        self.assertAlmostEqual(node1.RxnMZ['Combo 1'], -1810.0745, 4)
        self.assertAlmostEqual(node6.RxnFX['Combo 1'], -11.6877, 4)
        self.assertAlmostEqual(node6.RxnFY['Combo 1'], 30, 4)
        self.assertAlmostEqual(node6.RxnMZ['Combo 1'], 1810.0745, 4)
