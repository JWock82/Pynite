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
from numpy import allclose

class Test_Plates(unittest.TestCase):
    ''' Tests of analyzing plate elements. '''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_plate_displacement(self):
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Example 12.1
        # Units for this model are pounds and inches
        plModel= FEModel3D()

        plModel.add_node('N1', 0, 0, 0)
        plModel.add_node('N2', 10, 0, 0)
        plModel.add_node('N3', 20, 0, 0)
        plModel.add_node('N4', 0, 10, 0)
        plModel.add_node('N5', 10, 10, 0)
        plModel.add_node('N6', 20, 10, 0)
        plModel.add_node('N7', 0, 20, 0)
        plModel.add_node('N8', 10, 20, 0)
        plModel.add_node('N9', 20, 20, 0)

        plModel.add_plate('P1', 'N1', 'N2', 'N5', 'N4', 0.1, 30000000, 0.3)
        plModel.add_plate('P2', 'N2', 'N3', 'N6', 'N5', 0.1, 30000000, 0.3)
        plModel.add_plate('P3', 'N4', 'N5', 'N8', 'N7', 0.1, 30000000, 0.3)
        plModel.add_plate('P4', 'N5', 'N6', 'N9', 'N8', 0.1, 30000000, 0.3)

        plModel.add_node_load('N5', 'FZ', -100)

        plModel.def_support('N1', True, True, True, True, True, True)
        plModel.def_support('N2', True, True, True, True, True, True)
        plModel.def_support('N3', True, True, True, True, True, True)
        plModel.def_support('N4', True, True, True, True, True, True)
        plModel.def_support('N6', True, True, True, True, True, True)
        plModel.def_support('N7', True, True, True, True, True, True)
        plModel.def_support('N8', True, True, True, True, True, True)
        plModel.def_support('N9', True, True, True, True, True, True)

        plModel.def_support('N5', True, True, False, False, False, True)

        # Check to see if the stiffness matrix for each plate is symmetric
        # print(allclose(plModel.Plates[0].K(), plModel.Plates[0].K().T))
        # print(allclose(plModel.Plates[1].K(), plModel.Plates[1].K().T))
        # print(allclose(plModel.Plates[2].K(), plModel.Plates[2].K().T))
        # print(allclose(plModel.Plates[3].K(), plModel.Plates[3].K().T))

        # Check to see if the global stiffness matrix is symmetric
        # print(allclose(plModel.K(Renumber=True), plModel.K(Renumber=False).T))

        plModel.analyze(check_statics=True)
        # Test: displacement of N5 in Z direction
        calculated_displacement = plModel.Nodes['N5'].DZ['Combo 1']
        expected_displacement = -0.0861742424242424
        self.assertAlmostEqual(calculated_displacement/expected_displacement, 1.0, 2)