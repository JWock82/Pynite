# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import sys
from io import StringIO

class Test_Spring_Elements(unittest.TestCase):
    ''' Tests of spring members.'''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
        
    def test_spring_elements(self): 
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Example 2.1
        # Units for this model are pounds and inches
        system = FEModel3D()
        system.add_node('1', 0, 0, 0)
        system.add_node('2', 30, 0, 0)
        system.add_node('3', 10, 0, 0)
        system.add_node('4', 20, 0, 0)
        # Add spring members
        system.add_spring('S1', '1', '3', 1000)
        system.add_spring('S2', '3', '4', 2000)
        system.add_spring('S3', '4', '2', 3000)
        # Define supports
        system.def_support('1', True, True, True, True, True, True)
        system.def_support('2', True, True, True, True, True, True)
        system.def_support('3', False, True, True, True, True, True)
        system.def_support('4', False, True, True, True, True, True)
        # Add node loads
        system.add_node_load('4', 'FX', 5000)
        system.analyze(True)
        # Check results
        # correct_values = [('3', 0.9090909090909092),
        #                   ('4', 1.3636363636363638),
        #                   ('1', -909.0909090909091),
        #                   ('2', -4090.9090909090914)]
        n3_DX = system.Nodes['3'].DX['Combo 1']
        self.assertAlmostEqual(n3_DX/ 0.9090909090909092, 1.0, 2)

        n4_DX = system.Nodes['4'].DX['Combo 1']
        self.assertAlmostEqual(n4_DX/1.3636363636363638, 1.0, 2)
        
        n1_rxn = system.Nodes['1'].RxnFX['Combo 1']
        self.assertAlmostEqual(n1_rxn/-909.0909090909091, 1.0, 2)
        
        n2_rxn = system.Nodes['2'].RxnFX['Combo 1']
        self.assertAlmostEqual(n2_rxn/-4090.9090909090914, 1.0, 2)
        