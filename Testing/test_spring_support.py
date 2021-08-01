# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2021 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import sys
from io import StringIO

class Test_Spring_Supports(unittest.TestCase):
    ''' Tests of spring supports.'''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
        
    def test_beam_on_elastic_foundation(self):
        '''        
        Matrix Structural Analysis, 2nd Edition
        William McGuire, Richard H. Gallagher, Ronald D. Ziemian
        Example 4.15
        Units for this model are kips and inches
        '''

        # Create a new model
        boef = FEModel3D()

        # Define nodes
        for i in range(17):

            # Add nodes spaced at 15"
            boef.AddNode('N' + str(i + 1), i*15, 0, 0)

            # Add supports to the nodes
            if i == 0 or i == 16:
                boef.DefineSupport('N' + str(i + 1), True, True, True, True, False, False)
            else:
                boef.DefineSupport('N' + str(i + 1), False, 22.5, False, False, False, False)

        # Define member material properties
        E = 29000   # ksi
        G = 11200   # ksi
        A = 10.3    # in^2
        Iz = 128.5  # in^4 (strong axis)
        Iy = 42.6   # in^4 (weak axis)
        J = 0.769   # in^4

        # Define members
        for i in range(16):

            # Add the members
            boef.AddMember('M' + str(i + 1), 'N' + str(i + 1), 'N' + str(i + 2), E, G, Iy, Iz, J, A)
        
        # Add a point load at midspan
        boef.AddNodeLoad('N9', 'FY', -40)

        # Analyze the model
        boef.Analyze()

        print(boef.Members['M8'].MinMoment('Mz'))
        print(boef.Members['M8'].MaxMoment('Mz'))

        # Check that results are within 5% of the expected answer
        self.assertLess(abs(boef.Nodes['N9'].DY['Combo 1']/(-0.238)) - 1, 0.05, 'Failed beam on elastic foundation test.')
        self.assertLess(abs(boef.Members['M8'].MinMoment('Mz')/547) - 1, 0.05, 'Failed beam on elastic foundation test.')