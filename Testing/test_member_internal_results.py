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

class TestMemberInternalResults(unittest.TestCase):
    """
    Tests for member internal results
    """

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_beam_shear(self):
        """
        Units for this model are kips and feet
        """

        # Define a new beam
        beam = FEModel3D()

        # Define the nodes
        beam.add_node('N1', 0, 0, 0)
        beam.add_node('N2', 10, 0, 0)

        # Define the supports
        beam.def_support('N1', True, True, True, True, False, False)
        beam.def_support('N2', True, True, True, False, False, False)

        # Define beam section proerties
        J = 400/12**4
        Iy = 200/12**4
        Iz = 200/12**4
        E = 29000*144
        G = 11200*144
        A = 12/12**2

        # Create the beam
        beam.add_member('M1', 'N1', 'N2', E, G, Iy, Iz, J, A)

        # Add a mid-span node
        beam.add_node('N3', 5, 0, 0)

        # Add a member distributed load
        beam.add_member_dist_load('M1', 'FY', -0.5, -0.5, case='D')
        beam.add_load_combo('D', {'D':1.0})

        # Analyze the model
        beam.analyze_linear()

        # from PyNite.Visualization import Renderer
        # renderer = Renderer(beam)
        # renderer.combo_name = 'D'
        # renderer.annotation_size = 1
        # renderer.render_model()

        # Check the shear at midspan - it should be zero
        self.assertAlmostEqual(beam.Members['M1'].shear('Fy', 5, 'D'), 0, 2)
