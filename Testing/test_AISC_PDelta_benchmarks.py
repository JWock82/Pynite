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

class Test_AISC_Benchmark(unittest.TestCase):
    """
    AISC's benchmark problems used to determine if a second order analysis 
    procedure is rigorous enough.
    """

    def setUp(self):

        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_AISC_Benchmark(self):

        # Create the cantilever model
        cantilever = FEModel3D()

        # Define the column and its properties
        L = 20 # ft
        H = 5 # kips lateral load
        P = 100 # kips axial load
        G = 11200*12**2 # shear modulus (ksf)
        E = 29000*12**2 # modulus of elasticity (ksf)
        I = 100/12**4 # moment of inertia (ft^4)

        # Add nodes along the length of the column to capture P-little-delta effects
        num_nodes = 6
        for i in range(num_nodes):
            # Add nodes
            cantilever.add_node('N' + str(i+1), 0, i*L/6, 0)
        
        # Add the member
        cantilever.add_member('M1', 'N1', 'N6', E, G, I, I, 200/12**4, 10/12**2)

        # Add a fixed support at the base of the column
        cantilever.def_support('N1', True, True, True, True, True, True)

        # Add a -10 kip axial load to the top of the column
        cantilever.add_node_load('N6', 'FY', -P)

        # Add a 5 kip lateral load to the top of the column
        cantilever.add_node_load('N6', 'FX', H)

        # Perform 2nd order analysis
        cantilever.analyze_PDelta()

        from PyNite.Visualization import Renderer
        renderer = Renderer(cantilever)
        renderer.annotation_size = 0.5
        renderer.window_width = 750
        renderer.window_height = 750
        renderer.deformed_shape = True
        renderer.render_model()

        # The moment at the base of the column
        calculated_moment = cantilever.Nodes['N1'].RxnMZ['Combo 1']

        # the deflection at the top of the column
        calculated_displacement = cantilever.Nodes['N6'].DX['Combo 1']*12

        # Calculate the AISC benchmark problem solution:
        alpha = (P*L**2/(E*I))**0.5
        Mmax = H*L*(math.tan(alpha)/alpha)
        ymax = H*L**3/(3*E*I)*(3*(math.tan(alpha)-alpha)/alpha**3)
        
        # Compare the calculation results
        self.assertAlmostEqual(calculated_moment/Mmax, 1.0, 1)
        self.assertAlmostEqual(calculated_displacement/(ymax*12), 1.0, 1)
