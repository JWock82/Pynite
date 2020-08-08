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
    ''' Subclass of TestCase, which is executed by unittest.main

    Each test is a method with a name starting with 'test'.
    The setUp method is run prior to each test.
    The tearDown method is run after each test.
    '''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_AISC_Benckmark(self):
        # Create the cantilever model
        cantilever = FEModel3D()

        # Define the column and its properties
        L = 20 # ft
        H = 5 # kips lateral load
        P = 100 # kips axial load
        G = 11200*12**2 # shear modulus (ksf)
        E = 29000*12**2 # modulus of elasticity (ksf)
        I = 100/12**4 # moment of inertia (ft^4)

        # Break the column into several segments in order to capture P-little-delta effects
        num_segs = 5
        num_nodes = num_segs + 1
        for i in range(num_nodes):
            # Add nodes
            cantilever.AddNode(str(i+1), 0, i*L/(num_segs), 0)

        for i in range(num_segs):
            # Add members between nodes
            cantilever.AddMember(str(i+1), str(i+1), str(i+2), E, G, I, I, 200/12**4, 10/12**2)

        # Add a fixed support at the base of the column
        cantilever.DefineSupport('1', True, True, True, True, True, True)

        # Add a -10 kip axial load to the top of the column
        cantilever.AddNodeLoad(str(num_nodes), 'FY', -P)

        # Add a 5 kip lateral load to the top of the column
        cantilever.AddNodeLoad(str(num_nodes), 'FX', H)

        # Perform 2nd order analysis
        cantilever.Analyze_PDelta()

        # The moment at the base of the column
        calculated_moment = cantilever.GetNode('1').RxnMZ['Combo 1']

        # the deflection at the top of the column
        calculated_displacement = cantilever.GetNode(str(num_nodes)).DX['Combo 1']*12

        # Calculate the AISC benchmark problem solution:
        alpha = (P*L**2/(E*I))**0.5
        Mmax = H*L*(math.tan(alpha)/alpha)
        ymax = H*L**3/(3*E*I)*(3*(math.tan(alpha)-alpha)/alpha**3)
        
        # Compare the calculation results
        self.assertAlmostEqual(calculated_moment/Mmax, 1.0, 1)
        self.assertAlmostEqual(calculated_displacement/(ymax*12), 1.0, 1)

class Test_2D_Frame(unittest.TestCase):
    ''' Subclass of TestCase, which is executed by unittest.main

    Each test is a method with a name starting with 'test'.
    The setUp method is run prior to each test.
    The tearDown method is run after each test.
    '''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_AISC_Benckmark(self):
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
        print('Calculated displacements: ', frame.GetNode('N3').DY, frame.GetNode('N4').DY, frame.GetNode('N3').RZ, frame.GetNode('N4').RZ)
