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
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Problem 5.30
        # Units for this model are kips and inches
        frame = FEModel3D()
        # Define the nodes
        frame.add_node('N1', 0, 0, 0)
        frame.add_node('N2', 0, 30*12, 0)
        frame.add_node('N3', 15*12, 40*12, 0)
        frame.add_node('N4', 35*12, 40*12, 0)
        frame.add_node('N5', 50*12, 30*12, 0)
        frame.add_node('N6', 50*12, 0, 0)
        # Define the supports
        frame.def_support('N1', True, True, True, True, True, True)
        frame.def_support('N6', True, True, True, True, True, True)
        # Create members (all members will have the same properties in this example)
        J = 250
        Iy = 250
        Iz = 200
        E = 30000
        G = 250
        A = 12
        frame.add_member('M1', 'N1', 'N2', E, G, Iy, Iz, J, A)
        frame.add_member('M2', 'N2', 'N3', E, G, Iy, Iz, J, A)
        frame.add_member('M3', 'N3', 'N4', E, G, Iy, Iz, J, A)
        frame.add_member('M4', 'N4', 'N5', E, G, Iy, Iz, J, A)
        frame.add_member('M5', 'N5', 'N6', E, G, Iy, Iz, J, A)
        # Add nodal loads
        frame.add_node_load('N3', 'FY', -30)
        frame.add_node_load('N4', 'FY', -30)
        # Analyze the model
        frame.analyze()
        # subTest context manager prints which portion fails, if any
        correct_values = [('N1', {'RxnFX': 11.6877,
                                  'RxnFY': 30,
                                  'RxnMZ': -1810.0745}),
                          ('N6', {'RxnFX': -11.6877,
                                  'RxnFY': 30,
                                  'RxnMZ': 1810.0745})]
        for name, values in correct_values:
            with self.subTest(node=name):
                node = frame.Nodes[name]
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(node.RxnFX['Combo 1']/values['RxnFX'], 1.0, 2)
                self.assertAlmostEqual(node.RxnFY['Combo 1']/values['RxnFY'], 1.0, 2)
                self.assertAlmostEqual(node.RxnMZ['Combo 1']/values['RxnMZ'], 1.0, 2)

    def test_XY_member_ptload(self):
        frame = FEModel3D()
        # Add nodes
        frame.add_node('N1', 0, 0, 0)        # ft
        frame.add_node('N2', 0, 7.667, 0)    # ft
        frame.add_node('N3', 7.75, 7.667, 0) # ft
        frame.add_node('N4', 7.75, 0, 0)     # ft
        # Add supports
        frame.def_support('N1', True, True, True, True, True, False)
        frame.def_support('N4', True, True, True, True, True, False)
        # Define material and section properties for a W8x24
        E = 29000*12**2     # ksf
        G = 1111200*12**2   # ksf
        Iy = 18.3/12**4     # ft^4
        Iz = 82.7/12**4     # ft^4
        J = 0.346/12**4     # ft^4
        A = 5.26/12**2      # in^2
        # Define members
        frame.add_member('M1', 'N1', 'N2', E, G, Iy, Iz, J, A)
        frame.add_member('M2', 'N2', 'N3', E, G, Iy, Iz, J, A)
        frame.add_member('M3', 'N4', 'N3', E, G, Iy, Iz, J, A)
        # Add loads to the frame
        frame.add_member_pt_load('M2', 'Fy', -5, 7.75/2)       # 5 kips @ midspan
        frame.add_member_dist_load('M2', 'Fy', -0.024, -0.024) # W8x24 self-weight
        # Analyze the frame
        frame.analyze()
        calculated_RZ = frame.Nodes['N1'].RZ['Combo 1']
        # Update the expected value to an appropriate precision
        expected_RZ = 0.00022794540510395617
        self.assertAlmostEqual(calculated_RZ/expected_RZ, 1.0, 2)
    
    def test_YZ_gravity_load(self):
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Problem 5.30
        # Units for this model are kips and inches
        frame = FEModel3D()
        # Define the nodes
        frame.add_node('N1', 0, 0, 0)
        frame.add_node('N2', 0, 30*12, 0)
        frame.add_node('N3', 0, 40*12, 15*12)
        frame.add_node('N4', 0, 40*12, 35*12)
        frame.add_node('N5', 0, 30*12, 50*12)
        frame.add_node('N6', 0, 0, 50*12)
        # Define the supports
        frame.def_support('N1', True, True, True, True, True, True)
        frame.def_support('N6', True, True, True, True, True, True)
        # Create members (all members will have the same properties in this example)
        J = 250
        Iy = 250
        Iz = 200
        E = 30000
        G = 250
        A = 12
        frame.add_member('M1', 'N1', 'N2', E, G, Iz, Iy, J, A)
        frame.add_member('M2', 'N2', 'N3', E, G, Iy, Iz, J, A)
        frame.add_member('M3', 'N3', 'N4', E, G, Iy, Iz, J, A)
        frame.add_member('M4', 'N4', 'N5', E, G, Iy, Iz, J, A)
        frame.add_member('M5', 'N5', 'N6', E, G, Iz, Iy, J, A)
        # Add nodal loads
        frame.add_node_load('N3', 'FY', -30)
        frame.add_node_load('N4', 'FY', -30)
        # Analyze the model
        frame.analyze()
        # subTest context manager prints which portion fails, if any
        # Check reactions at N1 and N6
        correct_reactions = [('N1', {'RxnFZ': 11.6877,
                                     'RxnFY': 30,
                                     'RxnMX': 1810.0745}),
                             ('N6', {'RxnFZ': -11.6877,
                                     'RxnFY': 30,
                                     'RxnMX': -1810.0745})]
        for name, values in correct_reactions:
            with self.subTest(node=name):
                node = frame.Nodes[name]
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(node.RxnFZ['Combo 1']/values['RxnFZ'], 1.0, 2)
                self.assertAlmostEqual(node.RxnFY['Combo 1']/values['RxnFY'], 1.0, 2)
                self.assertAlmostEqual(node.RxnMX['Combo 1']/values['RxnMX'], 1.0, 2)
        # Check displacements at N3 and N4
        correct_displacements = [('N3', {'DY': -6.666757,
                                         'RX':  0.032}),
                                 ('N4', {'DY': -6.666757,
                                         'RX': -0.032})]
        for name, values in correct_displacements:
            with self.subTest(node=name):
                node = frame.Nodes[name]
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(node.DY['Combo 1']/values['DY'], 1.0, 2)
                self.assertAlmostEqual(node.RX['Combo 1']/values['RX'], 1.0, 2)

    def test_XZ_ptload(self):
        # A simply supported beam with a point load.
        # Units used in this example are inches, and kips
        SimpleBeam = FEModel3D()
        # Add nodes (14 ft = 168 in apart)
        SimpleBeam.add_node("N1", 0, 0, 0)
        SimpleBeam.add_node("N2", 0, 0, 168)
        # Add a beam with the following properties:
        A = 20
        E = 29000
        G = 11400
        Iy = 100
        Iz = 150
        J = 250
        SimpleBeam.add_member("M1", "N1", "N2", E, G, Iy, Iz, J, A)
        # Provide simple supports
        SimpleBeam.def_support("N1", True, True, True, False, False, True)
        SimpleBeam.def_support("N2", True, True, True, False, False, False)
        # Add a point load of 5 kips at the midspan of the beam
        SimpleBeam.add_member_pt_load("M1", "Fy", 5, 7 * 12)
        # Analyze the beam
        SimpleBeam.analyze(False)
        # Print reactions at each end of the beam
        correct_reactions = [('N1', -2.5),
                             ('N2', -2.5)]
        for node_name, rxn in correct_reactions:
            with self.subTest(node=node_name):
                calculated_reaction = SimpleBeam.Nodes[node_name].RxnFY['Combo 1']
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(calculated_reaction/rxn, 1.0, 2)
                