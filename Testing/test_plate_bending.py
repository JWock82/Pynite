# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
from PyNite.Mesh import RectangleMesh
import math
import sys
from io import StringIO
from numpy import allclose

class Test_Plates(unittest.TestCase):
    """
    Tests for plate bending
    """

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_plate_displacement(self):
        """
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Example 12.1
        # Units for this model are pounds and inches
        """

        plate_model= FEModel3D()

        plate_model.add_node('N1', 0, 0, 0)
        plate_model.add_node('N2', 10, 0, 0)
        plate_model.add_node('N3', 20, 0, 0)
        plate_model.add_node('N4', 0, 10, 0)
        plate_model.add_node('N5', 10, 10, 0)
        plate_model.add_node('N6', 20, 10, 0)
        plate_model.add_node('N7', 0, 20, 0)
        plate_model.add_node('N8', 10, 20, 0)
        plate_model.add_node('N9', 20, 20, 0)

        plate_model.add_material('Steel', 30000000, 11200000, 0.3, 0.284)

        plate_model.add_plate('P1', 'N1', 'N2', 'N5', 'N4', 0.1, 'Steel')
        plate_model.add_plate('P2', 'N2', 'N3', 'N6', 'N5', 0.1, 'Steel')
        plate_model.add_plate('P3', 'N4', 'N5', 'N8', 'N7', 0.1, 'Steel')
        plate_model.add_plate('P4', 'N5', 'N6', 'N9', 'N8', 0.1, 'Steel')

        plate_model.add_node_load('N5', 'FZ', -100)

        plate_model.def_support('N1', True, True, True, True, True, True)
        plate_model.def_support('N2', True, True, True, True, True, True)
        plate_model.def_support('N3', True, True, True, True, True, True)
        plate_model.def_support('N4', True, True, True, True, True, True)
        plate_model.def_support('N6', True, True, True, True, True, True)
        plate_model.def_support('N7', True, True, True, True, True, True)
        plate_model.def_support('N8', True, True, True, True, True, True)
        plate_model.def_support('N9', True, True, True, True, True, True)

        plate_model.def_support('N5', True, True, False, False, False, True)

        # Check to see if the stiffness matrix for each plate is symmetric
        # print(allclose(plate_model.Plates[0].K(), plate_model.Plates[0].K().T))
        # print(allclose(plate_model.Plates[1].K(), plate_model.Plates[1].K().T))
        # print(allclose(plate_model.Plates[2].K(), plate_model.Plates[2].K().T))
        # print(allclose(plate_model.Plates[3].K(), plate_model.Plates[3].K().T))

        # Check to see if the global stiffness matrix is symmetric
        # print(allclose(plate_model.K(Renumber=True), plate_model.K(Renumber=False).T))

        plate_model.analyze(check_statics=True, sparse=False)
        # Test: displacement of N5 in Z direction
        calculated_displacement = plate_model.Nodes['N5'].DZ['Combo 1']
        expected_displacement = -0.0861742424242424
        self.assertAlmostEqual(calculated_displacement/expected_displacement, 1.0, 2)

    def test_hydrostatic_plate(self):

        # Create the model
        plate_model = FEModel3D()

        # Define geometry
        t = 1  # ft
        mesh_size = 1  # ft
        a = 10  # ft
        b = 15  # ft

        # Define a material
        E = 57000*math.sqrt(4500)*12**2  # psf
        G = 0.4*E  # psf
        nu = 1/6
        rho = 150  # psf
        plate_model.add_material('Concrete', E, G, nu, rho)
        
        # Generate the mesh of plates
        plate_model.add_rectangle_mesh('MSH1', mesh_size, a, b, t, 'Concrete', kx_mod=1, ky_mod=1,
                                       element_type='Rect')
        
        # Generate the mesh
        plate_model.Meshes['MSH1'].generate()

        # Add supports to the sides and base of the wall
        for node in plate_model.Nodes.values():
            if node.X == 0 or node.X == a or node.Y == 0:
                plate_model.def_support(node.name, True, True, True, True, True, True)
        
        # Add hydrostatic loads to the elements
        for element in plate_model.Plates.values():
            Yavg = (element.i_node.Y + element.j_node.Y + element.m_node.Y + element.n_node.Y)/4
            p = 62.4*(b - Yavg)
            plate_model.add_plate_surface_pressure(element.name, p, 'Hydrostatic')
        
        # Add a load combination to the model
        plate_model.add_load_combo('F', {'Hydrostatic': 1.0})
        
        # Analyze the model
        plate_model.analyze()

        # Get the maximum deflection in the model at the top of the wall
        DZ_calcd = max([node.DZ['F'] for node in plate_model.Nodes.values() if node.Y == b])
        
        # Find the maximum deflection at the top of the wall from Timoshenko's Table 45
        q = 62.4*b
        D = E*t**3/(12*(1 - nu**2))
        DZ_expected = 0.00042*q*a**4/D

        # Check that the PyNite calculated values are within 15% of the Timoshenko calculated
        # values.
        self.assertLess(abs(DZ_calcd/DZ_expected - 1), 0.15, 'Failed Timoshenko rectangle hydrostatic test.')
    
    def test_hydrostatic_quad(self):

        # Create the model
        quad_model = FEModel3D()

        # Define geometry
        t = 1  # ft
        mesh_size = 1  # ft
        a = 10  # ft
        b = 15  # ft

        # Define a material
        E = 57000*math.sqrt(4500)*12**2  # psf
        G = 0.4*E  # psf
        nu = 1/6
        rho = 150  # psf
        quad_model.add_material('Concrete', E, G, nu, rho)
        
        # Generate the mesh of quads
        quad_model.add_rectangle_mesh('MSH1', mesh_size, a, b, t, 'Concrete', kx_mod=1, ky_mod=1,
                                      element_type='Quad')
        quad_model.Meshes['MSH1'].generate()

        # Add supports to the sides and base of the wall
        for node in quad_model.Nodes.values():
            if node.X == 0 or node.X == a or node.Y == 0:
                quad_model.def_support(node.name, True, True, True, True, True, True)
        
        # Add hydrostatic loads to the elements
        for element in quad_model.Quads.values():
            Yavg = (element.i_node.Y + element.j_node.Y + element.m_node.Y + element.n_node.Y)/4
            p = 62.4*(b - Yavg)
            quad_model.add_quad_surface_pressure(element.name, p, 'Hydrostatic')
        
        # Add a load combination to the model
        quad_model.add_load_combo('F', {'Hydrostatic': 1.0})
        
        # Analyze the model
        quad_model.analyze()

        # Get the maximum deflection in the model at the top of the wall
        DZ_calcd = max([node.DZ['F'] for node in quad_model.Nodes.values() if node.Y == b])
        
        # Find the maximum deflection at the top of the wall from Timoshenko's Table 45
        q = 62.4*b
        D = E*t**3/(12*(1 - nu**2))
        DZ_expected = 0.00042*q*a**4/D

        # Check that the PyNite calculated values are within 15% of the Timoshenko calculated
        # values.
        self.assertLess(abs(DZ_calcd/DZ_expected - 1), 0.15, 'Failed Timoshenko quadrilateral hydrostatic test.') 