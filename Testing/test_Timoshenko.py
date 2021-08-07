import unittest
from PyNite import FEModel3D
from PyNite.Mesh import FrustrumMesh, CylinderMesh, RectangleMesh
import math
import sys
from io import StringIO

class Test_Timoshenko(unittest.TestCase):
    # Tests against problems with known solutions from `Theory of Plates & Shells` by Timoshenko

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal behavior
        sys.stdout = sys.__stdout__

    def test_hydrostatic_plate(self):

        # Establish problem parameters
        t = 1  # ft
        E = 57000*math.sqrt(4500)*12**2  # psf
        nu = 1/6
        mesh_size = 1  # ft
        a = 10  # ft
        b = 15  # ft
        
        # Generate the mesh of plates
        plate_mesh = RectangleMesh(t, E, nu, mesh_size, a, b, element_type='Rect')

        # Create the model and add the plates
        plate_model = FEModel3D()
        plate_model.add_mesh(plate_mesh)

        # Add supports to the sides and base of the wall
        for node in plate_model.Nodes.values():
            if node.X == 0 or node.X == a or node.Y == 0:
                plate_model.def_support(node.Name, True, True, True, True, True, True)
        
        # Add hydrostatic loads to the elements
        for element in plate_model.Plates.values():
            Yavg = (element.iNode.Y + element.jNode.Y + element.mNode.Y + element.nNode.Y)/4
            p = 62.4*(b - Yavg)
            plate_model.add_plate_surface_pressure(element.Name, p, 'Hydrostatic')
        
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

        # Establish problem parameters
        t = 1  # ft
        E = 57000*math.sqrt(4500)*12**2  # psf
        nu = 1/6
        mesh_size = 1  # ft
        a = 10  # ft
        b = 15  # ft
        
        # Generate the mesh of plates
        plate_mesh = RectangleMesh(t, E, nu, mesh_size, a, b, element_type='Quad')

        # Create the model and add the plates
        plate_model = FEModel3D()
        plate_model.add_mesh(plate_mesh)

        # Add supports to the sides and base of the wall
        for node in plate_model.Nodes.values():
            if node.X == 0 or node.X == a or node.Y == 0:
                plate_model.def_support(node.Name, True, True, True, True, True, True)
        
        # Add hydrostatic loads to the elements
        for element in plate_model.Quads.values():
            Yavg = (element.iNode.Y + element.jNode.Y + element.mNode.Y + element.nNode.Y)/4
            p = 62.4*(b - Yavg)
            plate_model.add_quad_surface_pressure(element.Name, p, 'Hydrostatic')
        
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
        self.assertLess(abs(DZ_calcd/DZ_expected - 1), 0.15, 'Failed Timoshenko quadrilateral hydrostatic test.')    
