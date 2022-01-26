# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2022 D. Craig Brinck, SE; tamalone1
"""

from msilib.schema import tables
import unittest
from PyNite import FEModel3D
from PyNite.Mesh import CylinderMesh, RectangleMesh
import sys
from io import StringIO

class TestShearWalls(unittest.TestCase):
    
    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
    
    def test_quad_shear_wall(self):

        sw = FEModel3D()

        E = 57000*(4000)**0.5/1000*12**2
        nu = 0.17
        G = E/(2*(1 + nu))
        t = 1
        L = 10
        H = 20
        A = L*t
        I = t*L**3/12

        mesh_size = 1
        mesh = RectangleMesh(mesh_size, L, H, t, E, nu, element_type='Quad')
        mesh.generate()

        sw.add_mesh(mesh)

        V = 1000
        for node in sw.Nodes.values():
            if node.Y == 0:
                sw.def_support(node.Name, True, True, True, True, True, True)
            elif node.Y == H:
                sw.add_node_load(node.Name, 'FX', V/11)
        
        sw.analyze()

        # Calculated solution
        delta1 = max([node.DX['Combo 1'] for node in sw.Nodes.values()])

        # Theoretical solution
        delta2 = V*H**3/(3*E*I) + 1.2*V*H/(G*A)
        
        # Check that the solution matches the theoretical solution within 0.1%
        self.assertLess(abs(1 - delta1/delta2), 0.001, 'Failed quad shear wall test.')

    def test_rect_shear_wall(self):

        sw = FEModel3D()

        E = 57000*(4000)**0.5/1000*12**2
        nu = 0.17
        G = E/(2*(1 + nu))
        t = 1
        L = 10
        H = 20
        A = L*t
        I = t*L**3/12

        mesh_size = 1
        mesh = RectangleMesh(mesh_size, L, H, t, E, nu, element_type='Rect')
        mesh.generate()

        sw.add_mesh(mesh)

        V = 1000
        for node in sw.Nodes.values():
            if node.Y == 0:
                sw.def_support(node.Name, True, True, True, True, True, True)
            elif node.Y == H:
                sw.add_node_load(node.Name, 'FX', V/11)
        
        sw.analyze()

        # Calculated solution
        delta1 = max([node.DX['Combo 1'] for node in sw.Nodes.values()])

        # Theoretical solution
        delta2 = V*H**3/(3*E*I) + 1.2*V*H/(G*A)
        
        # Check that the solution matches the theoretical solution within 0.1%
        self.assertLess(abs(1 - delta1/delta2), 0.001, 'Failed rect plate shear wall test.')

    def test_cracked_rect_shear_wall(self):

        sw = FEModel3D()

        E = 57000*(4000)**0.5/1000*12**2
        nu = 0.17
        G = E/(2*(1 + nu))
        t = 1
        L = 10
        H = 20
        A = L*t
        I = 0.35*t*L**3/12

        mesh_size = 1
        mesh = RectangleMesh(mesh_size, L, H, t, E, nu, ky_mod=0.35, element_type='Rect')
        mesh.generate()

        sw.add_mesh(mesh)

        V = 1000
        for node in sw.Nodes.values():
            if node.Y == 0:
                sw.def_support(node.Name, True, True, True, True, True, True)
            elif node.Y == H:
                sw.add_node_load(node.Name, 'FX', V/11)
        
        sw.analyze()

        # Calculated solution
        delta1 = max([node.DX['Combo 1'] for node in sw.Nodes.values()])

        # Theoretical solution
        delta2 = V*H**3/(3*E*I) + 1.2*V*H/(G*A)
        
        # Check that the solution matches the theoretical solution within 2%
        self.assertLess(abs(1 - delta1/delta2), 0.02, 'Failed cracked rect plate shear wall test.')