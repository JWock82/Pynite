"""
Test the needs_update flag functionality for Meshes, ShearWalls, and MatFoundations.
"""

import unittest
import sys
sys.path.insert(0, '..')

from Pynite.FEModel3D import FEModel3D


class TestNeedsUpdateFlag(unittest.TestCase):
    """Test cases for the needs_update flag on Mesh, ShearWall, and MatFoundation objects."""

    def test_mesh_needs_update_after_adding_opening(self):
        """Test that needs_update is set to True when adding opening to generated mesh."""
        # Create a model and add a rectangular mesh
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_rectangle_mesh('M1', 1, 10, 10, 1, 'Concrete', 0.5, 1)
        
        # Initially, mesh should not be generated
        mesh = model.meshes['M1']
        self.assertFalse(mesh.is_generated)
        self.assertFalse(mesh.needs_update)
        
        # Generate the mesh
        mesh.generate()
        self.assertTrue(mesh.is_generated)
        self.assertFalse(mesh.needs_update)
        
        # Add an opening to the already-generated mesh
        mesh.add_rect_opening('O1', 2, 2, 2, 2)
        self.assertTrue(mesh.is_generated)
        self.assertTrue(mesh.needs_update)
        
    def test_shear_wall_needs_update_after_adding_opening(self):
        """Test that needs_update is set to True when adding opening to generated shear wall."""
        # Create a model and add a shear wall
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        
        model.add_shear_wall('SW1', 1, 20, 10, 0.25, 'Concrete')
        wall = model.shear_walls['SW1']
        
        # Initially, wall should not be generated
        self.assertFalse(wall.is_generated)
        self.assertFalse(wall.needs_update)
        
        # Generate the wall
        wall.generate()
        self.assertTrue(wall.is_generated)
        self.assertFalse(wall.needs_update)
        
        # Add an opening to the already-generated wall
        wall.add_opening('D1', 5, 2, 3, 4)
        self.assertTrue(wall.needs_update)
        
    def test_shear_wall_needs_update_after_adding_support(self):
        """Test that needs_update is set to True when adding support to generated wall."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        
        model.add_shear_wall('SW1', 1, 20, 10, 0.25, 'Concrete')
        wall = model.shear_walls['SW1']
        wall.generate()
        
        self.assertTrue(wall.is_generated)
        self.assertFalse(wall.needs_update)
        
        # Add support to generated wall
        wall.add_support(0, 0, 20)
        self.assertTrue(wall.needs_update)
        
    def test_shear_wall_needs_update_after_adding_story(self):
        """Test that needs_update is set to True when adding story to generated wall."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        
        model.add_shear_wall('SW1', 1, 20, 10, 0.25, 'Concrete')
        wall = model.shear_walls['SW1']
        wall.generate()
        
        self.assertTrue(wall.is_generated)
        self.assertFalse(wall.needs_update)
        
        # Add story to generated wall
        wall.add_story('Level 2', 10, 0, 20)
        self.assertTrue(wall.needs_update)
        
    def test_shear_wall_needs_update_after_adding_flange(self):
        """Test that needs_update is set to True when adding flange to generated wall."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        
        model.add_shear_wall('SW1', 1, 20, 10, 0.25, 'Concrete')
        wall = model.shear_walls['SW1']
        wall.generate()
        
        self.assertTrue(wall.is_generated)
        self.assertFalse(wall.needs_update)
        
        # Add flange to generated wall
        wall.add_flange(0.25, 2, 5, 0, 10, 'Concrete', '+z')
        self.assertTrue(wall.needs_update)
        
    def test_shear_wall_needs_update_after_assigning_material(self):
        """Test that needs_update is set to True when assigning material to generated wall."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_material('Concrete2', 4000, 1440, 0.2, 0.15)
        
        model.add_shear_wall('SW1', 1, 20, 10, 0.25, 'Concrete')
        wall = model.shear_walls['SW1']
        wall.generate()
        
        self.assertTrue(wall.is_generated)
        self.assertFalse(wall.needs_update)
        
        # Assign new material to generated wall
        wall.asign_material('Concrete2', 0.3, 5, 15, 0, 10)
        self.assertTrue(wall.needs_update)
        
    def test_mat_foundation_needs_update_after_adding_opening(self):
        """Test that needs_update is set to True when adding opening to generated mat foundation."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        
        model.add_mat_foundation('Mat1', 1, 20, 30, 1, 'Concrete', 100)
        mat = model.mats['Mat1']
        
        # Initially, mat should not be generated
        self.assertFalse(mat.is_generated)
        self.assertFalse(mat.needs_update)
        
        # Generate the mat
        mat.generate()
        self.assertTrue(mat.is_generated)
        self.assertFalse(mat.needs_update)
        
        # Add an opening to the already-generated mat
        mat.add_rect_opening('O1', 5, 5, 10, 10)
        self.assertTrue(mat.needs_update)
        
    def test_mat_foundation_needs_update_after_adding_point_load(self):
        """Test that needs_update is set to True when adding point load to generated mat."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        
        model.add_mat_foundation('Mat1', 1, 20, 30, 1, 'Concrete', 100)
        mat = model.mats['Mat1']
        mat.generate()
        
        self.assertTrue(mat.is_generated)
        self.assertFalse(mat.needs_update)
        
        # Add point load to generated mat
        mat.add_mat_pt_load([10, 15], 'FY', -100)
        self.assertTrue(mat.needs_update)
        
    def test_mesh_regeneration_clears_needs_update(self):
        """Test that regenerating a mesh clears the needs_update flag."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_rectangle_mesh('M1', 1, 10, 10, 1, 'Concrete', 0.5, 1)
        
        mesh = model.meshes['M1']
        mesh.generate()
        mesh.add_rect_opening('O1', 2, 2, 2, 2)
        
        # After adding opening, needs_update should be True
        self.assertTrue(mesh.needs_update)
        
        # Regenerate the mesh
        mesh.generate()
        
        # After regeneration, needs_update should be False
        self.assertFalse(mesh.needs_update)
        self.assertTrue(mesh.is_generated)


if __name__ == '__main__':
    unittest.main()
