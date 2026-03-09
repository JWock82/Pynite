"""
Test that analysis only regenerates meshes that need it.
"""

import unittest
import sys
sys.path.insert(0, '..')

from Pynite.FEModel3D import FEModel3D


class TestAnalysisRegeneration(unittest.TestCase):
    """Test cases for selective mesh regeneration during analysis."""

    def test_analysis_regenerates_mesh_needing_update(self):
        """Test that analysis regenerates a mesh when needs_update is True."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_rectangle_mesh('M1', 1, 10, 10, 1, 'Concrete', 0.5, 1)
        
        mesh = model.meshes['M1']
        mesh.generate()
        
        # Record the number of elements after first generation
        initial_element_count = len(mesh.elements)
        
        # Add an opening - this sets needs_update to True
        mesh.add_rect_opening('O1', 2, 2, 2, 2)
        self.assertTrue(mesh.needs_update)
        
        # Analyze the model - should regenerate the mesh
        model.analyze()
        
        # After analysis, needs_update should be False and mesh should have fewer elements
        self.assertFalse(mesh.needs_update)
        self.assertTrue(mesh.is_generated)
        self.assertLess(len(mesh.elements), initial_element_count)
        
    def test_analysis_skips_mesh_not_needing_update(self):
        """Test that analysis skips regeneration when mesh doesn't need it."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_rectangle_mesh('M1', 1, 10, 10, 1, 'Concrete', 0.5, 1)
        
        mesh = model.meshes['M1']
        mesh.generate()
        
        # Store element IDs to verify they're not regenerated
        original_element_ids = set(mesh.elements.keys())
        
        # Analyze without making changes
        model.analyze()
        
        # Element IDs should be unchanged (no regeneration occurred)
        self.assertEqual(set(mesh.elements.keys()), original_element_ids)
        
    def test_analysis_generates_ungenerated_mesh(self):
        """Test that analysis generates meshes that haven't been generated yet."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_rectangle_mesh('M1', 1, 10, 10, 1, 'Concrete', 0.5, 1)
        
        mesh = model.meshes['M1']
        self.assertFalse(mesh.is_generated)
        
        # Analyze without manually generating
        model.analyze()
        
        # Mesh should now be generated
        self.assertTrue(mesh.is_generated)
        self.assertFalse(mesh.needs_update)
        self.assertGreater(len(mesh.elements), 0)
        
    def test_analysis_regenerates_shear_wall_needing_update(self):
        """Test that analysis regenerates a shear wall when needs_update is True."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_shear_wall('SW1', 1, 20, 10, 0.25, 'Concrete')
        
        wall = model.shear_walls['SW1']
        wall.generate()
        
        initial_element_count = len(model.meshes['SW1'].elements)
        
        # Add an opening - sets needs_update to True
        wall.add_opening('D1', 5, 2, 3, 4)
        self.assertTrue(wall.needs_update)
        
        # Analyze the model
        model.analyze()
        
        # After analysis, needs_update should be False
        self.assertFalse(wall.needs_update)
        self.assertTrue(wall.is_generated)
        # Should have fewer elements due to opening
        self.assertLess(len(model.meshes['SW1'].elements), initial_element_count)
        
    def test_analysis_regenerates_mat_foundation_needing_update(self):
        """Test that analysis regenerates a mat foundation when needs_update is True."""
        model = FEModel3D()
        model.add_material('Concrete', 3605, 1440, 0.17, 0.145)
        model.add_mat_foundation('Mat1', 2, 20, 30, 1, 'Concrete', 100)
        
        mat = model.mats['Mat1']
        mat.generate()
        
        initial_element_count = len(mat.elements)
        initial_node_count = len(mat.nodes)
        
        # Add an opening - sets needs_update to True
        mat.add_rect_opening('O1', 5, 5, 10, 10)
        self.assertTrue(mat.needs_update)
        
        # Analyze the model
        model.analyze()
        
        # After analysis, needs_update should be False
        self.assertFalse(mat.needs_update)
        self.assertTrue(mat.is_generated)
        # Element count might increase or decrease depending on mesh refinement,
        # but it should be different after regeneration with opening
        self.assertNotEqual(len(mat.elements), initial_element_count)


if __name__ == '__main__':
    unittest.main()
