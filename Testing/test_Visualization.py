import unittest
from PyNite import FEModel3D
from PyNite.Visualization import Renderer, render_model
from PIL import Image

class TestRenderer(unittest.TestCase):

    def setUp(self):

        # Mock the model with necessary attributes
        self.model = FEModel3D()

        self.model.add_node('N1', 0, 0, 0)
        self.model.add_node('N2', 20, 0, 0)

        self.model.def_support('N1', True, True, True, True, True, True)
        self.model.def_support('N2', False, True, True, False, False, False)

        self.model.add_material('Steel', 29000/144, 11200/144, 0.3, 0.49, 36/144)
        
        self.model.add_section('Custom', 20, 100, 200, 150)
        
        self.model.add_member('M1', 'N1', 'N2', 'Steel', 'Custom')

        self.model.add_member_dist_load('M1', 'Fy', -1, -1, case='D')

        self.model.add_load_combo('1.4D', {'D': 1.4})

        self.model.analyze_linear()

        # Create the Renderer instance
        self.renderer = Renderer(self.model)

        # Set the renderer window to render offscreen
        self.renderer.window.SetOffScreenRendering(1)

        # Set up the load combo to be visualized
        self.renderer.combo_name = '1.4D'
        self.render_loads = True

    def test_render_model_old(self):
        self.render_model(self.model, combo_name='1.4D', render_loads=True, interact=False)
    
    def test_set_annotation_size(self):
        self.renderer.set_annotation_size(10)
        self.assertEqual(self.renderer.annotation_size, 10)

    def test_set_deformed_shape(self):
        self.renderer.set_deformed_shape(True)
        self.assertTrue(self.renderer.deformed_shape)

    def test_set_deformed_scale(self):
        self.renderer.set_deformed_scale(50)
        self.assertEqual(self.renderer.deformed_scale, 50)

    def test_set_render_nodes(self):
        self.renderer.set_render_nodes(False)
        self.assertFalse(self.renderer.render_nodes)

    def test_set_render_loads(self):
        self.renderer.set_render_loads(False)
        self.assertFalse(self.renderer.render_loads)

    def test_set_color_map(self):
        self.renderer.set_color_map('Mx')
        self.assertEqual(self.renderer.color_map, 'Mx')

    def test_set_combo_name(self):
        self.renderer.set_combo_name('Combo 2')
        self.assertEqual(self.renderer.combo_name, 'Combo 2')
        self.assertIsNone(self.renderer.case)

    def test_set_case(self):
        self.renderer.set_case('Case 2')
        self.assertEqual(self.renderer.case, 'Case 2')
        self.assertIsNone(self.renderer.combo_name)

    def test_set_show_labels(self):
        self.renderer.set_show_labels(False)
        self.assertFalse(self.renderer.labels)

    def test_set_scalar_bar(self):
        self.renderer.set_scalar_bar(True)
        self.assertTrue(self.renderer.scalar_bar)

    def test_set_scalar_bar_text_size(self):
        self.renderer.set_scalar_bar_text_size(30)
        self.assertEqual(self.renderer.scalar_bar_text_size, 30)

    def test_window_size_properties(self):
        self.renderer.window.SetSize(800, 600)
        self.assertEqual(self.renderer.window_width, 800)
        self.assertEqual(self.renderer.window_height, 600)

        self.renderer.window_width = 1024
        self.assertEqual(self.renderer.window.GetSize()[0], 1024)

        self.renderer.window_height = 768
        self.assertEqual(self.renderer.window.GetSize()[1], 768)

    def test_render_model(self):
        self.renderer.render_model(interact=False)
        self.renderer.update.assert_called_once()
        self.renderer.window.Render.assert_called_once()

    def test_screenshot(self):
        result = self.renderer.screenshot(filepath='console', interact=False)
        self.assertIsInstance(result, Image)
