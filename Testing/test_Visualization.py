from io import BytesIO
import unittest
from unittest.mock import MagicMock
import vtk
from Pynite import FEModel3D
from Pynite.Visualization import Renderer
from IPython.display import Image

class TestRenderer(unittest.TestCase):

    def setUp(self):

        # Create a model to visualize
        self.beam_model = FEModel3D()

        # Add nodes
        self.beam_model.add_node('N1', 0, 0, 0)
        self.beam_model.add_node('N2', 10, 0, 0)
        self.beam_model.add_node('N3', 20, 0, 0)
        self.beam_model.add_node('N4', 30, 0, 0)
        self.beam_model.add_node('N5', 40, 0, 0)

        # Add various support conditions for rendering
        self.beam_model.def_support('N1', True, True, True, True, True, True)      # Fully Fixed
        self.beam_model.def_support('N2', True, True, True, False, False, False)   # Pinned
        self.beam_model.def_support('N3', True, True, False, False, False, False)  # Fixed in XY
        self.beam_model.def_support('N4', False, True, True, False, False, False)  # Fixed in YZ
        self.beam_model.def_support('N5', False, False, False, True, True, True)   # Rotational Fixity Only

        self.beam_model.add_material('Steel', 29000/144, 11200/144, 0.3, 0.49, 36/144)
        
        self.beam_model.add_section('Custom', 20, 100, 200, 150)
        
        self.beam_model.add_member('M1', 'N1', 'N5', 'Steel', 'Custom')

        self.beam_model.add_member_dist_load('M1', 'Fy', -1, -1, case='D')

        self.beam_model.add_load_combo('1.4D', {'D': 1.4})

        self.beam_model.analyze_linear()

        # Create the Renderer instance
        self.renderer = Renderer(self.beam_model)

        # Set the renderer window to render offscreen
        self.renderer.window.SetOffScreenRendering(1)

        # Set up the load combo to be visualized
        self.renderer.combo_name = '1.4D'
        self.renderer.render_loads = True

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

        # Mock the update and render methods
        # self.renderer.update = MagicMock()
        # self.renderer.window.Render = MagicMock()

        # Call the render_model method
        self.renderer.render_model(interact=False)

        # Assert that the update and render methods were called
        # self.renderer.update.assert_called_once_with(True)
        # self.renderer.window.Render.assert_called_once_with()

    def test_screenshot_console(self):
        
        # Mock the render_model method
        # self.renderer.render_model = MagicMock(return_value=self.renderer.window)

        # Mock the vtkWindowToImageFilter and vtkPNGWriter
        # vtk.vtkWindowToImageFilter = MagicMock()
        # vtk.vtkPNGWriter = MagicMock()

        # Call the screenshot method
        result = self.renderer.screenshot(filepath='console', interact=False, reset_camera=True)

        # Assert that the render_model method was called
        # self.renderer.render_model.assert_called_once_with(False, True)

        # Assert that the result is an instance of IPython.display.Image
        self.assertIsInstance(result, Image)

    def test_screenshot_bytesio(self):

        # Call the screenshot method
        result = self.renderer.screenshot(filepath='BytesIO', interact=False, reset_camera=True)

        # Assert that the result is an instance of BytesIO
        self.assertIsInstance(result, BytesIO)

    def test_screenshot_file(self):

        # Call the screenshot method
        self.renderer.screenshot(filepath='test.png', interact=False, reset_camera=True)
