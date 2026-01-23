"""
Test coverage improvements for Material, Section, and LoadCombo classes.
"""

import unittest
import numpy as np
from Pynite import FEModel3D
from Pynite.Material import Material
from Pynite.Section import Section, SteelSection
from Pynite.LoadCombo import LoadCombo


class TestMaterial(unittest.TestCase):
    """Tests for Material class."""

    def setUp(self):
        self.model = FEModel3D()

    def test_material_creation(self):
        """Test basic material creation."""
        mat = Material(self.model, 'Steel', 29000, 11200, 0.3, 0.490)
        self.assertEqual(mat.name, 'Steel')
        self.assertEqual(mat.E, 29000)
        self.assertEqual(mat.G, 11200)
        self.assertEqual(mat.nu, 0.3)
        self.assertEqual(mat.rho, 0.490)
        self.assertIsNone(mat.fy)

    def test_material_with_yield_strength(self):
        """Test material creation with yield strength."""
        mat = Material(self.model, 'Steel', 29000, 11200, 0.3, 0.490, fy=50)
        self.assertEqual(mat.fy, 50)

    def test_material_model_reference(self):
        """Test that material holds reference to model."""
        mat = Material(self.model, 'Steel', 29000, 11200, 0.3, 0.490)
        self.assertIs(mat.model, self.model)


class TestLoadCombo(unittest.TestCase):
    """Tests for LoadCombo class."""

    def test_loadcombo_creation_empty(self):
        """Test creating an empty load combination."""
        combo = LoadCombo('D', factors={})
        self.assertEqual(combo.name, 'D')
        self.assertEqual(combo.factors, {})
        self.assertIsNone(combo.combo_tags)

    def test_loadcombo_creation_with_factors(self):
        """Test creating load combination with initial factors."""
        combo = LoadCombo('1.2D+1.6L', factors={'D': 1.2, 'L': 1.6})
        self.assertEqual(combo.factors['D'], 1.2)
        self.assertEqual(combo.factors['L'], 1.6)

    def test_loadcombo_with_tags(self):
        """Test load combination with tags."""
        combo = LoadCombo('Strength', combo_tags=['strength', 'ultimate'])
        self.assertEqual(combo.combo_tags, ['strength', 'ultimate'])

    def test_add_load_case(self):
        """Test adding load cases to combination."""
        combo = LoadCombo('My Combo')
        combo.AddLoadCase('D', 1.2)
        self.assertEqual(combo.factors['D'], 1.2)
        combo.AddLoadCase('L', 1.6)
        self.assertEqual(combo.factors['L'], 1.6)

    def test_delete_load_case(self):
        """Test deleting load cases from combination."""
        combo = LoadCombo('My Combo', factors={'D': 1.2, 'L': 1.6})
        combo.DeleteLoadCase('L')
        self.assertNotIn('L', combo.factors)
        self.assertIn('D', combo.factors)

    def test_overwrite_load_case(self):
        """Test overwriting an existing load case."""
        combo = LoadCombo('My Combo', factors={'D': 1.2})
        combo.AddLoadCase('D', 1.5)
        self.assertEqual(combo.factors['D'], 1.5)


class TestSection(unittest.TestCase):
    """Tests for Section class."""

    def setUp(self):
        self.model = FEModel3D()

    def test_section_creation(self):
        """Test basic section creation."""
        section = Section(self.model, 'W8x31', 9.13, 37.1, 110, 0.536)
        self.assertEqual(section.name, 'W8x31')
        self.assertEqual(section.A, 9.13)
        self.assertEqual(section.Iy, 37.1)
        self.assertEqual(section.Iz, 110)
        self.assertEqual(section.J, 0.536)

    def test_section_model_reference(self):
        """Test that section holds reference to model."""
        section = Section(self.model, 'W8x31', 9.13, 37.1, 110, 0.536)
        self.assertIs(section.model, self.model)

    def test_section_phi_not_implemented(self):
        """Test that base Section.Phi() raises NotImplementedError."""
        section = Section(self.model, 'W8x31', 9.13, 37.1, 110, 0.536)
        with self.assertRaises(NotImplementedError):
            section.Phi(100, 50, 50)

    def test_section_gradient_calculation(self):
        """Test gradient calculation for base Section class."""
        section = Section(self.model, 'W8x31', 9.13, 37.1, 110, 0.536)
        
        # Create a subclass that implements Phi for testing
        class TestSection(Section):
            def Phi(self, fx=0, my=0, mz=0):
                return (fx**2 + my**2 + mz**2) / 1000
        
        test_section = TestSection(self.model, 'Test', 9.13, 37.1, 110, 0.536)
        gradient = test_section.G(100, 50, 25)
        
        # Gradient should be a 3x1 array
        self.assertEqual(gradient.shape, (3, 1))
        # Values should be numerical approximations of derivatives
        self.assertGreater(gradient[0, 0], 0)  # Positive derivative for fx


class TestSteelSection(unittest.TestCase):
    """Tests for SteelSection class."""

    def setUp(self):
        self.model = FEModel3D()
        self.model.add_material('Steel', 29000, 11200, 0.3, 0.490, fy=50)

    def test_steel_section_creation(self):
        """Test steel section creation."""
        section = SteelSection(
            self.model, 'W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 'Steel'
        )
        self.assertEqual(section.name, 'W8x31')
        self.assertEqual(section.A, 9.13)
        self.assertEqual(section.Zy, 14.1)
        self.assertEqual(section.Zz, 30.4)

    def test_steel_section_radii_of_gyration(self):
        """Test that radii of gyration are calculated correctly."""
        section = SteelSection(
            self.model, 'W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 'Steel'
        )
        # ry = sqrt(Iy/A), rz = sqrt(Iz/A)
        expected_ry = np.sqrt(37.1 / 9.13)
        expected_rz = np.sqrt(110 / 9.13)
        self.assertAlmostEqual(section.ry, expected_ry, places=5)
        self.assertAlmostEqual(section.rz, expected_rz, places=5)

    def test_steel_section_phi_elastic(self):
        """Test Phi calculation when section is elastic."""
        section = SteelSection(
            self.model, 'W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 'Steel'
        )
        # Small loads should result in Phi < 1
        phi = section.Phi(fx=10, my=5, mz=5)
        self.assertLess(phi, 1.0)

    def test_steel_section_phi_plastic(self):
        """Test Phi calculation when section is plastic."""
        section = SteelSection(
            self.model, 'W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 'Steel'
        )
        # Large loads should result in Phi > 1
        phi = section.Phi(fx=200, my=500, mz=1000)
        self.assertGreater(phi, 1.0)

    def test_steel_section_gradient_elastic(self):
        """Test gradient calculation when section is elastic."""
        section = SteelSection(
            self.model, 'W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 'Steel'
        )
        # Small loads should result in zero gradient
        gradient = section.G(fx=10, my=5, mz=5)
        self.assertEqual(gradient.shape, (6, 1))
        np.testing.assert_array_almost_equal(gradient, np.zeros((6, 1)))

    def test_steel_section_gradient_plastic(self):
        """Test gradient calculation when section is plastic."""
        section = SteelSection(
            self.model, 'W8x31', 9.13, 37.1, 110, 0.536, 14.1, 30.4, 'Steel'
        )
        # Large loads should result in non-zero gradient
        gradient = section.G(fx=200, my=500, mz=1000)
        self.assertEqual(gradient.shape, (6, 1))
        # At least some components should be non-zero
        self.assertGreater(np.abs(gradient).max(), 0)


if __name__ == '__main__':
    unittest.main()
