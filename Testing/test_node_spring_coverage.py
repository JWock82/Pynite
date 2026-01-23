"""
Test coverage improvements for Node3D and Spring3D classes.
"""

import unittest
from Pynite import FEModel3D


class TestNode3D(unittest.TestCase):
    """Tests for Node3D class."""

    def setUp(self):
        self.model = FEModel3D()
        self.model.add_node('N1', 0, 0, 0)
        self.node = self.model.nodes['N1']

    def test_node_creation(self):
        """Test basic node creation."""
        self.assertEqual(self.node.name, 'N1')
        self.assertEqual(self.node.X, 0)
        self.assertEqual(self.node.Y, 0)
        self.assertEqual(self.node.Z, 0)

    def test_node_coordinates(self):
        """Test node with different coordinates."""
        self.model.add_node('N2', 10, 20, 30)
        node = self.model.nodes['N2']
        self.assertEqual(node.X, 10)
        self.assertEqual(node.Y, 20)
        self.assertEqual(node.Z, 30)

    def test_node_model_reference(self):
        """Test that node holds reference to model."""
        self.assertIs(self.node.model, self.model)

    def test_node_supports_initialized_false(self):
        """Test that all supports start as False."""
        self.assertFalse(self.node.support_DX)
        self.assertFalse(self.node.support_DY)
        self.assertFalse(self.node.support_DZ)
        self.assertFalse(self.node.support_RX)
        self.assertFalse(self.node.support_RY)
        self.assertFalse(self.node.support_RZ)

    def test_node_springs_initialized_none(self):
        """Test that all springs start as None."""
        self.assertIsNone(self.node.spring_DX[0])
        self.assertIsNone(self.node.spring_DY[0])
        self.assertIsNone(self.node.spring_DZ[0])
        self.assertIsNone(self.node.spring_RX[0])
        self.assertIsNone(self.node.spring_RY[0])
        self.assertIsNone(self.node.spring_RZ[0])

    def test_node_enforced_displacements_none(self):
        """Test that enforced displacements start as None."""
        self.assertIsNone(self.node.EnforcedDX)
        self.assertIsNone(self.node.EnforcedDY)
        self.assertIsNone(self.node.EnforcedDZ)
        self.assertIsNone(self.node.EnforcedRX)
        self.assertIsNone(self.node.EnforcedRY)
        self.assertIsNone(self.node.EnforcedRZ)

    def test_node_distance_zero(self):
        """Test distance calculation between same node."""
        distance = self.node.distance(self.node)
        self.assertEqual(distance, 0)

    def test_node_distance_simple(self):
        """Test distance calculation between two nodes."""
        self.model.add_node('N2', 3, 4, 0)
        node2 = self.model.nodes['N2']
        distance = self.node.distance(node2)
        self.assertEqual(distance, 5)  # 3-4-5 triangle

    def test_node_distance_3d(self):
        """Test distance calculation in 3D."""
        self.model.add_node('N2', 1, 2, 2)
        node2 = self.model.nodes['N2']
        distance = self.node.distance(node2)
        expected = (1**2 + 2**2 + 2**2)**0.5
        self.assertAlmostEqual(distance, expected, places=5)

    def test_node_load_lists_initialized(self):
        """Test that node load lists are initialized."""
        self.assertEqual(self.node.NodeLoads, [])

    def test_node_displacement_dicts_initialized(self):
        """Test that displacement dictionaries are initialized empty."""
        self.assertEqual(self.node.DX, {})
        self.assertEqual(self.node.DY, {})
        self.assertEqual(self.node.DZ, {})
        self.assertEqual(self.node.RX, {})
        self.assertEqual(self.node.RY, {})
        self.assertEqual(self.node.RZ, {})

    def test_node_reaction_dicts_initialized(self):
        """Test that reaction dictionaries are initialized empty."""
        self.assertEqual(self.node.RxnFX, {})
        self.assertEqual(self.node.RxnFY, {})
        self.assertEqual(self.node.RxnFZ, {})
        self.assertEqual(self.node.RxnMX, {})
        self.assertEqual(self.node.RxnMY, {})
        self.assertEqual(self.node.RxnMZ, {})

    def test_node_contour_list(self):
        """Test that contour list is initialized."""
        self.assertEqual(self.node.contour, [])


class TestSpring3D(unittest.TestCase):
    """Tests for Spring3D class."""

    def setUp(self):
        self.model = FEModel3D()
        self.model.add_node('N1', 0, 0, 0)
        self.model.add_node('N2', 10, 0, 0)
        self.node_i = self.model.nodes['N1']
        self.node_j = self.model.nodes['N2']

    def test_spring_creation(self):
        """Test basic spring creation."""
        self.model.add_spring('S1', 'N1', 'N2', 1000)
        spring = self.model.springs['S1']
        self.assertEqual(spring.name, 'S1')
        self.assertEqual(spring.ks, 1000)
        self.assertIs(spring.i_node, self.node_i)
        self.assertIs(spring.j_node, self.node_j)

    def test_two_way_spring(self):
        """Test two-way spring (default)."""
        self.model.add_spring('S1', 'N1', 'N2', 1000, tension_only=False, comp_only=False)
        spring = self.model.springs['S1']
        self.assertFalse(spring.tension_only)
        self.assertFalse(spring.comp_only)

    def test_tension_only_spring(self):
        """Test tension-only spring."""
        self.model.add_spring('S1', 'N1', 'N2', 1000, tension_only=True)
        spring = self.model.springs['S1']
        self.assertTrue(spring.tension_only)

    def test_compression_only_spring(self):
        """Test compression-only spring."""
        self.model.add_spring('S1', 'N1', 'N2', 1000, comp_only=True)
        spring = self.model.springs['S1']
        self.assertTrue(spring.comp_only)

    def test_spring_length(self):
        """Test spring length calculation."""
        self.model.add_node('N3', 3, 4, 0)
        self.model.add_node('N4', 0, 0, 0)
        self.model.add_spring('S2', 'N4', 'N3', 1000)
        spring = self.model.springs['S2']
        # Distance should be 5 (3-4-5 triangle)
        self.assertEqual(spring.L(), 5)

    def test_spring_active_dict(self):
        """Test spring active dictionary initialization."""
        self.model.add_spring('S1', 'N1', 'N2', 1000)
        spring = self.model.springs['S1']
        # Active dict should be initialized
        self.assertEqual(spring.active, {})


if __name__ == '__main__':
    unittest.main()
