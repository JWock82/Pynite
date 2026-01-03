"""
Tests for mesh regeneration behavior.

These tests verify that when a mesh is regenerated:
1. Old nodes and elements are removed from the model
2. Nodes shared with elements outside the mesh are preserved
3. The mesh can be regenerated multiple times with different properties
"""

import pytest
from Pynite import FEModel3D


def test_rectangle_mesh_regeneration():
    """Test that regenerating a RectangleMesh removes old nodes/elements."""
    
    # Create a simple model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 4, 3, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    
    # Generate the mesh
    mesh.generate()
    initial_node_count = len(model.nodes)
    initial_element_count = len(model.quads)
    
    # Store the initial node and element names
    initial_nodes = set(model.nodes.keys())
    initial_elements = set(model.quads.keys())
    
    # Regenerate the mesh with the same parameters
    mesh.generate()
    
    # Check that node and element counts haven't changed
    assert len(model.nodes) == initial_node_count, \
        f"Node count changed after regeneration: {initial_node_count} -> {len(model.nodes)}"
    assert len(model.quads) == initial_element_count, \
        f"Element count changed after regeneration: {initial_element_count} -> {len(model.quads)}"
    
    # Verify the same nodes and elements exist
    assert set(model.nodes.keys()) == initial_nodes, \
        "Node names changed after regeneration"
    assert set(model.quads.keys()) == initial_elements, \
        "Element names changed after regeneration"


def test_mesh_regeneration_with_shared_nodes():
    """Test that nodes shared with other elements are preserved during regeneration."""
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    # Get a node from the mesh that we'll attach to an external element
    shared_node = list(mesh.nodes.values())[0]
    
    # Create a new node outside the mesh
    external_node_name = model.add_node('ExtNode', 10, 10, 10)
    
    # Add a quad element that shares a node with the mesh
    # First, create three more nodes for the quad
    n2 = model.add_node('ExtNode2', 11, 10, 10)
    n3 = model.add_node('ExtNode3', 11, 11, 10)
    n4 = model.add_node('ExtNode4', 10, 11, 10)
    
    # Add the external quad with one shared node
    model.add_quad('ExtQuad', shared_node.name, n2, n3, n4, 0.25, 'Steel')
    
    initial_nodes = set(model.nodes.keys())
    
    # Regenerate the mesh
    mesh.generate()
    
    # The shared node should still exist
    assert shared_node.name in model.nodes, \
        f"Shared node {shared_node.name} was incorrectly removed during regeneration"
    
    # The external quad should still exist
    assert 'ExtQuad' in model.quads, \
        "External quad was removed during mesh regeneration"
    
    # The external nodes should still exist
    assert external_node_name in model.nodes
    assert n2 in model.nodes
    assert n3 in model.nodes
    assert n4 in model.nodes


def test_mesh_parameter_change_regeneration():
    """Test that a mesh can be regenerated after changing parameters."""
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 4, 3, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    initial_element_count = len(model.quads)
    
    # Change the mesh size to create more elements
    mesh.mesh_size = 0.5
    
    # Regenerate the mesh
    mesh.generate()
    
    # Should have more elements with a smaller mesh size
    assert len(model.quads) > initial_element_count, \
        "Mesh did not regenerate with new parameters"


def test_cylinder_mesh_regeneration():
    """Test that CylinderMesh regeneration works correctly."""
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a cylinder mesh
    mesh_name = model.add_cylinder_mesh('Cyl1', 1, 5, 10, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    initial_node_count = len(model.nodes)
    initial_element_count = len(model.quads)
    
    # Regenerate the mesh
    mesh.generate()
    
    # Check that counts haven't changed
    assert len(model.nodes) == initial_node_count
    assert len(model.quads) == initial_element_count


def test_annulus_mesh_regeneration():
    """Test that AnnulusMesh regeneration works correctly."""
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add an annulus mesh
    mesh_name = model.add_annulus_mesh('Annulus1', 1, 10, 5, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    initial_node_count = len(model.nodes)
    initial_element_count = len(model.quads)
    
    # Regenerate the mesh
    mesh.generate()
    
    # Check that counts haven't changed
    assert len(model.nodes) == initial_node_count
    assert len(model.quads) == initial_element_count


def test_multiple_meshes_independent_regeneration():
    """Test that regenerating one mesh doesn't affect other meshes."""
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add two separate meshes
    mesh1_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel', origin=[0, 0, 0])
    mesh2_name = model.add_rectangle_mesh('Mesh2', 1, 2, 2, 0.25, 'Steel', origin=[10, 0, 0])
    
    mesh1 = model.meshes[mesh1_name]
    mesh2 = model.meshes[mesh2_name]
    
    mesh1.generate()
    mesh2.generate()
    
    mesh2_nodes = set(mesh2.nodes.keys())
    mesh2_elements = set(mesh2.elements.keys())
    
    # Regenerate mesh1
    mesh1.generate()
    
    # Mesh2's nodes and elements should be unchanged
    assert set(mesh2.nodes.keys()) == mesh2_nodes, \
        "Mesh2 nodes were affected by Mesh1 regeneration"
    assert set(mesh2.elements.keys()) == mesh2_elements, \
        "Mesh2 elements were affected by Mesh1 regeneration"
    
    # All mesh2 nodes should still be in the model
    for node_name in mesh2_nodes:
        assert node_name in model.nodes, \
            f"Node {node_name} from Mesh2 was removed when Mesh1 was regenerated"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
