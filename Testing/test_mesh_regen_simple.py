"""
Simple test for mesh regeneration behavior (no pytest required).
"""

import sys
sys.path.insert(0, r'c:\Users\craig\Documents\Python\Pynite')

from Pynite import FEModel3D


def test_rectangle_mesh_regeneration():
    """Test that regenerating a RectangleMesh removes old nodes/elements."""
    
    print("\n=== Test 1: Rectangle Mesh Regeneration ===")
    
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
    
    print(f"Initial: {initial_node_count} nodes, {initial_element_count} elements")
    
    # Store the initial node and element names
    initial_nodes = set(model.nodes.keys())
    initial_elements = set(model.quads.keys())
    
    # Regenerate the mesh with the same parameters
    mesh.generate()
    
    final_node_count = len(model.nodes)
    final_element_count = len(model.quads)
    
    print(f"After regeneration: {final_node_count} nodes, {final_element_count} elements")
    
    # Check that node and element counts haven't changed
    assert final_node_count == initial_node_count, \
        f"Node count changed: {initial_node_count} -> {final_node_count}"
    assert final_element_count == initial_element_count, \
        f"Element count changed: {initial_element_count} -> {final_element_count}"
    
    # Verify the same nodes and elements exist
    assert set(model.nodes.keys()) == initial_nodes, "Node names changed"
    assert set(model.quads.keys()) == initial_elements, "Element names changed"
    
    print("✓ PASSED: Mesh regenerated correctly with no duplicates")


def test_mesh_regeneration_with_shared_nodes():
    """Test that nodes shared with other elements are preserved during regeneration."""
    
    print("\n=== Test 2: Shared Nodes Preservation ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    print(f"Mesh generated with {len(mesh.nodes)} nodes")
    
    # Get a node from the mesh that we'll attach to an external element
    shared_node = list(mesh.nodes.values())[0]
    print(f"Shared node: {shared_node.name}")
    
    # Create nodes for an external quad
    n2 = model.add_node('ExtNode2', 10, 10, 10)
    n3 = model.add_node('ExtNode3', 11, 10, 10)
    n4 = model.add_node('ExtNode4', 10, 11, 10)
    
    # Add the external quad with one shared node
    model.add_quad('ExtQuad', shared_node.name, n2, n3, n4, 0.25, 'Steel')
    print(f"Added external quad sharing node {shared_node.name}")
    
    nodes_before = len(model.nodes)
    
    # Regenerate the mesh
    print("Regenerating mesh...")
    mesh.generate()
    
    nodes_after = len(model.nodes)
    print(f"Nodes before: {nodes_before}, after: {nodes_after}")
    
    # The shared node should still exist
    assert shared_node.name in model.nodes, \
        f"Shared node {shared_node.name} was incorrectly removed"
    
    # The external quad should still exist
    assert 'ExtQuad' in model.quads, "External quad was removed"
    
    # The external nodes should still exist
    assert n2 in model.nodes
    assert n3 in model.nodes
    assert n4 in model.nodes
    
    print(f"✓ PASSED: Shared node {shared_node.name} preserved correctly")


def test_mesh_parameter_change_regeneration():
    """Test that a mesh can be regenerated after changing parameters."""
    
    print("\n=== Test 3: Parameter Change Regeneration ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 4, 3, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    initial_element_count = len(model.quads)
    print(f"Initial mesh size = 1, elements = {initial_element_count}")
    
    # Change the mesh size to create more elements
    mesh.mesh_size = 0.5
    print("Changed mesh size to 0.5")
    
    # Regenerate the mesh
    mesh.generate()
    
    final_element_count = len(model.quads)
    print(f"After regeneration: elements = {final_element_count}")
    
    # Should have more elements with a smaller mesh size
    assert final_element_count > initial_element_count, \
        f"Mesh did not regenerate properly: {initial_element_count} -> {final_element_count}"
    
    print(f"✓ PASSED: Mesh regenerated with new parameters ({initial_element_count} -> {final_element_count} elements)")


def test_multiple_meshes_independent_regeneration():
    """Test that regenerating one mesh doesn't affect other meshes."""
    
    print("\n=== Test 4: Independent Mesh Regeneration ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add two separate meshes at different locations
    mesh1_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel', origin=[0, 0, 0])
    mesh2_name = model.add_rectangle_mesh('Mesh2', 1, 2, 2, 0.25, 'Steel', origin=[10, 0, 0])
    
    mesh1 = model.meshes[mesh1_name]
    mesh2 = model.meshes[mesh2_name]
    
    mesh1.generate()
    mesh2.generate()
    
    mesh2_node_count = len(mesh2.nodes)
    mesh2_element_count = len(mesh2.elements)
    mesh2_nodes = set(mesh2.nodes.keys())
    
    print(f"Mesh1 nodes: {len(mesh1.nodes)}, Mesh2 nodes: {mesh2_node_count}")
    
    # Regenerate mesh1
    print("Regenerating Mesh1...")
    mesh1.generate()
    
    # Mesh2's nodes and elements should be unchanged
    assert len(mesh2.nodes) == mesh2_node_count, \
        f"Mesh2 node count changed: {mesh2_node_count} -> {len(mesh2.nodes)}"
    assert set(mesh2.nodes.keys()) == mesh2_nodes, "Mesh2 node names changed"
    
    # All mesh2 nodes should still be in the model
    for node_name in mesh2_nodes:
        assert node_name in model.nodes, \
            f"Node {node_name} from Mesh2 was removed"
    
    print("✓ PASSED: Mesh2 unaffected by Mesh1 regeneration")


if __name__ == '__main__':
    print("=" * 60)
    print("Testing Mesh Regeneration Behavior")
    print("=" * 60)
    
    try:
        test_rectangle_mesh_regeneration()
        test_mesh_regeneration_with_shared_nodes()
        test_mesh_parameter_change_regeneration()
        test_multiple_meshes_independent_regeneration()
        
        print("\n" + "=" * 60)
        print("ALL TESTS PASSED ✓")
        print("=" * 60)
        
    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
