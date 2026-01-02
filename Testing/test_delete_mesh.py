"""
Tests for mesh deletion functionality.

These tests verify that the delete_mesh() method:
1. Removes the mesh from the model
2. Removes the mesh's elements from the model
3. Removes the mesh's nodes that aren't shared with other elements
4. Preserves nodes and elements that are attached to the mesh
5. Handles errors appropriately
"""

import sys
sys.path.insert(0, r'c:\Users\craig\Documents\Python\Pynite')

from Pynite import FEModel3D


def test_delete_mesh_basic():
    """Test basic mesh deletion."""
    
    print("\n=== Test 1: Basic Mesh Deletion ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 4, 3, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    nodes_before = len(model.nodes)
    elements_before = len(model.quads)
    
    print(f"Before deletion: {nodes_before} nodes, {elements_before} quads")
    print(f"Meshes in model: {list(model.meshes.keys())}")
    
    # Delete the mesh
    model.delete_mesh('Mesh1')
    
    nodes_after = len(model.nodes)
    elements_after = len(model.quads)
    
    print(f"After deletion: {nodes_after} nodes, {elements_after} quads")
    print(f"Meshes in model: {list(model.meshes.keys())}")
    
    # Verify the mesh is gone
    assert 'Mesh1' not in model.meshes, "Mesh was not removed from model.meshes"
    
    # Verify all nodes and elements are gone
    assert nodes_after == 0, f"Nodes not fully deleted: {nodes_after} remaining"
    assert elements_after == 0, f"Elements not fully deleted: {elements_after} remaining"
    
    print("✓ PASSED: Mesh deleted successfully")


def test_delete_mesh_with_attached_elements():
    """Test that attached elements are preserved during mesh deletion."""
    
    print("\n=== Test 2: Mesh Deletion with Attached Elements ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    # Get a node from the mesh to attach an external element
    shared_node = list(mesh.nodes.values())[0]
    shared_node_name = shared_node.name
    
    print(f"Shared node: {shared_node_name}")
    
    # Create external nodes and attach a quad to the mesh's node
    n2 = model.add_node('ExtNode2', 10, 10, 10)
    n3 = model.add_node('ExtNode3', 11, 10, 10)
    n4 = model.add_node('ExtNode4', 10, 11, 10)
    
    # Add external quad sharing a node with the mesh
    model.add_quad('ExtQuad', shared_node_name, n2, n3, n4, 0.25, 'Steel')
    
    nodes_before = len(model.nodes)
    quads_before = len(model.quads)
    
    print(f"Before deletion: {nodes_before} nodes, {quads_before} quads")
    
    # Delete the mesh
    model.delete_mesh('Mesh1')
    
    nodes_after = len(model.nodes)
    quads_after = len(model.quads)
    
    print(f"After deletion: {nodes_after} nodes, {quads_after} quads")
    
    # Verify the shared node is still there
    assert shared_node_name in model.nodes, \
        f"Shared node {shared_node_name} was deleted"
    
    # Verify the external quad is still there
    assert 'ExtQuad' in model.quads, \
        "External quad was deleted"
    
    # Verify we still have the external nodes plus the shared node
    assert nodes_after == 4, f"Expected 4 nodes (1 shared + 3 external), got {nodes_after}"
    assert quads_after == 1, f"Expected 1 quad (external), got {quads_after}"
    
    print(f"✓ PASSED: Attached elements preserved ({shared_node_name}, ExtQuad)")


def test_delete_mesh_with_attached_members():
    """Test that member elements attached to mesh nodes are preserved."""
    
    print("\n=== Test 3: Mesh Deletion with Attached Members ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    model.add_section('Section1', 1.0, 1.0, 1.0, 0.1)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    # Get a node from the mesh to attach a member
    mesh_node = list(mesh.nodes.values())[0]
    mesh_node_name = mesh_node.name
    
    # Create an external node and a member connecting to the mesh node
    external_node = model.add_node('ExtNode', 10, 10, 10)
    
    # Add a member from external node to mesh node
    model.add_member('Member1', mesh_node_name, external_node, 'Steel', 'Section1')
    
    members_before = len(model.members)
    
    print(f"Before deletion: {members_before} members")
    
    # Delete the mesh
    model.delete_mesh('Mesh1')
    
    members_after = len(model.members)
    
    print(f"After deletion: {members_after} members")
    
    # Verify the member is still there
    assert 'Member1' in model.members, "Member was deleted"
    
    # Verify the mesh node is still there
    assert mesh_node_name in model.nodes, \
        f"Mesh node {mesh_node_name} was deleted (needed by member)"
    
    # Verify the external node is still there
    assert external_node in model.nodes, "External node was deleted"
    
    print(f"✓ PASSED: Member and mesh node preserved")


def test_delete_mesh_with_attached_springs():
    """Test that spring elements attached to mesh nodes are preserved."""
    
    print("\n=== Test 3b: Mesh Deletion with Attached Springs ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a rectangle mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel')
    mesh = model.meshes[mesh_name]
    mesh.generate()
    
    # Get a node from the mesh to attach a spring
    mesh_node = list(mesh.nodes.values())[0]
    mesh_node_name = mesh_node.name
    
    # Create an external node and a spring connecting to the mesh node
    external_node = model.add_node('ExtNode', 10, 10, 10)
    
    # Add a spring from external node to mesh node
    model.add_spring('Spring1', mesh_node_name, external_node, 1000.0)
    
    springs_before = len(model.springs)
    
    print(f"Before deletion: {springs_before} springs")
    
    # Delete the mesh
    model.delete_mesh('Mesh1')
    
    springs_after = len(model.springs)
    
    print(f"After deletion: {springs_after} springs")
    
    # Verify the spring is still there
    assert 'Spring1' in model.springs, "Spring was deleted"
    
    # Verify the mesh node is still there
    assert mesh_node_name in model.nodes, \
        f"Mesh node {mesh_node_name} was deleted (needed by spring)"
    
    # Verify the external node is still there
    assert external_node in model.nodes, "External node was deleted"
    
    print(f"✓ PASSED: Spring and mesh node preserved")



def test_delete_nonexistent_mesh():
    """Test that deleting a nonexistent mesh raises an error."""
    
    print("\n=== Test 4: Delete Nonexistent Mesh ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Try to delete a mesh that doesn't exist
    try:
        model.delete_mesh('NonexistentMesh')
        assert False, "Should have raised KeyError"
    except KeyError as e:
        print(f"✓ PASSED: Correctly raised KeyError: {e}")


def test_delete_multiple_meshes():
    """Test deleting multiple meshes independently."""
    
    print("\n=== Test 5: Delete Multiple Meshes ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add two meshes at different locations
    mesh1_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel', origin=[0, 0, 0])
    mesh2_name = model.add_rectangle_mesh('Mesh2', 1, 2, 2, 0.25, 'Steel', origin=[10, 0, 0])
    
    model.meshes[mesh1_name].generate()
    model.meshes[mesh2_name].generate()
    
    mesh1_quad_count = len(model.meshes[mesh1_name].elements)
    mesh2_quad_count = len(model.meshes[mesh2_name].elements)
    
    print(f"Before deletion: {list(model.meshes.keys())}")
    print(f"Total nodes: {len(model.nodes)}, Total quads: {len(model.quads)}")
    print(f"Mesh1 has {mesh1_quad_count} quads, Mesh2 has {mesh2_quad_count} quads")
    
    # Delete mesh1
    model.delete_mesh('Mesh1')
    
    mesh_names_after_1 = list(model.meshes.keys())
    nodes_after_1 = len(model.nodes)
    quads_after_1 = len(model.quads)
    
    print(f"After deleting Mesh1: {mesh_names_after_1}")
    print(f"Nodes: {nodes_after_1}, Quads: {quads_after_1}")
    
    # Verify Mesh1 is gone but Mesh2 is still there
    assert 'Mesh1' not in model.meshes
    assert 'Mesh2' in model.meshes
    assert quads_after_1 == mesh2_quad_count, f"Mesh2's quads should remain ({mesh2_quad_count}), got {quads_after_1}"
    
    # Delete mesh2
    model.delete_mesh('Mesh2')
    
    mesh_names_final = list(model.meshes.keys())
    
    print(f"After deleting Mesh2: {mesh_names_final}")
    
    # Verify both are gone
    assert len(model.meshes) == 0
    
    print("✓ PASSED: Multiple meshes deleted independently")


def test_delete_mesh_flags_unsolved():
    """Test that deleting a mesh flags the model as unsolved."""
    
    print("\n=== Test 6: Delete Mesh Flags Model as Unsolved ===")
    
    # Create a model
    model = FEModel3D()
    model.add_material('Steel', 29000, 11200, 0.3, 490)
    
    # Add a mesh
    mesh_name = model.add_rectangle_mesh('Mesh1', 1, 2, 2, 0.25, 'Steel')
    model.meshes[mesh_name].generate()
    
    # Simulate a solved state
    model.solution = 'Linear'
    print(f"Model solution state before deletion: {model.solution}")
    
    # Delete the mesh
    model.delete_mesh('Mesh1')
    
    print(f"Model solution state after deletion: {model.solution}")
    
    # Verify the model is flagged as unsolved
    assert model.solution is None, \
        f"Model should be unsolved, but solution = {model.solution}"
    
    print("✓ PASSED: Model flagged as unsolved")


if __name__ == '__main__':
    print("=" * 60)
    print("Testing Mesh Deletion Functionality")
    print("=" * 60)
    
    try:
        test_delete_mesh_basic()
        test_delete_mesh_with_attached_elements()
        test_delete_mesh_with_attached_members()
        test_delete_mesh_with_attached_springs()
        test_delete_nonexistent_mesh()
        test_delete_multiple_meshes()
        test_delete_mesh_flags_unsolved()
        
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
