"""
Pytest tests for FEModel3D.merge_duplicate_nodes.

This module focuses on correctness of:
- No-op when nodes are spaced beyond tolerance
- Pair merge behavior and canonical node selection
- Rewiring of elements (springs) to the canonical node
- Merging of support flags and spring supports
- Return value format: list of (deleted_node, merged_into_node) tuples
"""

from __future__ import annotations

import pytest

from Pynite.FEModel3D import FEModel3D


def test_no_merge_unique_nodes():
    model = FEModel3D()
    model.add_node("N1", 0.0, 0.0, 0.0)
    model.add_node("N2", 1.0, 0.0, 0.0)

    removed = model.merge_duplicate_nodes(tolerance=1e-6)

    assert removed == []
    assert len(model.nodes) == 2
    assert set(model.nodes.keys()) == {"N1", "N2"}


def test_simple_pair_merge_canonical_first_added():
    model = FEModel3D()
    # Add A first, then B very close to A (within tolerance)
    model.add_node("A", 0.0, 0.0, 0.0)
    model.add_node("B", 0.0, 0.0, 0.0005)

    merge_mappings = model.merge_duplicate_nodes(tolerance=0.001)

    # Expect B to be merged into A (A is canonical as first encountered)
    assert len(merge_mappings) == 1
    assert merge_mappings[0] == ("B", "A")
    assert len(model.nodes) == 1
    assert "A" in model.nodes


def test_pair_not_merged_when_beyond_tolerance():
    model = FEModel3D()
    model.add_node("A", 0.0, 0.0, 0.0)
    # Distance is slightly greater than tolerance
    model.add_node("B", 0.0, 0.0, 0.0011)

    removed = model.merge_duplicate_nodes(tolerance=0.001)

    assert removed == []
    assert set(model.nodes.keys()) == {"A", "B"}


def test_cluster_merge_and_element_rewire():
    model = FEModel3D()
    # Tight cluster around origin; all should collapse to first added node
    n1 = model.add_node("N1", 0.0, 0.0, 0.0)
    n2 = model.add_node("N2", 0.0004, 0.0, 0.0)
    n3 = model.add_node("N3", 0.0, 0.0004, 0.0)
    n4 = model.add_node("N4", 0.0, 0.0, 0.0004)

    # Add a couple of springs to ensure element node references are updated
    model.add_spring("S1", n2, n3, ks=1000.0)
    model.add_spring("S2", n3, n4, ks=2000.0)

    merge_mappings = model.merge_duplicate_nodes(tolerance=0.001)

    # Only N1 should remain
    assert set(model.nodes.keys()) == {"N1"}

    # All others should be merged into N1
    deleted_nodes = {m[0] for m in merge_mappings}
    target_nodes = {m[1] for m in merge_mappings}
    assert deleted_nodes == {"N2", "N3", "N4"}
    assert target_nodes == {"N1"}  # All merged into N1

    # Springs should be rewired to N1 at both ends
    s1 = model.springs["S1"]
    s2 = model.springs["S2"]
    assert s1.i_node.name == "N1" and s1.j_node.name == "N1"
    assert s2.i_node.name == "N1" and s2.j_node.name == "N1"


def test_support_and_spring_support_merge():
    model = FEModel3D()
    # F will be canonical; G will be merged into F
    model.add_node("F", 0.0, 0.0, 0.0)
    model.add_node("G", 0.0, 0.0, 0.0002)

    # Apply support/spring to G, which should transfer to F
    model.def_support("G", support_DX=True, support_DZ=True)
    model.def_support_spring("G", dof="DY", stiffness=5000.0)

    merge_mappings = model.merge_duplicate_nodes(tolerance=0.001)
    assert merge_mappings == [("G", "F")]

    f = model.nodes["F"]
    assert f.support_DX is True and f.support_DZ is True
    assert f.spring_DY[0] == 5000.0


def test_merge_mapping_returns_correct_tuple_format():
    """Verify that merge_duplicate_nodes returns list of (deleted, merged_into) tuples."""
    model = FEModel3D()
    model.add_node("X", 0.0, 0.0, 0.0)
    model.add_node("Y", 0.0, 0.0, 0.0001)

    merge_mappings = model.merge_duplicate_nodes(tolerance=0.001)

    # Should return a list of tuples
    assert isinstance(merge_mappings, list)
    assert len(merge_mappings) == 1
    assert isinstance(merge_mappings[0], tuple)
    assert len(merge_mappings[0]) == 2

    # First element is deleted node, second is the node it was merged into
    deleted_node, merged_into_node = merge_mappings[0]
    assert deleted_node == "Y"
    assert merged_into_node == "X"


def test_merge_mapping_multiple_separate_clusters():
    """Test merge mappings when there are multiple separate clusters of duplicates."""
    model = FEModel3D()
    # Cluster 1: A and B at origin
    model.add_node("A", 0.0, 0.0, 0.0)
    model.add_node("B", 0.0, 0.0, 0.0001)
    # Cluster 2: C and D far from origin
    model.add_node("C", 100.0, 0.0, 0.0)
    model.add_node("D", 100.0, 0.0, 0.0001)

    merge_mappings = model.merge_duplicate_nodes(tolerance=0.001)

    # Should have two merge operations
    assert len(merge_mappings) == 2

    # Convert to dict for easier checking
    merge_dict = dict(merge_mappings)
    assert merge_dict["B"] == "A"
    assert merge_dict["D"] == "C"

    # Only A and C should remain
    assert set(model.nodes.keys()) == {"A", "C"}


def test_merge_mapping_chain_collapses_to_first():
    """When multiple nodes are close, they should all merge into the first one."""
    model = FEModel3D()
    model.add_node("First", 0.0, 0.0, 0.0)
    model.add_node("Second", 0.0001, 0.0, 0.0)
    model.add_node("Third", 0.0002, 0.0, 0.0)

    merge_mappings = model.merge_duplicate_nodes(tolerance=0.001)

    # All should merge into "First"
    assert len(merge_mappings) == 2
    for deleted, merged_into in merge_mappings:
        assert merged_into == "First"

    deleted_nodes = [m[0] for m in merge_mappings]
    assert set(deleted_nodes) == {"Second", "Third"}


def test_merge_mapping_can_reconstruct_original_references():
    """Verify that merge mappings can be used to track where deleted nodes went."""
    model = FEModel3D()
    model.add_node("Keep", 0.0, 0.0, 0.0)
    model.add_node("Delete1", 0.0, 0.0, 0.0001)
    model.add_node("Delete2", 0.0, 0.0, 0.0002)

    merge_mappings = model.merge_duplicate_nodes(tolerance=0.001)

    # Build a lookup dictionary from deleted nodes to their new locations
    redirect_map = {deleted: target for deleted, target in merge_mappings}

    # Verify we can use this map to redirect references
    assert redirect_map.get("Delete1") == "Keep"
    assert redirect_map.get("Delete2") == "Keep"
    assert redirect_map.get("Keep") is None  # Keep was not deleted


if __name__ == "__main__":
    import sys
    import pytest
    # Allow running this file directly for convenience
    sys.exit(pytest.main([__file__, "-q"]))

