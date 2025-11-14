"""
Pytest tests for FEModel3D.merge_duplicate_nodes.

This module focuses on correctness of:
- No-op when nodes are spaced beyond tolerance
- Pair merge behavior and canonical node selection
- Rewiring of elements (springs) to the canonical node
- Merging of support flags and spring supports

Note: File name is `node_merge.py` per request. A lightweight wrapper
`test_node_merge.py` re-exports these tests for pytest auto-discovery.
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

    removed = model.merge_duplicate_nodes(tolerance=0.001)

    # Expect B to be merged into A (A is canonical as first encountered)
    assert set(removed) == {"B"}
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

    removed = model.merge_duplicate_nodes(tolerance=0.001)

    # Only N1 should remain
    assert set(model.nodes.keys()) == {"N1"}
    # All others should be removed
    assert set(removed) == {"N2", "N3", "N4"}

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

    removed = model.merge_duplicate_nodes(tolerance=0.001)
    assert set(removed) == {"G"}

    f = model.nodes["F"]
    assert f.support_DX is True and f.support_DZ is True
    assert f.spring_DY[0] == 5000.0


if __name__ == "__main__":
    import sys
    import pytest
    # Allow running this file directly for convenience
    sys.exit(pytest.main([__file__, "-q"]))

