"""Tests for reusing one factorization of K11 across load combinations in `analyze_linear`.

`analyze_linear` assembles the stiffness matrix once. With the sparse solver it now also factorizes
the partitioned matrix K11 once (`Analysis._factorize_K11`) and reuses that factorization for every
load combination, instead of calling `spsolve` (a new factorization) for each one. These tests
check that the results are the same as before, that the factorization is done only once, and that
unstable models are still reported.
"""

# Standard library imports
import math
from unittest.mock import patch

# Third party imports
import pytest

# Pynite imports
from Pynite import FEModel3D
from Pynite import Analysis


def _plate_wall():
    """A 10 ft x 15 ft concrete wall of quads, fixed on three sides, with two load cases and three
    load combinations (the same wall as `test_hydrostatic_quad`, plus a uniform pressure case)."""

    model = FEModel3D()

    # Geometry (ft) and concrete properties (psf)
    t, mesh_size, a, b = 1, 1, 10, 15
    E = 57000*math.sqrt(4500)*12**2
    model.add_material('Concrete', E, 0.4*E, 1/6, 150)

    # Mesh of quads
    model.add_rectangle_mesh('MSH1', mesh_size, a, b, t, 'Concrete', element_type='Quad')
    model.meshes['MSH1'].generate()

    # Fix the sides and the base of the wall
    for node in model.nodes.values():
        if node.X == 0 or node.X == a or node.Y == 0:
            model.def_support(node.name, True, True, True, True, True, True)

    # Hydrostatic pressure (case 'F') and a uniform pressure (case 'W') on every quad
    for quad in model.quads.values():
        Yavg = (quad.i_node.Y + quad.j_node.Y + quad.m_node.Y + quad.n_node.Y)/4
        model.add_quad_surface_pressure(quad.name, 62.4*(b - Yavg), 'F')
        model.add_quad_surface_pressure(quad.name, 20.0, 'W')

    # Three load combinations, so the factorization is reused at least twice
    model.add_load_combo('1.0F', {'F': 1.0})
    model.add_load_combo('1.4F', {'F': 1.4})
    model.add_load_combo('1.2F+0.6W', {'F': 1.2, 'W': 0.6})

    return model


def _displacements(model):
    """All nodal displacements of every load combination, in a fixed order."""

    results = []
    for combo in model.load_combos:
        for node in model.nodes.values():
            results += [node.DX[combo], node.DY[combo], node.DZ[combo],
                        node.RX[combo], node.RY[combo], node.RZ[combo]]
    return results


def test_same_results_as_spsolve():
    """Reusing the factorization gives the same displacements as calling `spsolve` per combination."""

    # New behaviour: factorize once and reuse
    model = _plate_wall()
    model.analyze_linear(sparse=True)
    reused = _displacements(model)

    # Previous behaviour: no reusable factorization, so `spsolve` runs for every combination
    model = _plate_wall()
    with patch.object(Analysis, '_factorize_K11', return_value=None):
        model.analyze_linear(sparse=True)
    per_combo = _displacements(model)

    # The two solutions must agree to round-off
    assert reused == pytest.approx(per_combo, rel=1e-9, abs=1e-15)


def test_same_results_as_dense_solver():
    """The sparse path with the reused factorization agrees with the dense solver."""

    model = _plate_wall()
    model.analyze_linear(sparse=True)
    sparse_results = _displacements(model)

    model = _plate_wall()
    model.analyze_linear(sparse=False)
    dense_results = _displacements(model)

    assert sparse_results == pytest.approx(dense_results, rel=1e-8, abs=1e-15)


def test_factorized_only_once():
    """K11 is factorized once for all three load combinations."""

    model = _plate_wall()
    with patch.object(Analysis, '_factorize_K11', wraps=Analysis._factorize_K11) as factorize:
        model.analyze_linear(sparse=True)
    assert factorize.call_count == 1


def test_unstable_model_still_raises():
    """A mechanism (beam with an internal hinge and no support at the hinge) still raises."""

    model = FEModel3D()
    model.add_section('W10x33', 0.006264504, 0.000015234, 0.0000712, 0.000000241)
    model.add_material('steel', 200000000, 77000000, 0.3, 77)
    model.add_node('n1', 0, 0, 0)
    model.add_node('n2', 1, 0, 0)
    model.add_node('n3', 2, 0, 0)
    model.add_member('m1', 'n1', 'n2', 'steel', 'W10x33')
    model.add_member('m2', 'n2', 'n3', 'steel', 'W10x33')
    model.def_support('n1', 1, 1, 1, 1, 1, 0)
    model.def_support('n3', 0, 1, 1, 1, 1, 0)
    model.add_member_dist_load('m1', 'Fy', 10, 10)
    model.add_member_dist_load('m2', 'Fy', 10, 10)
    model.def_releases('m2', Rzi=True)

    with pytest.raises(Exception):
        model.analyze_linear(sparse=True)
