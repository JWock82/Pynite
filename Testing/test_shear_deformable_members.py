"""Closed-form checks for shear-deformable (Timoshenko-Ehrenfest) members.

Each expected value is the Euler-Bernoulli bending term plus the shear term
from integrating V/(k*G*A) along the member. A timber-like section
(E/G = 16) is used so the shear term is a meaningful share of the result.
"""

import math

import pytest

from Pynite import FEModel3D

E = 11e9       # Pa
G = 0.69e9     # Pa
NU = 0.3
RHO = 0.0
B = 0.2        # m
H = 0.6        # m
A = B*H
IZ = B*H**3/12  # strong axis, bending under local y loads
IY = H*B**3/12
J = 3.33e-3
K = 5/6         # shear correction factor for a solid rectangle
L = 3.0         # m
P = 10e3        # N
W = 5e3         # N/m


def build_model(nodes, members, shear_deformable=True):
    model = FEModel3D()
    model.add_material('Wood', E, G, NU, RHO)
    model.add_section('Rect', A, IY, IZ, J, ksy=K, ksz=K)
    for name, x in nodes:
        model.add_node(name, x, 0, 0)
    for name, i_node, j_node in members:
        model.add_member(name, i_node, j_node, 'Wood', 'Rect', shear_deformable=shear_deformable)
    return model


def cantilever(shear_deformable=True):
    model = build_model([('N1', 0), ('N2', L)], [('M1', 'N1', 'N2')], shear_deformable)
    model.def_support('N1', True, True, True, True, True, True)
    return model


def test_cantilever_tip_load_includes_shear_deflection():
    model = cantilever()
    model.add_node_load('N2', 'FY', -P)
    model.analyze()

    expected = P*L**3/(3*E*IZ) + P*L/(K*G*A)
    assert math.isclose(-model.nodes['N2'].DY['Combo 1'], expected, rel_tol=1e-9)
    assert math.isclose(-model.members['M1'].deflection('dy', L, 'Combo 1'), expected, rel_tol=1e-9)


def test_cantilever_without_shear_flag_is_euler_bernoulli():
    model = cantilever(shear_deformable=False)
    model.add_node_load('N2', 'FY', -P)
    model.analyze()

    assert math.isclose(-model.nodes['N2'].DY['Combo 1'], P*L**3/(3*E*IZ), rel_tol=1e-9)


def test_minor_axis_cantilever_uses_iy_and_ksz():
    model = cantilever()
    model.add_node_load('N2', 'FZ', -P)
    model.analyze()

    expected = P*L**3/(3*E*IY) + P*L/(K*G*A)
    assert math.isclose(-model.nodes['N2'].DZ['Combo 1'], expected, rel_tol=1e-9)


def two_span_model(fixed_ends):
    """A single span split at midspan so the midspan result is a nodal displacement."""
    model = build_model([('N1', 0), ('N2', L/2), ('N3', L)], [('M1', 'N1', 'N2'), ('M2', 'N2', 'N3')])
    if fixed_ends:
        model.def_support('N1', True, True, True, True, True, True)
        model.def_support('N3', True, True, True, True, True, True)
    else:
        model.def_support('N1', True, True, True, True, False, False)
        model.def_support('N3', False, True, True, False, False, False)
    model.add_member_dist_load('M1', 'Fy', -W, -W)
    model.add_member_dist_load('M2', 'Fy', -W, -W)
    model.analyze()
    return model


def test_simply_supported_uniform_load_midspan_deflection():
    model = two_span_model(fixed_ends=False)

    expected = 5*W*L**4/(384*E*IZ) + W*L**2/(8*K*G*A)
    assert math.isclose(-model.nodes['N2'].DY['Combo 1'], expected, rel_tol=1e-9)


def test_fixed_fixed_uniform_load_deflection_and_end_moments():
    model = two_span_model(fixed_ends=True)

    # For a symmetric fixed-fixed span the end moments do not depend on shear flexibility.
    expected_deflection = W*L**4/(384*E*IZ) + W*L**2/(8*K*G*A)
    assert math.isclose(-model.nodes['N2'].DY['Combo 1'], expected_deflection, rel_tol=1e-9)
    assert math.isclose(abs(model.nodes['N1'].RxnMZ['Combo 1']), W*L**2/12, rel_tol=1e-9)


def test_propped_cantilever_prop_reaction_depends_on_shear():
    model = cantilever()
    model.def_support('N2', False, True, True, False, False, False)
    model.add_member_dist_load('M1', 'Fy', -W, -W)
    model.analyze()

    # Compatibility at the prop: the released deflection under the load equals
    # the deflection caused by the prop force, both including shear.
    released = W*L**4/(8*E*IZ) + W*L**2/(2*K*G*A)
    per_unit_force = L**3/(3*E*IZ) + L/(K*G*A)
    expected = released/per_unit_force
    reaction = model.nodes['N2'].RxnFY['Combo 1']
    assert math.isclose(reaction, expected, rel_tol=1e-9)
    assert reaction > 3*W*L/8


def single_simply_supported_member():
    model = build_model([('N1', 0), ('N2', L)], [('M1', 'N1', 'N2')])
    model.def_support('N1', True, True, True, True, False, False)
    model.def_support('N2', False, True, True, False, False, False)
    return model


@pytest.mark.xfail(strict=True, reason='member.deflection() between nodes uses Euler-Bernoulli integration and omits the shear term')
def test_interior_deflection_includes_shear_under_uniform_load():
    model = single_simply_supported_member()
    model.add_member_dist_load('M1', 'Fy', -W, -W)
    model.analyze()

    expected = 5*W*L**4/(384*E*IZ) + W*L**2/(8*K*G*A)
    assert math.isclose(-model.members['M1'].deflection('dy', L/2, 'Combo 1'), expected, rel_tol=1e-6)


@pytest.mark.xfail(strict=True, reason='member.deflection() between nodes uses Euler-Bernoulli integration and omits the shear term')
def test_interior_deflection_includes_shear_under_point_load():
    model = single_simply_supported_member()
    model.add_member_pt_load('M1', 'Fy', -P, L/2)
    model.analyze()

    expected = P*L**3/(48*E*IZ) + P*L/(4*K*G*A)
    assert math.isclose(-model.members['M1'].deflection('dy', L/2, 'Combo 1'), expected, rel_tol=1e-6)
