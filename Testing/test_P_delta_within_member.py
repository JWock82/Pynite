"""Tests P-little-delta results between the nodes of a member against closed-form beam-column solutions.

The closed forms are the standard ones for a prismatic member (Timoshenko & Gere, Theory of Elastic
Stability, Ch. 1). Compression is positive, k = sqrt(|P|/EI) and u = k*L/2. None of these models have
a node inside the member, so all of the second-order effect between the nodes comes from the segments.
"""
import math

import pytest

from Pynite import FEModel3D

E, I_strong, I_weak, L = 10000.0, 1.0e6, 4.0e5, 3000.0

# Each bending axis, with the load direction, moment and deflection names, and its moment of inertia
AXES = [('Fy', 'Mz', 'dy', I_strong), ('Fz', 'My', 'dz', I_weak)]


def _column(P, load_direction, fixed_ends, H=0.0, w=0.0):
    """Returns a column with no internal nodes, axial force `P`, and a transverse load in case 'H'."""

    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', L, 0, 0)
    model.add_material('Material', E, E/15, 0.3, 0)
    # The area is large so the closed forms, which ignore axial strain, apply to high accuracy
    model.add_section('Section', 5e6, I_weak, I_strong, 1e5)
    model.add_member('M1', 'N1', 'N2', 'Material', 'Section')
    model.def_support('N1', True, True, True, True, fixed_ends, fixed_ends)
    model.def_support('N2', False, True, True, False, fixed_ends, fixed_ends)
    model.add_node_load('N2', 'FX', -P, case='P')
    if H:
        model.add_member_pt_load('M1', load_direction, H, L/2, case='H')
    if w:
        model.add_member_dist_load('M1', load_direction, w, w, case='H')
    model.add_load_combo('Combo', {'P': 1.0, 'H': 1.0})
    model.analyze_PDelta()
    return model.members['M1']


@pytest.mark.parametrize('load_direction, M_name, d_name, I', AXES)
@pytest.mark.parametrize('P_ratio', [0.1, 0.3, 0.6, 0.9, -0.1, -0.6, -2.0])
def test_pinned_column_with_point_load(load_direction, M_name, d_name, I, P_ratio):
    """A point load at mid-span. Issue #350 is this case with no node at the load."""

    EI = E*I
    P = P_ratio*math.pi**2*EI/L**2
    H = 100.0
    k = math.sqrt(abs(P)/EI)
    u = k*L/2

    if P > 0:
        M_exact = H/(2*k)*math.tan(u)
        d_exact = H/(2*P*k)*(math.tan(u) - u)
    else:
        M_exact = H/(2*k)*math.tanh(u)
        d_exact = H/(2*-P*k)*(u - math.tanh(u))

    member = _column(P, load_direction, False, H=H)

    assert abs(member.moment(M_name, L/2, 'Combo')) == pytest.approx(M_exact, rel=1e-5)
    assert abs(member.deflection(d_name, L/2, 'Combo')) == pytest.approx(d_exact, rel=1e-5)


@pytest.mark.parametrize('load_direction, M_name, d_name, I', AXES)
@pytest.mark.parametrize('P_ratio', [0.1, 0.3, 0.6, 0.9, -0.6, -2.0])
def test_fixed_column_with_uniform_load(load_direction, M_name, d_name, I, P_ratio):
    """The end moments need the known end slopes, as the global solution's are first-order."""

    EI = E*I
    P = P_ratio*4*math.pi**2*EI/L**2
    w = 0.2
    k = math.sqrt(abs(P)/EI)
    u = k*L/2

    if P > 0:
        M_end = w/k**2*(u/math.tan(u) - 1)
        M_mid = w/k**2*(u/math.sin(u) - 1)
    else:
        M_end = w/k**2*(1 - u/math.tanh(u))
        M_mid = w/k**2*(1 - u/math.sinh(u))

    member = _column(P, load_direction, True, w=w)

    assert abs(member.moment(M_name, 0, 'Combo')) == pytest.approx(abs(M_end), rel=1e-5)
    assert abs(member.moment(M_name, L/2, 'Combo')) == pytest.approx(abs(M_mid), rel=1e-5)


@pytest.mark.parametrize('load_direction, M_name, d_name, I', AXES)
def test_deflected_shape_reaches_the_end_node(load_direction, M_name, d_name, I):
    """The chained segments must close on the end of the member, on both sides of a point load."""

    P = 0.6*math.pi**2*E*I/L**2
    member = _column(P, load_direction, False, H=100.0)

    assert member.deflection(d_name, L, 'Combo') == pytest.approx(0, abs=1e-8)
    assert member.deflection(d_name, L/2 - 1e-6, 'Combo') == pytest.approx(member.deflection(d_name, L/2 + 1e-6, 'Combo'), abs=1e-6)


def test_fixed_column_past_the_pinned_buckling_load():
    """A fixed column is stable past pi^2*EI/L^2, where a start-slope-only solution would be singular."""

    EI = E*I_strong
    P = math.pi**2*EI/L**2
    w = 0.2
    k = math.sqrt(P/EI)
    u = k*L/2

    member = _column(P, 'Fy', True, w=w)

    assert abs(member.moment('Mz', L/2, 'Combo')) == pytest.approx(w/k**2*(u/math.sin(u) - 1), rel=1e-5)


def test_taut_member_stays_finite():
    """Very high tension gives first-order results, not numbers that have overflowed or lost accuracy."""

    EI = E*I_strong
    H = 100.0
    P = -(100/L)**2*EI  # k*L = 100

    member = _column(P, 'Fy', False, H=H)

    # The global solution leaves a small residual end moment at this tension, hence the tolerance
    assert abs(member.moment('Mz', L/2, 'Combo')) == pytest.approx(H*L/4, rel=1e-3)
    assert abs(member.deflection('dy', L/2, 'Combo')) == pytest.approx(H*L**3/(48*EI), rel=1e-3)


def test_no_axial_force_matches_first_order():
    """With no axial force the beam-column solution has to reduce to the first-order equations."""

    EI = E*I_strong
    H = 100.0

    member = _column(0.0, 'Fy', False, H=H)

    assert abs(member.moment('Mz', L/2, 'Combo')) == pytest.approx(H*L/4, rel=1e-9)
    assert abs(member.deflection('dy', L/2, 'Combo')) == pytest.approx(H*L**3/(48*EI), rel=1e-9)
