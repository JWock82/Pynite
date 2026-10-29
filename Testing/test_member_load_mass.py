"""Verify load-to-mass projection using independently known member axes."""

import numpy as np
import pytest

from Pynite import FEModel3D


def make_cantilever(end=(8, 6, 0), rotation=0):
    """Create a ten-metre member with a fixed base and SI material properties."""
    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', *end)
    model.def_support('N1', True, True, True, True, True, True)
    # Pynite's rho is weight per volume; gravity converts it to mass per volume.
    model.add_material('Steel', 200e9, 80e9, 0.3, 7800*9.81)
    model.add_section('Section', 0.01, 8.33e-6, 8.33e-6, 1.67e-5)
    model.add_member('M', 'N1', 'N2', 'Steel', 'Section', rotation=rotation)
    model.add_load_combo('Mass', {'D': 1.5})
    return model


# Rows are the local unit axes expressed in global coordinates, from geometry.
# No expected load vector is calculated using Member3D.T().
GEOMETRIES = [
    pytest.param((10, 0, 0), 0, ((1, 0, 0), (0, 1, 0), (0, 0, 1)), id='horizontal'),
    pytest.param((8, 6, 0), 0, ((0.8, 0.6, 0), (-0.6, 0.8, 0), (0, 0, 1)), id='inclined'),
    pytest.param((0, 10, 0), 0, ((0, 1, 0), (-1, 0, 0), (0, 0, 1)), id='vertical'),
    pytest.param((10, 0, 0), 90, ((1, 0, 0), (0, 0, 1), (0, -1, 0)), id='rolled'),
]


@pytest.mark.parametrize('end,rotation,axes', GEOMETRIES)
@pytest.mark.parametrize('load_axis', range(3))
@pytest.mark.parametrize('mass_direction', ['X', 'Y', 'Z'])
@pytest.mark.parametrize('kind', ['point', 'distributed'])
def test_local_load_mass_projection(end, rotation, axes, load_axis, mass_direction, kind):
    """A force component divided by g gives mass in the selected GLOBAL direction."""
    model = make_cantilever(end, rotation)
    direction = ('Fx', 'Fy', 'Fz')[load_axis]
    if kind == 'point':
        model.add_member_pt_load('M', direction, -1000, 5, case='D')
    else:
        # The full-span uniform load has resultant -1000 N at midspan.
        model.add_member_dist_load('M', direction, -100, -100, case='D')

    # The unrelated load case must not be converted into mass.
    model.add_member_pt_load('M', 'Fy', -99999, 1, case='Excluded')
    mass, x = model.members['M']._calc_load_mass('Mass', mass_direction, gravity=9.81)
    component = axes[load_axis]['XYZ'.index(mass_direction)]
    assert mass == pytest.approx(abs(1.5*1000*component)/9.81, abs=1e-12)
    assert x == pytest.approx(5)


@pytest.mark.parametrize('mass_direction', ['X', 'Y', 'Z'])
def test_global_point_load_mass(mass_direction):
    """Global force inputs retain their direct component interpretation."""
    model = make_cantilever()
    for direction, force in [('FX', -200), ('FY', -300), ('FZ', -400)]:
        model.add_member_pt_load('M', direction, force, 4, case='D')
    mass, x = model.members['M']._calc_load_mass('Mass', mass_direction, 9.81)
    force = {'X': 200, 'Y': 300, 'Z': 400}[mass_direction]
    assert mass == pytest.approx(1.5*force/9.81)
    assert x == pytest.approx(4)


def test_point_load_mass_combination_and_location():
    """Load factors apply before finding the force resultant and its location."""
    model = make_cantilever()
    model.add_member_pt_load('M', 'Fy', -1000, 2, case='D')
    model.add_member_pt_load('M', 'Fy', -2000, 8, case='L')
    model.add_member_pt_load('M', 'Mz', 500, 3, case='D')
    model.add_load_combo('Mass', {'D': 1.5, 'L': 0.5})
    mass, x = model.members['M']._calc_load_mass('Mass', 'Y', 9.81)
    assert mass == pytest.approx(0.8*(1500 + 1000)/9.81)
    assert x == pytest.approx((1500*2 + 1000*8)/(1500 + 1000))


@pytest.mark.parametrize('direction,component', [('Fx', 0.6), ('Fy', 0.8), ('Fz', 0)])
def test_tapered_local_load_total_mass(direction, component):
    """The projected trapezoid area gives total mass, independently of its lumping."""
    model = make_cantilever()
    model.add_member_dist_load('M', direction, -100, -200, 2, 8, case='D')
    mass, _ = model.members['M']._calc_load_mass('Mass', 'Y', 9.81)
    # This checks force integration only. Distributed-load centroid is separate behavior.
    assert mass == pytest.approx(abs(1.5*component*(100 + 200)/2*(8 - 2))/9.81)


@pytest.mark.parametrize('kind', ['point', 'distributed'])
def test_modal_local_and_global_loads_are_equivalent(kind):
    """Representing the same physical load in another basis preserves frequencies."""
    frequencies = []
    mass_matrices = []
    for local in [True, False]:
        model = make_cantilever()
        # e_y = (-0.6, 0.8, 0), so -1000 e_y = (600, -800, 0) N.
        loads = [('Fy', -1000)] if local else [('FX', 600), ('FY', -800)]
        for direction, force in loads:
            if kind == 'point':
                model.add_member_pt_load('M', direction, force, 5, case='D')
            else:
                model.add_member_dist_load('M', direction, force/10, force/10, case='D')

        # Add the same consistent self-weight mass to both models.
        model.add_member_self_weight('FY', -1, case='Self weight')
        model.add_load_combo('Mass', {'D': 1.5, 'Self weight': 1})
        model.analyze_modal(num_modes=2, mass_combo_name='Mass', gravity=9.81)
        frequencies.append(model.frequencies.copy())
        mass_matrices.append(model.M('Mass', 'Y', 9.81, sparse=False))

    assert np.all(np.isfinite(frequencies[0]))
    assert np.all(frequencies[0] > 0)
    np.testing.assert_allclose(mass_matrices[0], mass_matrices[1], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(frequencies[0], frequencies[1], rtol=1e-9, atol=1e-12)
