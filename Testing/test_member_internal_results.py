from Pynite import FEModel3D
import math


def test_beam_internal_forces():
    """
    Units for this model are kips and feet
    """

    # Define a new beam
    beam = FEModel3D()

    # Define the nodes
    beam.add_node('N1', 0, 0, 0)
    beam.add_node('N2', 10, 0, 0)

    # Define the supports (simply supported)
    beam.def_support('N1', True, True, True, True, False, False)
    beam.def_support('N2', True, True, True, True, False, False)

    # Define beam section proerties
    J = 400/12**4
    Iy = 200/12**4
    Iz = 200/12**4
    A = 12/12**2
    beam.add_section('Section', A, Iy, Iz, J)

    # Define a material
    E = 29000*144  # ksf
    G = 11200*144  # ksf
    nu = 0.3
    rho = 0.490  # pcf
    beam.add_material('Steel', E, G, nu, rho)

    # Create the beam
    beam.add_member('M1', 'N1', 'N2', 'Steel', 'Section')

    # Add a mid-span node to force the model to split a physical member for this test
    beam.add_node('N3', 5, 0, 0)

    # Add a member distributed load along the strong axis
    beam.add_member_dist_load('M1', 'FY', -0.5, -0.5, case='D')
    beam.add_member_dist_load('M1', 'FY', -0.75, -0.75, case='L')

    # Add a member distributed laod along the weak axis
    beam.add_member_dist_load('M1', 'FZ', -0.5, -0.5, case='D')
    beam.add_member_dist_load('M1', 'FZ', -0.75, -0.75, case='L')

    # Add some load combinations
    beam.add_load_combo('D', {'D': 1.0}, 'blc')
    beam.add_load_combo('L', {'L': 1.0}, 'blc')
    beam.add_load_combo('1.2D + 1.6L', {'D': 1.2, 'L': 1.6}, 'strength')

    # Analyze the model
    beam.analyze_linear()

    # from Pynite.Visualization import Renderer
    # renderer = Renderer(beam)
    # renderer.combo_name = 'D'
    # renderer.annotation_size = 1
    # renderer.render_model()

    # Check the shear diagram
    assert math.isclose(beam.members['M1'].shear('Fy', 0, 'D'), 2.5, abs_tol=0.01), 'Fy internal shear test failed at start of member.'
    assert math.isclose(beam.members['M1'].shear('Fz', 0, 'D'), 2.5, abs_tol=0.01), 'Fz internal shear test failed at start of member.'
    assert math.isclose(beam.members['M1'].shear('Fy', 10, 'D'), -2.5, abs_tol=0.01), 'Fy internal shear test failed at end of member.'
    assert math.isclose(beam.members['M1'].shear('Fz', 10, 'D'), -2.5, abs_tol=0.01), 'Fz internal shear test failed at end of member.'
    assert math.isclose(beam.members['M1'].shear('Fy', 5, 'D'), 0, abs_tol=0.01), 'Fy internal shear test failed at midpoint of member.'
    assert math.isclose(beam.members['M1'].shear('Fz', 5, 'D'), 0, abs_tol=0.01), 'Fz internal shear test failed at midpoint of member.'

    assert math.isclose(beam.members['M1'].max_shear('Fy', 'D'), 2.5, abs_tol=0.01), 'Fy internal max shear test failed.'
    assert math.isclose(beam.members['M1'].max_shear('Fz', 'D'), 2.5, abs_tol=0.01), 'Fz internal max shear test failed.'
    assert math.isclose(beam.members['M1'].min_shear('Fy', 'D'), -2.5, abs_tol=0.01), 'Fy internal min shear test failed.'
    assert math.isclose(beam.members['M1'].min_shear('Fz', 'D'), -2.5, abs_tol=0.01), 'Fz internal min shear test failed.'
    max_shear_env, max_shear_combo = beam.members['M1'].max_shear('Fy', ['strength'])
    assert math.isclose(max_shear_env, 9.0, abs_tol=0.01), 'Failed envelope shear test.'
    assert max_shear_combo == '1.2D + 1.6L', 'Failed envelope shear governing combo test.'
    min_shear_env, min_shear_combo = beam.members['M1'].min_shear('Fz', ['strength'])
    assert math.isclose(min_shear_env, -9.0, abs_tol=0.01), 'Failed envelope shear test.'
    assert min_shear_combo == '1.2D + 1.6L', 'Failed envelope shear governing combo test.'

    # Check the moment diagram
    assert math.isclose(beam.members['M1'].moment('Mz', 0, 'D'), 0, abs_tol=2), 'Mz internal moment test failed at start of member.'
    assert math.isclose(beam.members['M1'].moment('My', 0, 'D'), 0, abs_tol=2), 'My internal moment test failed at start of member.'
    assert math.isclose(beam.members['M1'].moment('Mz', 5, 'D'), -6.25, abs_tol=2), 'Mz internal moment test failed at midpoint of member.'
    assert math.isclose(beam.members['M1'].moment('My', 5, 'D'), -6.25, abs_tol=2), 'My internal moment test failed at midpoint of member.'
    assert math.isclose(beam.members['M1'].moment('Mz', 10, 'D'), 0, abs_tol=2), 'Mz internal moment test failed at end of member.'
    assert math.isclose(beam.members['M1'].moment('My', 10, 'D'), 0, abs_tol=2), 'My internal moment test failed at end of member.'

    assert math.isclose(beam.members['M1'].min_moment('Mz', 'D'), -6.25, abs_tol=2), 'Mz internal min moment test failed.'
    assert math.isclose(beam.members['M1'].min_moment('My', 'D'), -6.25, abs_tol=2), 'My internal min moment test failed.'
    assert math.isclose(beam.members['M1'].max_moment('Mz', 'D'), 0, abs_tol=2), 'Mz internal max moment test failed.'
    assert math.isclose(beam.members['M1'].max_moment('My', 'D'), 0, abs_tol=2), 'My internal max moment test failed.'
    min_moment_env, min_moment_combo = beam.members['M1'].min_moment('Mz', ['strength'])
    assert math.isclose(min_moment_env, -22.5, abs_tol=2), 'Failed member Mz envelope results test.'
    assert min_moment_combo == '1.2D + 1.6L', 'Failed member Mz envelope governing combo test.'
    max_moment_env, max_moment_combo = beam.members['M1'].max_moment('My', ['strength'])
    assert math.isclose(max_moment_env, 0, abs_tol=2), 'Failed member My envelope results test.'

    # Check the deflected shape
    assert math.isclose(beam.members['M1'].deflection('dy', 0, 'D')*12, 0, abs_tol=0.00001), 'dy internal deflection test failed at start of member.'
    assert math.isclose(beam.members['M1'].deflection('dz', 0, 'D')*12, 0, abs_tol=0.00001), 'dz internal deflection test failed at start of member.'
    assert math.isclose(beam.members['M1'].deflection('dz', 5, 'D')*12, 5*(-0.5)*10**4/(384*E*Iz)*12, abs_tol=0.00001), 'dz internal deflection test failed at midpoint of member.'
    assert math.isclose(beam.members['M1'].deflection('dy', 5, 'D')*12, 5*(-0.5)*10**4/(384*E*Iy)*12, abs_tol=0.00001), 'dy internal deflection test failed at midpoint of member.'
    assert math.isclose(beam.members['M1'].deflection('dy', 10, 'D')*12, 0, abs_tol=0.00001), 'dy internal deflection test failed at end of member.'
    assert math.isclose(beam.members['M1'].deflection('dz', 10, 'D')*12, 0, abs_tol=0.00001), 'dz internal deflection test failed at end of member.'

    assert math.isclose(beam.members['M1'].max_deflection('dy', 'D')*12, 0, abs_tol=0.00001), 'dy internal max deflection test failed.'
    assert math.isclose(beam.members['M1'].max_deflection('dz', 'D')*12, 0, abs_tol=0.00001), 'dz internal max deflection test failed.'
    assert math.isclose(beam.members['M1'].min_deflection('dz', 'D')*12, 5*(-0.5)*10**4/(384*E*Iz)*12, abs_tol=0.00001), 'dz internal min deflection test failed.'
    assert math.isclose(beam.members['M1'].min_deflection('dy', 'D')*12, 5*(-0.5)*10**4/(384*E*Iy)*12, abs_tol=0.00001), 'dy internal min deflection test failed.'
    min_defl_dz, min_defl_dz_combo = beam.members['M1'].min_deflection('dz', ['strength'])
    assert math.isclose(min_defl_dz*12, 5*(-1.8)*10**4/(384*E*Iz)*12, abs_tol=0.00001), 'Failed member dz envelope test.'
    assert min_defl_dz_combo == '1.2D + 1.6L', 'Failed member dz envelope governing combo test.'
    min_defl_dy, min_defl_dy_combo = beam.members['M1'].min_deflection('dy', ['strength'])
    assert math.isclose(min_defl_dy*12, 5*(-1.8)*10**4/(384*E*Iy)*12, abs_tol=0.00001), 'Failed member dy envelope test.'
    assert min_defl_dy_combo == '1.2D + 1.6L', 'Failed member dy envelope governing combo test.'


def test_max_moment_PDelta_weak_axis():
    """
    `max_moment('My')` should include P-little-delta effects after a P-Delta analysis, the same
    way `max_moment('Mz')` and `moment('My', x)` do.
    """

    # Simply supported beam with an axial compression of 200 and a uniform load of 1 about the
    # weak axis
    beam = FEModel3D()
    beam.add_node('N1', 0, 0, 0)
    beam.add_node('N2', 10, 0, 0)
    beam.def_support('N1', True, True, True, True, False, False)
    beam.def_support('N2', False, True, True, False, False, False)
    beam.add_material('Steel', 29000, 11200, 0.3, 0.49)
    beam.add_section('Section', 20, 100, 150, 250)
    beam.add_member('M1', 'N1', 'N2', 'Steel', 'Section')
    beam.add_node_load('N2', 'FX', -200)
    beam.add_member_dist_load('M1', 'FZ', 1, 1)

    # Run a P-Delta analysis
    beam.analyze_PDelta()
    member = beam.members['M1']

    # The first-order moment at midspan is wL^2/8 = 12.5. The axial load amplifies it.
    M_mid = member.moment('My', 5)
    assert M_mid > 12.505, 'P-Delta did not amplify the weak axis moment.'

    # The maximum moment should match the amplified moment at midspan
    assert math.isclose(member.max_moment('My'), M_mid, rel_tol=1e-6), \
        'My max moment ignores P-Delta.'


if __name__ == '__main__':
    test_beam_internal_forces()
    test_max_moment_PDelta_weak_axis()
