from math import isclose
from Pynite import FEModel3D


def test_plastic_beam():
    """
    Matrix Structural Analysis, 2nd Edition - Examples 8.6 (p. 228) & 10.4 (p. 282)
    """

    # Create the model
    plastic_beam = FEModel3D()

    # Define a material
    E = 29000  # ksi
    G = 11200  # ksi
    nu = 0.3
    rho = 0.490/12**3  # kci
    fy = 50  # ksi
    plastic_beam.add_material('Stl_A992', E, G, nu, rho, fy)

    # Define a cross-section
    plastic_beam.add_steel_section('W12x65', 19.1, 20, 533, 1, 15, 96.8, 'Stl_A992')

    # Add nodes
    plastic_beam.add_node('N1', 0, 0, 0)
    plastic_beam.add_node('N2', 8*12, 0, 0)
    plastic_beam.add_node('N3', 24*12, 0, 0)

    # Add supports
    plastic_beam.def_support('N1', True, True, True, True, True, True)
    plastic_beam.def_support('N3', False, True, True, False, False, False)

    # Add a member
    plastic_beam.add_member('M1', 'N1', 'N2', 'Stl_A992', 'W12x65')
    plastic_beam.add_member('M2', 'N2', 'N3', 'Stl_A992', 'W12x65')

    # Add a load
    plastic_beam.add_node_load('N2', 'FY', 0.3, 'Push')
    plastic_beam.add_node_load('N3', 'FX', -1.0, 'Push')

    # Add a load combination
    plastic_beam.add_load_combo('Pushover', {'Push': 1.0})

    # Analysis the model
    plastic_beam._not_ready_yet_analyze_pushover(log=True, check_stability=False, push_combo='Pushover', max_iter=30, tol=0.01, sparse=False)

    # Get the resulting moments
    calculated_moment = plastic_beam.members['M1'].min_moment('Mz')
    expected_moment = -0.5*(10*12)**2/8

    assert isclose(calculated_moment, expected_moment), f'Calculated Moment: {calculated_moment} \nExpected Moment: {expected_moment}'


if __name__ == '__main__':

    test_plastic_beam()