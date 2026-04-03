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
    plastic_beam.add_node('N2', 24*12, 0, 0)

    # Add supports
    plastic_beam.def_support('N1', True, True, True, True, True, True)
    plastic_beam.def_support('N2', False, True, True, False, False, False)

    # Add a member
    plastic_beam.add_member('M1', 'N1', 'N2', 'Stl_A992', 'W12x65')

    # Add a load
    plastic_beam.add_member_pt_load('M1', 'Fy', -1.0, 8*12, 'Push')
    plastic_beam.add_node_load('N2', 'FX', -1.0, 'Push')

    # Add a load combination
    plastic_beam.add_load_combo('Pushover', {'Push': 1.0})

    # Analysis the model
    plastic_beam._not_ready_yet_analyze_pushover(log=True, check_stability=False, push_combo='Pushover', max_iter=30, tol=0.01, sparse=False)

    # Get the resulting moments
    M_a = plastic_beam.members['M1'].moment('Mz', x=0.0, combo_name='Pushover')
    M_b = plastic_beam.members['M1'].moment('Mz', x=8.0*12.0, combo_name='Pushover')
    M_ae = 3752
    M_be = -3752

    assert isclose(M_a, M_ae), f'Calculated Moment: {M_a} \nExpected Moment: {M_ae}'
    assert isclose(M_b, M_be), f'Calculated Moment: {M_b} \nExpected Moment: {M_be}'


if __name__ == '__main__':

    test_plastic_beam()