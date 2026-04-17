from math import isclose
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from Pynite import FEModel3D  # noqa: E402


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
    plastic_beam.add_member('M1', 'N1', 'N3', 'Stl_A992', 'W12x65')

    # Add a load
    # P = 259.3
    # P = 325.7
    P = 346.2
    plastic_beam.add_node_load('N2', 'FY', -0.3*P, 'Push')
    plastic_beam.add_node_load('N3', 'FX', -P, 'Push')

    # Add load combinations
    # Primary combo is required for pushover analysis to have a base case to apply pushover loads to
    plastic_beam.add_load_combo('Primary', {})
    plastic_beam.add_load_combo('Pushover', {'Push': 0.01})

    traces = {
        'Fixed End Moment': lambda combo_name: plastic_beam.members['M1'].moment('Mz', x=0.0, combo_name=combo_name),
        'Load Point Moment': lambda combo_name: plastic_beam.members['M1'].moment('Mz', x=96.0, combo_name=combo_name),
        'Phi N1': lambda combo_name: plastic_beam.members['M1'].section.Phi(
            plastic_beam.members['M1']._fxi.get(combo_name, 0.0),
            plastic_beam.members['M1']._myi.get(combo_name, 0.0),
            plastic_beam.members['M1']._mzi.get(combo_name, 0.0)
        ),
        'Load Point Deflection': lambda combo_name: plastic_beam.nodes['N2'].DY[combo_name],
    }

    # Analysis the model
    plastic_beam.analyze_pushover(log=True, check_stability=False, push_combo='Pushover', max_iter=30, tol=0.01, sparse=False, traces=traces)

    # Plot the traces one by one
    plastic_beam.plot_pushover_trace('Fixed End Moment', combo_name='Primary')
    plastic_beam.plot_pushover_trace('Load Point Moment', combo_name='Primary')
    plastic_beam.plot_pushover_trace('Phi N1', combo_name='Primary')
    plastic_beam.plot_pushover_trace('Load Point Deflection', combo_name='Primary')

    # Get the resulting moments from the Primary combo (where pushover results are stored)
    M_N1 = plastic_beam.members['M1'].moment('Mz', x=0.0, combo_name='Primary')
    M_N2 = plastic_beam.members['M1'].moment('Mz', x=96.0, combo_name='Primary')
    print(f'Final Mz at x=0   : {M_N1}')
    print(f'Final Mz at x=96  : {M_N2}')
    M_N1e = 3752
    M_N2e = -3752
    # M_hinge = 4148

    assert isclose(M_N1, M_N1e, rel_tol=0.10), f'Calculated Moment: {M_N1} \nExpected Moment: {M_N1e}'
    assert isclose(M_N2, M_N2e, rel_tol=0.10), f'Calculated Moment: {M_N2} \nExpected Moment: {M_N2e}'


if __name__ == '__main__':

    test_plastic_beam()
