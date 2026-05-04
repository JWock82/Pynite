from math import isclose

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt  # noqa: E402
from numpy import array  # noqa: E402

from Pynite.FEModel3D import FEModel3D  # noqa: E402


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
    P = 334.99
    plastic_beam.add_node_load('N2', 'FY', -0.3*P, 'Push')
    plastic_beam.add_node_load('N3', 'FX', -P, 'Push')

    # Add load combinations
    # Primary combo is required for pushover analysis to have a base case to apply pushover loads to
    plastic_beam.add_load_combo('Primary', {})
    plastic_beam.add_load_combo('Pushover', {'Push': 0.01})
    
    # Set up traces for the moment at `b`, both looking back toward `a` and looking forward toward `c`.
    traces = {
        'M_ba': lambda combo_name: plastic_beam.members['M1'].sub_members['M1a'].moment('Mz', x=96.0, combo_name=combo_name),
        'M_bc': lambda combo_name: plastic_beam.members['M1'].sub_members['M1b'].moment('Mz', x=0.0, combo_name=combo_name),   
    }

    plastic_beam.analyze_pushover(
        log=True,
        check_stability=False,
        push_combo='Pushover',
        max_iter=30,
        tol=0.01,
        sparse=False,
        traces=traces,
        P_Delta=False,
    )

    # Plot the two traces on the same graph for each combo
    # for combo_name, trace_data in plastic_beam._pushover_traces.items():
    #     steps = list(range(len(trace_data['M_ba'])))
    #     plt.figure()
    #     plt.plot(steps, trace_data['M_ba'], label='M_ba (M1a end)')
    #     plt.plot(steps, trace_data['M_bc'], label='M_bc (M1b start)')
    #     plt.xlabel('Pushover Step')
    #     plt.ylabel('Moment at b (kip-in)')
    #     plt.title(f'Moment Traces at Node b — Combo: {combo_name}')
    #     plt.legend()
    #     plt.tight_layout()
    # plt.show()

    # Compare the traces to each other
    for combo_name, trace_data in plastic_beam._pushover_traces.items():
        for step, (M_ba, M_bc) in enumerate(zip(trace_data['M_ba'], trace_data['M_bc'])):
            assert isclose(M_ba, M_bc, rel_tol=0.01), f"Traces do not match for combo {combo_name} step {step}: M_ba={M_ba}, M_bc={M_bc}"

if __name__ == '__main__':

    test_plastic_beam()
