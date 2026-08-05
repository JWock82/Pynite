from math import isclose

import matplotlib
matplotlib.use('TkAgg')
import numpy as np  # noqa: E402

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
    )

    # Plot the two traces on the same graph for each combo
    # import matplotlib.pyplot as plt
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


def test_moment_frame_pushover_example():
    """
    Matrix Analysis of Structures, 2nd Edition - Example 10.5 (Moment Frame Pushover)
    """

    # Create a new finite element model
    frame = FEModel3D()

    # Define the nodes of the frame
    frame.add_node('a', 0, 0, 0)
    frame.add_node('b', 0, 24, 0)
    frame.add_node('c', 20, 24, 0)
    frame.add_node('d', 60, 24, 0)
    frame.add_node('e', 60, 0, 0)

    # Intermediate nodes from the published example geometry
    frame.add_node('f', 0, 12, 0)
    frame.add_node('g', 60, 12, 0)
    frame.add_node('h', 10, 24, 0)
    frame.add_node('i', 40, 24, 0)

    # Define steel material properties
    E = 29000*12**2  # ksf
    G = 11200*12**2  # ksf
    nu = 0.3
    rho = 490/1000  # kips/ft^3
    fy = 36*12**2  # ksf
    frame.add_material('Steel', E, G, nu, rho, fy)

    # Define section properties
    frame.add_steel_section('W10x45', 13.3/12**2, 53.4/12**4, 248/12**4, 1.51/12**4, 20.3/12**3, 54.9/12**3, 'Steel')
    frame.add_steel_section('W27x84', 24.8/12**2, 106/12**4, 2850/12**4, 2.81/12**4, 33.2/12**3, 244/12**3, 'Steel')

    # Define the members of the frame
    frame.add_member('ab', 'a', 'b', 'Steel', 'W10x45')
    frame.add_member('bc', 'b', 'c', 'Steel', 'W27x84')
    frame.add_member('cd', 'c', 'd', 'Steel', 'W27x84')
    frame.add_member('ed', 'e', 'd', 'Steel', 'W10x45')

    # Define the supports
    frame.def_support('a', True, True, True, True, True, False)
    frame.def_support('e', True, True, True, True, True, False)
    frame.def_support('b', False, False, True, False, False, False)
    frame.def_support('d', False, False, True, False, False, False)

    # Define pushover loading pattern
    frame.add_node_load('b', 'FX', 0.99*6, 'Push')
    frame.add_node_load('c', 'FY', 0.99*-60, 'Push')
    frame.add_node_load('d', 'FY', 0.99*-120, 'Push')

    # Define load combinations
    frame.add_load_combo('Primary', {})
    frame.add_load_combo('Pushover', {'Push': 0.01})

    traces = {
        'Node d Drift': lambda combo_name: frame.nodes['d'].DX[combo_name]*12,
    }

    frame.analyze_pushover(
        log=True,
        check_stability=False,
        push_combo='Pushover',
        traces=traces,
    )

    # 1) Check the moment at point c (take magnitude because sign is convention-dependent).
    moment_at_c = abs(frame.members['bc'].moment('Mz', x=20.0, combo_name='Primary'))
    assert isclose(moment_at_c, 731.9, rel_tol=0.02)

    drift_trace = np.array(frame._pushover_traces['Primary']['Node d Drift'], dtype=float)
    assert len(drift_trace) > 2

    # 2) Confirm bilinear behavior and identify slope change near step 85.
    
    # First, compute incremental slope (drift increase per pushover step). If the response is
    # piecewise linear, this slope should be nearly constant in each segment and then jump once
    # at the transition point.
    slopes = np.diff(drift_trace)

    # The largest absolute change in consecutive slopes marks the most likely "kink" location.
    # `np.diff(slopes)` lives between step indices, so we add 1 to map back to the trace index.
    slope_change_step = int(np.argmax(np.abs(np.diff(slopes))) + 1)

    # Example 10.5 indicates this transition should occur near step 85; allow a small tolerance
    # band so tiny numerical differences between platforms do not produce false test failures.
    assert 82 <= slope_change_step <= 88

    pre = drift_trace[:slope_change_step + 1]
    post = drift_trace[slope_change_step:]

    # Fit a straight line to each branch and verify both branches are strongly linear.
    # High R^2 on each side plus different slopes is our automated bilinear check.
    x_pre = np.arange(len(pre), dtype=float)
    fit_pre = np.polyfit(x_pre, pre, 1)
    pred_pre = fit_pre[0]*x_pre + fit_pre[1]
    r2_pre = 1 - np.sum((pre - pred_pre)**2)/np.sum((pre - np.mean(pre))**2)

    x_post = np.arange(len(post), dtype=float)
    fit_post = np.polyfit(x_post, post, 1)
    pred_post = fit_post[0]*x_post + fit_post[1]
    r2_post = 1 - np.sum((post - pred_post)**2)/np.sum((post - np.mean(post))**2)

    # Both sides should be near-linear, and post-kink slope should be steeper than pre-kink slope.
    assert r2_pre > 0.99
    assert r2_post > 0.99
    assert fit_post[0] > fit_pre[0]

    # 3) Check displacement just before slope change.
    assert 3.0 < drift_trace[slope_change_step] < 4.0

    # 4) Check final displacement at end of the trace.
    assert 5.0 < drift_trace[-1] < 6.0

if __name__ == '__main__':

    test_plastic_beam()
    test_moment_frame_pushover_example()
