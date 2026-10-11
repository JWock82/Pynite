import math

from Pynite import FEModel3D

def test_tc_braced_frame():

    # Create a new 3D frame model
    frame = FEModel3D()

    # Add nodes to the frame
    frame.add_node('N1', 0, 0, 0)
    frame.add_node('N2', 15, 0, 0)
    frame.add_node('N3', 0, 15, 0)
    frame.add_node('N4', 15, 15, 0)

    # Define material properties
    E = 29000/144
    G = 0.4*E
    nu = 0.17
    rho = 150/1000
    frame.add_material('Steel', E, G, nu, rho)

    # Define beam member section
    Iz = 204/12**4
    Iy = 17.3/12**4
    A = 7.62/144
    J = 0.300/12**4
    frame.add_section('W12x26', A, Iy, Iz, J)

    # Define column member section
    Iz = 171/12**4
    Iy = 36.6/12**4
    A = 9.71/144
    J = 0.583/12**4
    frame.add_section('W10x33', A, Iy, Iz, J)

    # Define brace member section
    Iz = 3.67/12**4
    Iy = 3.67/12**4
    A = 2.40/144
    J = 0.0832/12**4
    frame.add_section('L4x4x5/16', A, Iy, Iz, J)

    # Add members to the frame
    frame.add_member('C1', 'N1', 'N3', 'Steel', 'W10x33')
    frame.add_member('C2', 'N2', 'N4', 'Steel', 'W10x33')
    frame.add_member('B1', 'N3', 'N4', 'Steel', 'W12x26')
    frame.add_member('Br1', 'N1', 'N4', 'Steel', 'L4x4x5/16', tension_only=True)
    frame.add_member('Br2', 'N2', 'N3', 'Steel', 'L4x4x5/16', tension_only=True)

    # Release strong & weak axis moments at the i-ends; release both beam ends to keep the frame symmetric.
    frame.def_releases('C1', False, False, False, False, True, True,
                       False, False, False, False, False, False)
    frame.def_releases('C2', False, False, False, False, True, True,
                       False, False, False, False, False, False)
    frame.def_releases('B1', False, False, False, False, True, True,
                       False, False, False, False, True, True)
    frame.def_releases('Br1', False, False, False, False, True, True,
                       False, False, False, False, True, True)
    frame.def_releases('Br2', False, False, False, False, True, True,
                       False, False, False, False, True, True)

    # Fully support the base nodes
    frame.def_support('N1', True, True, True, True, True, True)
    frame.def_support('N2', True, True, True, True, True, True)

    # Support the top nodes out-of-plane
    frame.def_support('N3', False, False, True, False, False, False)
    frame.def_support('N4', False, False, True, False, False, False)

    # Add vertical distributed loads to the beam to place both braces into compression
    frame.add_member_dist_load('B1', 'Fy', -0.6, -0.6, case='D')
    frame.add_member_dist_load('B1', 'Fy', -1.5, -1.5, case='L')

    # Add a lateral load to the frame
    frame.add_node_load('N3', 'FX', 20, 'E')

    # Set up load combinations
    frame.add_load_combo('1.4D', {'D': 1.4})
    frame.add_load_combo('1.2D + 1.6L', {'D': 1.2, 'L': 1.6})
    frame.add_load_combo('1.2D + 1.0E + 1.0L', {'D': 1.2, 'E': 1.0, 'L': 1.0})

    # Render the model if this script is run directly
    if __name__ == "__main__":
        from Pynite.Visualization import Renderer
        rndr = Renderer(frame)
        rndr.combo_name = '1.2D + 1.0E + 1.0L'
        rndr.render_loads = True
        rndr.annotation_size = 1
        rndr.render_model()

    # Perform the analysis
    frame.analyze(log=True)

    # Check that both braces are removed from the model for the gravity load combinations
    for combo in ['1.4D', '1.2D + 1.6L']:
        assert frame.members['Br1'].active[combo] == False, "Br1 should be inactive for gravity load combinations"
        assert frame.members['Br2'].active[combo] == False, "Br2 should be inactive for gravity load combinations"
                                
    # Next, run checks on the lateral load combo
    for combo in ['1.2D + 1.0E + 1.0L']:

        # Check that only brace 1 is active for the lateral load combination
        assert frame.members['Br1'].active[combo] == True, "Br1 should be active for the lateral load combination"
        assert frame.members['Br2'].active[combo] == False, "Br2 should be inactive for the lateral load combination"

        # Check that the deflections along the length of the inactive brace is the linear interpolation of the member end deflections
        dy = frame.members['Br1'].deflection_array('dy', 20, combo)

        # Local transverse (y) displacements of the brace end nodes
        d = frame.members['Br1'].d(combo)
        dy_N1 = float(d[1, 0])
        dy_N4 = float(d[7, 0])
        for i in range(20):
            # Linear interpolation of the deflection along the length of the brace
            dy_interp = dy_N1 + (dy_N4 - dy_N1) * i / 19
            assert math.isclose(float(dy[1][i]), dy_interp, rel_tol=1e-5), f"Deflection at point {i} along Br1 does not match linear interpolation"

if __name__ == "__main__":
    test_tc_braced_frame()
