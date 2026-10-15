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

    # Add the self weight of all members, so the inactive braces carry member loads
    frame.add_member_self_weight('FY', -1, 'D')

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

    # Check vertical equilibrium for each combo. This includes the self weight of the inactive braces.
    total_self_weight = sum(m.material.rho*m.section.A*m.L() for m in frame.members.values())
    dead = 0.6*15 + total_self_weight
    live = 1.5*15
    expected_FY = {'1.4D': 1.4*dead,
                   '1.2D + 1.6L': 1.2*dead + 1.6*live,
                   '1.2D + 1.0E + 1.0L': 1.2*dead + 1.0*live}
    for combo, expected in expected_FY.items():
        reaction_FY = frame.nodes['N1'].RxnFY[combo] + frame.nodes['N2'].RxnFY[combo]
        assert math.isclose(reaction_FY, expected, rel_tol=1e-6), f"Vertical reactions do not balance the applied loads for {combo}"

    # Check that only brace 1 is active for the lateral load combination (the inactive braces are checked below)
    assert frame.members['Br1'].active['1.2D + 1.0E + 1.0L'] == True, "Br1 should be active for the lateral load combination"

    def check_inactive_brace(name, combo):

        brace = frame.members[name]
        assert brace.active[combo] == False, f"{name} should be inactive for {combo}"

        # An inactive brace has no stiffness, so its end forces are just its fixed end reactions from self weight
        fer = brace.fer(combo)
        shear = brace.shear_array('Fy', 5, combo)
        moment = brace.moment_array('Mz', 5, combo)
        axial = brace.axial_array(5, combo)
        assert math.isclose(shear[1][0], fer[1, 0], rel_tol=1e-6), f"{name} i-end shear should equal its fixed end reaction for {combo}"
        assert math.isclose(moment[1][0], fer[5, 0], abs_tol=1e-9), f"{name} i-end moment should equal its fixed end reaction for {combo}"
        assert math.isclose(axial[1][0], fer[0, 0], rel_tol=1e-6), f"{name} i-end axial force should equal its fixed end reaction for {combo}"

        # The braces are pinned at both ends, so the moment diagram must close to zero at both ends
        assert math.isclose(moment[1][-1], 0, abs_tol=1e-9), f"{name} j-end moment should be zero for {combo}"

        # Deflections are a straight line between the end node displacements
        d = brace.d(combo)
        for direction, i_dof, j_dof in [('dx', 0, 6), ('dy', 1, 7), ('dz', 2, 8)]:
            d_i = float(d[i_dof, 0])
            d_j = float(d[j_dof, 0])
            d_arr = brace.deflection_array(direction, 5, combo)
            for k in range(5):
                assert math.isclose(d_arr[1][k], d_i + (d_j - d_i)*k/4, rel_tol=1e-6, abs_tol=1e-12), f"{name} {direction} deflection is not linear for {combo}"

    # Check the internal results of the inactive braces in each combo they are inactive for
    check_inactive_brace('Br1', '1.4D')
    check_inactive_brace('Br2', '1.4D')
    check_inactive_brace('Br1', '1.2D + 1.6L')
    check_inactive_brace('Br2', '1.2D + 1.6L')
    check_inactive_brace('Br2', '1.2D + 1.0E + 1.0L')

if __name__ == "__main__":
    test_tc_braced_frame()
