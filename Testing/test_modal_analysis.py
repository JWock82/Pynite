"""
Tests for modal analysis functionality in Pynite
"""

import pytest
import numpy as np
from Pynite import FEModel3D


def analytical_cantilever_frequency(mode_num, L, E, I, rho, A, g):
    """
    Calculates analytical natural frequencies for a cantilever beam.

    Formula: f_n = (β_n² / (2π L²)) * sqrt(EI / (ρA))
    where β_n are the roots of cos(β)cosh(β) + 1 = 0
    """
    # First 5 roots of cos(β)cosh(β) + 1 = 0
    beta_values = [1.875104, 4.694091, 7.854757, 10.995541, 14.137168]

    if mode_num > len(beta_values):
        return None

    beta_n = beta_values[mode_num - 1]
    frequency = (beta_n**2 / (2 * np.pi * L**2)) * np.sqrt(E * I / (rho / g * A))

    return frequency


def test_cantilever_frequency():
    """Test natural frequencies of a cantilever beam against theoretical values"""

    # Create a simple cantilever beam model
    model = FEModel3D()

    # Add nodes for a cantilever beam
    # Some mass is lost to the supports. Intermediate nodes will help mitigate this issue.
    L = 20  # meters
    for i in range(11):
        model.add_node(f'N{i+1}', L*i/10, 0, 0)
        model.def_support(f'N{i+1}', True, False, True, True, True, False)

    # Fixed support at left end
    model.def_support('N1', True, True, True, True, True, True)

    # Add material and section
    E = 200e9  # Steel in Pa
    G = 80e9
    rho = 7800  # kg/m³
    A = 0.01    # m²
    I = 8.33e-6 # m⁴
    J = 1.67e-6 # m⁴

    model.add_material('Steel', E, G, 0.3, rho)
    model.add_section('BeamSection', A, I, I, J)

    # Add members
    model.add_member('M1', 'N1', 'N11', 'Steel', 'BeamSection')

    # Add member self weight
    model.add_member_self_weight('FY', -1.0, 'D1')
    for node_name in model.nodes.keys():
        if node_name == 'N1' or node_name == 'N11':
            model.add_node_load(node_name, 'FY', -A*rho*L/20, 'D2')
        else:
            model.add_node_load(node_name, 'FY', -A*rho*L/10, 'D2')

    # Add a mass load combination
    model.add_load_combo('Consistent Mass Combo', {'D1': 1.0})
    model.add_load_combo('Lumped Mass Combo', {'D2': 1.0})

    # Perform modal analysis
    model.analyze_modal(num_modes=3, mass_combo_name='Consistent Mass Combo', mass_direction='Y', gravity=9.81, log=False)
    f_consistent = model.frequencies

    model.analyze_modal(num_modes=3, mass_combo_name='Lumped Mass Combo', mass_direction='Y', gravity=9.81, log=False)
    f_lumped = model.frequencies

    # Theoretical validation
    f_analytical = [analytical_cantilever_frequency(i, L, E, I, rho, A, 9.81) for i in range(1, 4)]

    assert all(
               abs(f_num - f_theory) / f_theory < 0.05
               for f_num, f_theory in zip(f_consistent, f_analytical)
               ), "Calculated frequencies should match theoretical values within 5%"

    assert all(
               abs(f_num - f_theory) / f_theory < 0.05
               for f_num, f_theory in zip(f_lumped, f_analytical)
               ), "Calculated frequencies should match theoretical values within 5%"

def test_simple_frame_modes():
    """Test modal analysis on a simple 2D frame"""

    model = FEModel3D()

    # Create a simple portal frame
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', 5, 0, 0)  
    model.add_node('N3', 0, 3, 0)
    model.add_node('N4', 5, 3, 0)

    # Fixed bases
    model.def_support('N1', True, True, True, True, True, True)
    model.def_support('N2', True, True, True, True, True, True)

    # Add material with density
    model.add_material('Steel', 200e9, 80e9, 0.3, 7800)
    model.add_section('Column', 0.1, 8.33e-5, 8.33e-5, 1.67e-5)
    model.add_section('Beam', 0.08, 4.17e-5, 4.17e-5, 8.33e-6)

    # Add members
    model.add_member('Col1', 'N1', 'N3', 'Steel', 'Column')
    model.add_member('Col2', 'N2', 'N4', 'Steel', 'Column')
    model.add_member('Beam1', 'N3', 'N4', 'Steel', 'Beam')

    # Add self-weight
    model.add_member_self_weight('FY', 1.0)

    # Run modal analysis
    model.analyze_modal(num_modes=4, log=False)
    frequencies = model.frequencies

    assert all(freq > 0 for freq in frequencies)


def test_mass_increase():
    """Test that frequencies drop as mass increases"""

    model = FEModel3D()

    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', 5, 0, 0)

    model.def_support('N1', True, True, True, True, True, True)

    # Add material with density
    model.add_material('Steel', 200e9, 80e9, 0.3, 7800)
    model.add_section('Beam', 0.1, 8.33e-5, 8.33e-5, 1.67e-5)
    model.add_member('M1', 'N1', 'N2', 'Steel', 'Beam')

    # Add loads to the model that will be converted to mass
    model.add_member_self_weight(global_direction='FY', factor=1.0, case='Mass 1')
    model.add_node_load('N2', 'FY', -1000, case='Mass 2')  # 1000 N downward

    # Add a load combinations for mass
    model.add_load_combo('MassCombo1', {'Mass 1': 1.0})
    model.add_load_combo('MassCombo2', {'Mass 1': 1.0, 'Mass 2': 1.0})

    # Run the analysis
    model.analyze_modal(num_modes=1, mass_combo_name="MassCombo1", mass_direction='Y', log=False)
    freq1 = model.frequencies
    model.analyze_modal(num_modes=1, mass_combo_name="MassCombo2", mass_direction='Y', log=False)
    freq2 = model.frequencies

    # Frequency should be lower when additional mass is included
    assert freq2[0] < freq1[0], 'Frequencies did not drop as mass increased.'


def test_lumped_vs_consistent_mass():
    """Compare lumped and consistent mass formulations"""

    model = FEModel3D()

    # Create 11 nodes spaced 10 units apart
    num = 10
    for n in range(11):
        model.add_node(f'N{n}', (n)*10/num, 0, 0)

    # Add a fixed support to the first node
    model.def_support('N1', True, True, True, True, True, True)

    # Define material and section properties
    A = 0.1
    rho = 7800  # Unit weight
    model.add_material('Steel', 200e9, 80e9, 0.3, rho)
    model.add_section('Beam', A, 8.33e-5, 8.33e-5, 1.67e-5)

    # Create 10 members
    for n in range(10):

        # Define the member
        model.add_member(f'M{n + 1}', f'N{n}', f'N{n + 1}', 'Steel', 'Beam')

        # Calculate and apply the self-weight manually as distributed loads
        # Pynite will form a lumped mass matrix from these loads
        sw = rho*A
        model.add_member_dist_load(f'M{n + 1}', 'FY', -sw, -sw, case='Lumped')

    # Create a mass load combination
    model.add_load_combo('Lumped Mass Combo', {'Lumped': 1.0})

    # Perform modal analysis
    model.analyze_modal(num_modes=2, mass_combo_name='Lumped Mass Combo', log=False)
    freq_lumped = model.frequencies[0]

    # Apply the self-weight loads in a separate load case using the self weight command
    # Pynite creates a consistent mass matrix with the self-weight command
    model.add_member_self_weight('FY', 1, 'Consistent')

    # Create another mass load combination
    model.add_load_combo('Consistent Mass Combo', {'Consistent': 1.0})

    # Perform modal analysis
    model.analyze_modal(num_modes=2, mass_combo_name='Consistent Mass Combo', log=False)
    freq_consistent = model.frequencies[0]

    # Both should give reasonable results
    assert freq_lumped > 0
    assert freq_consistent > 0

    # They won't be identical but should be close
    ratio = freq_lumped/freq_consistent

    # Typically within 10-20% for fundamental frequency
    assert 0.8 < ratio < 1.2

@pytest.mark.parametrize("num_modes", [1, 3, 5])
def test_different_mode_counts(num_modes):
    """Test that different numbers of modes can be requested"""

    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', 5, 0, 0)
    model.def_support('N1', True, True, True, True, True, True)

    model.add_material('Steel', 200e9, 80e9, 0.3, 7800)
    model.add_section('Beam', 0.1, 8.33e-5, 8.33e-5, 1.67e-5)
    model.add_member('M1', 'N1', 'N2', 'Steel', 'Beam')
    model.add_member_self_weight('FY', 1)

    model.analyze_modal(num_modes=num_modes, log=False)

    assert len(model.frequencies) == num_modes


if __name__ == '__main__':

    test_mass_increase()
    test_cantilever_frequency()
    test_different_mode_counts(3)
    test_lumped_vs_consistent_mass()
    test_simple_frame_modes()
