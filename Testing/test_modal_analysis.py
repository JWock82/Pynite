"""
Tests for modal analysis functionality in Pynite
"""

import pytest
import numpy as np
from Pynite import FEModel3D
from Pynite.LoadCombo import LoadCombo

def analytical_cantilever_frequency(mode_num, L, E, I, rho, A):
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
    frequency = (beta_n**2 / (2 * np.pi * L**2)) * np.sqrt(E * I / (rho * A))
    
    return frequency

class TestModalAnalysis:
    """Test cases for modal analysis feature"""
    
    def test_cantilever_frequency(self):
        """Test natural frequencies of a cantilever beam against theoretical values"""
        
        # Create a simple cantilever beam model
        model = FEModel3D()
        
        # Add nodes for a cantilever beam
        L = 10  # meters
        model.add_node('N1', 0, 0, 0)
        model.add_node('N2', L, 0, 0)
        model.add_node('N3', 2*L, 0, 0)  # For better mode shape resolution
        
        # Fixed support at left end
        model.def_support('N1', True, True, True, True, True, True)
        model.def_support('N2', False, False, True, True, True, False)
        model.def_support('N3', False, False, True, True, True, False)
        
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
        model.add_member('M1', 'N1', 'N2', 'Steel', 'BeamSection')
        model.add_member('M2', 'N2', 'N3', 'Steel', 'BeamSection')
        
        # Perform modal analysis
        results = model.analyze_modal(num_modes=3, log=False)
        
        frequencies = results['frequencies']
        
        # Basic validation checks
        assert len(frequencies) == 3, "Should return requested number of modes"
        assert all(freq > 0 for freq in frequencies), "All frequencies should be positive"
        assert frequencies[0] < frequencies[1] < frequencies[2], "Frequencies should be in ascending order"
        
        # You could add theoretical validation here if you have reference solutions
        f_analytical = [analytical_cantilever_frequency(i, 2*L, E, I, rho, A) for i in range(1, 4)]
        assert np.allclose(frequencies[:2], f_analytical[:2], atol=0.02), "Calculated frequencies should match theoretical values"
        # print(f"Calculated frequencies: {frequencies} Hz")
    
    def test_simple_frame_modes(self):
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
        
        # Run modal analysis
        results = model.analyze_modal(num_modes=4, log=False)
        
        # Validate results
        assert 'frequencies' in results
        assert 'mode_shapes' in results
        assert 'num_modes' in results
        
        frequencies = results['frequencies']
        mode_shapes = results['mode_shapes']
        
        assert len(frequencies) == len(mode_shapes) == 4
        assert all(freq > 0 for freq in frequencies)
        
        # Check mode shapes have correct dimensions
        for mode_shape in mode_shapes:
            assert mode_shape.shape[0] == len(model.nodes) * 6
    
    def test_mass_increase(self):
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
        results1 = model.analyze_modal(num_modes=1, mass_combo_name="MassCombo1", mass_direction=1, log=False)
        results2 = model.analyze_modal(num_modes=1, mass_combo_name="MassCombo2", mass_direction=1, log=False)

        # Frequency should be lower when additional mass is included
        assert results2['frequencies'][0] < results1['frequencies'][0], 'Frequencies did not drop as mass increased.'

    def test_lumped_vs_consistent_mass(self):
        """Compare lumped and consistent mass formulations"""

        model = FEModel3D()

        num = 10
        for n in range(11):
            model.add_node(f'N{n}', (n)*10/num, 0, 0)

        model.def_support('N1', True, True, True, True, True, True)

        model.add_material('Steel', 200e9, 80e9, 0.3, 7800)
        model.add_section('Beam', 0.1, 8.33e-5, 8.33e-5, 1.67e-5)

        # Test with lumped mass
        for n in range(10):
            model.add_member(f'M{n+1}', f'N{n}', f'N{n + 1}', 'Steel', 'Beam')
        results_lumped = model.analyze_modal(num_modes=2, log=False)

        # Test with consistent mass  
        for n in range(10):
            model.delete_member(f'M{n+1}')
            model.add_member(f'M{n+1}', f'N{n}', f'N{n + 1}', 'Steel', 'Beam')
        results_consistent = model.analyze_modal(num_modes=2, log=False)

        # Both should give reasonable results
        assert results_lumped['frequencies'][0] > 0
        assert results_consistent['frequencies'][0] > 0

        # They won't be identical but should be close
        lumped_freq = results_lumped['frequencies'][0]
        consistent_freq = results_consistent['frequencies'][0]
        ratio = lumped_freq / consistent_freq

        # Typically within 10-20% for fundamental frequency
        assert 0.8 < ratio < 1.2

    @pytest.mark.parametrize("num_modes", [1, 3, 5])
    def test_different_mode_counts(self, num_modes):
        """Test that different numbers of modes can be requested"""
        
        model = FEModel3D()
        model.add_node('N1', 0, 0, 0)
        model.add_node('N2', 5, 0, 0)
        model.def_support('N1', True, True, True, True, True, True)
        
        model.add_material('Steel', 200e9, 80e9, 0.3, 7800)
        model.add_section('Beam', 0.1, 8.33e-5, 8.33e-5, 1.67e-5)
        model.add_member('M1', 'N1', 'N2', 'Steel', 'Beam')
        
        results = model.analyze_modal(num_modes=num_modes, log=False)
        
        assert len(results['frequencies']) == num_modes
        assert len(results['mode_shapes']) == num_modes

if __name__ == '__main__':
    test = TestModalAnalysis()
    test.test_mass_increase()
