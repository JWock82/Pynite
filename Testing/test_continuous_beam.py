from Pynite import FEModel3D
import math
# import matplotlib, matplotlib.pyplot as plt

def create_model():

    continuous_beam = FEModel3D()

    continuous_beam.add_node('N1', 0, 0, 0)
    continuous_beam.add_node('N2', 10*12, 0, 0)
    continuous_beam.add_node('N3', 20*12, 0, 0)

    continuous_beam.add_material("Material", 29000, 7700, 0.3, 0.001)
    continuous_beam.add_section("Section", A=20, Iy=100, Iz=150, J=250)

    continuous_beam.add_member('M1', 'N1', 'N3', "Material", "Section")

    continuous_beam.def_support('N1', True, True, True, True, True, False)  # Constrained for torsion at 'N1'
    continuous_beam.def_support('N2', False, True, False, False, False, False) # Not constrained for torsion at 'N2'
    continuous_beam.def_support('N3', False, True, False, False, False, False) # Same for 'N3'

    continuous_beam.add_member_dist_load("M1", "Fy", -10, -10, 0, 20*12, case="D")

    continuous_beam.add_load_combo('1.0D', {'D':1.0})
    continuous_beam.analyze()

    return continuous_beam

def test_continuous_beam_moments():
    model = create_model()
    member = model.members['M1']
    assert math.isclose(member.max_moment('Mz', '1.0D'), 18000.0), 'Incorrect max moment during continuous beam test'
    assert math.isclose(max(member.moment_array('Mz', 11, '1.0D')[1]), 18000.00), 'Incorrect max moment in moment array during continuous beam test.'

# def test_plots():

#     model = create_model()

#     assert isinstance(model.members['M1'].plot_shear('Fy', '1.0D', 100), matplotlib.figure.Figure)
#     plt.close()

#     assert isinstance(model.members['M1'].plot_moment('Mz', '1.0D', 100), matplotlib.figure.Figure)
#     plt.close()
    
#     assert isinstance(model.members['M1'].plot_torque('1.0D', 100), matplotlib.figure.Figure)
#     plt.close()
    
#     assert isinstance(model.members['M1'].plot_axial('1.0D', 100), matplotlib.figure.Figure)
#     plt.close()
    
#     assert isinstance(model.members['M1'].plot_deflection('dy', '1.0D', 100), matplotlib.figure.Figure)
#     plt.close()

if __name__ == '__main__':
    test_continuous_beam_moments()
    # test_plots()