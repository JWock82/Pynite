from Pynite import FEModel3D
import math
# import matplotlib, matplotlib.pyplot as plt

def create_3_support_model():
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

def test_3_support_beam_moments():
    
    model = create_3_support_model()
    member = model.members['M1']
    assert math.isclose(member.max_moment('Mz', '1.0D'), 18000.0), 'Incorrect max moment during continuous beam test'
    assert math.isclose(max(member.moment_array('Mz', 11, '1.0D')[1]), 18000.00), 'Incorrect max moment in moment array during continuous beam test.'

def test_2_support_beam_moments():

    model = FEModel3D()

    model.add_material("default", 1, 1, 1, 1, 1)
    model.add_section("default", 1, 1, 1, 1)

    model.add_node("0", 0, 0, 0)
    model.add_node("1", 10, 0, 0)
    model.add_node("2", 13, 0, 0)

    model.def_support("0", True, True, True, True, True, False)
    model.def_support("1", False, True, False, False, False, False)

    model.add_member("M0", "0", "2", "default", "default")
    model.add_member_dist_load("M0", "Fy", -10, -10, 0, 13, case='load')
    model.add_load_combo("combo", {"load": 1.0})

    model.analyze(log=True, check_statics=True)

    member = model.members['M0']
    fe_max_moment = member.max_moment("Mz", "combo")
    fe_min_moment = member.min_moment("Mz", "combo")
    # member.plot_moment("Mz", "combo", 2000)

    # Analytic solution
    w = -10
    l = 10
    a = 3

    # Beam formula source: Handbook of Steel Construction, CISC, 11th Ed. Pg. 5-138 diag. 24
    analytic_min_moment = (w * (l + a)**2 * (l - a)**2 ) / (8 * l**2)
    analytic_max_moment = -w * a**2 / 2

    try:
        assert math.isclose(fe_max_moment, analytic_max_moment)
    except AssertionError:
        print(fe_max_moment, analytic_max_moment)
    assert math.isclose(fe_min_moment, analytic_min_moment)

if __name__ == '__main__':
    pass
    # test_2_support_beam_moments()
    # test_continuous_beam_moments()
    # test_plots()
