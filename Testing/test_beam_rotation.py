from PyNite import FEModel3D
import pytest as pt

# Tests
def test_beam_rotation():

    beam = FEModel3D()

    beam.add_node('N1', 0, 0, 0)
    beam.add_node('N2', 10, 0, 0)
    beam.def_support('N1', True, True, True, True, False, False)
    beam.def_support('N2', False, True, True, False, False, False)
    beam.add_material('Steel', 29000/12**2, 11200/12**2, 0.3, 0.490, 60)
    beam.add_section('W12x26', 7.65/12**2, 17.3/12**4, 204/12**4, 0.3/12**4)
    beam.add_member('M1', 'N1', 'N2', 'Steel', 'W12x26', 45)
    beam.add_member_dist_load('M1', 'FY', -2, -2)

    # Analyze the beam
    beam.analyze()

    # Obtain the max/min moment about each axis
    Mz = beam.members['M1'].min_moment('Mz')
    My = beam.members['M1'].max_moment('My')

    # The moment about each axis should be identical
    assert Mz == pt.approx(-17.677, rel=1e-2)
    assert My == pt.approx(17.677, rel=1e-2)
