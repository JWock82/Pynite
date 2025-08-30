from Pynite import FEModel3D
import pytest as pt

# Tests
def test_beam_rotation():

    # Create a model
    beam = FEModel3D()

    # Add nodes to the model
    beam.add_node('N1', 0, 0, 0)
    beam.add_node('N2', 10, 0, 0)

    # Define supports
    beam.def_support('N1', True, True, True, True, False, False)
    beam.def_support('N2', False, True, True, False, False, False)

    # Add a material and section properties
    beam.add_material('Steel', 29000/12**2, 11200/12**2, 0.3, 0.490, 60)
    beam.add_section('W12x26', 7.65/12**2, 17.3/12**4, 204/12**4, 0.3/12**4)

    # Create the member and add a distributed load
    # The beam will be rotated 45 degrees about its x-axis
    beam.add_member('M1', 'N1', 'N2', 'Steel', 'W12x26', 45)
    beam.add_member_dist_load('M1', 'FY', -2, -2)

    # Analyze the beam
    beam.analyze()

    # Obtain the max/min moment about each axis
    Mz_max = beam.members['M1'].max_moment('Mz')
    Mz_min = beam.members['M1'].min_moment('Mz')
    My_max = beam.members['M1'].max_moment('My')
    My_min = beam.members['M1'].min_moment('My')

    # Check the maximum and minimum moments
    # Note that the strong and weak axis have different positive moment sign conventions due to the right hand rule
    assert Mz_max == pt.approx(0.0, rel=1e-2)
    assert Mz_min == pt.approx(-17.677, rel=1e-2)
    assert My_max == pt.approx(17.677, rel=1e-2)
    assert My_min == pt.approx(0.0, rel=1e-2)

if __name__ == '__main__':
    test_beam_rotation()
