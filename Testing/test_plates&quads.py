# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

from Pynite import FEModel3D
import math

def test_plate_displacement():
    """
    # A First Course in the Finite Element Method, 4th Edition
    # Daryl L. Logan
    # Example 12.1
    # Units for this model are pounds and inches
    """

    plate_model= FEModel3D()

    plate_model.add_node('N1', 0, 0, 0)
    plate_model.add_node('N2', 10, 0, 0)
    plate_model.add_node('N3', 20, 0, 0)
    plate_model.add_node('N4', 0, 10, 0)
    plate_model.add_node('N5', 10, 10, 0)
    plate_model.add_node('N6', 20, 10, 0)
    plate_model.add_node('N7', 0, 20, 0)
    plate_model.add_node('N8', 10, 20, 0)
    plate_model.add_node('N9', 20, 20, 0)

    plate_model.add_material('Steel', 30000000, 11200000, 0.3, 0.284)

    plate_model.add_plate('P1', 'N1', 'N2', 'N5', 'N4', 0.1, 'Steel')
    plate_model.add_plate('P2', 'N2', 'N3', 'N6', 'N5', 0.1, 'Steel')
    plate_model.add_plate('P3', 'N4', 'N5', 'N8', 'N7', 0.1, 'Steel')
    plate_model.add_plate('P4', 'N5', 'N6', 'N9', 'N8', 0.1, 'Steel')

    plate_model.add_node_load('N5', 'FZ', -100)

    plate_model.def_support('N1', True, True, True, True, True, True)
    plate_model.def_support('N2', True, True, True, True, True, True)
    plate_model.def_support('N3', True, True, True, True, True, True)
    plate_model.def_support('N4', True, True, True, True, True, True)
    plate_model.def_support('N6', True, True, True, True, True, True)
    plate_model.def_support('N7', True, True, True, True, True, True)
    plate_model.def_support('N8', True, True, True, True, True, True)
    plate_model.def_support('N9', True, True, True, True, True, True)

    plate_model.def_support('N5', True, True, False, False, False, True)

    # Check to see if the stiffness matrix for each plate is symmetric
    # print(allclose(plate_model.plates[0].K(), plate_model.plates[0].K().T))
    # print(allclose(plate_model.plates[1].K(), plate_model.plates[1].K().T))
    # print(allclose(plate_model.plates[2].K(), plate_model.plates[2].K().T))
    # print(allclose(plate_model.plates[3].K(), plate_model.plates[3].K().T))

    # Check to see if the global stiffness matrix is symmetric
    # print(allclose(plate_model.K(Renumber=True), plate_model.K(Renumber=False).T))

    plate_model.analyze(check_statics=True, sparse=False)
    # Test: displacement of N5 in Z direction
    calculated_displacement = plate_model.nodes['N5'].DZ['Combo 1']
    expected_displacement = -0.0861742424242424

    # Check that results are withing 1% of the expected value
    assert abs(calculated_displacement/expected_displacement - 1.0) < 0.01

def test_hydrostatic_plate():

    # Create the model
    plate_model = FEModel3D()

    # Define geometry
    t = 1  # ft
    mesh_size = 1  # ft
    a = 10  # ft
    b = 15  # ft

    # Define a material
    E = 57000*math.sqrt(4500)*12**2  # psf
    G = 0.4*E  # psf
    nu = 1/6
    rho = 150  # psf
    plate_model.add_material('Concrete', E, G, nu, rho)
    
    # Generate the mesh of plates
    plate_model.add_rectangle_mesh('MSH1', mesh_size, a, b, t, 'Concrete', kx_mod=1, ky_mod=1,
                                    element_type='Rect')
    
    # Generate the mesh
    plate_model.meshes['MSH1'].generate()

    # Add supports to the sides and base of the wall
    for node in plate_model.nodes.values():
        if node.X == 0 or node.X == a or node.Y == 0:
            plate_model.def_support(node.name, True, True, True, True, True, True)
    
    # Add hydrostatic loads to the elements
    for element in plate_model.plates.values():
        Yavg = (element.i_node.Y + element.j_node.Y + element.m_node.Y + element.n_node.Y)/4
        p = 62.4*(b - Yavg)
        plate_model.add_plate_surface_pressure(element.name, p, 'Hydrostatic')
    
    # Add a load combination to the model
    plate_model.add_load_combo('F', {'Hydrostatic': 1.0})
    
    # Analyze the model
    plate_model.analyze()

    # Get the maximum deflection in the model at the top of the wall
    DZ_calcd = max([node.DZ['F'] for node in plate_model.nodes.values() if node.Y == b])
    
    # Find the maximum deflection at the top of the wall from Timoshenko's Table 45
    q = 62.4*b
    D = E*t**3/(12*(1 - nu**2))
    DZ_expected = 0.00042*q*a**4/D

    # Check that the Pynite calculated values are within 15% of the Timoshenko calculated
    # values.
    assert abs(DZ_calcd/DZ_expected - 1) < 0.15, 'Failed Timoshenko rectangle hydrostatic test.'

def test_hydrostatic_quad():

    # Create the model
    quad_model = FEModel3D()

    # Define geometry
    t = 1  # ft
    mesh_size = 1  # ft
    a = 10  # ft
    b = 15  # ft

    # Define a material
    E = 57000*math.sqrt(4500)*12**2  # psf
    G = 0.4*E  # psf
    nu = 1/6
    rho = 150  # psf
    quad_model.add_material('Concrete', E, G, nu, rho)
    
    # Generate the mesh of quads
    quad_model.add_rectangle_mesh('MSH1', mesh_size, a, b, t, 'Concrete', kx_mod=1, ky_mod=1,
                                    element_type='Quad')
    quad_model.meshes['MSH1'].generate()

    # Add supports to the sides and base of the wall
    for node in quad_model.nodes.values():
        if node.X == 0 or node.X == a or node.Y == 0:
            quad_model.def_support(node.name, True, True, True, True, True, True)
    
    # Add hydrostatic loads to the elements
    for element in quad_model.quads.values():
        Yavg = (element.i_node.Y + element.j_node.Y + element.m_node.Y + element.n_node.Y)/4
        p = 62.4*(b - Yavg)
        quad_model.add_quad_surface_pressure(element.name, p, 'Hydrostatic')
    
    # Add a load combination to the model
    quad_model.add_load_combo('F', {'Hydrostatic': 1.0})
    
    # Analyze the model
    quad_model.analyze()

    # Get the maximum deflection in the model at the top of the wall
    DZ_calcd = max([node.DZ['F'] for node in quad_model.nodes.values() if node.Y == b])
    
    # Find the maximum deflection at the top of the wall from Timoshenko's Table 45
    q = 62.4*b
    D = E*t**3/(12*(1 - nu**2))
    DZ_expected = 0.00042*q*a**4/D

    # Check that the Pynite calculated values are within 15% of the Timoshenko calculated
    # values.
    assert abs(DZ_calcd/DZ_expected - 1) < 0.15, 'Failed Timoshenko quadrilateral hydrostatic test.'

def test_circular_hopper():
    """This test currently only tests that the code runs with no exceptions. It doesn't necessarily mean it's checking stresses correctly."""

    model = FEModel3D()
    model.add_material('Steel', 29000*144, 11200*144, 0.3, 0.490)
    model.add_frustrum_mesh('Hopper', 1, 10, 2, 10, 0.25, 'Steel')
    model.meshes['Hopper'].generate()

    # Add supports at the top/springline of the hopper
    for node in model.nodes.values():
        if math.isclose(node.Y, 0, abs_tol=0.01):
            model.def_support(node.name, True, True, True, False, False, False)
    
    # Add loads to each quad
    for quad in model.quads.values():
        model.add_quad_surface_pressure(quad.name, -0.2)

    # Solve the model
    model.analyze_linear()

    from Pynite.Visualization import Renderer
    rndr = Renderer(model)
    rndr.annotation_size = 0.25
    rndr.color_map = 'Sy'
    rndr.render_loads = False
    rndr.render_model(interact=False) # interact=False prevents the test from stalling in a headless environment
    
    # from Pynite.Rendering import Renderer
    # rndr = Renderer(model)
    # rndr.annotation_size = 0.25
    # rndr.color_map = 'Sy'
    # rndr.render_loads = False
    # rndr.render_model()

if __name__ == '__main__':
    test_hydrostatic_quad()