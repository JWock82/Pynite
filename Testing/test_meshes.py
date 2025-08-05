from Pynite import FEModel3D, ShearWall
from math import isclose


def test_rect_mesh():

    model = FEModel3D()

    E = 57*(5000)**0.5*12**2
    G = 0.4*E
    nu = 0.15
    rho = 0.150
    model.add_material('Concrete', E, G, nu, rho)

    model.add_rectangle_mesh('MSH1', 1, 10, 15, 0.667, 'Concrete', 1, 0.7, (0, 0, 0), 'XY', element_type='Quad')

    model.meshes['MSH1'].generate()

    bott_nodes = [node for node in model.nodes.values() if isclose(node.Y, 0.0)]
    top_nodes = [node for node in model.nodes.values() if isclose(node.Y, 15.0)]
    num_top_nodes = len(top_nodes)
    node_shear = 1000/num_top_nodes

    for node in top_nodes:
        model.add_node_load(node.name, 'FX', node_shear, 'E')

    for node in bott_nodes:
        model.def_support(node.name, True, True, True, True, True, True)

    model.add_load_combo('1.0E', {'E': 1.0})

    model.analyze()

    delta_max = 0.0
    for node in top_nodes:

        if node.DX['1.0E'] > delta_max:

            delta_max = node.DX['1.0E']

    k = 1000/(delta_max*12)

    print(f'Stiffness, k = {k} kip/in')

    # Check that the stiffness is as expected
    assert round(k) == 1369, 'Failed rectangular mesh stiffness test.'


def test_shear_wall():

    shear_wall = ShearWall()

    shear_wall.L = 10
    shear_wall.H = 15
    shear_wall.ky_mod = 0.70

    E = 57*(5000)**0.5*12**2
    G = 0.4*E
    nu = 0.15
    rho = 0.150
    t = 0.667
    shear_wall.add_material('Concrete', E, G, nu, rho, t)
    shear_wall.add_support(0.0)
    shear_wall.add_story('Roof', 15.0)
    shear_wall.add_shear('Roof', 1000, 'E')

    shear_wall.add_load_combo('1.0E', {'E': 1.0})

    shear_wall.generate()
    shear_wall.model.analyze()

    k = shear_wall.stiffness('Roof')/12

    print(f'Stiffness, k = {k} kip/in')

    # shear_wall.render(combo_name='Stiffness: Roof')

    # Check that the stiffness is as expected
    assert round(k) == 1369, 'Failed shear wall stiffness test.'


def test_PCA_7_quad():
    """
    Tests against the example from Section 7 of "Circular Concrete Tanks
    without Prestressing" by PCA.
    """

    # Create a new finite element model
    tank_model = FEModel3D()

    H = 20     # Tank wall height (ft)
    D = 54     # Tank inside diameter (ft)
    R = D/2    # Tank inside radius (ft)
    t = 10/12  # Tank wall thickness (ft)

    w = 62.5   # Liquid unit weight (pcf)

    fc = 4000                    # Concrete compressive strength (psi)
    E = 57000*(fc)**0.5*(12**2)  # Concrete modulus of elasticity (psf)
    nu = 0.25  # 0.17            # Poisson's ratio for concrete
    tank_model.add_material('Concrete', E, 0.4*E, nu, 150)

    mesh_size = 1       # Desired mesh size (ft)
    center = [0, 0, 0]  # Origin (X, Y, Z)
    axis = 'Y'          # Axis of revolution

    tank_model.add_cylinder_mesh('MSH1', mesh_size, R, H, t, 'Concrete', 1, 1, center, axis, element_type='Quad')
    tank_model.meshes['MSH1'].generate()

    # Add hydrostatic loads to the elements
    for element in tank_model.quads.values():

        avg_Y = (element.i_node.Y + element.j_node.Y
                + element.m_node.Y + element.n_node.Y)/4

        p = (H - avg_Y)*w

        tank_model.add_quad_surface_pressure(element.name, p)

    # Add fixed supports to the base
    for node in tank_model.nodes.values():
        if node.Y == 0:
            tank_model.def_support(node.name, True, True, True, True, True, True)

    # Analyze the model
    tank_model.analyze()

    # Max/min moment and max hoop tension as determined by PCA.
    My_max_PCA = 14804/1.3/1.7
    My_min_PCA = -3756/1.3/1.7
    Sx_PCA = 55945/1.3/1.7

    # From Timoshenko Section 117 (p. 485)
    # The Timoshenko solution yields similar results to the PCA solution, but with a slightly larger margin of error
    beta = (3*(1 - nu**2)/(R**2*t**2))**0.25  # Equation 275
    My_max_Tim = (1 - 1/(beta*H))*w*R*H*t/(12*(1 - nu**2))**0.5
    Qy_max_Tim = -(w*R*H*t)/(12*(1 - nu**2))**0.5*(2*beta - 1/H)

    # Find the max/min moments at the top of any element
    My_min = -max([element.moment(0, 1)[1, 0] for element in tank_model.quads.values()])
    My_max = -min([element.moment(0, 1)[1, 0] for element in tank_model.quads.values()])

    # Find the max hoop tension at the center of any element
    Sx = max([element.membrane(0, 0)[0, 0] for element in tank_model.quads.values()])*t

    # Find the maximum reaction at the base of the tank
    RMy = max([node.RxnMX['Combo 1'] for node in tank_model.nodes.values()])/mesh_size

    # Check that the Pynite calculated values are within 3% of the calculated PCA values.
    assert abs(1 - My_max/My_max_PCA) < 0.03, 'Failed quad cylinder flexure test.'
    assert abs(1 - RMy/My_max_PCA) < 0.03, 'Failed quad cylinder flexure test.'
    assert abs(1 - My_min/My_min_PCA) < 0.03, 'Failed quad cylinder flexure test.'

    # Check the sign convention for local y-axis bending
    assert My_max > 0, 'Failed quad cylinder sign convention test'

    # Check the expected hoop tension
    assert abs(1 - Sx/20000) < 0.03, 'Failed quad cylinder hoop tension test.'

    # Render the model
    # from Pynite.Rendering import Renderer
    # rndr = Renderer(tank_model)
    # rndr.annotation_size = 0.25
    # rndr.color_map = 'Sx'
    # rndr.render_loads = False
    # rndr.render_model()


def test_PCA_7_rect():
    """
    Tests against the example from Section 7 of "Circular Concrete Tanks without Prestressing" by PCA.
    """

    # Create a new finite element model
    tank_model = FEModel3D()

    H = 20     # Tank wall height (ft)
    D = 54     # Tank inside diameter (ft)
    R = D/2    # Tank inside radius (ft)
    t = 10/12  # Tank wall thickness (ft)

    w = 62.5   # Liquid unit weight (pcf)

    fc = 4000                    # Concrete compressive strength (psi)
    E = 57000*(fc)**0.5*(12**2)  # Concrete modulus of elasticity (psf)
    nu = 0.25  # 0.17            # Poisson's ratio for concrete
    tank_model.add_material('Concrete', E, 0.4*E, nu, 150)

    mesh_size = 2       # Desired mesh size (ft)
    center = [0, 0, 0]  # Origin (X, Y, Z)
    axis = 'Y'          # Axis of revolution

    # Add a cylinder mesh to the model
    tank_model.add_cylinder_mesh('MSH1', mesh_size, R, H, t, 'Concrete', 1, 1, center, axis, element_type='Rect')

    # Generate the mesh prior to running so we can work with it
    tank_model.meshes['MSH1'].generate()

    # Add hydrostatic loads to the elements
    for element in tank_model.plates.values():

        avg_Y = (element.i_node.Y + element.j_node.Y
                + element.m_node.Y + element.n_node.Y)/4

        p = (H - avg_Y)*w

        tank_model.add_plate_surface_pressure(element.name, p)

    # Add fixed supports to the base
    for node in tank_model.nodes.values():
        if node.Y == 0:
            tank_model.def_support(node.name, True, True, True, True, True, True)

    # Analyze the model
    tank_model.analyze()

    # Max/min moment and max hoop tension as determined by PCA.
    My_max_PCA = 14804/1.3/1.7
    My_min_PCA = -3756/1.3/1.7
    Sx_PCA = 55945/1.3/1.7

    # From Timoshenko Section 117 (p. 485)
    # The Timoshenko solution yields similar results to the PCA solution
    beta = (3*(1 - nu**2)/(R**2*t**2))**0.25  # Equation 275
    My_max_Tim = (1 - 1/(beta*H))*w*R*H*t/(12*(1 - nu**2))**0.5
    Qy_max_Tim = -(w*R*H*t)/(12*(1 - nu**2))**0.5*(2*beta - 1/H)

    My_max = tank_model.meshes['MSH1'].max_moment('My')
    My_min = tank_model.meshes['MSH1'].min_moment('My')
    Sx = max([element.membrane(element.width()/2, element.height()/2)[0, 0] for element in tank_model.plates.values()])*t

    # Check that the Pynite calculated values are within 8% of expected values. With a finer mesh the results are known to converge even closer, but 8% allows the model to run faster.
    assert abs(1 - My_max/My_max_PCA) < 0.08, 'Failed plate cylinder flexure test.'
    assert abs(1 - My_min/My_min_PCA) < 0.08, 'Failed plate cylinder flexure test.'
    assert My_max > 0, 'Failed plate cylinder sign convention test'
    assert abs(1 - Sx/20000) < 0.05, 'Failed plate cylinder hoop tension test.'


if __name__ == '__main__':

    # test_rect_mesh()
    test_shear_wall()
    # test_PCA_7_quad()
