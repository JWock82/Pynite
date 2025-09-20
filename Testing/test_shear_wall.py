from Pynite import FEModel3D
from Pynite.Rendering import Renderer
from math import isclose


def test_shear_wall():

    # Create a new finite element model
    model = FEModel3D()

    # Define a material for our shear wall
    fm = 2000/1000*144  # ksf
    Em = 900*fm  # ksf
    Gm = 0.4*Em  # ksf
    nu = 0.17
    rho_m = 0.140  # kcf
    model.add_material('CMU', Em, Gm, nu, rho_m)

    # Create a wall in the 'XY' plane and the 'YZ' plane
    for wall_name in ['W1', 'W2']:

        if wall_name == 'W1':
            plane = 'XY'
            origin = [0, 0, 0]
        else:
            plane = 'YZ'
            origin = [40, 0, 0]

        # Add a new shear wall to the model
        model.add_shear_wall(wall_name, mesh_size=1, length=26, height=16, thickness=0.667, material_name='CMU', ky_mod=1, plane=plane, origin=origin)

        # Add openings to the wall
        # Pynite gives the option to add a tie's stiffness over the opening to simulate collectors/drag-struts tying the wall segments together. Since none of the openings defined below rise to the top of the wall, ties are not needed.
        model.shear_walls[wall_name].add_opening('Door 1', x_start=2, y_start=0, width=4, height=12, tie=None)
        model.shear_walls[wall_name].add_opening('Window 1', x_start=8, y_start=8, width=4, height=4, tie=None)
        model.shear_walls[wall_name].add_opening('Window 2', x_start=14, y_start=8, width=4, height=4, tie=None)
        model.shear_walls[wall_name].add_opening('Door 2', x_start=20, y_start=0, width=4, height=12, tie=None)

        # Add support across the entire base of the wall
        model.shear_walls[wall_name].add_support(elevation=0, x_start=0, x_end=26)

        # Add a story to the shear wall where loads will be applied
        # Note that `x_start` and `x_end` are optional arguments. If they are omitted, the story's length will default to the full length of the wall. `x_start` and `x_end` can be used to simulate loading from a partial depth diaphragm.
        model.shear_walls[wall_name].add_story('Roof', elevation=16, x_start=0, x_end=26)

        # Add a seismic shear force of 100 kips to the roof
        model.shear_walls[wall_name].add_shear(story_name='Roof', force=100, case='E')

        model.shear_walls[wall_name].add_flange(thickness=0.667, width=0.75*16, x=0, y_start=0, y_end=16, material='CMU', side='-z')
        model.shear_walls[wall_name].add_flange(thickness=0.667, width=0.75*16, x=7, y_start=0, y_end=16, material='CMU', side='+z')
        model.shear_walls[wall_name].add_flange(thickness=0.667, width=0.75*16, x=26, y_start=0, y_end=16, material='CMU', side='-z')

    # Add a load combination to the model
    model.add_load_combo('1.0E', {'E': 1.0}, combo_tags='strength')

    # Analyze the model. Use the linear solver for greater speed
    model.analyze_linear(log=True, check_statics=True)

    # Find the stiffness of the shear wall when a shear is applied to the roof (usefull for rigid diaphragm analysis)
    k = model.shear_walls[wall_name].stiffness('Roof')

    # Check that the calculated stiffness is within 5% of the expected value
    # assert abs(1 - k/k_expected) <= 0.05

    # rndr = Renderer(model)
    # rndr.annotation_size = 0.25
    # rndr.combo_name = '1.0E'
    # rndr.color_map = 'Txy'
    # rndr.scalar_bar = True
    # rndr.render_loads = True
    # rndr.deformed_shape = True
    # rndr.deformed_scale = 300
    # rndr.render_model()

def test_quad_shear_wall():

    sw = FEModel3D()

    # Add a material
    E = 57000*(4000)**0.5/1000*12**2
    nu = 0.17
    G = E/(2*(1 + nu))
    sw.add_material('Concrete', E, G, nu, 0.150)

    # Define section properties
    t = 1
    L = 10
    H = 20
    A = L*t
    I = t*L**3/12

    mesh_size = 1
    sw.add_rectangle_mesh('MSH1', mesh_size, L, H, t, 'Concrete', element_type='Quad')
    sw.meshes['MSH1'].generate()

    V = 1000
    for node in sw.nodes.values():
        if node.Y == 0:
            sw.def_support(node.name, True, True, True, True, True, True)
        elif node.Y == H:
            sw.add_node_load(node.name, 'FX', V/11)

    sw.analyze()

    # Calculated solution
    delta1 = max([node.DX['Combo 1'] for node in sw.nodes.values()])

    # Theoretical solution
    delta2 = V*H**3/(3*E*I) + 1.2*V*H/(G*A)

    # Check that the solution matches the theoretical solution within 0.1%
    assert abs(1 - delta1/delta2) < 0.001, 'Failed quad shear wall test.'

    rndr = Renderer(sw)
    rndr.annotation_size = 5
    rndr.window_width = 700
    rndr.window_height = 700
    rndr.plotter.off_screen = True


def test_rect_shear_wall():

    sw = FEModel3D()

    # Add a material
    E = 57000*(4000)**0.5/1000*12**2
    nu = 0.17
    G = E/(2*(1 + nu))
    sw.add_material('Concrete', E, G, nu, 0.150)

    # Define section properties
    t = 1
    L = 10
    H = 20
    A = L*t
    I = t*L**3/12

    mesh_size = 1
    sw.add_rectangle_mesh('MSH1', mesh_size, L, H, t, 'Concrete', element_type='Rect')
    sw.meshes['MSH1'].generate()

    V = 1000
    for node in sw.nodes.values():
        if node.Y == 0:
            sw.def_support(node.name, True, True, True, True, True, True)
        elif node.Y == H:
            sw.add_node_load(node.name, 'FX', V/11)

    sw.analyze()

    # Calculated solution
    delta1 = max([node.DX['Combo 1'] for node in sw.nodes.values()])

    # Theoretical solution
    delta2 = V*H**3/(3*E*I) + 1.2*V*H/(G*A)

    # Check that the solution matches the theoretical solution within 0.1%
    assert abs(1 - delta1/delta2) < 0.001, 'Failed rect plate shear wall test.'


def test_cracked_rect_shear_wall():

    sw = FEModel3D()

    # Define a material
    E = 57000*(4000)**0.5/1000*12**2
    nu = 0.17
    G = E/(2*(1 + nu))
    sw.add_material('Concrete', E, G, nu, 0.150)

    # Define geometry and section properties
    t = 1
    L = 10
    H = 20
    A = L*t
    I = 0.35*t*L**3/12

    mesh_size = 1
    sw.add_rectangle_mesh('MSH1', mesh_size, L, H, t, 'Concrete', ky_mod=0.35, element_type='Rect')
    sw.meshes['MSH1'].generate()

    V = 1000
    for node in sw.nodes.values():
        if node.Y == 0:
            sw.def_support(node.name, True, True, True, True, True, True)
        elif node.Y == H:
            sw.add_node_load(node.name, 'FX', V/11)

    sw.analyze()

    # Calculated solution
    delta1 = max([node.DX['Combo 1'] for node in sw.nodes.values()])

    # Theoretical solution
    delta2 = V*H**3/(3*E*I) + 1.2*V*H/(G*A)

    # Check that the solution matches the theoretical solution within 2%
    assert abs(1 - delta1/delta2) < 0.02, 'Failed cracked rect plate shear wall test.'


def test_shear_wall_openings():
    # This example demonstrates how to analyze a shear wall with openings. It
    # follows Section 10.5.3 of "Masonry Structures - Behavior and Design, 2nd
    # Edition" by Robert G. Drysdale, Ahmad A. Hamid, and Lawrie R. Baker. The
    # solution given in that text is obtained using an approximation method that
    # isn't nearly as accurate as the finite element method, so some differences
    # in the final results are expected.

    # Create a finite element model
    model = FEModel3D()

    # Set material properties for the wall (2 ksi masonry)
    f_m = 2000        # Masonry compressive strength (psi)
    E = 900*f_m/1000  # Masonry modulus of elasticity (ksi)
    nu = 0.17         # Poisson's ratio for masonry
    rho = 78/1000/12**2/7.625  # Masonry unit weight (kip/in^3)
    model.add_material('Masonry', E, 0.4*E, nu, rho)

    # Choose a desired mesh size. The program will try to stick to this as best as it can.
    mesh_size = 6  # in

    # Set the wall's dimensions
    width = 26*12   # Wall overall width (in)
    height = 16*12  # Wall overall height (in)
    t = 8           # Masonry thickness (in)

    # Generate the rectangular mesh
    # The effects of cracked masonry can be modeled by adjusting the `ky_mod` factor. For this example
    # uncracked masonry will be used to match the textbook problem.
    model.add_rectangle_mesh('MSH1', mesh_size, width, height, t, 'Masonry', kx_mod=1, ky_mod=1,
                                origin=[0, 0, 0], plane='XY', element_type='Rect')

    # Add a 4' wide x 12' tall door opening to the mesh
    model.meshes['MSH1'].add_rect_opening(name='Door 1', x_left=2*12, y_bott=0*12, width=4*12, height=12*12)

    # Add a 4' wide x 4' tall window opening to the mesh
    model.meshes['MSH1'].add_rect_opening(name='Window 1', x_left=8*12, y_bott=8*12, width=4*12, height=4*12)

    # Add another 4' wide x 4' tall window opening to the mesh
    model.meshes['MSH1'].add_rect_opening(name='Window 2', x_left=14*12, y_bott=8*12, width=4*12, height=4*12)

    # Add another 4' wide x 12' tall door opening to the mesh
    model.meshes['MSH1'].add_rect_opening(name='Door 2', x_left=20*12, y_bott=0*12, width=4*12, height=12*12)

    # Generate the mesh now that we've defined all the openings
    model.meshes['MSH1'].generate()

    # Shear at the top of the wall
    V = 100  # kip

    # The shear at the top of the wall will be distributed evently to all the
    # nodes at the top of the wall. Determine how many nodes are at the top of the
    # wall.
    n = len([node for node in model.nodes.values() if isclose(node.Y, height)])
    v = V/n

    # Add supports and loads to the nodes
    for node in model.nodes.values():

        # Determine if the node is at the base of the wall
        if isclose(node.Y, 0):
            # Fix the base of the wall
            model.def_support(node.name, True, True, True, True, True, True)
        # Determine if the node is at the top of the wall
        elif isclose(node.Y, height):
            # Add out-of-plane support (provided by the diaphragm)
            model.def_support(node.name, False, False, True, False, False, False)
            # Add a seismic shear load to the top of the wall (applied by the diaphragm)
            model.add_node_load(node.name, 'FX', v, case='E')

    # Add a load combination named 'Seismic' with a factor of 1.0 applied to any loads designated as
    # 'E'.
    model.add_load_combo('Seismic', {'E': 1.0})

    # Analyze the model
    model.analyze(log=True, check_statics=True)

    # # Render the model and plot the `Txy` shears.
    # # window = render_model(model, text_height=1, render_loads=True, deformed_shape=True,
    # #                       deformed_scale=200, color_map='Txy', scalar_bar=False,
    # #                       combo_name='Seismic', labels=False, screenshot='console')
    # from Pynite.Visualization import Renderer
    # renderer = Renderer(model)
    # renderer.combo_name = 'Seismic'
    # renderer.color_map = 'Txy'
    # renderer.annotation_size = 1
    # renderer.deformed_shape = True
    # renderer.deformed_scale = 200
    # renderer.scalar_bar = True
    # # renderer.render_model()
    # renderer.screenshot()

    # Print the maximum displacement
    # d_max = max([node.DX['Seismic'] for node in model.nodes.values()])
    # print('Max displacement: ', d_max, 'in')
    # print('Expected displacement from reference text: ', 7.623/E*t, 'in')
    # print('Wall rigidity: ', V/d_max, 'kips/in')
    # print('Expected rigidity from reference text: ', 1/7.623*E*t, 'kips/in')


if __name__ == '__main__':
    test_shear_wall()
