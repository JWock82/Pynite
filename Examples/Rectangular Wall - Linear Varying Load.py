# This example demonstrates how to analyze a tank wall for hydrostatic loads.

# Import a few libraries from Pynite that we'll need
from Pynite.FEModel3D import FEModel3D

# Create a finite element model
model = FEModel3D()

# Set material properties for the wall (4500 psi concrete)
E = 57000*(4500)**0.5*12**2  # psf
G = 0.4*E
nu = 1/6
rho = 150  # pcf
model.add_material('Concrete', E, G, nu, rho)

# Choose a mesh size
mesh_size = 1

# Set the wall's dimensions
width = 25   # ft
height = 15  # ft
t = 1        # ft

# Set the liquid height
HL = 15  # 12.5  # ft

# Create a rectangular mesh. A few things worth noting:
# 1. We will build it from rectangular plate elements. Alternatively we could build it from
#    quadrilateral elements by setting 'element_type' equal to 'Quad' instead. Since the geometry
#    is undistorted, and the wall panel is relatively thin (insignificant transverse shear
#    deformations), rectangular plate elements will produce better plate corner stress results.
# 2. It'd be nice to have the mesh hit the liquid level, so we'll add a control point at the
#    liquid level to the list of control points along the mesh's local y-axis.
model.add_rectangle_mesh('MSH1', mesh_size, width, height, t, 'Concrete', kx_mod=1, ky_mod=1,
                         origin=[0, 0, 0], plane='XY', y_control=[HL], element_type='Quad')

# Generate the mesh prior to analysis. This step is optional, but it allows you to work with it to
# it before analyzing.
model.meshes['MSH1'].generate()

# Step through each quadrilateral and rectangular plate in the model
for element in list(model.quads.values()) + list(model.plates.values()):

    # Calculate the average elevation of the element based on the positions of its nodes
    Yavg = (element.i_node.Y + element.j_node.Y + element.m_node.Y + element.n_node.Y)/4

    # Determine if the element falls below the liquid level
    if Yavg < HL:

        # Add hydrostatic loads to the element
        if model.meshes['MSH1'].element_type == 'Rect':
            model.add_plate_surface_pressure(element.name, 62.4*(HL - Yavg), case='F')
        else:
            model.add_quad_surface_pressure(element.name, 62.4*(HL - Yavg), case='F')

# Add fully fixed supports at left, right, and bottom of the wall
for node in model.nodes.values():
    if round(node.Y, 10) == 0 or round(node.X, 10) == 0 or round(node.X, 10) == width:
        model.def_support(node.name, True, True, True, True, True, True)

# Add a load combination named '1.4F' with a factor of 1.4 applied to any loads designated as 'Hydrostatic'.
model.add_load_combo('1.4F', {'F': 1.4})

# Analyze the model
model.analyze(log=True, check_statics=True)

# Render the model and plot the `Mx` moments.
from Pynite.Rendering import Renderer
renderer = Renderer(model)
renderer.annotation_size = 0.2
renderer.render_loads = True
renderer.deformed_shape = True
renderer.deformed_scale = 1000
renderer.color_map = 'Qy'
renderer.combo_name = '1.4F'
renderer.show_labels = True
renderer.scalar_bar = True
renderer.scalar_bar_text_size = 12
renderer.render_model()

# Timoshenko solution
b = HL
a = width
qo = HL*62.4*1.4
Qx = 0.136*qo*a
Qy = 0.248*qo*a
Mx = -0.0131*qo*a**2
My = -0.0242*qo*a**2

# Pynite solution
Qx_pn = model.quads['Q176'].shear(-1, 0, True, '1.4F')[0, 0]
Qy_pn = model.meshes['MSH1'].max_shear('Qy', '1.4F')
Mx_pn = model.quads['Q176'].moment(-1, 0, True, '1.4F')[0, 0]
My_pn = model.meshes['MSH1'].min_moment('My', '1.4F')

# Comparison of solutions
print('Max Moment at Side Mid-Height of Wall, Mx | Pynite: ', Mx_pn, '| Timoshenko: ', Mx)
print('Max Moment at Base of Wall, My | Pynite: ', My_pn, '| Timoshenko: ', My)
print('Max Shear at Side Mid-Height of Wall, Qx | Pynite: ', Qx_pn, '| Timoshenko: ', Qx)
print('Max Shear at Base of Wall, My | Pynite: ', Qy_pn, '| Timoshenko: ', Qy)

# Print the maximum and minumum displacements
D = E*t**3/(12*(1-nu**2))
d = 0.00069*qo*a**4/D
print('Max displacement: ', max([node.DZ['1.4F'] for node in model.nodes.values()]))
print('Timoshenko Solution for displacement, d: ', d)
