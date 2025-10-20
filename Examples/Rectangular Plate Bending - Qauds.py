# This is an example of quadrilateral elements used to model a rectangular wall under uniform load.
# This example demonstrates some of the nuances of using quadrilaterals in Pynite. Pynite's
# quadrilaterals are based on an DKMQ formulation, which produces very accurate results for thick
# and thin plates for the most part.

# Pynite has another element which is a rectangular plate bending element. See the "Rectangular
# Tank Wall - Hydrostatic Loads" example using that element.
t = 1  # ft
width = 10  # ft
height = 20  # ft
nu = 0.17
mesh_size = 1  # ft
load = 250  # psf

# Create a finite element model
from Pynite import FEModel3D
model = FEModel3D()

# Define concrete in the model
E = 57000*(4000)**0.5*12**2  # psf
model.add_material('Concrete', E, 0.4*E, 0.17, 150)

# Add the mesh to the model
model.add_rectangle_mesh('MSH1', mesh_size, width, height, t, 'Concrete', 1, 1, [0, 0, 0], 'XY', element_type='Quad')

# Pynite automatically generates each mesh when it analyzes. Sometimes it's convenient to generate
# the nodes and elements in the mesh in advance so that we can manipulate the mesh prior to
# analysis.
model.meshes['MSH1'].generate()

# Step through each quadrilateral in the model
for element in model.quads.values():
    # Add loads to the element
    model.add_quad_surface_pressure(element.name, load, case='W')

# Add fully fixed supports on all side of the wall
for node in model.nodes.values():
    if (round(node.Y, 10) == 0 or round(node.Y, 10) == height or round(node.X, 10) == 0 or
        round(node.X, 10) == width):
        model.def_support(node.name, True, True, True, True, True, True)

# Add a load combination named '1.0W' with a factor of 1.0 applied to any loads designated as 'W'
model.add_load_combo('1.0W', {'W': 1.0})

# Analyze the model
model.analyze(check_statics=True)

# +-----------------------+
# | Discussion of Results |
# +-----------------------+

# Set up a renderer for the wall. The quad mesh will be set to show 'Mx' results.
from Pynite.Rendering import Renderer
renderer = Renderer(model)
renderer.annotation_size = mesh_size/6
renderer.deformed_shape = False
renderer.combo_name = '1.0W'
renderer.color_map = 'Mx'
renderer.render_loads = True

# Render the model
renderer.render_model()

# The it should be noted that the rendered contours are smoothed. Smoothing averages the corner
# stresses from every quad framing into each node. This leads to a much more accurate contour.
# An unsmoothed plot would essentially show quad center stresses at the quad element corners.

# The mesh object in Pynite has built-in methods for finding the extreme values. Find the max
# moments for the mesh.
Mx_max = model.meshes['MSH1'].max_moment('Mx', '1.0W')
Mx_min = model.meshes['MSH1'].min_moment('Mx', '1.0W')
My_max = model.meshes['MSH1'].max_moment('My', '1.0W')
My_min = model.meshes['MSH1'].min_moment('My', '1.0W')

# Here are the expected results from Timoshenko's "Theory of Plates and Shells" Table 35, p. 202.
# Note that the deflection values for the Pynite solution are slightly larger, due to transverse
# shear deformations being accounted for.
D = E*t**3/(12*(1-nu**2))
print('Solution from Timoshenko Table 35 for b/a = 2.0:')
print('Expected displacement: ', 0.00254*load*width**4/D)
print('Expected Mx at Center:', -0.0412*load*width**2)
print('Expected Mx at Edges:', 0.0829*load*width**2)
print('Expected My at Center:', 0.0158*load*width**2)
print('Expected My at Top & Bottom:', -0.0571*load*width**2)

# Print Pynite's calculated maximum/minimum bending moments
print(f'Calculated Maximum Mx (at Center): {Mx_max}')
print(f'Calculated Minimum Mx (at Edges): {Mx_min}')
print(f'Calculated Maximum My (at Center): {My_max}')
print(f'Calculated Minimum My (at Top & Bottom): {My_min}')

# Note that Pynite and Timoshenko have opposite sign conventions, but the magnitudes are very close
# to Timoshenko's published solutions.
