# This is an example of quadrilateral elements used to model a rectangular wall under uniform load.
# This example demonstrates some of the nuances of using quadrilaterals in PyNite. PyNite's
# quadrilaterals are based on an MITC4 formulation, which produces very accurate results for thick
# and thin plates for the most part. One problem however is that stress results at quadrilateral
# corners can be difficult to determine. This example highlights the issue and shows how to work
# around it.

# PyNite has another element which is a rectangular plate bending element. Usually for this type of
# problem that element is better suited, as long as the wall does not have significant transverse
# shear deformations, and the elements are not skewed. See the "Rectangular Tank Wall - 
# Hydrostatic Loads" example using that element.
t = 1  # ft 
width = 10  # ft
height = 20  # ft 
nu = 0.17
mesh_size = 1  # ft
load = 250  # psf

# Create a finite element model
from PyNite import FEModel3D
model = FEModel3D()

# Define concrete in the model
E = 57000*(4000)**0.5*12**2  # psf
model.add_material('Concrete', E, 0.4*E, 0.17, 150)

# Add the mesh to the model
model.add_rectangle_mesh('MSH1', mesh_size, width, height, t, 'Concrete', 1, 1, [0, 0, 0], 'XY', element_type='Quad')

# PyNite automatically generates each mesh when it analyzes. Sometimes it's convenient to generate
# the nodes and elements in the mesh in advance so that we can manipulate the mesh prior to
# analysis.
model.Meshes['MSH1'].generate()

# Step through each quadrilateral in the model
for element in model.Quads.values():
    # Add loads to the element
    model.add_quad_surface_pressure(element.name, load, case='W')

# Add fully fixed supports on all side of the wall
for node in model.Nodes.values():
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
from PyNite.Rendering import Renderer
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

# Here are the expected results from Timoshenko's "Theory of Plates and Shells" Table 35, p. 202.
# Note that the deflection values for the PyNite solution are slightly larger, due to transverse
# shear deformations being accounted for.
D = E*t**3/(12*(1-nu**2))
print('Solution from Timoshenko Table 35 for b/a = 2.0:')
print('Expected displacement: ', 0.00254*load*width**4/D)
print('Expected Mx at Center:', -0.0412*load*width**2)
print('Expected Mx at Edges:', 0.0829*load*width**2)
print('Expected My at Center:', 0.0158*load*width**2)
print('Expected My at Top & Bottom:', -0.0571*load*width**2)

# It should be noted that even the smoothed Mx contours are off by nearly 30% from the theoretical
# solution at the wall boundaries. Because there are no adjacent quads at the boundaries, PyNite
# cannot smooth the results there, and center stresses are being reported at the boundaries
# instead of corner stresses. So what's really going on here?

# MITC4 elements are very accurate, but the internal bending stress results at the corners are not.
# They are secondary values extrapolated from the derivatives of primary values. The corner bending
# stresses appear to more accurately represent the bending stresses at the center of the element.
# While corner STRESSES have this problem, the corner FORCES do not. To find the bending stresses
# at the quad nodes we'll get the moment at the corner and divide it by 1/2 of the plate width to
# convert it to a stress result. Note that when we talk about quad "stresses" we're really talking
# about forces per unit length of the element.

# With 200 quads in the model, quad 101 is on the left edge of the wall at mid-height, and carries
# the maximum moment.

# Get the force vector for quad 101 for load combination '1.0W'
f_vector = model.Quads['Q101'].f('1.0W')

# Although the MITC4 element nodes are defined by the user in the order [i, j, m, n], internally
# PyNite formulates the plate in the order [m, n, i, j] to be consistent with the literature used
# to derive this element. Therefore the MITC4 element's force vector is arranged as follows:
# ************************************************************************************************************************************************************
# f vector: [[fx_m, fy_m, fz_m, mx_m, my_m, mz_m, fx_n, fy_n, fz_n, mx_n, my_n, mz_n, fx_i, fy_i, fz_i, mx_i, my_i, mz_i, fx_j, fy_j, fz_j, mx_j, my_j, mz_j]]
# index:    [[  0 ,   1 ,   2 ,   3 ,   4 ,   5 ,   6 ,   7 ,   8 ,   9 ,  10 ,  11 ,  12 ,  13 ,  14 ,  15 ,  16 ,  17 ,  18 ,  19 ,  20 ,  21 ,  22 ,  23 ]]
# ************************************************************************************************************************************************************
# nomenclature: fx_n = force in the local x-direction at the n-node
#               my_j = moment about the local y-axis at the j-node

# We are interested in the moment applied to the i-node about the quad's local y-axis (my_i). This
# is at index 16 in the f vector. Note that here Mx is the moment about the local y-axis, rather
# than about it's local x-axis. This can be confusing, but is a commonly used plate nomenclature.
Mx = f_vector[16, 0]

# We can find the max My moment similarly from quad 6
f_vector_11 = model.Quads['Q6'].f('1.0W')
My = f_vector_11[15, 0]

# Now we'll convert these values to a force per unit length. The height and width of the plate is
# `meshsize`.
Mx = Mx/(mesh_size/2)
My = My/(mesh_size/2)

# Print the correct maximum bending moment:
print('Calculated maximum Mx (back-calculated from qaud nodal forces): ', Mx)
print('Calculated maximum My (back-calculated from qaud nodal forces): ', My)

# These values are much closer to the Timoshenko solution than the direct stress results. The Mx
# solution is within 1%, and the My solution is within about 6%. That is very good convergence for
# the 1'x1' mesh size. If the mesh size is reduced to 0.5' the My solution is within about 4% of
# the Timoshenko solution. It should also be remembered that even the Timoshenko solution is an
# estimate.
