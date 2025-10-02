# This example uses Pynite's new `ShearWall` object instead of its `FEModel3D` object. The `ShearWall` object has a built-in `FEModel3D` and includes additional features built in to simplify analysis of complex shear wall shapes and loadings. It can automatically subdivide the wall into piers and coupling beams, help you quickly identify loads acting on individual piers or coupling beams (handy for seismic design), and report wall stiffness (used for rigid diaphragm analysis). I built this because many of the provisions in the seismic building codes are difficult to apply without a good tool, and many of the wall geometries used in modern structures are much more complex than simple rectangles without openings or wall intersections. Enjoy!

# Model units will be kips and feet   
from Pynite import FEModel3D

# Create a new finite element model
model = FEModel3D()

# Define material properties that will be used for the wall
E = 57000*(4500)**0.5/1000*12**2  # Modulus of elasticity of concrete (ksf)
G = 0.4*E                         # Shear modulus of concrete (ksf)
nu = 0.17                         # Poisson's ratio for concrete
rho = 0.150                       # Concrete unit weight (kcf)

# Add the material to the model
model.add_material('Concrete', E, G, nu, rho)

# Create the shear wall
mesh_size = 1       # Desired mesh size (ft)
length = 50         # Wall length (ft)
height = 25         # Wall height (ft)
mesh_size = 1       # Desired mesh size (ft)
t = 1               # Wall thickness (ft)
ky_mod = 0.35       # Wall stiffness adjustment (in the vertical direction) to account for cracking
plane = 'XY'        # Plane to place the wall parallel to
origin = [0, 0, 0]  # Origin in 3D space for the wall
model.add_shear_wall('Wall 1', mesh_size, length, height, t, 'Concrete', ky_mod, 'XY', [0, 0, 0])

# Add an opening to the shear wall
# The `tie` parameter can be used to add a tie above the top of the opening, forcing piers on either side to work together. This is only necessary for openings at the top of the wall where a drag strut or collector is used. The `tie` parameter accepts a tie's stiffness.
model.shear_walls['Wall 1'].add_opening('Opng 1', x_start=5, y_start=6, width=4, height=2, tie=None)

# Flanges (intersecting wall segments) can affect the behavior of shear walls greatly by increasing the stiffness of the wall or individual wall piers. Building codes typically require flanges to be considered in shear wall analysis. We can add as many flanges as we like to the wall here. The `side` parameter defines which side the flage is on ('+z' or '-z' in the wall's local coordinate system). For a tee-shaped intersection we would specify two flanges, one '+z' and one '-z'.
model.shear_walls['Wall 1'].add_flange(t, width=3, x=0, y_start=0, y_end=height, material='Concrete', side='+z')

# Add foundations to the model
# The shear wall needs support. Multiple supports can be specified at multiple elevations to allow for stepped walls to be modeled.
model.shear_walls['Wall 1'].add_support(elevation=0, x_start=0, x_end=length)

# Add stories/floors to the model
# `Pynite` will calculate and store the wall's stiffness at every story. Floors can be attached partially across the story, or all the way across. Shears applied to a story will be distributed evenly between the nodes that are part of the floor.
model.shear_walls['Wall 1'].add_story(story_name='Roof', elevation=height, x_start=0, x_end=length)
model.shear_walls['Wall 1'].add_shear(story_name='Roof', force=150, case='E')

# Add a general load combination
model.add_load_combo('1.2D+1.0E', {'D': 1.2, 'E': 1.0})

# Analyze the wall - note that the wall itself is not a model, but has a `FEModel3D` object built into it as an attribute. The `ShearWall` class is essentially a class that builds the model for you.
model.analyze_linear(log=True, check_statics=False, check_stability=False)

# Print out the wall stiffness at the roof (k/in)
print('Roof level stiffness: ', model.shear_walls['Wall 1'].stiffness('Roof')/12)

# Let's render the wall shear stresses
from Pynite.Rendering import Renderer
rndr = Renderer(model)
rndr.annotation_size = 2
rndr.combo_name = '1.2D+1.0E'
rndr.deformed_shape = True
rndr.deformed_scale = 1000
rndr.color_map = 'Txy'
rndr.render_model()

# Show the wall piers for which results are available
model.shear_walls['Wall 1'].draw_piers(show=True)

# Note: press `q` to close the window that pops up and move on. If you hit the `X` in the corner the program will kill the window before we can reuse it below on the coupling beams (sorry - that's a `pyvista` nuance I haven't figures out a workaround for yet)

# Print the pier results
model.shear_walls['Wall 1'].print_piers(combo_name='1.2D+1.0E')

# Show the wall coupling beams for which results are available
model.shear_walls['Wall 1'].draw_coupling_beams(show=True)

# Print the coupling beam results
model.shear_walls['Wall 1'].print_coupling_beams(combo_name='1.2D+1.0E')

# Note: At this point all the pier and coupling beam results are stored in the `wall.piers` and `wall.coupling_beams` dictionaries for further use if desired.
