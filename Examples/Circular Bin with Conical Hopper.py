from math import pi
from PyNite.FEModel3D import FEModel3D
from PyNite.Mesh import FrustrumMesh, CylinderMesh
from PyNite.Visualization import render_model

t = 0.25/12
E = 29000*1000*12**2
nu = 0.3
kx_mod=1
ky_mod=1
mesh_size = 0.5
r_shell = 7.5
r_hopper = 1.5
h_shell = 20
h_hopper = 10
center = [0, 0, 0]
axis = 'Y'
start_node = 'N1'
start_element = 'Q1'

# Generate the conical hopper mesh
hopper_mesh = FrustrumMesh(mesh_size, r_shell, r_hopper, h_hopper, t, E, nu, kx_mod, ky_mod,
                           center, axis, start_node, start_element)

# Determine how many elements make up the circumference at the top of the hopper. This will
# determine the number of elements making up the circumference of the cylindrical shell.
n_hopper = hopper_mesh.num_quads_outer

# Find an unused node and element name to start the cylinder off with
first_node = 'N' + str(len(hopper_mesh.nodes) + 1)
first_element = 'Q' + str(len(hopper_mesh.elements) + 1)

# Generate the cylindrical shell
shell_mesh = CylinderMesh(mesh_size, r_shell, h_shell, t, E, nu, kx_mod, ky_mod, center, axis,
                          first_node, first_element, n_hopper)

# Create a finite element model
model = FEModel3D()

# Add the two meshes to the model
model.add_mesh(hopper_mesh)
model.add_mesh(shell_mesh)

# The two meshes have overlapping duplicate nodes that need to be merged
import cProfile
# cProfile.run('model.merge_duplicate_nodes()', sort='cumtime')
model.merge_duplicate_nodes()

# Add hydrostatic pressure to each element in the model
for element in model.Quads.values():
    Yavg = (element.i_node.Y + element.j_node.Y + element.m_node.Y + element.n_node.Y)/4
    model.add_quad_surface_pressure(element.Name, 62.4*(h_shell - Yavg), case='Hydrostatic')

# Add supports at the springline
for node in model.Nodes.values():
    if round(node.Y, 10) == 0:
        model.def_support(node.Name, True, True, True, False, False, False)

# Add a load combination
model.add_load_combo('1.4F', {'Hydrostatic': 1.4})

# Analyze the model
# cProfile.run('model.analyze()', sort='cumtime')
model.analyze()

# Render the model. Labels and loads will be turned off to speed up interaction.
render_model(model, 0.1, render_loads=True, color_map='dz', combo_name='1.4F', labels=False)