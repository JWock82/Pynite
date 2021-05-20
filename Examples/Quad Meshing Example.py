from PyNite.FEModel3D import FEModel3D
from PyNite.Mesh import FrustrumMesh, CylinderMesh
from PyNite.Visualization import RenderModel

t = 0.25/12
E = 29000*1000*12**2
nu = 0.3
mesh_size = 0.5
r_shell = 5
r_hopper = 1
h_shell = 15
h_hopper = 10
center = [0, 0, 0]
start_node = 'N1'
start_element = 'Q1'

# Generate the conical hopper mesh
hopper_mesh = FrustrumMesh(t, E, nu, mesh_size, r_shell, r_hopper, h_hopper, center, start_node, start_element)

# Determine how many elements make up the circumference at the top of the hopper. This will
# determine the number of elements making up the circumference of the cylindrical shell.
n_hopper = hopper_mesh.num_quads_outer

# Find an unused node and element name to start the cylinder off with
first_node = 'N' + str(len(hopper_mesh.nodes) + 1)
first_element = 'Q' + str(len(hopper_mesh.elements) + 1)

# Generate the cylindrical shell
shell_mesh = CylinderMesh(t, E, nu, mesh_size, r_shell, h_shell, center, first_node, first_element, n_hopper)

# Create a finite element model
model = FEModel3D()

# Add the two meshes to the model
model.AddMesh(hopper_mesh)
model.AddMesh(shell_mesh)

# The two meshes have overlapping duplicate nodes that need to be merged
model.MergeDuplicateNodes()

# Add hydrostatic pressure to each element in the model
for element in model.Quads.values():
    Zavg = (element.iNode.Z + element.jNode.Z + element.mNode.Z + element.nNode.Z)/4
    model.AddQuadSurfacePressure(element.Name, -62.4*(h_shell - Zavg), case='Hydrostatic')

# Add 4 supports at the springline at quarter points
for node in model.Nodes.values():
    if round(node.Z, 10) == 0 and (round(node.X, 10) == 0 or round(node.Y, 10) == 0):
        model.DefineSupport(node.Name, True, True, True, False, False, False)

# Add a load combination
model.AddLoadCombo('1.4F', {'Hydrostatic': 1.4})

# Analyze the model
model.Analyze()

# Render the model. Labels and loads will be turned off to speed up interaction.
RenderModel(model, 0.1, render_loads=False, color_map='dz', combo_name='1.4F', labels=False)
