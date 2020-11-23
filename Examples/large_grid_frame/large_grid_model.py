from PyNite import FEModel3D, Visualization
import itertools
import os, csv
import inputfiles
import testing_io

# Constants for the member properties
E, G, Iy, Iz, J, A = 30*10**6, 12*10**6, 0.00100, 0.00100, 50, 0.010

# Initialize the model
model = FEModel3D()

# Import nodes from file
path = os.path.join(os.path.dirname(__file__), 'gridnodes.csv')
nodes_list = inputfiles.nodes_from_csv(path)

# Add the nodes to the model
for node in nodes_list:
    model.AddNode(*node)

# Import members from file
path = os.path.join(os.path.dirname(__file__), 'gridmembers.csv')
member_list = inputfiles.read_csv(path)

# Add the members to the model
for member in member_list:
    Name, iNode, jNode = member
    model.AddMember(Name, iNode, jNode, E, G, Iy, Iz, J, A)

# Add supports
for node in ('N0', 'N4', 'N220', 'N224'):
    model.DefineSupport(node, True, True, True, True, True, True)

# Add node loads
for node in ('N50','N54','N270','N274'):
    model.AddNodeLoad(node, 'FY', -10)

# Analyze the model
model.Analyze(check_statics=True)

# Render the model
Visualization.RenderModel(model,
                          deformed_shape=True,
                          text_height=0.25,
                          render_loads=True)
