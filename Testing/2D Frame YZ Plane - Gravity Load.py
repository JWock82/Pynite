# A First Course in the Finite Element Method, 4th Edition
# Daryl L. Logan
# Problem 5.30
# Units for this model are kips and inches

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new model
frame = FEModel3D()

# Define the nodes
frame.AddNode('N1', 0, 0, 0)
frame.AddNode('N2', 0, 30*12, 0)
frame.AddNode('N3', 0, 40*12, 15*12)
frame.AddNode('N4', 0, 40*12, 35*12)
frame.AddNode('N5', 0, 30*12, 50*12)
frame.AddNode('N6', 0, 0, 50*12)

# Define the supports
frame.DefineSupport('N1', True, True, True, True, True, True)
frame.DefineSupport('N6', True, True, True, True, True, True)

# Create members (all members will have the same properties in this example)
J = 250
Iy = 250
Iz = 200
E = 30000
G = 250
A = 12

frame.AddMember('M1', 'N1', 'N2', E, G, Iz, Iy, J, A)
frame.AddMember('M2', 'N2', 'N3', E, G, Iy, Iz, J, A)
frame.AddMember('M3', 'N3', 'N4', E, G, Iy, Iz, J, A)
frame.AddMember('M4', 'N4', 'N5', E, G, Iy, Iz, J, A)
frame.AddMember('M5', 'N5', 'N6', E, G, Iz, Iy, J, A)

# Add nodal loads
frame.AddNodeLoad('N3', 'FY', -30)
frame.AddNodeLoad('N4', 'FY', -30)

# Analyze the model
frame.Analyze()

# Render the model for viewing
Visualization.RenderModel(frame, text_height=15, deformed_shape=True, deformed_scale=10, render_loads=True)

node1 = frame.GetNode('N1')
node6 = frame.GetNode('N6')
print('Calculated reactions: ', node1.RxnFZ, node1.RxnFY, node1.RxnMX, node6.RxnFZ, node6.RxnFY, node6.RxnMX)
print('Expected reactions: ', 11.69, 30, 1810, -11.69, 30, -1810)
print('Calculated displacements: ', frame.GetNode('N3').DY, frame.GetNode('N4').DY, frame.GetNode('N3').RX, frame.GetNode('N4').RX)