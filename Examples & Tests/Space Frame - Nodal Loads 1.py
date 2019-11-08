# A First Course in the Finite Element Method, 4th Edition
# Daryl L. Logan
# Example 5.8
# Units for this model are kips and inches

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new model
frame = FEModel3D()

# Define the nodes
frame.AddNode('N1', 0, 0, 0)
frame.AddNode('N2', -100, 0, 0)
frame.AddNode('N3', 0, 0, -100)
frame.AddNode('N4', 0, -100, 0)

# Define the supports
frame.DefineSupport('N2', True, True, True, True, True, True)
frame.DefineSupport('N3', True, True, True, True, True, True)
frame.DefineSupport('N4', True, True, True, True, True, True)

# Create members (all members will have the same properties in this example)
J = 50
Iy = 100
Iz = 100
E = 30000
G = 10000
A = 10

frame.AddMember('M1', 'N1', 'N2', E, G, Iy, Iz, J, A)
frame.AddMember('M2', 'N1', 'N3', E, G, Iy, Iz, J, A)
frame.AddMember('M3', 'N1', 'N4', E, G, Iy, Iz, J, A)

# Add nodal loads
frame.AddNodeLoad('N1', 'FY', -50)
frame.AddNodeLoad('N1', 'MX', -1000)

frame.Analyze()
Visualization.RenderModel(frame)
print('Node 1 deformations:')
print('Calculated values: ', [frame.GetNode('N1').DX, frame.GetNode('N1').DY, frame.GetNode('N1').DZ, frame.GetNode('N1').RX, frame.GetNode('N1').RY, frame.GetNode('N1').RZ])
print('Expected values: ', [7.098e-5, -0.014, -2.352e-3, -3.996e-3, 1.78e-5, -1.033e-4])