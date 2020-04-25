# A First Course in the Finite Element Method, 4th Edition
# Daryl L. Logan
# Problem 5.58
# Units for this model are kips and inches

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new model
frame = FEModel3D()

# Define the nodes
frame.AddNode('N1', 0, 0, 0)
frame.AddNode('N2', 10*12, 0, 0)
frame.AddNode('N3', 10*12, 0, -10*12)
frame.AddNode('N4', 10*12, -20*12, -10*12)

# Define the supports
frame.DefineSupport('N1', True, True, True, True, True, True)
frame.DefineSupport('N4', True, True, True, True, True, True)

# Create members (all members will have the same properties in this example)
J = 100
Iy = 200
Iz = 1000
E = 30000
G = 10000
A = 100

frame.AddMember('M12', 'N1', 'N2', E, G, Iy, Iz, J, A)
frame.AddMember('M23', 'N2', 'N3', E, G, Iy, Iz, J, A)
frame.AddMember('M34', 'N3', 'N4', E, G, Iy, Iz, J, A)

# Add nodal loads
frame.AddNodeLoad('N2', 'FY', -5)
frame.AddNodeLoad('N2', 'MX', -100*12)
frame.AddNodeLoad('N3', 'FZ', 40)

# Analyze the frame
frame.Analyze()

print('Calculated results: ', frame.GetNode('N2').DY, frame.GetNode('N3').DZ)
print('Expected results: ', -0.063, 1.825)

# Render the model for viewing
Visualization.RenderModel(frame, text_height=5, deformed_shape=True, deformed_scale=40, render_loads=True)