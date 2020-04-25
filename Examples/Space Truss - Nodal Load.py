# Engineering Mechanics: Statics, 4th Edition
# Bedford and Fowler
# Problem 6.64
# Units for this model are meters and kilonewtons

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new model
truss = FEModel3D()

# Define the nodes
truss.AddNode('A', 1.1, -0.4, 0)
truss.AddNode('B', 1, 0, 0)
truss.AddNode('C', 0, 0, 0.6)
truss.AddNode('D', 0, 0, -0.4)
truss.AddNode('E', 0, 0.8, 0)

# Define the supports
truss.DefineSupport('C', True, True, True, True, True, True)
truss.DefineSupport('D', True, True, True, True, True, True)
truss.DefineSupport('E', True, True, True, True, True, True)

# Create members
# Member properties were not given for this problem, so assumed values will be used
# To make all the members act rigid, the modulus of elasticity will be set to a very large value
E = 99999999
truss.AddMember('AB', 'A', 'B', E, 100, 100, 100, 100, 100)
truss.AddMember('AC', 'A', 'C', E, 100, 100, 100, 100, 100)
truss.AddMember('AD', 'A', 'D', E, 100, 100, 100, 100, 100)
truss.AddMember('BC', 'B', 'C', E, 100, 100, 100, 100, 100)
truss.AddMember('BD', 'B', 'D', E, 100, 100, 100, 100, 100)
truss.AddMember('BE', 'B', 'E', E, 100, 100, 100, 100, 100)

# Release the moments at the ends of the members to make truss members
truss.DefineReleases('AC', False, False, False, False, True, True, \
                           False, False, False, False, True, True)
truss.DefineReleases('AD', False, False, False, False, True, True, \
                           False, False, False, False, True, True)
truss.DefineReleases('BC', False, False, False, False, True, True, \
                           False, False, False, False, True, True)
truss.DefineReleases('BD', False, False, False, False, True, True, \
                           False, False, False, False, True, True)
truss.DefineReleases('BE', False, False, False, False, True, True, \
                           False, False, False, False, True, True)

# Add nodal loads
truss.AddNodeLoad('A', 'FX', 10)
truss.AddNodeLoad('A', 'FY', 60)
truss.AddNodeLoad('A', 'FZ', 20)

# Analyze the model
truss.Analyze()

# Print results
print('Member BC calculated axial force: ' + str(truss.GetMember('BC').MaxAxial()))
print('Member BC expected axial force: 32.7 Tension')
print('Member BD calculated axial force: ' + str(truss.GetMember('BD').MaxAxial()))
print('Member BD expected axial force: 45.2 Tension')
print('Member BE calculated axial force: ' + str(truss.GetMember('BE').MaxAxial()))
print('Member BE expected axial force: 112.1 Compression')

# Render the model for viewing. The text height will be set to 50 mm.
# Because the members in this example are nearly rigid, there will be virtually no deformation. The deformed shape won't be rendered.
# The program has created a default load case 'Case 1' and a default load combo 'Combo 1' since we didn't specify any. We'll display 'Case 1'.
Visualization.RenderModel(truss, text_height=0.05, render_loads=True, case='Case 1')
