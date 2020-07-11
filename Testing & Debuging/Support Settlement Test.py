# This file tests the program's ability to model support settlements
# Units used in this test are inches and kips

# Structural Analysis, 3rd Edition
# Aslam Kassimali
# Example 13.14

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new finite element model for the beam
beam = FEModel3D()

# Add nodes
beam.AddNode('A', 0, 0, 0)
beam.AddNode('B', 20*12, 0, 0)
beam.AddNode('C', 40*12, 0, 0)
beam.AddNode('D', 60*12, 0, 0)

# Add members
A = 20
E = 29000
G = 11400
Iy = 1000
Iz = 7800
J = 8800
beam.AddMember('AB', 'A', 'B', E, G, Iy, Iz, J, A)
beam.AddMember('BC', 'B', 'C', E, G, Iy, Iz, J, A)
beam.AddMember('CD', 'C', 'D', E, G, Iy, Iz, J, A)

# Provide supports
beam.DefineSupport('A', True, True, True, True, False, False)
beam.DefineSupport('B', False, True, True, False, False, False)
beam.DefineSupport('C', False, True, True, False, False, False)
beam.DefineSupport('D', False, True, True, False, False, False)

# Add a uniform load to the beam
beam.AddMemberDistLoad('AB', 'Fy', -2/12, -2/12)
beam.AddMemberDistLoad('BC', 'Fy', -2/12, -2/12)
beam.AddMemberDistLoad('CD', 'Fy', -2/12, -2/12)

# Add support settlements
beam.AddNodeDisplacement('B', 'DY', -5/8)
beam.AddNodeDisplacement('C', 'DY', -1.5)
beam.AddNodeDisplacement('D', 'DY', -0.75)

# Analyze the beam
beam.Analyze()

# Render the beam's deformed shape
Visualization.RenderModel(beam, text_height=12, deformed_shape=True, deformed_scale=40, render_loads=True)

# Print reactions
print('Calculated Reaction at A:', beam.GetNode('A').RxnFY, 'k')
print('Calculated Reaction at B:', beam.GetNode('B').RxnFY, 'k')
print('Calculated Reaction at C:', beam.GetNode('C').RxnFY, 'k')
print('Calculated Reaction at D:', beam.GetNode('D').RxnFY, 'k')
print('Expected Reaction at A: -1.098 k')
print('Expected Reaction at B: 122.373 k')
print('Expected Reaction at C: -61.451 k')
print('Expected Reaction at D: 60.176 k')

# Print the shear diagrams
beam.GetMember('AB').PlotShear('Fy')
beam.GetMember('BC').PlotShear('Fy')
beam.GetMember('CD').PlotShear('Fy')

# Print the moment diagrams
beam.GetMember('AB').PlotMoment('Mz')
beam.GetMember('BC').PlotMoment('Mz')
beam.GetMember('CD').PlotMoment('Mz')