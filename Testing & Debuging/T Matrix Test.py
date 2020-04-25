# Matrix Structural Analysis, 2nd Edition
# William McGuire, Richard H. Gallagher, Ronald D. Ziemian
# Example 5.4
# Units for this model are kilonewtons and meters

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new model
truss = FEModel3D()

# Define the nodes
truss.AddNode('a', 2, 4, 8)
truss.AddNode('b', 0, 0, 0)
truss.AddNode('c', 8, 0, 0)
truss.AddNode('d', 8, 6, 0)
truss.AddNode('e', 0, 6, 0)

# Define the supports - All fixed - End releases will be used to create pins
truss.DefineSupport('b', True, True, True, True, True, True)
truss.DefineSupport('c', True, True, True, True, True, True)
truss.DefineSupport('d', True, True, True, True, True, True)
truss.DefineSupport('e', True, True, True, True, True, True)

# Define properties common to all members
E = 200000000  # kN/m^2
G = 0.4*E
J = 100
Iy = 100
Iz = 100

# Define members
truss.AddMember('ab', 'a', 'b', E, G, Iy, Iz, J, 20*10**3/1000**2)
truss.AddMember('ac', 'a', 'c', E, G, Iy, Iz, J, 30*10**3/1000**2)
truss.AddMember('ad', 'a', 'd', E, G, Iy, Iz, J, 40*10**3/1000**2)
truss.AddMember('ae', 'a', 'e', E, G, Iy, Iz, J, 30*10**3/1000**2)

# Define member end releases
truss.DefineReleases('ab', False, False, False, False, True, True,
                           False, False, False, False, True, True)
truss.DefineReleases('ac', False, False, False, False, True, True,
                           False, False, False, False, True, True)
truss.DefineReleases('ad', False, False, False, False, True, True,
                           False, False, False, False, True, True)
truss.DefineReleases('ae', False, False, False, False, True, True,
                           False, False, False, False, True, True)

# Add nodal loads
truss.AddNodeLoad('a', 'FY', 600)
truss.AddNodeLoad('a', 'FX', 200)
truss.AddNodeLoad('a', 'FZ', -800)

# Analyze the model
truss.Analyze(False)

# Render the model
Visualization.RenderModel(truss, text_height=0.2, deformed_shape=True, deformed_scale=1000, render_loads=True)

# Print results
a = truss.GetNode('a')
b = truss.GetNode('b')
c = truss.GetNode('c')
d = truss.GetNode('d')
e = truss.GetNode('e')

print('Node b Reactions: ', b.RxnFX, b.RxnFY, b.RxnFZ, b.RxnMX, b.RxnMY, b.RxnMZ)
print('Node c Reactions: ', c.RxnFX, c.RxnFY, c.RxnFZ, c.RxnMX, c.RxnMY, c.RxnMZ)
print('Node e Reactions: ', d.RxnFX, d.RxnFY, d.RxnFZ, d.RxnMX, d.RxnMY, d.RxnMZ)
print('Node e Reactions: ', e.RxnFX, e.RxnFY, e.RxnFZ, e.RxnMX, e.RxnMY, e.RxnMZ)
print('Node a Displacements: ', a.DX, a.DY, a.DZ)

# Results printed above matched the textbook problem results perfectly