# A First Course in the Finite Element Method, 4th Edition
# Daryl L. Logan
# Example 2.1
# Units for this model are pounds and inches

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

system = FEModel3D()

system.AddNode('1', 0, 0, 0)
system.AddNode('2', 30, 0, 0)
system.AddNode('3', 10, 0, 0)
system.AddNode('4', 20, 0, 0)

system.AddSpring('S1', '1', '3', 1000)
system.AddSpring('S2', '3', '4', 2000)
system.AddSpring('S3', '4', '2', 3000)

system.DefineSupport('1', True, True, True, True, True, True)
system.DefineSupport('2', True, True, True, True, True, True)

system.DefineSupport('3', False, True, True, True, True, True)
system.DefineSupport('4', False, True, True, True, True, True)

system.AddNodeLoad('4', 'FX', 5000)

system.Analyze(True)

print(system.GetNode('3').DX['Combo 1'])
print(system.GetNode('4').DX['Combo 1'])
print('')
print(system.GetNode('1').RxnFX['Combo 1'])
print(system.GetNode('2').RxnFX['Combo 1'])

from PyNite import Visualization
Visualization.RenderModel(system, text_height=0.5, deformed_shape=True,
                          deformed_scale=1)