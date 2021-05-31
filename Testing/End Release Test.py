# Simple beam created by using an end release at fixed supports
# Units used in this example are inches, and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D
from PyNite import Visualization

myModel = FEModel3D()

myModel.AddNode('N1', 0, 0, 0)
myModel.AddNode('N2', 10*12, 0, 0)

myModel.DefineSupport('N1', True, True, True, True, True, True)
myModel.DefineSupport('N2', True, True, True, True, True, True)

myModel.AddMember('M1', 'N1', 'N2', 29000, 11400, 100, 150, 250, 10)

myModel.DefineReleases('M1', False, False, False, False, False, True, \
                             False, False, False, False, False, True)

myModel.AddMemberDistLoad('M1', 'Fy', -0.5, -0.5)

myModel.Analyze()

Visualization.RenderModel(myModel, text_height=3, deformed_shape=True, deformed_scale=100, render_loads=True)

myModel.GetMember('M1').PlotShear('Fy')
myModel.GetMember('M1').PlotMoment('Mz')
print('Calculated moment: ' + str(myModel.GetMember('M1').MinMoment('Mz')))
print('Expected moment: ' + str(-0.5*(10*12)**2/8))