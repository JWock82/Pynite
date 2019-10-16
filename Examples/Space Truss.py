# Import `FEModel3D` and `Visualization` from `PyNite`
from PyNite import FEModel3D
from PyNite import Visualization

truss = FEModel3D()

truss.AddNode('A', 1.1, -0.4, 0)
truss.AddNode('B', 1, 0, 0)
truss.AddNode('C', 0, 0, 0.6)
truss.AddNode('D', 0, 0, -0.4)
truss.AddNode('E', 0, 0.8, 0)

truss.DefineSupport('C', True, True, True, True, True, True)
truss.DefineSupport('D', True, True, True, True, True, True)
truss.DefineSupport('E', True, True, True, True, True, True)

truss.AddMember('AB', 'A', 'B', 100, 100, 100, 100, 100, 100)
truss.AddMember('AC', 'A', 'C', 100, 100, 100, 100, 100, 100)
truss.AddMember('AD', 'A', 'D', 100, 100, 100, 100, 100, 100)
truss.AddMember('BC', 'B', 'C', 100, 100, 100, 100, 100, 100)
truss.AddMember('BD', 'B', 'D', 100, 100, 100, 100, 100, 100)
truss.AddMember('BE', 'B', 'E', 100, 100, 100, 100, 100, 100)

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

truss.AddNodeLoad('A', 'FX', 10)
truss.AddNodeLoad('A', 'FY', 60)
truss.AddNodeLoad('A', 'FZ', 20)

truss.Analyze()

print(truss.GetMember('BC').MinAxial())
print(truss.GetMember('BD').MinAxial())
print(truss.GetMember('BE').MaxAxial())

Visualization.RenderModel(truss, 0.05)
