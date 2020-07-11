#Units N e m

###Import FEModel3D from PyNite
from PyNite import FEModel3D
from PyNite import Visualization

###Create a new finite element model
Beam = FEModel3D()
L = 5 #m

###Nodes
Beam.AddNode("N1", 0, 0, 0)
Beam.AddNode("N2", L, 0, 0)

###Beams (30x50 cm)
E = 2.1e11 #N/m^2
G = 1
Iy = 0.001125 #m^4
Iz = 0.003125 #m^4
J = 1
A = 0.15 #m^2

Beam.AddMember("M1", "N1", "N2", E, G, Iy, Iz, J, A)

###Supports
Beam.DefineSupport("N1", True, True, True, True, True, True)
Beam.DefineSupport("N2", True, True, True, False, True, True)

###Load
Beam.AddMemberDistLoad("M1", "Fx", 10, 10, 0, 5)

# Analyze
Beam.Analyze()

# Diagrams V,M and displacements
Beam.GetMember("M1").PlotAxial()

# Member fixed end reaction vector
print('M1 Displacement Vector: ', Beam.GetMember('M1').d())
print('M1 Fixed End Reaction Vector: ', Beam.GetMember('M1').fer())

#Reactions
print("R1=",Beam.GetNode("N1").RxnFX,"N")
print("R2=",Beam.GetNode("N2").RxnFX,"N")