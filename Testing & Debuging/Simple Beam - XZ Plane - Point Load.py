# Example of a simply supported beam with a point load.
# Units used in this example are inches, and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new finite element model
SimpleBeam = FEModel3D()

# Add nodes (14 ft = 168 in apart)
SimpleBeam.AddNode("N1", 0, 0, 0)
SimpleBeam.AddNode("N2", 0, 0, 168)

# Add a beam with the following properties:
A = 20
E = 29000
G = 11400
Iy = 100
Iz = 150
J = 250
SimpleBeam.AddMember("M1", "N1", "N2", E, G, Iy, Iz, J, A)

# Provide simple supports
SimpleBeam.DefineSupport("N1", True, True, True, False, False, True)
SimpleBeam.DefineSupport("N2", True, True, True, False, False, False)

# Add a point load of 5 kips at the midspan of the beam
SimpleBeam.AddMemberPtLoad("M1", "Fy", 5, 7 * 12)

# Analyze the beam
SimpleBeam.Analyze(False)

# Print the shear, moment, and deflection diagrams
SimpleBeam.GetMember("M1").PlotShear("Fy")
SimpleBeam.GetMember("M1").PlotMoment("Mz")
SimpleBeam.GetMember("M1").PlotDeflection("dy")

# Print reactions at each end of the beam
print('Left Support Reaction:', SimpleBeam.GetNode("N1").RxnFY['Combo 1'])
print('Right Support Reacton:', SimpleBeam.GetNode("N2").RxnFY['Combo 1'])