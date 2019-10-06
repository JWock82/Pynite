# Example of a moment frame with a lateral load
# Units used in this example are inches, and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new finite element model
MomentFrame = FEModel3D()

# Add nodes (frame is 15 ft wide x 12 ft tall)
MomentFrame.AddNode("N1", 0, 0, 0)
MomentFrame.AddNode("N2", 0, 12*12, 0)
MomentFrame.AddNode("N3", 15*12, 12*12, 0)
MomentFrame.AddNode("N4", 15*12, 0*12, 0)

# Add columns with the following properties:
# E = 29000 ksi, G = 11400 ksi, Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 10 in^2
MomentFrame.AddMember("M1", "N1", "N2", 29000, 11400, 100, 150, 250, 10)
MomentFrame.AddMember("M2", "N4", "N3", 29000, 11400, 100, 150, 250, 10)

# Add a beam with the following properties:
# E = 29000 ksi, G = 11400 ksi, Iy = 100 in^4, Iz = 250 in^4, J = 250 in^4, A = 15 in^2
MomentFrame.AddMember("M3", "N2", "N3", 29000, 11400, 100, 250, 250, 15)

# Provide fixed supports at the bases of the columns
MomentFrame.DefineSupport("N1", True, True, True, True, True, True)
MomentFrame.DefineSupport("N4", True, True, True, True, True, True)

# Add a nodal lateral load of 50 kips at the left side of the frame
MomentFrame.AddNodeLoad("N2", "FX", 50)

# Analyze the frame
MomentFrame.Analyze()

# Display the model for viewing
Visualization.RenderModel(MomentFrame)

# Print the shear, moment, and deflection diagrams for member 'M3'
MomentFrame.GetMember("M3").PlotShear("Fy")
MomentFrame.GetMember("M3").PlotMoment("Mz")
MomentFrame.GetMember("M3").PlotDeflection("dy")

# Print reactions at the bottom of the frame
print("Left Support X Reaction: {Rxn:.2f} kip".format(Rxn = MomentFrame.GetNode("N1").RxnFX))
print("Left Support Y Reaction: {Rxn:.2f} kip".format(Rxn = MomentFrame.GetNode("N1").RxnFY))
print("Left Support Moment Reaction: {Rxn:.2f} kip-ft".format(Rxn = MomentFrame.GetNode("N1").RxnMZ))

print("Right Support X Reacton: {Rxn:.2f} kip".format(Rxn = MomentFrame.GetNode("N4").RxnFX))
print("Right Support Y Reacton: {Rxn:.2f} kip".format(Rxn = MomentFrame.GetNode("N4").RxnFY))
print("Right Support Moment Reacton: {Rxn:.2f} kip-ft".format(Rxn = MomentFrame.GetNode("N4").RxnMZ))