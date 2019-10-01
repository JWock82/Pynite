# -*- coding: utf-8 -*-
# Example of a simply supported beam with a point load.
# Units used in this example are metric (SI)

import numpy as np

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# 1: No beam-reference point
# 2: Ref @ [0,0,-10] - should match Case 1
# 3: Ref @ [0,-10,0] - twists the beam through 90 degrees

# *** Puzzle: Why do the node reactions flip upside down in Case 3 only? ***
# (Actually, they go to zero for [0,-10,10].)

### Commenting out the next 3 lines to isolate if 'beamref' is causing the problem
# case = 3

# if case == 1:
#     beamref = None
# if case == 2:
#     beamref = np.asarray([0,0,-10])
# if case == 3:
#     beamref = np.asarray([0,-10,0])

# Create a new finite element model
SimpleBeam = FEModel3D()

# Add nodes (3 metres apart on x-axis)
L = 3
SimpleBeam.AddNode("N1", -L/2, 0, 0)
SimpleBeam.AddNode("N2",  0,   0, 0)
SimpleBeam.AddNode("N3",  L/2, 0, 0)

# Add a beam with a solid rectangular section:
E = 209E9  # Young's elastic modulus             ### 209,000,000,000 kgf/m^2
G = 79.3E9 # shear modulus
b = 30E-2  # breadth (y-thickness)               ### 0.3 m
d = 1E-2   # depth   (z-thickness)               ### 0.01 m
Iyy = b**3 * d / 12 # 2nd moment of area
Izz = b * d**3 / 12 # 2nd moment of area
J = Iyy + Izz       # 2nd polar moment of area 
A = b * d           # cross-sectional area
SimpleBeam.AddMember("M1", "N1", "N2", E, G, Iyy, Izz, J, A) ### I removed the 'beamref' argument to isolate whether this is causing the problem
SimpleBeam.AddMember("M2", "N2", "N3", E, G, Iyy, Izz, J, A) ### I removed the 'beamref' argument to isolate whether this is causing the problem

# Provide simple supports - can't displace or twist at either end
SimpleBeam.DefineSupport("N1", True, True, True, True, False, False)
SimpleBeam.DefineSupport("N3", True, True, True, True, False, False)

# Add a point load at the midspan of the beam acting on the strong and weak axes
F = -1E3
SimpleBeam.AddNodeLoad("N2", "FY", F)
SimpleBeam.AddNodeLoad("N2", "FZ", F)

# Expected moment
Mz = -F * L / 4
My = -F * L / 4

# Calculate and print expected results
DY = F * L**3 / (48 * E * Izz)
DZ = F * L**3 / (48 * E * Iyy)
print('Node N2 expected Y-displacement = {d:.2f} mm; Expected max. z-z bending moment = {m:.2f}'.format(d=(DY * 1E3), m=Mz))
print('Node N2 expected Z-displacement = {d:.2f} mm; Expected max. y-y bending moment = {m:.2f}'.format(d=(DZ * 1E3), m=My))

# Analyze the beam
SimpleBeam.Analyze()

# Print calculated node displacements
print('Node N2 Y-displacement from analysis = {d:.2f} mm'.format(d=(SimpleBeam.GetNode("N2").DY * 1E3)))
print('Node N2 Z-displacement from analysis = {d:.2f} mm'.format(d=(SimpleBeam.GetNode("N2").DZ * 1E3)))

# Print the shear, moment, and deflection diagrams
SimpleBeam.GetMember("M1").PlotShear("Fz")
SimpleBeam.GetMember("M1").PlotMoment("My")
SimpleBeam.GetMember("M1").PlotDeflection("dz")

# Print reactions at each end of the beam
print("Left Support Reaction: {Rxn:.2f} N".format(Rxn = SimpleBeam.GetNode("N1").RxnFZ))
print("Right Support Reacton: {Rxn:.2f} N".format(Rxn = SimpleBeam.GetNode("N3").RxnFZ))
