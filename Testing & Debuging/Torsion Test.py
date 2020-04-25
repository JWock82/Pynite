# Testing torsional point loads
# Units used in this example are inches, and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Create a new finite element model
TorqueBeam = FEModel3D()

# Add nodes (14 ft = 168 in apart)
TorqueBeam.AddNode('N1', 0, 0, 0)
TorqueBeam.AddNode('N2', 168, 0, 0)

# Add a beam with the following properties:
# E = 29000 ksi, G = 11400 ksi, Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 20 in^2
TorqueBeam.AddMember('M1', 'N1', 'N2', 29000, 11400, 100, 150, 250, 20)

# Provide fixed supports
TorqueBeam.DefineSupport('N1', False, True, True, True, True, True)
TorqueBeam.DefineSupport('N2', True, True, True, True, True, True)

# Add a point load of 5 kip-ft and 10 kip-ft at 3ft and 11 ft along the beam respectively
TorqueBeam.AddMemberPtLoad('M1', 'Mx', 5, 3*12)
TorqueBeam.AddMemberPtLoad('M1', 'Mx', 10, 11*12)

# Analyze the beam and perform a statics check
TorqueBeam.Analyze(check_statics=True)

from PyNite import Visualization
Visualization.RenderModel(TorqueBeam, text_height=5, deformed_shape=True, deformed_scale=30, render_loads=True)

# Print the torsion diagram
TorqueBeam.GetMember('M1').PlotTorsion()

# Print reactions at each end of the beam
print('**Individual Support Reactions**')
print('Left Support Reaction: {Rxn:.2f} kip'.format(Rxn = TorqueBeam.GetNode('N1').RxnMX['Combo 1']))
print('Right Support Reacton: {Rxn:.2f} kip'.format(Rxn = TorqueBeam.GetNode('N2').RxnMX['Combo 1']))

# Print the max/min torques on the beam
print('**Member Internal Torsion**')
print('Maximum Torque: {Rxn:.2f} kip-ft'.format(Rxn = TorqueBeam.GetMember('M1').MaxTorsion()))
print('Minimum Torque: {Rxn:.2f} kip-ft'.format(Rxn = TorqueBeam.GetMember('M1').MinTorsion()))
