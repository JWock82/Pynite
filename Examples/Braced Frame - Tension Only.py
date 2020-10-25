# Example of a basic 2D tension-only braced frame with gravity and lateral
# loads. Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Create a new finite element model
BracedFrame = FEModel3D()

# Add nodes (frame is 15 ft wide x 12 ft tall)
BracedFrame.AddNode('N1', 0, 0, 0)
BracedFrame.AddNode('N2', 0, 12*12, 0)
BracedFrame.AddNode('N3', 15*12, 12*12, 0)
BracedFrame.AddNode('N4', 15*12, 0*12, 0)

# Define column properties (use W10x33 from the AISC Manual):
E = 29000 # ksi
G = 11400 # ksi
Iy = 36.6 # in^4
Iz = 171 # in^4
J = 0.58 # in^4
A = 9.71 # in^2

# Define the columns
BracedFrame.AddMember('Col1', 'N1', 'N2', E, G, Iy, Iz, J, A)
BracedFrame.AddMember('Col2', 'N4', 'N3', E, G, Iy, Iz, J, A)

# Define beam properties (Use W8x24)
Iy = 18.3 # in^4
Iz = 82.7 # in^4
J = 0.346 # in^4
A = 7.08 # in^2

# Define the beams
BracedFrame.AddMember('Beam', 'N2', 'N3', E, G, Iy, Iz, J, A)
BracedFrame.DefineReleases('Beam', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Define the brace properties
# We'll use a section with L/r <= 300 which is a common rule of thumb for
# tension members. We'll use L4x4x1/4.
Iy = 3 # in^4
Iz = 3 # in^4
J = 0.0438 # in^4
A = 1.94 # in^2

# Define the braces
BracedFrame.AddMember('Brace1', 'N1', 'N3', E, G, Iy, Iz, J, A,
                      tension_only=True)
BracedFrame.AddMember('Brace2', 'N4', 'N2', E, G, Iy, Iz, J, A,
                      tension_only=True)

# Release the brace ends to form an axial member
BracedFrame.DefineReleases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
BracedFrame.DefineReleases('Brace2', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Provide pinned supports at the bases of the columns (also restrained about
# the Y-axis for stability)
BracedFrame.DefineSupport('N1', SupportDX=True, SupportDY=True, SupportDZ=True,
                          SupportRY=True)
BracedFrame.DefineSupport('N4', SupportDX=True, SupportDY=True, SupportDZ=True,
                          SupportRY=True)

# Stabilize the frame in the global Z-direction so it doesn't tip over
# out-of-plane.
BracedFrame.DefineSupport('N2', SupportDZ=True)
BracedFrame.DefineSupport('N3', SupportDZ=True)

# Add self weight dead loads to the frame.
# Note that we could leave 'x1' and 'x2' undefined below and it would default
# to the full member length. Note also that the direction uses lowercase
# notations to indicate member local coordinate systems. Brace loads have been
# neglected.
BracedFrame.AddMemberDistLoad('Beam', Direction='Fy', w1=-0.024/12,
                              w2=-0.024/12, x1=0, x2=15*12, case='D')
BracedFrame.AddMemberDistLoad('Col1', Direction='Fx', w1=-0.033/12,
                              w2=-0.033/12, x1=0, x2=12*12, case='D')
BracedFrame.AddMemberDistLoad('Col2', Direction='Fx', w1=-0.033/12,
                              w2=-0.033/12, x1=0, x2=12*12, case='D')

# Add nodal wind loads of 25 kips to each side of the frame. Note that the
# direction uses uppercase notation to indicate model global coordinate
# system.
BracedFrame.AddNodeLoad('N2', Direction='FX', P=25, case='W')
BracedFrame.AddNodeLoad('N3', Direction='FX', P=25, case='W')

# Create load combinations
# Note that the load combination '1.4D' has no lateral load, but does have
# gravity load. The gravity load forces the braces to receive minor
# compression, which causes them to be deactivated on the first iteration.
# Once deactivated the model is unstable and an exception is thrown. This is
# normal and correct behavior. Some programs allow for slight compression in
# tension-only members to avoid this problem. Load combination '1.4D' ahs been
# commented out, but you can uncomment it to see for yourself what happens.

# BracedFrame.AddLoadCombo('1.4D', factors={'D':1.4})
BracedFrame.AddLoadCombo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
BracedFrame.AddLoadCombo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Analyze the braced frame
# P-Delta analysis could also be performed using BracedFrame.Analyze_PDelta().
# Generally, P-Delta analysis will have little effect on a model of a braced
# frame, as there is usually very little bending moment in the members.
BracedFrame.Analyze()

# Display the deformed shape of the structure magnified 50 times with the text
# height 5 model units (inches) high
from PyNite import Visualization
Visualization.RenderModel(BracedFrame, text_height=5, deformed_shape=True,
                          deformed_scale=50, combo_name='1.2D+1.0W')

# Plot the axial load diagrams for the braces. We should see no compression on
# 'Brace2' and 64 kips on 'Brace1' if the tension-only analysis worked
# correctly.
BracedFrame.GetMember('Brace1').PlotAxial(combo_name='1.2D+1.0W')
BracedFrame.GetMember('Brace2').PlotAxial(combo_name='1.2D+1.0W')

# Report the frame reactions for the load combination '1.2D+1.0W'. We should
# see a -50 kip horizontal reaction at node 'N1', and a zero kip reaction at
# node 'N4' if the tension-only analysis worked correctly. Similarly, we
print('1.2D+1.0W: N1 reaction FY =',
      '{:.3f}'.format(BracedFrame.GetNode('N1').RxnFY['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N1 reaction FX =',
      '{:.3f}'.format(BracedFrame.GetNode('N1').RxnFX['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N4 reaction FY =',
      '{:.3f}'.format(BracedFrame.GetNode('N4').RxnFY['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N4 reaction FX =',
      '{:.3f}'.format(BracedFrame.GetNode('N4').RxnFX['1.2D+1.0W']), 'kip')