# Example of a basic 2D tension-only braced frame with gravity and lateral
# loads. Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Create a new finite element model
BracedFrame = FEModel3D()

# Add nodes (frame is 15 ft wide x 12 ft tall)
BracedFrame.add_node('N1', 0, 0, 0)
BracedFrame.add_node('N2', 0, 12*12, 0)
BracedFrame.add_node('N3', 15*12, 12*12, 0)
BracedFrame.add_node('N4', 15*12, 0*12, 0)

# Define column properties (use W10x33 from the AISC Manual):
E = 29000 # ksi
G = 11400 # ksi
Iy = 36.6 # in^4
Iz = 171 # in^4
J = 0.58 # in^4
A = 9.71 # in^2

# Define the columns
BracedFrame.add_member('Col1', 'N1', 'N2', E, G, Iy, Iz, J, A)
BracedFrame.add_member('Col2', 'N4', 'N3', E, G, Iy, Iz, J, A)

# Define beam properties (Use W8x24)
Iy = 18.3 # in^4
Iz = 82.7 # in^4
J = 0.346 # in^4
A = 7.08 # in^2

# Define the beams
BracedFrame.add_member('Beam', 'N2', 'N3', E, G, Iy, Iz, J, A)
BracedFrame.def_releases('Beam', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Define the brace properties
# We'll use a section with L/r <= 300 which is a common rule of thumb for
# tension members. We'll use L4x4x1/4.
Iy = 3 # in^4
Iz = 3 # in^4
J = 0.0438 # in^4
A = 1.94 # in^2

# Define the braces
BracedFrame.add_member('Brace1', 'N1', 'N3', E, G, Iy, Iz, J, A,
                      tension_only=True)
BracedFrame.add_member('Brace2', 'N4', 'N2', E, G, Iy, Iz, J, A,
                      tension_only=True)

# Release the brace ends to form an axial member
BracedFrame.def_releases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
BracedFrame.def_releases('Brace2', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Provide pinned supports at the bases of the columns (also restrained about
# the Y-axis for stability)
BracedFrame.def_support('N1', support_DX=True, support_DY=True, support_DZ=True,
                          support_RY=True)
BracedFrame.def_support('N4', support_DX=True, support_DY=True, support_DZ=True,
                          support_RY=True)

# Stabilize the frame in the global Z-direction so it doesn't tip over
# out-of-plane.
BracedFrame.def_support('N2', support_DZ=True)
BracedFrame.def_support('N3', support_DZ=True)

# Add self weight dead loads to the frame.
# Note that we could leave 'x1' and 'x2' undefined below and it would default
# to the full member length. Note also that the direction uses lowercase
# notations to indicate member local coordinate systems. Brace loads have been
# neglected.
BracedFrame.add_member_dist_load('Beam', Direction='Fy', w1=-0.024/12,
                              w2=-0.024/12, x1=0, x2=15*12, case='D')
BracedFrame.add_member_dist_load('Col1', Direction='Fx', w1=-0.033/12,
                              w2=-0.033/12, x1=0, x2=12*12, case='D')
BracedFrame.add_member_dist_load('Col2', Direction='Fx', w1=-0.033/12,
                              w2=-0.033/12, x1=0, x2=12*12, case='D')

# Add nodal wind loads of 25 kips to each side of the frame. Note that the
# direction uses uppercase notation to indicate model global coordinate
# system.
BracedFrame.add_node_load('N2', Direction='FX', P=25, case='W')
BracedFrame.add_node_load('N3', Direction='FX', P=25, case='W')

# Create load combinations
# Note that the load combination '1.4D' has no lateral load, but does have
# gravity load. The gravity load forces the braces to receive minor
# compression, which causes them to be deactivated on the first iteration.
# Once deactivated the model is unstable and an exception is thrown. This is
# normal and correct behavior. Some programs allow for slight compression in
# tension-only members to avoid this problem. Load combination '1.4D' ahs been
# commented out, but you can uncomment it to see for yourself what happens.

# BracedFrame.add_load_combo('1.4D', factors={'D':1.4})
BracedFrame.add_load_combo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
BracedFrame.add_load_combo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Analyze the braced frame
# P-Delta analysis could also be performed using BracedFrame.analyze_PDelta().
# Generally, P-Delta analysis will have little effect on a model of a braced
# frame, as there is usually very little bending moment in the members.
BracedFrame.analyze()

# Display the deformed shape of the structure magnified 50 times with the text
# height 5 model units (inches) high
from PyNite import Visualization
Visualization.render_model(BracedFrame, annotation_size=5, deformed_shape=True,
                          deformed_scale=50, combo_name='1.2D+1.0W')

# Plot the axial load diagrams for the braces. We should see no compression on
# 'Brace2' and 64 kips on 'Brace1' if the tension-only analysis worked
# correctly.
BracedFrame.Members['Brace1'].plot_axial(combo_name='1.2D+1.0W')
BracedFrame.Members['Brace2'].plot_axial(combo_name='1.2D+1.0W')

# Report the frame reactions for the load combination '1.2D+1.0W'. We should
# see a -50 kip horizontal reaction at node 'N1', and a zero kip reaction at
# node 'N4' if the tension-only analysis worked correctly. Similarly, we
print('1.2D+1.0W: N1 reaction FY =',
      '{:.3f}'.format(BracedFrame.Nodes['N1'].RxnFY['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N1 reaction FX =',
      '{:.3f}'.format(BracedFrame.Nodes['N1'].RxnFX['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N4 reaction FY =',
      '{:.3f}'.format(BracedFrame.Nodes['N4'].RxnFY['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N4 reaction FX =',
      '{:.3f}'.format(BracedFrame.Nodes['N4'].RxnFX['1.2D+1.0W']), 'kip')