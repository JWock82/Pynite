# Example of a basic 2D tension-only braced frame with gravity and lateral
# loads. Units used for the model in this example are inches and kips.

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

# Define a brace (tension and compression - both ways)
BracedFrame.AddMember('Brace1', 'N1', 'N3', E, G, Iy, Iz, J, A)

# Let's add spring supports to the base of the structure. We'll add a couple of
# extra nodes at the base of the structure that will receive the springs. The
# lengths of these springs is irrelevant, since we're defining a spring constant
# that is independent of the spring's length. Only the direction of the spring
# matters as it defines the direction of the spring's stiffness. The nodes will
# be directly below N1 and N4 in the Y-direction.
BracedFrame.AddNode('N1s', 0, -2*12, 0)
BracedFrame.AddNode('N4s', 15*12, -2*12, 0)

BracedFrame.AddSpring('Spring1','N1', 'N1s', 10000, tension_only=True,
                      comp_only=False)
BracedFrame.AddSpring('Spring2','N4', 'N4s', 10000, tension_only=False,
                      comp_only=True) # The structure would be unstable if
                                      # this was tension only

# Release the brace ends to form an axial member
BracedFrame.DefineReleases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Springs only carry axial loads, nothing else, so we'll need to stabilize
# the column bases in the other directions. The column bases will be
# supported by the springs vertically. For the other directions (horizontally
# and about the Y-axis) we'll need to provide supports.
BracedFrame.DefineSupport('N1', SupportDX=True, SupportDZ=True, SupportRY=True)
BracedFrame.DefineSupport('N4', SupportDX=True, SupportDZ=True, SupportRY=True)

# Fix the nodes supporting the bottoms of the springs. Note that even though
# we're fixing these nodes, the only reactions the supports will carry will
# be in the Y-direction, due to the fact that the spring only has stiffness in
# that direction. We fix the node so that it's not free to spin or translate
# in the other directions however. If we didn't the node would be unstable and
# the model would crash. PyNite is unforgiving in this regard. Every degree of
# freedom (3 translations and 3 rotations) at every node must be stabilized so
# it's not free to move infinitely.
BracedFrame.DefineSupport('N1s', SupportDX=True, SupportDY=True, SupportDZ=True, SupportRX=True, SupportRY=True, SupportRZ=True)
BracedFrame.DefineSupport('N4s', SupportDX=True, SupportDY=True, SupportDZ=True, SupportRX=True, SupportRY=True, SupportRZ=True)

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
# gravity load. The gravity load forces the tension only spring to receive
# minor compression, which causes it to be deactivated on the first iteration.
# Once deactivated the model is unstable and an exception is thrown. This is
# normal and correct behavior. Load combination '1.4D' has been commented out,
# but you can uncomment it to see for yourself what happens.

# BracedFrame.AddLoadCombo('1.4D', factors={'D':1.4})
BracedFrame.AddLoadCombo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
BracedFrame.AddLoadCombo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Analyze the braced frame
# P-Delta analysis could also be performed using BracedFrame.Analyze_PDelta().
# Generally, P-Delta analysis will have little effect on a model of a braced
# frame, as there is usually very little bending moment in the members.
BracedFrame.Analyze()

# Display the deformed shape of the structure magnified 50 times with the text
# height 5 model units (inches) high.
from PyNite import Visualization
Visualization.RenderModel(BracedFrame, text_height=5, deformed_shape=True,
                          deformed_scale=50, combo_name='1.2D+1.0W')

# We should see upward displacement at N1 and downward displacement at N4 if
# our springs worked correctly
print('N1 displacement in Y =', BracedFrame.GetNode('N1').DY['1.2D+1.0W'])
print('N4 displacement in Y =', BracedFrame.GetNode('N4').DY['1.2D+1.0W'])
