# Example of a basic 2D tension-only braced frame with gravity and lateral loads
# Units used for the model in this example are inches and kips

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

# Define the braces
Iy = 10
Iz = 10
J = 0.3
A = 5
BracedFrame.AddMember('Brace1', 'N1', 'N3', E, G, Iy, Iz, J, A, tension_only=True)
BracedFrame.AddMember('Brace2', 'N4', 'N2', E, G, Iy, Iz, J, A, tension_only=True)

# Release the brace ends to form an axial member
BracedFrame.DefineReleases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
BracedFrame.DefineReleases('Brace2', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Provide pinned supports at the bases of the columns (also restrained about the Y-axis for stability)
BracedFrame.DefineSupport('N1', SupportDX=True, SupportDY=True, SupportDZ=True, SupportRY=True)
BracedFrame.DefineSupport('N4', SupportDX=True, SupportDY=True, SupportDZ=True, SupportRY=True)

# Stabilize the frame in the global Z-direction
BracedFrame.DefineSupport('N2', SupportDZ=True)
BracedFrame.DefineSupport('N3', SupportDZ=True)

# Add self weight dead loads to the frame
# Note that we could leave 'x1' and 'x2' undefined below and it would default to the full member length
# Note also that the direction uses lowercase notations to indicate member local coordinate systems
# BracedFrame.AddMemberDistLoad('Beam', Direction='Fy', w1=-0.024/12, w2=-0.024/12, x1=0, x2=15*12, case='D')
# BracedFrame.AddMemberDistLoad('Col1', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')
# BracedFrame.AddMemberDistLoad('Col2', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')

# Add a nodal wind load of 10 kips at the left side of the frame
# Note that the direction uses uppercase notation to indicate model global coordinate system
BracedFrame.AddNodeLoad('N2', Direction='FX', P=25, case='W')
BracedFrame.AddNodeLoad('N3', Direction='FX', P=25, case='W')

# Create two load combinations
BracedFrame.AddLoadCombo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
BracedFrame.AddLoadCombo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Display the deformed shape of the structure magnified 50 times with the text height 5 model units (inches) high
from PyNite import Visualization
Visualization.RenderModel(BracedFrame, text_height=5, combo_name='1.2D+1.0W')

# Analyze the braced frame
BracedFrame.Analyze_PDelta()

# Display the deformed shape of the structure magnified 50 times with the text height 5 model units (inches) high
from PyNite import Visualization
Visualization.RenderModel(BracedFrame, text_height=5, deformed_shape=True, deformed_scale=50, combo_name='1.2D+1.0W')

# Plot the axial load diagrams for the braces
BracedFrame.GetMember('Brace1').PlotAxial(combo_name='1.2D+1.0W')
BracedFrame.GetMember('Brace2').PlotAxial(combo_name='1.2D+1.0W')

# Report the frame reactions
print('1.2D+1.0W: N1 reaction FY = ', BracedFrame.GetNode('N1').RxnFY['1.2D+1.0W'])
print('1.2D+1.0W: N4 reaction FY = ', BracedFrame.GetNode('N4').RxnFY['1.2D+1.0W'])
print('0.9D+1.0W: N1 reaction FY = ', BracedFrame.GetNode('N1').RxnFY['0.9D+1.0W'])
print('0.9D+1.0W: N4 reaction FY = ', BracedFrame.GetNode('N4').RxnFY['0.9D+1.0W'])

# Report the frame deflection at Node 3
print('1.2D+1.0W: Frame drift = ', BracedFrame.GetNode('N3').DX['1.2D+1.0W'])
print('0.9D+1.0W: Frame drift = ', BracedFrame.GetNode('N3').DX['0.9D+1.0W'])