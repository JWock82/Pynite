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
BracedFrame.AddMember('Brace1', 'N4', 'N2', E, G, Iy, Iz, J, A, tension_only=True)
BracedFrame.DefineReleases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
BracedFrame.DefineReleases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)


# Provide pinned supports at the bases of the columns
BracedFrame.DefineSupport('N1', SupportDX=True, SupportDY=True, SupportDZ=True)
BracedFrame.DefineSupport('N4', SupportDX=True, SupportDY=True, SupportDZ=True)

# Stabilize the frame in the global Z-direction
BracedFrame.DefineSupport('N2', SupportDZ=True)
BracedFrame.DefineSupport('N3', SupportDZ=True)

# Add self weight dead loads to the frame
# Note that we could leave 'x1' and 'x2' undefined below and it would default to the full member length
# Note also that the direction uses lowercase notations to indicate member local coordinate systems
BracedFrame.AddMemberDistLoad('Beam', Direction='Fy', w1=-0.024/12, w2=-0.024/12, x1=0, x2=15*12, case='D')
BracedFrame.AddMemberDistLoad('Col1', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')
BracedFrame.AddMemberDistLoad('Col2', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')

# Add a nodal wind load of 10 kips at the left side of the frame
# Note that the direction uses uppercase notation to indicate model global coordinate system
BracedFrame.AddNodeLoad('N2', Direction='FX', P=50, case='W')

# Create two load combinations
BracedFrame.AddLoadCombo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
BracedFrame.AddLoadCombo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Analyze the braced frame
BracedFrame.Analyze()

# Display the deformed shape of the structure magnified 50 times with the text height 5 model units (inches) high
from PyNite import Visualization
Visualization.RenderModel(BracedFrame, text_height=5, deformed_shape=True, deformed_scale=50, combo_name='1.2D+1.0W')

# Plot the moment diagram for the beam
BracedFrame.GetMember('Beam').PlotMoment('Mz', combo_name='1.2D+1.0W')

# Plot the deflection of the column
BracedFrame.GetMember('Col1').PlotDeflection('dy', combo_name='1.2D+1.0W')

# Find the maximum shear in the first column
print('Column Shear Force:', BracedFrame.GetMember('Col1').MaxShear('Fy', combo_name='1.2D+1.0W'), 'kip')

# Find the frame lateral drift
# Note that the deflections are stored as a dictionary in the node object by load combination, so [] are used instead of () below
print('Frame Lateral Drift:', BracedFrame.GetNode('N2').DX['1.2D+1.0W'], 'in')

# Find the maximum uplift reaction
print('Left Support Y Reaction:', BracedFrame.GetNode('N1').RxnFY['0.9D+1.0W'], 'kip')