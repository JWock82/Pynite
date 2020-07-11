# Example of a basic 2D moment frame with gravity and lateral loads
# Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Create a new finite element model
MomentFrame = FEModel3D()

# Add nodes (frame is 15 ft wide x 12 ft tall)
MomentFrame.AddNode('N1', 0, 0, 0)
MomentFrame.AddNode('N2', 0, 12*12, 0)
MomentFrame.AddNode('N3', 15*12, 12*12, 0)
MomentFrame.AddNode('N4', 15*12, 0*12, 0)

# Define column properties (use W10x33 from the AISC Manual):
E = 29000 # ksi
G = 11400 # ksi
Iy = 36.6 # in^4
Iz = 171 # in^4
J = 0.58 # in^4
A = 9.71 # in^2

# Define the columns
MomentFrame.AddMember('Col1', 'N1', 'N2', E, G, Iy, Iz, J, A)
MomentFrame.AddMember('Col2', 'N4', 'N3', E, G, Iy, Iz, J, A)

# Define beam properties (Use W8x24)
Iy = 18.3 # in^4
Iz = 82.7 # in^4
J = 0.346 # in^4
A = 7.08 # in^2

# Define the beams
MomentFrame.AddMember('Beam', 'N2', 'N3', E, G, Iy, Iz, J, A)

# Provide fixed supports at the bases of the columns
MomentFrame.DefineSupport('N1', SupportDX=True, SupportDY=True, SupportDZ=True, SupportRX=True, SupportRY=True, SupportRZ=True)
MomentFrame.DefineSupport('N4', SupportDX=True, SupportDY=True, SupportDZ=True, SupportRX=True, SupportRY=True, SupportRZ=True)

# Add self weight dead loads to the frame
# Note that we could leave 'x1' and 'x2' undefined below and it would default to the full member length
# Note also that the direction uses lowercase notations to indicate member local coordinate systems
MomentFrame.AddMemberDistLoad('Beam', Direction='Fy', w1=-0.024/12, w2=-0.024/12, x1=0, x2=15*12, case='D')
MomentFrame.AddMemberDistLoad('Col1', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')
MomentFrame.AddMemberDistLoad('Col2', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')

# Add a nodal wind load of 10 kips at the left side of the frame
# Note that the direction uses uppercase notation to indicate model global coordinate system
MomentFrame.AddNodeLoad('N2', Direction='FX', P=10, case='W')

# Create two load combinations
MomentFrame.AddLoadCombo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
MomentFrame.AddLoadCombo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Perform a P-Big Delta analysis on the frame
# Note that to capture P-little delta effects the members should idealy be broken into three segments each
MomentFrame.Analyze_PDelta()

# A first-order analysis could be done with the following line instead
# MomentFrame.Analyze()

# Display the deformed shape of the structure magnified 50 times with the text height 5 model units (inches) high
from PyNite import Visualization
Visualization.RenderModel(MomentFrame, text_height=5, deformed_shape=True, deformed_scale=50, combo_name='1.2D+1.0W')

# Plot the moment diagram for the beam
MomentFrame.GetMember('Beam').PlotMoment('Mz', combo_name='1.2D+1.0W')

# Plot the deflection of the column
MomentFrame.GetMember('Col1').PlotDeflection('dy', combo_name='1.2D+1.0W')

# Find the maximum shear in the first column
print('Column Shear Force:', MomentFrame.GetMember('Col1').MaxShear('Fy', combo_name='1.2D+1.0W'), 'kip')

# Find the frame lateral drift
# Note that the deflections are stored as a dictionary in the node object by load combination, so [] are used instead of () below
print('Frame Lateral Drift:', MomentFrame.GetNode('N2').DX['1.2D+1.0W'], 'in')

# Find the maximum uplift reaction
print('Left Support Y Reaction:', MomentFrame.GetNode('N1').RxnFY['0.9D+1.0W'], 'kip')