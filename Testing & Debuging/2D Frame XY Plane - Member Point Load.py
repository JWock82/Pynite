from PyNite import FEModel3D
from PyNite.Visualization import RenderModel

# Create a new finite element model
frame = FEModel3D()

# Add nodes
frame.AddNode('N1', 0, 0, 0) # ft
frame.AddNode('N2', 0, 7.667, 0) # ft
frame.AddNode('N3', 7.75, 7.667, 0) # ft
frame.AddNode('N4', 7.75, 0, 0) # ft

# Add supports
frame.DefineSupport('N1', True, True, True, True, True, False)
frame.DefineSupport('N4', True, True, True, True, True, False)

# Define material and section properties for a W8x24
E = 29000*12**2 # ksf
G = 1111200*12**2 # ksf
Iy = 18.3/12**4 # ft^4
Iz = 82.7/12**4 # ft^4
J = 0.346/12**4 # ft^4
A = 5.26/12**2 # in^2

# Define members
frame.AddMember('M1', 'N1', 'N2', E, G, Iy, Iz, J, A)
frame.AddMember('M2', 'N2', 'N3', E, G, Iy, Iz, J, A)
frame.AddMember('M3', 'N4', 'N3', E, G, Iy, Iz, J, A)

# Add loads to the frame
frame.AddMemberPtLoad('M2', 'Fy', -5, 7.75/2) # 5 kips @ midspan
frame.AddMemberDistLoad('M2', 'Fy', -0.024, -0.024) # W8x24 self-weight

# Analyze the frame
frame.Analyze()

# Render the model
RenderModel(frame, text_height=0.2, deformed_shape=True, deformed_scale=100, render_loads=True)

frame.GetMember('M2').PlotMoment('Mz')
frame.GetMember('M3').PlotDeflection('dy')
print(frame.GetNode('N1').RZ)