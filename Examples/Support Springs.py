from PyNite import FEModel3D

# Create a new finite element model for a beam
beam = FEModel3D()

# Add nodes to the model
# All units will be in kips and inches
beam.add_node('N1', 0, 0, 0)
beam.add_node('N2', 10*12, 0, 0)

# Add a member to the model and a member point load
beam.add_member('M1', 'N1', 'N2', 29000, 11200, 100, 100, 200, 20)
beam.add_member_pt_load('M1', 'Fy', -10, 5*12)

# Define simple supports, except with a spring support at each end in the 'Y' direction.
# The stiffness on the left side will be 2.5 k/in. The stiffness on the right will be 5 k/in.
beam.def_support('N1', True, 2.5, True, True, False, False)
beam.def_support('N2', True, 5, True, False, False, False)

# Analyze the beam
beam.analyze(log=True)

# Print out support reactions
print('FY1 = ', beam.Nodes['N1'].RxnFY['Combo 1'])
print('FY2 = ', beam.Nodes['N2'].RxnFY['Combo 1'])

# Print out support displacements
print('DY1 = ', beam.Nodes['N1'].DY['Combo 1'])
print('DY2 = ', beam.Nodes['N2'].DY['Combo 1'])

# Render the model
from PyNite.Visualization import RenderModel
RenderModel(beam, text_height=3, deformed_shape=True)