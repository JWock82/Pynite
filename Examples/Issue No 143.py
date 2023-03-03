# Example of a simply supported beam with a uniform distributed load.
# Units used in this example are inches and kips
# This example does not use load combinations. The program will create a default load combindation called 'Combo 1'

# Import FEModel3D from PyNite
from PyNite import FEModel3D

# Create a new finite element model
SimpleBeam = FEModel3D()

# Add nodes (14 ft = 168 in apart)
SimpleBeam.add_node('N1', 0, 0, 0)
SimpleBeam.add_node('N2', 168, 0, 0)

# Add a beam with the following properties:
# E = 29000 ksi, G = 11400 ksi, Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 20 in^2
SimpleBeam.add_member('M1', 'N1', 'N2', 29000, 11400, 100, 150, 250, 20)

# Provide supports
SimpleBeam.def_support('N1', True, True, True, True, True, True)
SimpleBeam.def_support_spring('N2', 'RZ', 1, None)

# Add a uniform load of 200 lbs/ft to the beam
SimpleBeam.add_member_dist_load('M1', 'Fy', -200/1000/12, -200/1000/12, 0, 168)

# Analyze the beam
SimpleBeam.analyze()

# Print the shear, moment, and deflection diagrams
SimpleBeam.Members['M1'].plot_shear('Fy')
SimpleBeam.Members['M1'].plot_moment('Mz')
SimpleBeam.Members['M1'].plot_deflection('dy')

# Print reactions at each end of the beam
print('Left Support Reaction:', SimpleBeam.Nodes['N1'].RxnFY, 'kip')
print('Right Support Reacton:', SimpleBeam.Nodes['N2'].RxnFY, 'kip')
# Render the deformed shape of the beam magnified 100 times, with a text height of 5 inches

from PyNite import Visualization
Visualization.render_model(SimpleBeam, annotation_size=5, deformed_shape=True, deformed_scale=100, render_loads=True)