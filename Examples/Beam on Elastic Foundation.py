# -*- coding: utf-8 -*-
'''
This example shows how to generate a beam on elastic foundation by using spring supports. All units
in this model are expressed in terms of kips (force) and inches (length).
'''

# Create a new model for the beam
from PyNite import FEModel3D
boef = FEModel3D()

# Assign the length of the beam
L = 25*12  # in

# The beam will be modeled as a series of elements at a 12 inch spacing. A support spring will be
# added at the end of each segment.
num_segs = round(L / 12)
num_nodes = num_segs + 1
L_seg = L/num_segs

# Define the spring stiffness
ks = 22.5  # k/in

# Generate the nodes
for i in range(num_nodes):

    # Add nodes spaced at 15"
    boef.add_node('N' + str(i + 1), i*L_seg, 0, 0)

    # Add supports to the nodes
    if i == 0 or i == num_nodes - 1:
        # Pinned supports at beam ends
        boef.def_support('N' + str(i + 1), True, True, True, True, False, False)
    else:
        # Spring supports at all other locations
        boef.def_support_spring('N' + str(i + 1), 'DY', ks, '-')

# Define member material properties (W8x35)
E = 29000   # ksi
G = 11200   # ksi
A = 10.3    # in^2
Iz = 127    # in^4 (strong axis)
Iy = 42.6   # in^4 (weak axis)
J = 0.769   # in^4

# Define members
for i in range(num_segs):

    # Add the members
    boef.add_member('M' + str(i + 1), 'N' + str(i + 1), 'N' + str(i + 2), E, G, Iy, Iz, J, A)
        
# Add a point load at midspan
if num_nodes % 2 == 0:  # Checks if there is a node at midspan or a member between nodes
    mid_member = 'M' + str(int(num_segs/2))
    boef.add_member_pt_load(mid_member, 'Fy', -40, L_seg/2)
else:
    mid_node = 'N' + str(int(num_nodes/2))
    boef.add_node_load(mid_node, 'FY', -40)

# Analyze the model
boef.analyze()

# Render the mdoel with the deformed shape
from PyNite.Visualization import render_model
render_model(boef, text_height=1.5, deformed_shape=True)

# Find and print the maximum displacement
d_max = min([node.DY['Combo 1'] for node in boef.Nodes.values()])
print('Maximum displacement: ', round(d_max, 4), 'in')

# Find and print the minimum moment
M_max = min([segment.min_moment('Mz') for segment in boef.Members.values()])
print('Minimum moment: ', round(M_max), 'k-in')