'''
This example shows how to use the repair_model function to automatically detect
and insert intermediate member nodes when generating a 2d braced frame . All
units in this model are expressed in terms of kips (force) and inches (length).
'''

from PyNite import FEModel3D

# Create a new finite element model
f = FEModel3D()

# Define Shape
length_verticals = 12*12
length_horizontals = 12*12
count_horizontals = 4
count_verticals = 4
strong_params = {'E': 29000, 'G': 11200, 'Iy': 11.4, 'Iz': 118, 'J': 0.239, 'A': 6.49}  # W10X22

# Add Horizontals
releases_horizontals = {'Ryi': True, 'Rzi': True, 'Ryj': True, 'Rzj': True, }
names_horizontals = []
for i in range(count_horizontals):
    y = i / (count_horizontals - 1) * length_verticals
    name_node_left = f.add_node(None, 0, y, 0)
    name_node_right = f.add_node(None, length_horizontals, y, 0)
    name_member = f.add_member(f'H{i}', name_node_left, name_node_right, **strong_params)
    f.def_releases(name_member, **releases_horizontals)
    names_horizontals.append(name_member)

# Add Verticals
names_verticals = []
for i in range(count_verticals):
    x = i / (count_verticals - 1) * length_horizontals
    name_node_low = f.add_node(None, x, 0, 0)
    name_node_high = f.add_node(None, x, length_verticals, 0)
    name_member = f.add_member(f'V{i}', name_node_low, name_node_high, **strong_params)
    names_verticals.append(name_member)

# Add and define supports
support_node_left = f.find_node_by_coordinates((0, 0, 0))
support_node_right = f.find_node_by_coordinates((length_horizontals, 0, 0))

f.def_support(
    support_node_left,
    support_DX=True, support_DY=True, support_DZ=True,
    support_RX=True, support_RY=True, support_RZ=True)
f.def_support(
    support_node_right,
    support_DX=True, support_DY=True, support_DZ=True,
    support_RX=True, support_RY=True, support_RZ=True)

# Add a distributed load. Should be on the horizontal closest to y=0.
member = list(f.Members.values())[0]
f.add_member_dist_load(
    Member=member.name,
    Direction='Fz',
    w1=0.2, w2=0.8,
    x1=0.2, x2=0.8*member.L())

# Add a distributed load. Should be the vertical closest to x=0.
member = list(f.Members.values())[count_horizontals]
f.add_member_dist_load(
    Member=member.name,
    Direction='Fz',
    w1=0.3, w2=0.7,
    x1=0.3, x2=0.7*member.L())


# Add a member point load. Should be on the horizontal furthest from y=0.
member = list(f.Members.values())[count_horizontals-1]
f.add_member_pt_load(
    Member=member.name,
    Direction='Fz',
    P=3.14, x=0.25*member.L(),
    case='Case 1')

# Add a member point load. Should be on the vertical furthest from x=0.
member = list(f.Members.values())[-1]
f.add_member_pt_load(
    Member=member.name,
    Direction='Fz',
    P=1.57, x=0.25*member.L(),
    case='Case 1')

# Clean-up
f.repair(
    merge_duplicates=True,  # must be run to avoid duplicate nodes
    tolerance_node_distance=1e-3,  # maximum 3d distance between coordinates to determine duplicate node
    tolerance_intersection=1e-3,  # analagous to minimum intersection angle; equivalent to sin-1(angle)
    )

# Analyze
f.analyze(log=False, check_statics=True)
theoretical_members = 2 * count_horizontals * count_verticals - count_horizontals - count_verticals
actual_members = len(f.Members)

# Display the deformed shape of the structure magnified 1 times with the text
# height 2 model units (inches) high.
from PyNite import Visualization
Visualization.render_model(f, annotation_size=2, deformed_shape=True, deformed_scale=1)
