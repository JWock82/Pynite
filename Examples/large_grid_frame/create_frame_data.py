## -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 tamalone1
"""
import numpy as np
import itertools, csv, os

# Nodes coordinates in a 3D rectangle
x_values = np.linspace(0, 10, 5)
y_values = np.linspace(0, 50, 11)
z_values = np.linspace(0, 10, 5)

# Create a 3D grid of nodes (list of lists of lists)
grid = []
# A series number for each node
id_number = itertools.count()
# Take each x-coordinate
for i, x in enumerate(x_values):
    # Each element is a list of y-coordinates
    grid.append([])
    for j, y in enumerate(y_values):
        # Each element is a list of z-coordinates
        grid[i].append([])
        for k, z in enumerate(z_values):
            # Pack together the corresponding name and coordinates
            grid[i][j].append(('N'+str(next(id_number)), x, y, z))

# Write the nodes to a CSV file
filename = os.path.join(os.path.dirname(__file__), 'gridnodes.csv')
with open(filename, mode='w', newline='') as f:
    csv_writer = csv.writer(f)
    nodes_at_top_level = []
    # Write the header
    csv_writer.writerow(['Name', 'X', 'Y', 'Z'])
    # Write the node list to the file
    # chain() unpacks the first two layers (x- and y-coordinates)
    for nodes_xy in itertools.chain(*grid):
        # for-loop unpacks the third layer (z-coordinates)
        for node_z in nodes_xy:
            # Each item is a tuple, writerows() unpacks it into the file
            csv_writer.writerow(node_z)
            # if node y-value is the maximum
            if node_z[2] == max(y_values):
                # keep the node name for later
                nodes_at_top_level.append(node_z[0])

# Connect each point to the adjacent one
connections = {}
for row_idx, row in enumerate(grid):
    for col_idx, col in enumerate(row):
        for z_idx, value in enumerate(col):
            connections[value[0]] = []
            try:
                # Get value from next column
                connections[value[0]].append(grid[row_idx][col_idx + 1][z_idx][0])
            except IndexError:
                # Ignore IndexError caused by hitting the end of the line
                pass
            try:
                # Get value from next row
                connections[value[0]].append(grid[row_idx + 1][col_idx][z_idx][0])
            except IndexError:
                # Ignore IndexError caused by hitting the end of the line
                pass
            try:
                # Get value from next z-value
                connections[value[0]].append(grid[row_idx][col_idx][z_idx + 1][0])
            except IndexError:
                # Ignore IndexError caused by hitting the end of the line
                pass

member_list = []
member_id = itertools.count()
for first_node, connected_nodes in connections.items():
    for second_node in connected_nodes:
        member_list.append(('m' + str(next(member_id)), first_node, second_node))

filename = os.path.join(os.path.dirname(__file__), 'gridmembers.csv')
with open(filename, mode='w', newline='') as f:
    csv_writer = csv.writer(f)
    # Write the header
    csv_writer.writerow(['Name', 'iNode', 'jNode'])
    # Write the member list to the file
    csv_writer.writerows(member_list)


def filter_node_coordinates(nodes, coordinate, value):
    """ Return the nodes with coordinate == value. 
    
    node: Node3D, or object with coordinates
    coordinate: str, coordinate name
    value: value with appropriate type
    """
    # # Create list of matching results
    # matches = []
    # # Unpack nodes (iterable) and check each element
    # for node in nodes:
    #     if getattr(node, coordinate) == value:
    #         match.append(node)
    # return matches
    return getattr(node, coordinate) == value


# Create a file of supports
filename = os.path.join(os.path.dirname(__file__), 'gridsupports.csv')
with open(filename, mode='w', newline='') as f:
    support_types = ['Node', 'SupportDX', 'SupportDY', 'SupportDZ', 'SupportRX',
                     'SupportRY', 'SupportRZ']
    csv_writer = csv.DictWriter(f, fieldnames=support_types, restval=False)
    csv_writer.writeheader()

    # Take each x-coordinate
    for i in grid:
        # i = grid[i]
        for j in i:
            # j = grid[i][j]
            for k in j:
                # k = grid[i][j][k] = node tuple(name, X, Y, Z)
                if k[2] == 0:
                    # node is at Y=0, add base supports
                    name = k[0]
                    # Write the supports to the file
                    csv_writer.writerow({'Node': name,
                                         'SupportDX': True,
                                         'SupportDY': True,
                                         'SupportDZ': True})

# Create a file of loads
# Find all nodes at top level
filename = os.path.join(os.path.dirname(__file__), 'gridnodesloads.csv')
with open(filename, mode='w', newline='') as f:
    csv_writer = csv.writer(f)
    # Write the header
    csv_writer.writerow(['Name', 'Direction', 'P'])
    # Find the nodes at y_max and apply load to them
    for node in nodes_at_top_level:
        csv_writer.writerow([node, 'FY', -10.0])