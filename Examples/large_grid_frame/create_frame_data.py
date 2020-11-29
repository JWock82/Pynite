import numpy as np
import itertools, csv, os

x_values = np.linspace(0, 10, 5)
y_values = np.linspace(0, 50, 11)
z_values = np.linspace(0, 10, 5)

# 3D grid
grid = []
id_number = itertools.count()
for i, x in enumerate(x_values):
    grid.append([])
    for j, y in enumerate(y_values):
        grid[i].append([])
        for k, z in enumerate(z_values):
            grid[i][j].append(('N'+str(next(id_number)), x, y, z))

filename = os.path.join(os.path.dirname(__file__), 'gridnodes.csv')
with open(filename, mode='w', newline='') as f:
    csv_writer = csv.writer(f)
    # Write the header
    csv_writer.writerow(['Name', 'X', 'Y', 'Z'])
    # Write the member list to the file
    for nodes_xy in itertools.chain(grid):
        for node_z in nodes_xy:
            csv_writer.writerows(node_z)

# connect each point to the adjacent one
connections = {}
# all rows except last one
for row_idx, row in enumerate(grid):
    # all columns minus the last one
    for col_idx, col in enumerate(row):
        # all z-values except the last one
        for z_idx, value in enumerate(col):
            connections[value[0]] = []
            try:
                # Get value from next column
                connections[value[0]].append(grid[row_idx][col_idx + 1][z_idx][0])
            except IndexError:
                pass
            try:
                # Get value from next row
                connections[value[0]].append(grid[row_idx + 1][col_idx][z_idx][0])
            except IndexError:
                pass
            try:
                # Get value from next z-value
                connections[value[0]].append(grid[row_idx][col_idx][z_idx + 1][0])
            except IndexError:
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
