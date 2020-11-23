import os, csv

def read_csv(filename):
    # CSV import of members
    with open(filename, newline='') as f:
        csv_reader = csv.reader(f)
        # skip the first row in the file (the header)
        next(csv_reader)
        # Get the other rows
        csv_rows = [row for row in csv_reader]
    return csv_rows

def nodes_from_csv(filename):
    # Read in nodes from the file
    csv_node_list = read_csv(filename)
    for row in csv_node_list:
        # Convert the numerical fields to floats
        # default for csv is strings
        row[1:] = [float(elem) for elem in row[1:]]
    return csv_node_list