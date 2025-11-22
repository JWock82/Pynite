"""
This is an example from SAP2000, used by OpenSees:

https://openseespydoc.readthedocs.io/en/latest/src/PortalFrame2d.html

Two dimensional Frame: Eigenvalue & Static Loads


REFERENCES:

used in verification by SAP2000:
* SAP2000 Integrated Finite Element Analysis and Design of Structures, Verification Manual, Computers and Structures, 1997. Example 1.
* Seismo-struct (Example 10), SeismoStruct, Verification Report For Version 6, 2012. Example 11.

"""

from Pynite import FEModel3D
import math

# Create a new 3D model
model = FEModel3D()

sparse = True

# Set some properties
# Units: kip, in (consistent with the OpenSees model)
numBay = 2
numFloor = 7

bayWidth = 360.0  # in
storyHeights = [162.0, 162.0, 156.0, 156.0, 156.0, 156.0, 156.0]  # in

# Material properties
E = 29500.0  # ksi
G = 11346.15  # ksi (calculated from E and nu=0.3)
nu = 0.3
rho = 0.000000000001 # 0.490/12**3 # Density of steel

weight_Y = 190  # kip

# Section properties from AISC database with actual I_minor and J values
WSection = {
    'W14X176': [51.7, 2150., 838., 26.5],   # [A, I_major, I_minor, J]
    'W14X211': [62.1, 2670., 1030., 38.8],
    'W14X246': [72.3, 3230., 1200., 54.6],
    'W14X287': [84.4, 3910., 1390., 75.8],
    'W24X110': [32.5, 3330., 370., 4.19],
    'W24X130': [38.3, 4020., 473., 6.28],
    'W24X160': [47.1, 5120., 632., 11.8]
}

beams = ['W24X160', 'W24X160', 'W24X130', 'W24X130', 'W24X110', 'W24X110', 'W24X110']
eColumn = ['W14X246', 'W14X246', 'W14X246', 'W14X211', 'W14X211', 'W14X176', 'W14X176']
iColumn = ['W14X287', 'W14X287', 'W14X287', 'W14X246', 'W14X246', 'W14X211', 'W14X211']
columns = [eColumn, iColumn, eColumn]

# Add material
model.add_material('Steel', E, G, nu, rho)

# CREATING SECTION DATABASE FIRST - BEFORE ADDING MEMBERS
print("Creating section database...")
added_sections = set()
for section_name in set(WSection.keys()):
    if section_name not in added_sections:
        A = WSection[section_name][0]  # Area (in²)
        Iy = WSection[section_name][2]  # Weak axis moment of inertia (estimated, not critical for 2D)
        Iz = WSection[section_name][1]  # Strong axis moment of inertia (in⁴)
        J = WSection[section_name][3]  # Torsional constant (estimated)

        model.add_section(section_name, A, Iy, Iz, J)
        added_sections.add(section_name)
        print(f"Added section: {section_name}")

# Add nodes floor by floor
node_dict = {}  # To store node tags for element creation
node_tag = 1

y_loc = 0.0
for j in range(numFloor + 1):

    x_loc = 0.0

    for i in range(numBay + 1):

        # Add node at Z = 0 (2D frame in XY plane)
        model.add_node(f'N{node_tag}', x_loc, y_loc, 0.0)

        # Store node tag for element creation
        node_dict[(i, j)] = f'N{node_tag}'

        x_loc += bayWidth
        node_tag += 1

    # Move to next floor level
    if j < numFloor:

        y_loc += storyHeights[j]

# Support all nodes in the Z-direction (out-of-plane translation)
for node in model.nodes.values():
    model.def_support(node.name, False, False, True, False, False, False)

# Adjust first floor nodes to be fixed
model.def_support('N1', True, True, True, True, True, True)  # Fixed support
model.def_support('N2', True, True, True, True, True, True)  # Fixed support  
model.def_support('N3', True, True, True, True, True, True)  # Fixed support

# Add masses to master nodes (equivalent to equalDOF constraint in OpenSees)
# In OpenSees, nodes 4, 7, 10, 13, 16, 19, 22 are the master nodes for each floor
master_nodes = ['N4', 'N7', 'N10', 'N13', 'N16', 'N19', 'N22']
for node_name in master_nodes:
    model.add_node_load(node_name, 'FY', weight_Y, case='Mass')
    # Add small rotational mass for numerical stability
    model.add_node_load(node_name, 'MY', 1.0e-10, case='Mass')

# Add column elements
ele_tag = 1
for j in range(numBay + 1):  # For each column line
    thisColumn = columns[j]
    
    for i in range(numFloor):  # For each story

        node_tag1 = node_dict[(j, i)]      # Bottom node
        node_tag2 = node_dict[(j, i + 1)]  # Top node
        secType = thisColumn[i]
        
        # Add column member - USE PRE-DEFINED SECTION NAME
        model.add_member(f'M{ele_tag}', node_tag1, node_tag2, 'Steel', secType, lumped_mass=False)
        
        ele_tag += 1

# Add beam elements
for j in range(1, numFloor + 1):  # Start from first floor up
    secType = beams[j - 1]
    
    for i in range(numBay):  # For each bay

        node_tag1 = node_dict[(i, j)]      # Left node
        node_tag2 = node_dict[(i + 1, j)]  # Right node
        
        # Add beam member - USE PRE-DEFINED SECTION NAME
        model.add_member(f'M{ele_tag}', node_tag1, node_tag2, 'Steel', secType, lumped_mass=False)
        
        ele_tag += 1

# Create a load combination that includes only mass for modal analysis
model.add_load_combo('MassCombo', {'Mass': 1.0})

# Run modal analysis
print("Running modal analysis...")
model.analyze_modal(num_modes=7, mass_combo_name="MassCombo", mass_direction=1, gravity=386, sparse=sparse, log=False)  # X-direction

# Access results
frequencies = model.modal_results['frequencies']
mode_shapes = model.modal_results['mode_shapes']

print("\nNatural Frequencies (Hz) and Periods (s):")
print('{:>10}{:>15}{:>15}'.format('Mode', 'Frequency', 'Period'))
for i, freq in enumerate(frequencies):
    period = 1.0 / freq if freq > 0 else 0
    print(f"Mode {i+1}: {freq:.4f} Hz, {period:.4f} s")


# ==== Add static load case ('Static') =========================
# Apply static loads for comparison with OpenSees results
print("\nApplying static loads...")

# Apply the same loads as in OpenSees model
# Loads are applied to master nodes (equivalent to equalDOF constrained nodes)
static_loads = {
    'N22': 20.0,  # Top floor
    'N19': 15.0,  # 6th floor  
    'N16': 12.5,  # 5th floor
    'N13': 10.0,  # 4th floor
    'N10': 7.5,   # 3rd floor
    'N7': 5.0,    # 2nd floor
    'N4': 2.5     # 1st floor
}

for node_name, load_val in static_loads.items():
    model.add_node_load(node_name, 'FX', load_val, case='Static')

# Create static load combination
model.add_load_combo('StaticCombo', {'Static': 1.0})

# Analyze static loads
model.analyze()

# Get static analysis results for comparison
print("\nStatic Analysis Results:")

# 1. Displacement at top (node 22)
disp_top = model.nodes['N22'].DX['StaticCombo']
print(f"Displacement at top (N22): {disp_top:.5f} in")

# 2. Axial force in bottom left column (member 1)
# Get a dictionary of all internal forces at the member ends
mem_1 = model.members['M1']
axial_forces = mem_1.axial_array(5, 'StaticCombo')
axial_force = axial_forces[1][0] # 
print(axial_forces)
print(f"Axial force in bottom left column (M1): {axial_force:.2f} kip")

print(mem_1.moment_array('Mz',5,'StaticCombo'))
moments = mem_1.moment_array('Mz',5, 'StaticCombo')
moment = moments[1][0]
print(f"Moment in bottom left column (M1): {moment:.2f} kip-ft")

# ============= Comparisons with OpenSees Static Results ===============
# Expected results from OpenSees verification
expected_results = {
    'disp_top': 1.45076,
    'axial_force': 69.99, 
    'moment': 2324.68
}

print("\nComparison with OpenSees Static Results:")

axial_force = -axial_force  # needed for comparison
print("\nComparison with Expected Results:")
print(f"{'Parameter':<30}{'PyNite':>10}{'Expected':>10}{'Difference':>12}{'Percentage':>12}")
print(f"{'Displacement at top (ft)':<30}{disp_top:>10.5f}{expected_results['disp_top']:>10.2f}" + \
      f"{disp_top-expected_results['disp_top']:>12.5f}{100%(disp_top/expected_results['disp_top']-1):>11.2f}%")
print(f"{'Axial force (kip)':<30}{axial_force:>10.2f}{expected_results['axial_force']:>10.2f}" + \
      f"{axial_force-expected_results['axial_force']:>12.2f}{100*(axial_force/expected_results['axial_force']-1):>11.2f}%")
print(f"{'Moment (kip-ft)':<30}{moment:>10.2f}{expected_results['moment']:>10.2f}" + \
      f"{moment-expected_results['moment']:>12.2f}{100*(moment/expected_results['moment']-1):>11.2f}%")


# ============= Comparisons with OpenSees Modal Results ===============
print("\nComparison with OpenSees Modal Results:")

# Print period comparisons with SAP2000 and SeismoStruct
print("\n\nPeriod Comparisons:")
print('{:>10}{:>10}{:>15}{:>15}{:>15}{:>15}'.
      format('Mode', 'Diff', 'PyNite', 'OpenSees', 'SAP2000', 'SeismoStruct'))

comparisonResults = [[1.27321, 0.43128, 0.24204, 0.16018, 0.11899, 0.09506, 0.07951],
                     [1.2732, 0.4313, 0.2420, 0.1602, 0.1190, 0.0951, 0.0795],
                     [1.2732, 0.4313, 0.2420, 0.1602, 0.1190, 0.0951, 0.0795]]

periods = [((1.0 / freq) if freq > 0 else 0) for freq in frequencies]
diff = [f'{100*(p1/p2-1):6.2f}%' for p1, p2 in zip(periods, comparisonResults[0])]

for i in range(7):
    if i < len(frequencies):
        print('{:>10}{:>10}{:>15.4f}{:>15.4f}{:>15.4f}{:>15.4f}'.
              format(i + 1, diff[i], periods[i], comparisonResults[0][i], comparisonResults[1][i], comparisonResults[2][i]))