# A pinned column with an axial load and a lateral point load at mid-height. There is no node at the
# lateral load, so the second-order moment and deflection at mid-height are found by the member's
# internal P-little-delta solution. The results are compared to the closed-form beam-column solution.

from math import pi, sqrt, tan

from Pynite import FEModel3D

# Create a new finite element model
column = FEModel3D()

# Define the material and section properties (kips and inches)
E = 29000
G = 11200
Iz = 484
Iy = 51.4
A = 14.1
J = 1.45
column.add_material('Steel', E, G, 0.3, 0.490/12**3)
column.add_section('W14x48', A, Iy, Iz, J)

# Define the column. It has a single member with no nodes along its length.
L = 28*12
column.add_node('N1', 0, 0, 0)
column.add_node('N2', 0, L, 0)
column.add_member('M1', 'N1', 'N2', 'Steel', 'W14x48')

# Pin the base and the top about the global Z-axis. The top is free to move vertically. Twisting
# and out-of-plane rotation are restrained so that the column only bends in the XY plane.
column.def_support('N1', True, True, True, True, True, False)
column.def_support('N2', True, False, True, True, False, False)

# Apply an axial load of 60% of the Euler buckling load, and a lateral load at mid-height
Pe = pi**2*E*Iz/L**2
P = 0.6*Pe
H = 5
column.add_node_load('N2', 'FY', -P, 'Axial')
column.add_member_pt_load('M1', 'Fy', H, L/2, 'Lateral')
column.add_load_combo('Combo 1', {'Axial': 1.0, 'Lateral': 1.0})

# Analyze the model
column.analyze_PDelta()

# Calculate the closed-form results at mid-height
k = sqrt(P/(E*Iz))
u = k*L/2
M_exact = H/(2*k)*tan(u)
d_exact = H/(2*P*k)*(tan(u) - u)

# Compare the results
M_Pynite = abs(column.members['M1'].moment('Mz', L/2, 'Combo 1'))
d_Pynite = abs(column.members['M1'].deflection('dy', L/2, 'Combo 1'))
print(f'Moment at mid-height: Pynite {M_Pynite:.1f} k-in, closed-form {M_exact:.1f} k-in')
print(f'Deflection at mid-height: Pynite {d_Pynite:.4f} in, closed-form {d_exact:.4f} in')
