from Pynite.FEModel3D import FEModel3D

# Example 10.5 from Matrix Analysis of Structures, 2nd Edition by William McGuire, Richard H. Gallagher, and Ronald D. Ziemian
# Units in this example are in feet and kips

# Create a new finite element model
frame = FEModel3D()

# Define the nodes of the frame
frame.add_node('a', 0, 0, 0)  # Node a
frame.add_node('b', 0, 24, 0) # Node b
frame.add_node('c', 20, 24, 0) # Node c
frame.add_node('d', 60, 24, 0)  # Node d
frame.add_node('e', 60, 0, 0)  # Node e

# Intermediate nodes
frame.add_node('f', 0, 12, 0) # Node f
frame.add_node('g', 60, 12, 0) # Node g
frame.add_node('h', 10, 24, 0) # Node h
frame.add_node('i', 40, 24, 0) # Node i

# Define steel material properties
# Note that for pushover analysis the yield strength is required.
E = 29000*12**2  # Modulus of elasticity in ksf
G = 11200*12**2  # Shear modulus in ksf
nu = 0.3  # Poisson's ratio
rho = 490/1000  # Density in kips/ft^3
fy = 36*12**2  # Yield strength in ksf
frame.add_material('Steel', E, G, nu, rho, fy)

# Define the cross-sectional properties for a W10x45 steel section
A = 13.3/12**2  # Cross-sectional area in ft^2
Iz = 248/12**4  # Strong axis moment of inertia in ft^4
Iy = 53.4/12**4  # Weak axis moment of inertia in ft^4
J = 1.51/12**4  # Torsional constant in ft^4
Zz = 54.9/12**3  # Strong axis plastic section modulus in ft^3
Zy = 20.3/12**3  # Weak axis plastic section modulus in ft^3

# Note that we are using the (plastic) `add_steel_section` method instead of the generic (elastic)`add_section` method. This is because the pushover analysis requires plastic section properties.
frame.add_steel_section('W10x45', A, Iy, Iz, J, Zy, Zz, 'Steel')

# Define the cross-section properties for a W27x84 steel section
A = 24.8/12**2  # Cross-sectional area in ft^2
Iz = 2850/12**4  # Strong axis moment of inertia in ft^4
Iy = 106/12**4  # Weak axis moment of inertia in ft^4
J = 2.81/12**4  # Torsional constant in ft^4
Zz = 244/12**3  # Strong axis plastic section modulus in ft^3
Zy = 33.2/12**3  # Weak axis plastic section modulus in ft^3

frame.add_steel_section('W27x84', A, Iy, Iz, J, Zy, Zz, 'Steel')

# Define the members of the frame
frame.add_member('ab', 'a', 'b', 'Steel', 'W10x45')  # Column ab
frame.add_member('bc', 'b', 'c', 'Steel', 'W27x84')  # Beam bc
frame.add_member('cd', 'c', 'd', 'Steel', 'W27x84')  # Beam cd
frame.add_member('ed', 'e', 'd', 'Steel', 'W10x45')  # Column ed

# Define the supports
frame.def_support('a', True, True, True, True, True, False)  # Pinned support at node a
frame.def_support('e', True, True, True, True, True, False)  # Pinned support at node e
frame.def_support('b', False, False, True, False, False, False)  # Out-of-plane support at node b
frame.def_support('d', False, False, True, False, False, False)  # Out-of-plane support at node b

# Define the final pushover load pattern
frame.add_node_load('b', 'FX', 0.99*6, 'Push')    # 6 kips of horizontal load at node b
frame.add_node_load('c', 'FY', 0.99*-60, 'Push')  # 60 kips of vertical load at node c
frame.add_node_load('d', 'FY', 0.99*-120, 'Push') # 120 kips of vertical load at node d

# Define load combinations
# At least one primary elastic combo is required for pushover analysis to have a base case to apply pushover loads to. In this case, the Primary combo has no loads in it, so the pushover loads will be applied to an initially unloaded structure.
frame.add_load_combo('Primary', {})

# The pushover combo defines where the pushover loads are found, and what load factor to apply to
# them at each pushover step. The solver will solve any primary load combos first, and then start
# applying the pushover loads to the solutions incrementally until 100% of the pushover loads are
# applied. In this case, we apply 1% of the pushover loads at each step, so it will take 100 steps
# to apply 100% of the pushover loads.
frame.add_load_combo('Pushover', {'Push': 0.01})

# We can define traces to track values of interest at each step of the pushover analysis. We can
# trace just about anything by defining simple functions in a `traces` dictionary. `lambda`
# expressions are a convenient way to define simple trace functions that can access the pushover
# combo results for each load combo and load step. Traces are an optional feature for pushover
# analysis intended to help visualize the progression of key response quantities.
traces = {
    'Node d Drift': lambda combo_name: frame.nodes['d'].DX[combo_name]*12,
}

# Analyze the model
frame.analyze_pushover(log=True, push_combo='Pushover', traces=traces)

from Pynite.Visualization import Renderer
rndr = Renderer(frame)
rndr.combo_name = 'Primary'
rndr.member_diagrams = 'Mz'
rndr.render_model()

# Plot the trace of node d's drift throughout the pushover analysis
frame.plot_pushover_trace('Node d Drift', combo_name='Primary')