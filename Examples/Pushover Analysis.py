from math import isclose

from Pynite import FEModel3D


"""
Matrix Structural Analysis, 2nd Edition - Examples 8.6 (p. 228) & 10.4 (p. 282)
"""

# Create the model
plastic_beam = FEModel3D()

# Define a material. Yield stress is required for pushover analysis to determine when plastic hinges form.
E = 29000  # ksi
G = 11200  # ksi
nu = 0.3
rho = 0.490/12**3  # kci
fy = 50  # ksi
plastic_beam.add_material('Stl_A992', E, G, nu, rho, fy)

# Define a cross-section
plastic_beam.add_steel_section('W12x65', 19.1, 20, 533, 1, 15, 96.8, 'Stl_A992')

# Add nodes
plastic_beam.add_node('N1', 0, 0, 0)
plastic_beam.add_node('N2', 8*12, 0, 0)
plastic_beam.add_node('N3', 24*12, 0, 0)

# Add supports
plastic_beam.def_support('N1', True, True, True, True, True, True)
plastic_beam.def_support('N3', False, True, True, False, False, False)

# Add a member
plastic_beam.add_member('M1', 'N1', 'N3', 'Stl_A992', 'W12x65')

# Define the pushover loads.
# P = 259.3  # This load causes the first plastic hinge to form at the fixed end (N1)
P = 346.2  # This load causes a second plastic hinge to form at the load point (N2).
plastic_beam.add_node_load('N2', 'FY', -0.3*P, 'Push')
plastic_beam.add_node_load('N3', 'FX', -P, 'Push')

# Add load combinations

# At least one primary elastic combo is required for pushover analysis to have a base case to apply pushover loads to.
# In this case, the Primary combo has no loads in it, so the pushover loads will be applied to an initially unloaded structure.
plastic_beam.add_load_combo('Primary', {})

# The pushover combo defines where the pushover loads are found, and what load factor to apply to
# them at each pushover step. The solver will solve any primary load combos first, and then start
# applying the pushover loads to the solutions incrementally until 100% of the pushover loads are
# applied. In this case, we apply 1% of the pushover loads at each step, so it will take 100 steps
# to apply 100% of the pushover loads.
plastic_beam.add_load_combo('Pushover', {'Push': 0.01})

# We can also define a target displacement for the pushover analysis to stop at. The pushover
# analysis will continue until 100% of the pushover loads are applied, or until the target
# displacement is reached, whichever comes first. If a collapse mechanism forms without a target
# displacement defined, the analysis will continue indefinitely. Let's keep stop the analysis if
# 6" of displacement is ever reached at node 2. 
control_node = 'N2'
control_direction = 'DY'
control_limit = -6.0

# We can define traces to track values of interest at each step of the pushover analysis. We can
# trace just about anything by defining simple functions in a `traces` dictionary. `lambda`
# expressions are a convenient way to define simple trace functions that can access the pushover
# combo results for each load combo and load step. Traces are an optional feature for pushover
# analysis intended to help visualize the progression of key response quantities.
traces = {
    'Fixed End Moment': lambda combo_name: plastic_beam.members['M1'].moment('Mz', x=0.0, combo_name=combo_name),
    'Load Point Moment': lambda combo_name: plastic_beam.members['M1'].moment('Mz', x=96.0, combo_name=combo_name),
    'Load Point Deflection': lambda combo_name: plastic_beam.nodes['N2'].DY[combo_name],
}

# Analyze the model
plastic_beam.analyze_pushover(log=True, check_stability=False, push_combo='Pushover', control_node=control_node, control_direction=control_direction, control_limit=control_limit, traces=traces)

# Plot the traces one by one
# plastic_beam.plot_pushover_trace('Fixed End Moment', combo_name='Primary')
# plastic_beam.plot_pushover_trace('Load Point Moment', combo_name='Primary')
# plastic_beam.plot_pushover_trace('Load Point Deflection', combo_name='Primary')

plastic_beam.members['M1'].plot_moment('Mz', combo_name='Primary', n_points=100)
plastic_beam.members['M1'].plot_deflection('dy', combo_name='Primary', n_points=100)

# Get the final resulting moments from the Primary combo (where pushover results are stored)
M_N1 = plastic_beam.members['M1'].moment('Mz', x=0.0, combo_name='Primary')
M_N2 = plastic_beam.members['M1'].moment('Mz', x=96.0, combo_name='Primary')

print(f'Final Mz at x=0   : {M_N1}')
print(f'Final Mz at x=96  : {M_N2}')
