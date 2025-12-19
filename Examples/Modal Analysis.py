from Pynite import FEModel3D

model = FEModel3D()

# Create some nodes
model.add_node('N1', 0, 0, 0)
model.add_node('N2', 0, 10, 0)
model.add_node('N3', 0, 20, 0)
model.add_node('N4', 0, 30, 0)

# Add supports
model.def_support('N1', True, True, True, True, True, True)
model.def_support('N2', False, False, True, True, True, False)
model.def_support('N3', False, False, True, True, True, False)
model.def_support('N4', False, False, True, True, True, False)

# Define material properties
E = 29000*12**2  # ksf
G = 11200*12**2  # ksd
nu = 0.3
rho = 490/1000   # kcf

model.add_material('Steel', E, G, nu, rho)

# Define section properties (AISC W12x26)
A = 7.625/12**2  # ft^2
Iy = 17.3/12**4  # ft^4
Iz = 204/12**4   # ft^4
J = 0.300/12**4  # ft^4

model.add_section('W12x26', A, Iy, Iz, J)

# Add members to the model
model.add_member('M1', 'N1', 'N2', 'Steel', 'W12x26')
model.add_member('M2', 'N2', 'N3', 'Steel', 'W12x26')
model.add_member('M3', 'N3', 'N4', 'Steel', 'W12x26')

# Add some nodal dead loads (2 kips each)
model.add_node_load('N2', 'FY', -2, 'D')
model.add_node_load('N3', 'FY', -2, 'D')
model.add_node_load('N4', 'FY', -2, 'D')

# Define the dead load due to member self-weight
model.add_member_self_weight('FY', 1.0, 'D-SW')

# Create a load combination containing all loads to be converted to masses for modal analysis
model.add_load_combo('Mass', {'D': 1.0, 'D-SW': 1.0})

# Gather the parameters needed for modal analysis
num_modes = 3             # Since there are 3 stories in this structure, there are 3 dominant modes of vibration
mass_combo_name = 'Mass'  # The name of the mass combo we just created
gravity = 32.2            # ft/s^2
mass_direction = 'Y'      # Direction for load-to-mass conversion

# Run the modal analysis
model.analyze_modal(num_modes, mass_combo_name, mass_direction, gravity)

# Print the natural frequencies
for i, freq in enumerate(model.frequencies):
    print(f'Mode {i + 1} Frequency: {freq:3f} Hz')

# Pynite has stored each mode shape in its own load combination. We can render any mode shape as follows
from Pynite.Visualization import Renderer  # Using the VTK renderer for its speed, but we could also use the PyVista renderer the same way
rndr = Renderer(model)      # Create a new renderer based on the model
rndr.annotation_size = 1    # Text will display 1' high, and objects (e.g. supports) will scale accordingly
rndr.deformed_scale = 1     # Adjust as needed for your model
rndr.deformed_shape = True  # Make sure to turn on the deformed shape for mode shapes
rndr.combo_name = 'Mode 3'  # Pick any mode (e.g. 'Mode 1', 'Mode 2', 'Mode 3')
rndr.render_model()
