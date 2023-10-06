from PyNite import FEModel3D
from PyNite.Section import SteelSection

# Create the model
plastic_beam = FEModel3D()

# Define a material
E = 29000  # ksi
G = 11200  # ksi
nu = 0.3
rho = 0.490/12**3  # kci
fy = 50  # ksi
plastic_beam.add_material('Stl_A992', E, G, nu, rho, fy)

# Define a cross-section
W12 = SteelSection(plastic_beam, 'W12x65', 19.1, 20, 533, 1, 15, 96.8, 'Stl_A992')
plastic_beam.add_section('W12x65', W12)

# Add nodes
plastic_beam.add_node('N1', 0, 0, 0)
plastic_beam.add_node('N2', 8*12, 0, 0)
plastic_beam.add_node('N3', 24*12, 0, 0)

# Add supports
plastic_beam.def_support('N1', True, True, True, True, True, True)
plastic_beam.def_support('N2', False, True, True, False, False, False)

# Add a member
plastic_beam.add_member('M1', 'N1', 'N2', 'A992', section='W12x65')

# Add a load
plastic_beam.add_node_load('N2', 'FY', 0.3, 'Push')
plastic_beam.add_node_load('N3', 'FX', -1, 'Push')

# Add a load combination
plastic_beam.add_load_combo('Pushover', {'Push':1})

# Analysis the model
plastic_beam.analyze_pushover(True, True, True, 'Pushover', 1.02, 20, 30, True)