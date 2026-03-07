# Units used for the model in this example are inches and kips

# Import `FEModel3D` from `Pynite`
from Pynite import FEModel3D
from matplotlib import pyplot as plt
import numpy as np

# Import 'Visualization' for rendering the model
from Pynite import Visualization

# Create a new finite element model
simple_beam = FEModel3D()

# Add nodes (14 ft apart)
simple_beam.add_node('N1', 0, 0, 0)
simple_beam.add_node('N2', 14*12, 0, 0)

# Define a material
E = 29000       # Modulus of elasticity (ksi)
G = 11200       # Shear modulus of elasticity (ksi)
nu = 0.3        # Poisson's ratio
rho = 2.836e-4  # Density (kci)
simple_beam.add_material('Steel', E, G, nu, rho)

# Add a section with the following properties:
# Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 20 in^2
simple_beam.add_section('MySection', 20, 100, 150, 250)

#Add member
simple_beam.add_member('M1', 'N1', 'N2', 'Steel', 'MySection')

# Provide simple supports
simple_beam.def_support('N1', True, True, True, True, False, False)  # Constrained for torsion at 'N1'
simple_beam.def_support('N2', True, True, True, False, False, False)  # Not constrained for torsion at 'N2'

# Add a downward point load of 5 kips at the midspan of the beam
simple_beam.add_member_pt_load('M1', 'Fy', 5, 4*12, 'D')  # 5 kips Dead load
simple_beam.add_member_pt_load('M1', 'Fy', 8, 9*12, 'S')  # 8 kips Live load
simple_beam.add_member_pt_load('M1', 'Fy', -8, 12*12, 'L')  # 8 kips Live load

# Add load combinations
simple_beam.add_load_combo('1.4D', {'D': 1.4}, 'strength')
simple_beam.add_load_combo('1.2D+1.6L', {'D': 1.2, 'L': 1.6}, 'strength')
simple_beam.add_load_combo('1.2D+1.6S', {'D': 1.2, 'S': 1.6}, 'strength')

# Analyze the beam and perform a statics check
simple_beam.analyze(check_statics=True)

from Pynite.Visualization import Renderer
rndr = Renderer(simple_beam)
rndr.deformed_shape = True
rndr.deformed_scale = 300
rndr.render_loads = True
rndr.combo_name = '1.2D+1.6L'
rndr.render_model()

# Plot an envelope moment diagram for all the load combinations tagged as 'strength'
simple_beam.members['M1'].plot_moment('Mz', ['strength'], n_points=400)
