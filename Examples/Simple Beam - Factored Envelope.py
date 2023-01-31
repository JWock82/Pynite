# Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D
from matplotlib import pyplot as plt
import numpy as np

# Import 'Visualization' for rendering the model
from PyNite import Visualization

# Create a new finite element model
SimpleBeam = FEModel3D()

# Add nodes (14 ft apart)
SimpleBeam.add_node('N1', 0, 0, 0)
SimpleBeam.add_node('N2', 14*12, 0, 0)

# Add a beam with the following properties:
# E = 29000 ksi, G = 11400 ksi, Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 20 in^2
SimpleBeam.add_member('M1', 'N1', 'N2', 29000, 11400, 100, 150, 250, 20)

# Provide simple supports
SimpleBeam.def_support('N1', True, True, True, True, False, False)  # Constrained for torsion at 'N1'
SimpleBeam.def_support('N2', True, True, True, False, False, False) # Not constrained for torsion at 'N2'

# Add a downward point load of 5 kips at the midspan of the beam
SimpleBeam.add_member_pt_load('M1', 'Fy', 5, 4*12, 'D') # 5 kips Dead load
SimpleBeam.add_member_pt_load('M1', 'Fy', 8, 9*12, 'S') # 8 kips Live load
SimpleBeam.add_member_pt_load('M1', 'Fy', -8, 12*12, 'L') # 8 kips Live load

# Add load combinations
SimpleBeam.add_load_combo('1.4D', {'D':1.4})
SimpleBeam.add_load_combo('1.2D+1.6L', {'D':1.2, 'L':1.6})
SimpleBeam.add_load_combo('1.2D+1.6S', {'D':1.2, 'S':1.6})

# Analyze the beam and perform a statics check
SimpleBeam.analyze(check_statics=True)
# 
# Visualization.render_model(SimpleBeam, annotation_size=10, deformed_shape=True, deformed_scale=30, render_loads=True, combo_name='1.2D+1.6L')

# Plot the shear diagram with all load cases and max/min envelope
x, M1 = SimpleBeam.Members['M1'].moment_array("Mz", n_points=400, combo_name='1.4D')
_, M2 = SimpleBeam.Members['M1'].moment_array("Mz", n_points=400, combo_name='1.2D+1.6L')
_, M3 = SimpleBeam.Members['M1'].moment_array("Mz", n_points=400, combo_name='1.2D+1.6S')

max_envelope = np.maximum(np.maximum(M1, M2), M3)
min_envelope = np.minimum(np.minimum(M1, M2), M3)

plt.plot(x, np.zeros(len(x)), c="black", lw=3)
plt.plot(x, M1)
plt.plot(x, M2)
plt.plot(x, M3)
plt.plot(x, max_envelope, alpha=0.3, c="green", lw=5)
plt.plot(x, min_envelope, alpha=0.3, c="red", lw=5)