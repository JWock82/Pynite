# Example of a simply supported beam with a point load.
# Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

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
SimpleBeam.add_member_pt_load('M1', 'Fy', -5, 7*12, 'D') # 5 kips Dead load
SimpleBeam.add_member_pt_load('M1', 'Fy', -8, 7*12, 'L') # 8 kips Live load

# Add load combinations
SimpleBeam.add_load_combo('1.4D', {'D':1.4})
SimpleBeam.add_load_combo('1.2D+1.6L', {'D':1.2, 'L':1.6})

# Analyze the beam and perform a statics check
SimpleBeam.analyze(check_statics=True)

Visualization.RenderModel(SimpleBeam, text_height=10, deformed_shape=True, deformed_scale=30, render_loads=True, combo_name='1.2D+1.6L')

# Print the shear, moment, and deflection diagrams
SimpleBeam.GetMember('M1').plot_shear('Fy', '1.2D+1.6L')
SimpleBeam.GetMember('M1').plot_moment('Mz', '1.2D+1.6L')
SimpleBeam.GetMember('M1').plot_deflection('dy', '1.2D+1.6L')

# Print reactions at each end of the beam
print('Left Support Reaction:', SimpleBeam.GetNode('N1').RxnFY['1.2D+1.6L'], 'kip')
print('Right Support Reacton:', SimpleBeam.GetNode('N2').RxnFY['1.2D+1.6L'], 'kip')

# Print the max/min shears and moments in the beam
print('Maximum Shear:', SimpleBeam.GetMember('M1').max_shear('Fy', '1.2D+1.6L'), 'kip')
print('Minimum Shear:', SimpleBeam.GetMember('M1').min_shear('Fy', '1.2D+1.6L'), 'kip')
print('Maximum Moment:', SimpleBeam.GetMember('M1').max_moment('Mz', '1.2D+1.6L')/12, 'kip-ft')
print('Minimum Moment:', SimpleBeam.GetMember('M1').min_moment('Mz', '1.2D+1.6L')/12, 'kip-ft')

# Print the max/min deflections in the beam
print('Maximum Deflection:', SimpleBeam.GetMember('M1').max_deflection('dy', '1.2D+1.6L'), 'in')
print('Minimum Deflection:', SimpleBeam.GetMember('M1').min_deflection('dy', '1.2D+1.6L'), 'in')

# The following lines can be uncommented to create a PDF report. Follow the instructions on the
# wiki under "Generating PDF Reports" to prevent errors. The report will be output to the PyNite
# folder unless the 'output_path' variable below is modified.

# from PyNite import Reporting
# Reporting.CreateReport(SimpleBeam, output_filepath='.//PyNite Report.pdf', plates=False, plate_corner_forces=False, \
#                        plate_center_forces=False, plate_corner_membrane=False, plate_center_membrane=False)