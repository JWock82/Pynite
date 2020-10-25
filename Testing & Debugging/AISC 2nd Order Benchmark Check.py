# AISC's benchmark problem used to determine if a second order analysis procedure is rigorous enough
# Units used in this example are inches, and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D
from PyNite import Visualization
import math

# Create the cantilever model
cantilever = FEModel3D()

# Define the column and its properties
L = 20 # ft
H = 5 # kips lateral load
P = 100 # kips axial load
G = 11200*12**2 # shear modulus (ksf)
E = 29000*12**2 # modulus of elasticity (ksf)
I = 100/12**4 # moment of inertia (ft^4)

# Break the column into several segments in order to capture P-little-delta effects
num_segs = 5
num_nodes = num_segs + 1
for i in range(num_nodes):
    # Add nodes
    cantilever.AddNode(str(i+1), 0, i*L/(num_segs), 0)

for i in range(num_segs):
    # Add members between nodes
    cantilever.AddMember(str(i+1), str(i+1), str(i+2), E, G, I, I, 200/12**4, 10/12**2)

# Add a fixed support at the base of the column
cantilever.DefineSupport('1', True, True, True, True, True, True)

# Add a -10 kip axial load to the top of the column
cantilever.AddNodeLoad(str(num_nodes), 'FY', -P)

# Add a 5 kip lateral load to the top of the column
cantilever.AddNodeLoad(str(num_nodes), 'FX', H)

# Perform 2nd order analysis
cantilever.Analyze_PDelta() 

# Render the deformed shape
Visualization.RenderModel(cantilever, text_height=0.3, deformed_shape=True, deformed_scale=2, render_loads=True)

# Print the moment at the base of the column
print('PyNite Calculated Member Moment: ', cantilever.GetMember('1').Moment('Mz', 0.0, 'Combo 1'))
print('PyNite Calculated Reaction Moment: ', cantilever.GetNode('1').RxnMZ['Combo 1'])

# Print the deflection at the top of the column
print('PyNite Calculated Displacement: ', cantilever.GetNode(str(num_nodes)).DX['Combo 1']*12)

# Print the AISC benchmark problem solution:
alpha = (P*L**2/(E*I))**0.5
Mmax = H*L*(math.tan(alpha)/alpha)
ymax = H*L**3/(3*E*I)*(3*(math.tan(alpha)-alpha)/alpha**3)
print('AISC Benchmark Problem Moment: ', Mmax)
print('AISC Benchmark Problem Displacement: ', ymax*12)
