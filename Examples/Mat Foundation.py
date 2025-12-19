from Pynite import FEModel3D

# Create a new finite element model
model = FEModel3D()

# Create a new material (concrete)
E=0.35*57*(4000)**0.5*144
G=0.4*E
nu = 0.17
rho = 0.150
model.add_material('Concrete', E, G, nu, rho)

# Add a mat foundation to the model
ks = 100/1000*12**3  # Subgrade modulus of reaction (100 pci coverted to kcf)
model.add_mat_foundation('MAT1', mesh_size=1.0, length_X=50, length_Z=25, thickness=1.0, material_name='Concrete', ks=ks, origin=[0, 0, 0])

# Add point loads to the mat foundation
model.mats['MAT1'].add_mat_pt_load((10, 10), 'FY', -10, 'D')
model.mats['MAT1'].add_mat_pt_load((10, 10), 'MZ', 15, 'D')

# Add a rectangular cutout to the mat foundation
model.mats['MAT1'].add_rect_opening('OPNG1', 5, 5, 10, 15)

# Add a load combination to the model
model.add_load_combo('D', {'D': 1.0})

# Analyze the model
model.analyze(log=True, check_statics=True, num_steps=1)

from Pynite.Rendering import Renderer
rndr = Renderer(model)
rndr.annotation_size = 0.5
rndr.combo_name = 'D'
rndr.deformed_shape = True
rndr.deformed_scale = 100
rndr.color_map = 'Mx'
rndr.render_nodes = True
rndr.render_model()
