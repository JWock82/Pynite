from PyNite import FEModel3D
from PyNite.Visualization import render_model

beam_model = FEModel3D()
nodes = [
    (0, 0, 0),
    (1000, 0, 0),
    (3600, 0, 0),
    (4800, 0, 0),
]

for idx, node in enumerate(nodes):
    beam_model.add_node(f"N{idx}", *node)
    
beam_model.def_support("N1", 1, 1, 1, 1, 0, 0) # @ 1000
beam_model.def_support("N2", 0, 1, 1, 0, 0, 0) # @ 3600

beam_model.add_material('Mat1', 24500, 24500/3, 0.3, 0)

beam_model.add_member("M1", "N0", "N3", 'Mat1', Iy=1000, Iz=1200000000, J=120000, A=1000)

beam_model.add_member_dist_load("M1", "Fy", w1=-40, w2=-40, x1=0, x2=4800)
beam_model.add_member_pt_load("M1", "Fy", P=1000, x=4800)

beam_model.analyze(check_statics=True)

beam_model.Nodes['N1'].RxnFY
beam_model.Nodes["N2"].RxnFY

beam_model.Members["M1"].plot_shear(Direction="Fy", n_points=1000)
beam_model.Members["M1"].plot_moment(Direction="Mz", n_points=1000)
beam_model.Members["M1"].plot_deflection(Direction="dy", n_points=1000)