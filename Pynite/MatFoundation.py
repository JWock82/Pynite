import numpy as np

from Pynite.Mesh import RectangleMesh

class MatFoundation(RectangleMesh):

    def __init__(self, name, mesh_size, length_X, length_Z, thickness, material_name, model, ks, origin=[0, 0, 0]):

        super().__init__(mesh_size, length_X, length_Z, thickness, material_name, model, 1, 1, origin, 'XZ')

        self.name = name
        self.ks = ks
        self.pt_loads = []  # [XZ_coord, direction, magnitude, case]

    def add_rect_cutout(self, name, X_min, Z_min, X_max, Z_max):

        self.add_rect_opening(name, X_min, Z_min, X_max - X_min, Z_max - Z_min)

    def add_mat_pt_load(self, XZ_coord, direction, magnitude, case='Case 1'):

        self.x_control.append(XZ_coord[0])
        self.y_control.append(XZ_coord[1])
        self.pt_loads.append([XZ_coord, direction, magnitude, case])

    def generate(self):

        # Generate the mesh
        super().generate()

        # Add point loads to the model
        for node in self.nodes.values():
            for pt_load in self.pt_loads:
                if np.isclose(node.X, pt_load[0][0]) and np.isclose(node.Z, pt_load[0][1]):
                    self.model.add_node_load(node.name, pt_load[1], pt_load[2], pt_load[3])

        # Step through each node in the mat
        for node in self.nodes.values():

            # Initialize the tributary area to the node to zero
            trib = 0

            # Step through each plate in the model
            for plate in self.elements.values():

                # Determine if the plate is attached to the node
                if node.name in [plate.i_node.name, plate.j_node.name, plate.m_node.name, plate.n_node.name]:

                    # Add 1/4 the plate's area to the tributary area to the node
                    trib += abs(plate.j_node.X - plate.i_node.X)*abs(plate.m_node.Z - plate.j_node.Z)/4

            # Add a soil spring to the node
            self.model.def_support_spring(node.name, 'DY', self.ks*trib, '-')
