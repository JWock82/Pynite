from numpy import array, cross
from numpy.linalg import norm, dot
from math import isclose, acos
import Member3D

class PhysMember(Member3D):

    def __init__(self, name, i_node, j_node, E, G, Iy, Iz, J, A, model, aux_node=None,
                 tension_only=False, comp_only=False):
        
        super.__init__(name, i_node, j_node, E, G, Iy, Iz, J, A, model, aux_node, tension_only, comp_only)
        self.sub_members = {}
    
    def __subdivide(self):
        """
        Subdivides the physical member into sub-members at each node along the physical member
        """

        # Start a new list of nodes along the member
        int_nodes = []

        # Create a vector from the i-node to the j-node
        Xi, Yi, Zi = self.i_node.X, self.i_node.Y, self.i_node.Z
        Xj, Yj, Zj = self.j_node.X, self.j_node.Y, self.j_node.Z
        vector_ij = array([Xj-Xi, Yj-Yi, Zj-Zi])

        # Add the i-node and j-node to the list
        int_nodes.append([self.i_node, 0])
        int_nodes.append([self.j_node, norm(vector_ij)])

        # Step through each node in the model
        for node in self.model.Nodes.values():

            # Create a vector from the i-node to the current node
            X, Y, Z = node.X, node.Y, node.Z
            vector_in = array([X-Xi, Y-Yi, Z-Zi])

            # Calculate the angle between the two vectors
            angle = acos(dot(vector_in, vector_ij)/(norm(vector_in)*norm(vector_ij)))

            # Determine if the node is colinear with the member
            if isclose(angle, 0):

                # Determine if the node is on the member
                if norm(vector_in) < norm(vector_ij):

                    # Add the node to the list of intermediate nodes
                    int_nodes.append([node, norm(vector_in)])
        
        # Create a list of sorted intermediate nodes by distance from the i-node
        int_nodes = sorted(int_nodes, key=lambda x: x[1])

        # Break up the member into sub-members at each intermediate node
        for i, node in enumerate(int_nodes[1:]):

            # Generate the sub-member's name (physical member name + a, b, c, etc.)
            name = self.name + chr(i+96)

            i_node = int_nodes[i-1][0]
            j_node = node
            xi = int_nodes[i-1][1]
            xj = int_nodes[i][1]

            # Create a new sub-member
            self.model.add_member(name, int_nodes[i-1], node, self.E,
                                  self.G, self.Iy, self.Iz, self.J, self.A, self.aux_node,
                                  self.tension_only, self.comp_only)
            
            # Apply end releases if applicable
            if i == 1:
                self.model.Members[name].Releases[0:6] = self.Releases[0:6]
            if i == len(int_nodes) - 1:
                self.model.Members[name].Releases[6:12] = self.Releases[6:12]

            # Add distributed to the sub-member
            for dist_load in self.model.Members[name].DistLoads:
                
                # Find the start and end points of the distributed load in the physical member's
                # local coordinate system
                x1_load = dist_load[3]
                x2_load = dist_load[4]

                # Find the start and end points of the sub-member in the physical member's
                # local coordinate system
                x1_mem = ((int_nodes[i-1].X - Xi)**2 + (int_nodes[i-1].Y - Yi)**2 + (int_nodes[i-1].Z - Zi)**2)**0.5
                x2_mem = ((node.X - Xi)**2 + (node.Y - Yi)**2 + (node.Z - Zi)**2)**0.5

                # Determine if the distributed load should be applied to this segment
                if x1_load <= x2_mem and x2_load > x1_mem: 
                    
                    direction = dist_load[0]
                    w1 = dist_load[1]
                    w2 = dist_load[2]
                    case = dist_load[5]

                    # Equation describing the load as a function of x
                    w = lambda x: (w2 - w1)/(x2_load - x1_load)*(x - x1_load) + w1

                    # Chop up the distributed load for the sub-member
                    if x1_load > x1_mem:
                        x1 = x1_load - x1_mem
                    else:
                        x1 = 0
                        w1 = w(x1_mem)
                    
                    if x2_load < x2_mem:
                        x2 = x2_load - x1_mem
                    else:
                        x2 = x2_mem
                        w2 = w(x2_mem)

                    # Add the load to the sub-member
                    self.model.Members[name].DistLoads.append([direction, w1, w2, x1, x2, case])

            # Add point loads to the sub-member
            for pt_load in self.model.Members[name].PtLoads:
                
                direction = pt_load[0]
                P = pt_load[1]
                x = dist_load[2]
                case = pt_load[3]

                # Determine if the point load should be applied to this segment
                if x >= x1_mem and x < x2_mem or isclose(x, self.Length()):

                    x = x - x1_mem
                    
                    # Add the load to the sub-member
                    self.model.Members[name].PtLoads.append([direction, P, x, case])

