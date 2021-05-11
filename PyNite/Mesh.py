from PyNite.FEModel3D import FEModel3D
from PyNite.Node3D import Node3D
from PyNite.Quad3D import Quad3D
from PyNite.Visualization import RenderModel
from math import pi, sin, cos, isclose
from numpy import cross

class Mesh():

    def __init__(self, t, E, nu, start_node='N1', start_element='Q1'):

        self.t = t                          # Thickness
        self.E = E                          # Modulus of elasticity
        self.nu = nu                        # Poisson's ratio
        self.start_node = start_node        # The name of the first node in the mesh
        self.start_element = start_element  # The name of the first element in the mesh
        self.nodes = {}                     # A dictionary containing the nodes in the mesh
        self.elements = {}                  # A dictionary containing the elements in the mesh
            
class AnnulusMesh(Mesh):

    def __init__(self, t, E, nu, mesh_size, outer_radius, inner_radius, center=[0, 0, 0], start_node='N1', start_element='Q1'):

        super().__init__(t, E, nu, start_node, start_element)

        self.r1 = inner_radius
        self.r2 = outer_radius
        self.mesh_size = mesh_size
        self.center = center

        self.__mesh()
    
    def __mesh(self):
        
        t = self.t
        E = self.E
        nu = self.nu
        mesh_size = self.mesh_size
        r_outer = self.r2
        r_inner = self.r1
        n = int(self.start_node[1:])
        q = int(self.start_element[1:])

        circumf = 2*pi*r_inner           # Circumference of the ring at the inner radius
        n_circ = int(circumf/mesh_size)  # Number of times `mesh_size` fits in the circumferential direction

        # Mesh the annulus from the inside toward the outside
        while round(r_inner, 10) < round(r_outer, 10):
            
            radial = r_outer - r_inner                    # Remaining length in the radial direction to be meshed
            circumf = 2*pi*r_inner                        # Circumference of the ring at the inner radius
            b_circ = circumf/n_circ                       # Element width in the circumferential direction
            n_rad = int(radial/min(mesh_size, 3*b_circ))  # Number of times the plate width fits in the remaining unmeshed radial direction
            h_rad = radial/n_rad                          # Element height in the radial direction

            # Determine if the mesh is getting too big. If so the mesh will need to transition to a
            # finer mesh.
            if b_circ > 3*mesh_size:
                transition = True
            else:
                transition = False
        
            # Create a mesh of nodes for the ring
            if transition == True:
                ring = AnnulusTransRingMesh(t, E, nu, r_inner + h_rad, r_inner, n_circ, self.center, 'N' + str(n), 'Q' + str(q))
                n += 3*n_circ
                q += 4*n_circ
                n_circ *= 3
            else:
                ring = AnnulusRingMesh(t, E, nu, r_inner + h_rad, r_inner, n_circ, self.center, 'N' + str(n), 'Q' + str(q))
                n += n_circ
                q += n_circ
        
            # Add the newly generated nodes and elements to the overall mesh. Note that if duplicate
            # keys exist, the `.update()` method will overwrite them with the newly generated key value
            # pairs. This works in our favor by automatically eliminating duplicate nodes at the shared
            # boundaries between rings.
            self.nodes.update(ring.nodes)
            self.elements.update(ring.elements)

            # Prepare to move to the next ring
            r_inner += h_rad

class AnnulusRingMesh(Mesh):

    def __init__(self, t, E, nu, outer_radius, inner_radius, num_quads, center=[0, 0, 0], start_node='N1', start_element='Q1'):
        '''
        Parameters
        ----------
        direction : array
            A vector indicating the direction normal to the ring.
        '''

        super().__init__(t, E, nu, start_node=start_node, start_element=start_element)

        self.r1 = inner_radius
        self.r2 = outer_radius
        self.n = num_quads
        self.Xo = center[0]
        self.Yo = center[1]
        self.Zo = center[2]

        # Generate the nodes and elements
        self.__mesh()

    def __mesh(self):

        n = self.n  # Number of plates in the ring

        r1 = self.r1  # The inner radius of the ring
        r2 = self.r2  # The outer radius of the ring

        Xo = self.Xo  # Global X-coordinate of the center of the ring
        Yo = self.Yo  # Global Y-coordinate of the center of the ring
        Zo = self.Zo  # Global Z-coordinate of the center of the ring

        theta = 2*pi/self.n  # Angle between nodes in the ring

        # Each node number will be increased by the offset calculated below
        node_offset = int(self.start_node[1:]) - 1

        # Each element number will be increased by the offset calculated below
        element_offset = int(self.start_element[1:]) - 1

        # Generate the nodes that make up the ring, working from the inside to the outside
        angle = 0
        for i in range(1, 2*n + 1, 1):

            # Assign the node a name
            node_name = 'N' + str(i + node_offset)

            # Generate the inner radius of nodes
            if i <= n:
                angle = theta*(i - 1)
                x = Xo + r1*cos(angle)
                y = Yo + r1*sin(angle)
                z = Zo
            # Generate the outer radius of nodes
            else:
                angle = theta*((i - n) - 1)
                x = Xo + r2*cos(angle)
                y = Yo + r2*sin(angle)
                z = Zo
            
            self.nodes[node_name] = [node_name, x, y, z]

        # Generate the elements that make up the ring
        for i in range(1, n + 1, 1):

            # Assign the element a name
            element_name = 'Q' + str(i + element_offset)
            
            n_node = i
            i_node = i + n
            if i != n:
                m_node = i + 1
                j_node = i + 1 + n
            else:
                m_node = 1
                j_node = 1 + n

            self.elements[element_name] = [element_name, 'N' + str(i_node + node_offset),
                                                         'N' + str(j_node + node_offset),
                                                         'N' + str(m_node + node_offset),
                                                         'N' + str(n_node + node_offset),
                                                         self.t, self.E, self.nu]

class AnnulusTransRingMesh(Mesh):

    def __init__(self, t, E, nu, outer_radius, inner_radius, num_inner_quads, center=[0, 0, 0], start_node='N1', start_element='Q1'):
        '''
        Parameters
        ----------
        direction : array
            A vector indicating the direction normal to the ring.
        '''

        super().__init__(t, E, nu, start_node=start_node, start_element=start_element)

        self.r1 = inner_radius
        self.r2 = (inner_radius + outer_radius)/2
        self.r3 = outer_radius
        self.n = num_inner_quads
        self.Xo = center[0]
        self.Yo = center[1]
        self.Zo = center[2]

        # Create the mesh
        self.__mesh()

    def __mesh(self):

        n = self.n  # Number of plates in the outside of the ring (coarse mesh)

        r1 = self.r1  # The inner radius of the ring
        r2 = self.r2  # The center radius of the ring
        r3 = self.r3  # The outer radius of the ring

        Xo = self.Xo  # Global X-coordinate of the center of the ring
        Yo = self.Yo  # Global Y-coordinate of the center of the ring
        Zo = self.Zo  # Global Z-coordinate of the center of the ring

        theta1 = 2*pi/self.n      # Angle between nodes at the inner radius of the ring
        theta2 = 2*pi/(self.n*3)  # Angle between nodes at the center of the ring
        theta3 = 2*pi/(self.n*3)  # Angle between nodes at the outer radius of the ring

        # Each node number will be increased by the offset calculated below
        node_offset = int(self.start_node[1:]) - 1

        # Each element number will be increased by the offset calculated below
        element_offset = int(self.start_element[1:]) - 1

        # Generate the nodes that make up the ring, working from the inside to the outside
        angle = 0
        for i in range(1, 6*n + 1, 1):

            # Assign the node a name
            node_name = 'N' + str(i + node_offset)

            # Generate the inner radius of nodes
            if i <= n:
                angle = theta1*(i - 1)
                x = Xo + r1*cos(angle)
                y = Yo + r1*sin(angle)
                z = Zo
            # Generate the center radius of nodes
            elif i <= 3*n:
                if (i - n) == 1:
                    angle = theta2
                elif (i - n) % 2 == 0:
                    angle += theta2
                else:
                    angle += 2*theta2
                x = Xo + r2*cos(angle)
                y = Yo + r2*sin(angle)
                z = Zo
            # Generate the outer radius of nodes
            else:
                if (i - 3*n) == 1:
                    angle = 0
                else:
                    angle = theta3*((i - 3*n) - 1)
                x = Xo + r3*cos(angle)
                y = Yo + r3*sin(angle)
                z = Zo
            
            self.nodes[node_name] = [node_name, x, y, z]

        # Generate the elements that make up the ring
        for i in range(1, 4*n + 1, 1):

            # Assign the element a name
            element_name = 'Q' + str(i + element_offset)

            if i <= n:
                n_node = i
                j_node = 2*i + n
                i_node = 2*i + n - 1
                if i != n:
                    m_node = i + 1
                else:
                    m_node = 1
            elif (i - n) % 3 == 1:
                n_node = 1 + (i - (n + 1))//3
                m_node = i - (i - (n + 1))//3
                j_node = i + 2*n + 1
                i_node = i + 2*n
            elif (i - n) % 3 == 2:
                n_node = i - 1 - (i - (n + 1))//3
                m_node = i - (i - (n + 1))//3
                j_node = i + 2*n + 1
                i_node = i + 2*n            
            else:
                n_node = i - 1 - (i - (n + 1))//3
                i_node = i + 2*n
                if i != 4*n:
                    m_node = 2 + (i - (n + 1))//3
                    j_node = i + 2*n + 1
                else:
                    m_node = 1
                    j_node = 1 + 3*n

            self.elements[element_name] = [element_name, 'N' + str(i_node + node_offset),
                                                         'N' + str(j_node + node_offset),
                                                         'N' + str(m_node + node_offset),
                                                         'N' + str(n_node + node_offset),
                                                         self.t, self.E, self.nu]

class FrustrumMesh(AnnulusMesh):

    def __init__(self, t, E, nu, mesh_size, large_radius, small_radius, height, center=[0, 0, 0], start_node='N1', start_element='Q1'):
        
        # Create an annulus mesh
        super().__init__(t, E, nu, mesh_size, large_radius, small_radius, center, start_node, start_element)

        Xo = center[0]
        Yo = center[1]

        # Adjust the Z-cooridnate of each node
        for node in self.nodes.values():
            X = node[1]
            Y = node[2]
            Z = node[3]
            r = ((X - Xo)**2 + (Y - Yo)**2)**0.5
            Z += (r - large_radius)/(large_radius - small_radius)*height
            node[3] = Z

