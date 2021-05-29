from PyNite.Node3D import Node3D
from PyNite.Quad3D import Quad3D
from math import pi, sin, cos

#%%
class Mesh():
    """
    A parent class for meshes to inherit from.
    """

    def __init__(self, t, E, nu, start_node='N1', start_element='Q1'):

        self.t = t                          # Thickness
        self.E = E                          # Modulus of elasticity
        self.nu = nu                        # Poisson's ratio
        self.start_node = start_node        # The name of the first node in the mesh
        self.start_element = start_element  # The name of the first element in the mesh
        self.last_element = None            # The name of the last element in the mesh
        self.nodes = {}                     # A dictionary containing the nodes in the mesh
        self.elements = {}                  # A dictionary containing the elements in the mesh

#%%           
class AnnulusMesh(Mesh):
    """
    A mesh of quadrilaterals forming an annulus (a donut).
    """

    def __init__(self, t, E, nu, mesh_size, outer_radius, inner_radius, center=[0, 0, 0], start_node='N1', start_element='Q1'):

        super().__init__(t, E, nu, start_node, start_element)

        self.r1 = inner_radius
        self.r2 = outer_radius
        self.mesh_size = mesh_size
        self.center = center

        self.num_quads_inner = None
        self.num_quads_outer = None

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

        circumf = 2*pi*r_inner                                 # Circumference of the ring at the inner radius
        n_circ = int(circumf/mesh_size)  # Number of times `mesh_size` fits in the circumference
        self.num_quads_outer = n_circ

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
                self.num_quads_outer = n_circ
            else:
                ring = AnnulusRingMesh(t, E, nu, r_inner + h_rad, r_inner, n_circ, self.center, 'N' + str(n), 'Q' + str(q))
                n += n_circ
                q += n_circ
        
            # Add the newly generated nodes and elements to the overall mesh. Note that if duplicate
            # keys exist, the `.update()` method will overwrite them with the newly generated key value
            # pairs. This works in our favor by automatically eliminating duplicate nodes from the
            # dictionary.
            self.nodes.update(ring.nodes)
            self.elements.update(ring.elements)

            # Prepare to move to the next ring
            r_inner += h_rad

        # After calling the `.update()` method some elements are still attached to the duplicate
        # nodes that are no longer in the dictionary. Attach these plates to the nodes that are
        # still in the dictionary instead. 
        for element in self.elements.values():
            element.iNode = self.nodes[element.iNode.Name]
            element.jNode = self.nodes[element.jNode.Name]
            element.mNode = self.nodes[element.mNode.Name]
            element.nNode = self.nodes[element.nNode.Name]

#%%
class AnnulusRingMesh(Mesh):
    """
    A mesh of quadrilaterals forming an annular ring (a donut).
    """

    def __init__(self, t, E, nu, outer_radius, inner_radius, num_quads, center=[0, 0, 0], start_node='N1', start_element='Q1'):

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

        n = self.n  # Number of plates in the initial ring

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
                y = Yo
                z = Zo + r1*sin(angle)
            # Generate the outer radius of nodes
            else:
                angle = theta*((i - n) - 1)
                x = Xo + r2*cos(angle)
                y = Yo 
                z = Zo + r2*sin(angle)
            
            self.nodes[node_name] = Node3D(node_name, x, y, z)

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

            self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                               self.nodes['N' + str(j_node + node_offset)],
                                                               self.nodes['N' + str(m_node + node_offset)],
                                                               self.nodes['N' + str(n_node + node_offset)],
                                                               self.t, self.E, self.nu)

#%%
class AnnulusTransRingMesh(Mesh):
    """
    A mesh of quadrilaterals forming an annular ring (a donut) with the mesh getting finer on the outer
    edge.
    """

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
                y = Yo 
                z = Zo + r1*sin(angle)
            # Generate the center radius of nodes
            elif i <= 3*n:
                if (i - n) == 1:
                    angle = theta2
                elif (i - n) % 2 == 0:
                    angle += theta2
                else:
                    angle += 2*theta2
                x = Xo + r2*cos(angle)
                y = Yo
                z = Zo + r2*sin(angle)
            # Generate the outer radius of nodes
            else:
                if (i - 3*n) == 1:
                    angle = 0
                else:
                    angle = theta3*((i - 3*n) - 1)
                x = Xo + r3*cos(angle)
                y = Yo
                z = Zo + r3*sin(angle)
            
            self.nodes[node_name] = Node3D(node_name, x, y, z)

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

            self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                               self.nodes['N' + str(j_node + node_offset)],
                                                               self.nodes['N' + str(m_node + node_offset)],
                                                               self.nodes['N' + str(n_node + node_offset)],
                                                               self.t, self.E, self.nu)

#%%
class FrustrumMesh(AnnulusMesh):
    """
    A mesh of quadrilaterals forming a frustrum (a cone intersected by a horizontal plane at the top and bottom).
    """

    def __init__(self, t, E, nu, mesh_size, large_radius, small_radius, height, center=[0, 0, 0], start_node='N1', start_element='Q1'):
        
        # Create an annulus mesh
        super().__init__(t, E, nu, mesh_size, large_radius, small_radius, center, start_node, start_element)

        Xo = center[0]
        Zo = center[2]

        # Adjust the Z-cooridnate of each node
        for node in self.nodes.values():
            X = node.X
            Y = node.Y
            Z = node.Z
            r = ((X - Xo)**2 + (Z - Zo)**2)**0.5
            Y += (r - large_radius)/(large_radius - small_radius)*height
            node.Y = Y

#%%
class CylinderMesh(Mesh):
    """
    A mesh of quadrilaterals forming a cylinder.

    Parameters
    ----------
    mesh_size : number
        The desired mesh size. This value will only be used to mesh vertically if `num_quads` is
        specified. Otherwise it will be used to mesh the circumference too.
    start_node : string, optional
        The name of the first node in the mesh. The name must be formatted starting with a single
        letter followed by a number (e.g. 'N12'). The mesh will begin numbering nodes from this
        number. The default is 'N1'. 
    start_element : string, optional
        The name of the first element in the mesh. The name must be formatted starting with a
        single letter followed by a number (e.g. 'Q32'). The mesh will begin numbering elements
        from this number. The default is 'Q1'.
    num_quads : number, optional
        The number of quadrilaterals to divide the circumference into. If this value is omitted
        `mesh_size` will be used instead to calculate the number of quadrilaterals in the
        circumference. The default is `None`.
    """

    def __init__(self, t, E, nu, mesh_size, radius, height, center=[0, 0, 0], start_node='N1', start_element='Q1', num_quads=None):

        super().__init__(t, E, nu, start_node, start_element)

        self.radius = radius
        self.h = height
        self.mesh_size = mesh_size
        self.num_quads = num_quads
        self.center = center

        self.__mesh()
    
    def __mesh(self):
        
        t = self.t
        E = self.E
        nu = self.nu

        mesh_size = self.mesh_size     # Desired mesh size
        num_quads = self.num_quads  # Number of quadrilaterals in the ring
        n = self.num_quads

        radius = self.radius
        h = self.h
        y = self.center[1]
        n = int(self.start_node[1:])
        q = int(self.start_element[1:])

        # Determine the number of quads to mesh the circumference into
        if num_quads == None:
            num_quads = int(2*pi/mesh_size)

        # Mesh the cylinder from the bottom toward the top
        while round(y, 10) < round(h, 10):
            
            height = h - y                                # Remaining height to be meshed
            n_vert = int(height/mesh_size)  # Number of times the plate height fits in the remaining unmeshed height
            h_y = height/n_vert                         # Element height in the vertical direction
        
            # Create a mesh of nodes for the ring
            ring = CylinderRingMesh(t, E, nu, radius, h_y, num_quads, [0, y, 0], 'N' + str(n), 'Q' + str(q))
            n += num_quads
            q += num_quads
        
            # Add the newly generated nodes and elements to the overall mesh. Note that if duplicate
            # keys exist, the `.update()` method will overwrite them with the newly generated key value
            # pairs. This works in our favor by automatically eliminating duplicate nodes at the shared
            # boundaries between rings.
            self.nodes.update(ring.nodes)
            self.elements.update(ring.elements)

            # Prepare to move to the next ring
            y += h_y
        
        # After calling the `.update()` method some elements are still attached to the duplicate
        # nodes that are no longer in the dictionary. Attach these plates to the nodes that are
        # still in the dictionary instead. 
        for element in self.elements.values():
            element.iNode = self.nodes[element.iNode.Name]
            element.jNode = self.nodes[element.jNode.Name]
            element.mNode = self.nodes[element.mNode.Name]
            element.nNode = self.nodes[element.nNode.Name]

#%%
class CylinderRingMesh(Mesh):
    """
    A mesh of quadrilaterals forming a cylindrical ring.

    Parameters
    ----------
    start_node : string, optional
        The name of the first node in the mesh. The name must be formatted starting with a single
        letter followed by a number (e.g. 'N12'). The mesh will begin numbering nodes from this
        number. The default is 'N1'. 
    start_element : string, optional
        The name of the first element in the mesh. The name must be formatted starting with a
        single letter followed by a number (e.g. 'Q32'). The mesh will begin numbering elements
        from this number. The default is 'Q1'.
    num_quads : number
        The number of quadrilaterals to divide the circumference into.
    """

    def __init__(self, t, E, nu, radius, height, num_quads, center=[0, 0, 0], start_node='N1', start_element='Q1'):

        super().__init__(t, E, nu, start_node=start_node, start_element=start_element)

        self.radius = radius
        self.height = height

        self.num_quads = num_quads

        self.Xo = center[0]
        self.Yo = center[1]
        self.Zo = center[2]

        # Generate the nodes and elements
        self.__mesh()

    def __mesh(self):
        """
        Generates the nodes and elements in the mesh.
        """

        num_quads = self.num_quads  # Number of quadrilaterals in the ring
        n = self.num_quads

        radius = self.radius  # The radius of the ring
        height = self.height  # The height of the ring

        Xo = self.Xo  # Global X-coordinate of the center of the bottom of the ring
        Yo = self.Yo  # Global Y-coordinate of the center of the bottom of the ring
        Zo = self.Zo  # Global Z-coordinate of the center of the bottom of the ring
        
        # Calculate the angle between nodes in the circumference of the ring
        theta = 2*pi/num_quads

        # Each node number will be increased by the offset calculated below
        node_offset = int(self.start_node[1:]) - 1

        # Each element number will be increased by the offset calculated below
        element_offset = int(self.start_element[1:]) - 1

        # Generate the nodes that make up the ring
        angle = 0
        for i in range(1, 2*n + 1, 1):

            # Assign the node a name
            node_name = 'N' + str(i + node_offset)

            # Generate the bottom nodes of the ring
            if i <= n:
                angle = theta*(i - 1)
                x = Xo + radius*cos(angle)
                y = Yo
                z = Zo + radius*sin(angle)
            # Generate the top nodes of the ring
            else:
                angle = theta*((i - n) - 1)
                x = Xo + radius*cos(angle)
                y = Yo + height
                z = Zo + radius*sin(angle)
            
            self.nodes[node_name] = Node3D(node_name, x, y, z)

        # Generate the elements that make up the ring
        for i in range(1, n + 1, 1):

            # Assign the element a name
            element_name = 'Q' + str(i + element_offset)
            
            # Assign nodes to the element
            n_node = i
            i_node = i + n
            if i != n:
                m_node = i + 1
                j_node = i + 1 + n
            else:
                m_node = 1
                j_node = 1 + n

            # Create the element and add it to the `elements` dictionary
            self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                               self.nodes['N' + str(j_node + node_offset)],
                                                               self.nodes['N' + str(m_node + node_offset)],
                                                               self.nodes['N' + str(n_node + node_offset)],
                                                               self.t, self.E, self.nu)
