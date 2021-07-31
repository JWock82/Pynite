from PyNite.Node3D import Node3D
from PyNite.Quad3D import Quad3D
from PyNite.Plate3D import Plate3D
from math import pi, sin, cos, ceil
from bisect import bisect

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
    
    def max_shear(self, direction='Qx', combo=None):
        '''
        Returns the maximum shear in the mesh.
        
        Checks corner and center shears in all the elements in the mesh. The mesh must be part of
        a solved model prior to using this function.

        Parameters
        ----------
        direction : string, optional
            The direction to ge the maximum shear for. Options are 'Qx' or 'Qy'. Default
            is 'Qx'.
        combo : string, optional
            The name of the load combination to get the maximum shear for. If omitted, all load
            combinations will be evaluated.
        '''

        if direction == 'Qx':
            i = 0
        elif direction == 'Qy':
            i = 1
        else:
            raise Exception('Invalid direction specified for mesh shear results. Valid values are \'Qx\', or \'Qy\'')

        # Initialize the maximum value to None
        Q_max = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine whether the element is a rectangle or a quadrilateral
            if element.type == 'Rect':
                # Use the rectangle's local (x, y) coordinate system
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn, = 0, element.height()
            elif element.type == 'Quad':
                # Use the quad's natural (r, s) coordinate system
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1

            # Step through each load combination the element utilizes
            for load_combo in element.LoadCombos.values():

                # Determine if this load combination should be evaluated
                if combo == None or load_combo.name == combo:
                    
                    # Find the maximum shear in the element, checking each corner and the center
                    # of the element
                    Q_element = max([element.shear(xi, yi, load_combo.name)[i, 0],
                                     element.shear(xj, yj, load_combo.name)[i, 0],
                                     element.shear(xm, ym, load_combo.name)[i, 0],
                                     element.shear(xn, yn, load_combo.name)[i, 0],
                                     element.shear((xi + xj)/2, (yi + yn)/2, load_combo.name)[i, 0]])

                    # Determine if the maximum shear calculated is the largest encountered so far
                    if Q_max == None or Q_max < Q_element:
                        # Save this value if it's the largest
                        Q_max = Q_element
            
        # Return the largest value encountered from all the elements
        return Q_max
    
    def min_shear(self, direction='Qx', combo=None):
        '''
        Returns the minimum shear in the mesh.
        
        Checks corner and center shears in all the elements in the mesh. The mesh must be part of
        a solved model prior to using this function.

        Parameters
        ----------
        direction : string, optional
            The direction to ge the minimum shear for. Options are 'Qx' or 'Qy'. Default
            is 'Qx'.
        combo : string, optional
            The name of the load combination to get the minimum shear for. If omitted, all load
            combinations will be evaluated.
        '''

        if direction == 'Qx':
            i = 0
        elif direction == 'Qy':
            i = 1
        else:
            raise Exception('Invalid direction specified for mesh shear results. Valid values are \'Qx\', or \'Qy\'')

        # Initialize the minimum value to None
        Q_min = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine whether the element is a rectangle or a quadrilateral
            if element.type == 'Rect':
                # Use the rectangle's local (x, y) coordinate system
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn, = 0, element.height()
            elif element.type == 'Quad':
                # Use the quad's natural (r, s) coordinate system
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1

            # Step through each load combination the element utilizes
            for load_combo in element.LoadCombos.values():

                # Determine if this load combination should be evaluated
                if combo == None or load_combo.name == combo:
                    
                    # Find the minimum shear in the element, checking each corner and the center
                    # of the element
                    Q_element = min([element.shear(xi, yi, load_combo.name)[i, 0],
                                     element.shear(xj, yj, load_combo.name)[i, 0],
                                     element.shear(xm, ym, load_combo.name)[i, 0],
                                     element.shear(xn, yn, load_combo.name)[i, 0],
                                     element.shear((xi + xj)/2, (yi + yn)/2, load_combo.name)[i, 0]])

                    # Determine if the minimum shear calculated is the smallest encountered so far
                    if Q_min == None or Q_min > Q_element:
                        # Save this value if it's the smallest
                        Q_min = Q_element
            
        # Return the smallest value encountered from all the elements
        return Q_min

    def max_moment(self, direction='Mx', combo=None):
        '''
        Returns the maximum moment in the mesh.
        
        Checks corner and center moments in all the elements in the mesh. The mesh must be part of
        a solved model prior to using this function.

        Parameters
        ----------
        direction : string, optional
            The direction to ge the maximum moment for. Options are 'Mx', 'My', or 'Mxy'. Default
            is 'Mx'.
        combo : string, optional
            The name of the load combination to get the maximum moment for. If omitted, all load
            combinations will be evaluated.
        '''

        if direction == 'Mx':
            i = 0
        elif direction == 'My':
            i = 1
        elif direction == 'Mxy':
            i = 2
        else:
            raise Exception('Invalid direction specified for mesh moment results. Valid values are \'Mx\', \'My\', or \'Mxy\'')

        # Initialize the maximum value to None
        M_max = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine whether the element is a rectangle or a quadrilateral
            if element.type == 'Rect':
                # Use the rectangle's local (x, y) coordinate system
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn, = 0, element.height()
            elif element.type == 'Quad':
                # Use the quad's natural (r, s) coordinate system
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1

            # Step through each load combination the element utilizes
            for load_combo in element.LoadCombos.values():

                # Determine if this load combination should be evaluated
                if combo == None or load_combo.name == combo:
                    
                    # Find the maximum moment in the element, checking each corner and the center
                    # of the element
                    M_element = max([element.moment(xi, yi, load_combo.name)[i, 0],
                                     element.moment(xj, yj, load_combo.name)[i, 0],
                                     element.moment(xm, ym, load_combo.name)[i, 0],
                                     element.moment(xn, yn, load_combo.name)[i, 0],
                                     element.moment((xi + xj)/2, (yi + yn)/2, load_combo.name)[i, 0]])

                    # Determine if the maximum moment calculated is the largest encountered so far
                    if M_max == None or M_max < M_element:
                        # Save this value if it's the largest
                        M_max = M_element
            
        # Return the largest value encountered from all the elements
        return M_max
    
    def min_moment(self, direction='Mx', combo=None):
        '''
        Returns the minimum moment in the mesh.
        
        Checks corner and center moments in all the elements in the mesh. The mesh must be part of
        a solved model prior to using this function.

        Parameters
        ----------
        direction : string, optional
            The direction to ge the minimum moment for. Options are 'Mx', 'My', or 'Mxy'. Default
            is 'Mx'.
        combo : string, optional
            The name of the load combination to get the minimum moment for. If omitted, all load
            combinations will be evaluated.
        '''

        if direction == 'Mx':
            i = 0
        elif direction == 'My':
            i = 1
        elif direction == 'Mxy':
            i = 2
        else:
            raise Exception('Invalid direction specified for mesh moment results. Valid values are \'Mx\', \'My\', or \'Mxy\'')

        # Initialize the minimum value to None
        M_min = None

        # Step through each element in the mesh
        for element in self.elements.values():

            # Determine whether the element is a rectangle or a quadrilateral
            if element.type == 'Rect':
                # Use the rectangle's local (x, y) coordinate system
                xi, yi = 0, 0
                xj, yj = element.width(), 0
                xm, ym = element.width(), element.height()
                xn, yn, = 0, element.height()
            elif element.type == 'Quad':
                # Use the quad's natural (r, s) coordinate system
                xi, yi = -1, -1
                xj, yj = 1, -1
                xm, ym = 1, 1
                xn, yn = -1, 1

            # Step through each load combination the element utilizes
            for load_combo in element.LoadCombos.values():

                # Determine if this load combination should be evaluated
                if combo == None or load_combo.name == combo:
                    
                    # Find the minimum moment in the element, checking each corner and the center
                    # of the element
                    M_element = min([element.moment(xi, yi, load_combo.name)[i, 0],
                                     element.moment(xj, yj, load_combo.name)[i, 0],
                                     element.moment(xm, ym, load_combo.name)[i, 0],
                                     element.moment(xn, yn, load_combo.name)[i, 0],
                                     element.moment((xi + xj)/2, (yi + yn)/2, load_combo.name)[i, 0]])

                    # Determine if the minimum moment calculated is the smallest encountered so far
                    if M_min == None or M_min > M_element:
                        # Save this value if it's the smallest
                        M_min = M_element
            
        # Return the smallest value encountered from all the elements
        return M_min
    
#%%
class RectangleMesh(Mesh):

    def __init__(self, t, E, nu, mesh_size, width, height, origin=[0, 0, 0], plane='XY', 
                 x_control=[], y_control=[], start_node='N1', start_element='Q1',
                 element_type='Quad'):
        """
        A rectangular mesh of elements.

        Parameters
        ----------
        t : number
            Element thickness.
        E : number
            Element modulus of elasticity.
        nu : number
            Element poisson's ratio.
        mesh_size : number
            Desired mesh size.
        width : number
            The overall width of the mesh measured along its local x-axis.
        height : number
            The overall height of the mesh measured along its local y-axis.
        origin : list, optional
            The origin of the rectangular mesh's local coordinate system. The default is [0, 0, 0].
        plane : string, optional
            The plane the mesh will be parallel to. Options are 'XY', 'YZ', and 'XZ'. The default
            is 'XY'.
        x_control : list, optional
            A list of control points along the mesh's local x-axis work into the mesh.
        y_control : list, optional
            A list of control points along the mesh's local y-axis work into the mesh.
        start_node : string, optional
            A unique name for the first node in the mesh. The default is 'N1'.
        start_element : string, optional
            A unique name for the first element in the mesh. The default is 'Q1' or 'R1' depending
            on the type of element selected.
        element_type : string, optional
            The type of element to make the mesh out of. Either 'Quad' or 'Rect'. The default is
            'Quad'.

        Returns
        -------
        A new rectangular mesh object.

        """
        
        super().__init__(t, E, nu, start_node, start_element)
        self.mesh_size = mesh_size
        self.width = width
        self.height = height
        self.origin = origin
        self.plane = plane
        self.x_control = x_control
        self.y_control = y_control
        self.element_type = element_type
        self.__mesh()
    
    def __mesh(self):

        mesh_size = self.mesh_size
        width = self.width
        height = self.height
        Xo = self.origin[0]
        Yo = self.origin[1]
        Zo = self.origin[2]
        plane = self.plane
        x_control = self.x_control
        y_control = self.y_control
        element_type = self.element_type

        # Add the mesh's boundaries to the list of control points
        x_control.append(0)
        x_control.append(width)
        y_control.append(0)
        y_control.append(height)

        # Sort the control points and remove duplicate values
        x_control = sorted(set(x_control))
        y_control = sorted(set(y_control))
        
        # Each node number will be increased by the offset calculated below
        node_offset = int(self.start_node[1:]) - 1

        # Each element number will be increased by the offset calculated below
        element_offset = int(self.start_element[1:]) - 1

        # Determine which prefix to assign to new elements
        if element_type == 'Quad':
            element_prefix = 'Q'
        elif element_type == 'Rect':
            element_prefix = 'R'
        else:
            raise Exception('Invalid element type specified for RectangleMesh. Select \'Quad\' or \'Rect\'.')

        # Initialize node numbering
        node_num = 1

        # Step through each y control point (except the first one which is always zero)
        num_rows = 0
        num_cols = 0
        y, h = 0, None
        for j in range(1, len(y_control), 1):
            
            # If this is not the first iteration 'y' will be too high at this point.
            if j != 1:
                y -= h

            # Determine the mesh size between this y control point and the previous one
            ny = max(1, (y_control[j] - y_control[j - 1])/mesh_size)
            h = (y_control[j] - y_control[j - 1])/ceil(ny)

            # Adjust 'y' if this is not the first iteration.
            if j != 1:
                y += h

            # Generate nodes between the y control points
            while round(y, 10) <= round(y_control[j], 10):
                
                # Count the number of rows of plates as we go
                num_rows += 1

                # Step through each x control point (except the first one which is always zero)
                x, b = 0, None
                for i in range(1, len(x_control), 1):
                    
                    # 'x' needs to be adjusted for the same reasons 'y' needed to be adjusted
                    if i != 1:
                        x -= b

                    # Determine the mesh size between this x control point and the previous one
                    nx = max(1, (x_control[i] - x_control[i - 1])/mesh_size)
                    b = (x_control[i] - x_control[i - 1])/ceil(nx)

                    if i != 1:
                        x += b

                    # Generate nodes between the x control points
                    while round(x, 10) <= round(x_control[i], 10):
                        
                        # Count the number of columns of plates as we go
                        if y == 0:
                            num_cols += 1

                        # Assign the node a name
                        node_name = 'N' + str(node_num + node_offset)

                        # Calculate the node's coordinates
                        if plane == 'XY':
                            X = Xo + x
                            Y = Yo + y
                            Z = Zo + 0
                        elif plane == 'YZ':
                            X = Xo + 0
                            Y = Yo + y
                            Z = Zo + x
                        elif plane == 'XZ':
                            X = Xo + x
                            Y = Yo + 0
                            Z = Zo + y
                        else:
                            raise Exception('Invalid plane selected for RectangleMesh.')

                        # Add the node to the mesh
                        self.nodes[node_name] = Node3D(node_name, X, Y, Z)

                        # Move to the next x coordinate
                        x += b

                        # Move to the next node number
                        node_num += 1

                # Move to the next y coordinate
                y += h
        
        # At this point `num_cols` and `num_rows` represent the number of columns and rows of
        # nodes. We'll adjust these variables to be the number of columns and rows of elements
        # instead.
        num_cols -= 1
        num_rows -= 1
        
        # Create the elements
        r = 1
        n = 1
        for i in range(1, num_cols*num_rows + 1, 1):

            # Assign the element a name
            element_name = element_prefix + str(i + element_offset)

            # Find the attached nodes
            i_node = n + (r - 1)
            j_node = i_node + 1
            m_node = j_node + (num_cols + 1)
            n_node = m_node - 1

            if i % num_cols == 0:
                r += 1
            
            n += 1
            
            if element_type == 'Quad':
                self.elements[element_name] = Quad3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                                   self.nodes['N' + str(j_node + node_offset)],
                                                                   self.nodes['N' + str(m_node + node_offset)],
                                                                   self.nodes['N' + str(n_node + node_offset)],
                                                                   self.t, self.E, self.nu)
            else:
                self.elements[element_name] = Plate3D(element_name, self.nodes['N' + str(i_node + node_offset)],
                                                                    self.nodes['N' + str(j_node + node_offset)],
                                                                    self.nodes['N' + str(m_node + node_offset)],
                                                                    self.nodes['N' + str(n_node + node_offset)],
                                                                    self.t, self.E, self.nu)

#%%           
class AnnulusMesh(Mesh):
    """
    A mesh of quadrilaterals forming an annulus (a donut).
    """

    def __init__(self, t, E, nu, mesh_size, outer_radius, inner_radius, origin=[0, 0, 0], axis='Y',
                 start_node='N1', start_element='Q1'):

        super().__init__(t, E, nu, start_node, start_element)

        self.r1 = inner_radius
        self.r2 = outer_radius
        self.mesh_size = mesh_size
        self.origin = origin
        self.axis = axis

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
                ring = AnnulusTransRingMesh(t, E, nu, r_inner + h_rad, r_inner, n_circ,
                                            self.origin, self.axis, 'N' + str(n), 'Q' + str(q))
                n += 3*n_circ
                q += 4*n_circ
                n_circ *= 3
                self.num_quads_outer = n_circ
            else:
                ring = AnnulusRingMesh(t, E, nu, r_inner + h_rad, r_inner, n_circ, self.origin,
                                       self.axis, 'N' + str(n), 'Q' + str(q))
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

    def __init__(self, t, E, nu, outer_radius, inner_radius, num_quads, origin=[0, 0, 0], axis='Y',
                 start_node='N1', start_element='Q1'):

        super().__init__(t, E, nu, start_node=start_node, start_element=start_element)

        self.r1 = inner_radius
        self.r2 = outer_radius
        self.n = num_quads
        self.Xo = origin[0]
        self.Yo = origin[1]
        self.Zo = origin[2]

        self.axis = axis

        # Generate the nodes and elements
        self.__mesh()

    def __mesh(self):

        n = self.n  # Number of plates in the initial ring

        r1 = self.r1  # The inner radius of the ring
        r2 = self.r2  # The outer radius of the ring

        Xo = self.Xo  # Global X-coordinate of the center of the ring
        Yo = self.Yo  # Global Y-coordinate of the center of the ring
        Zo = self.Zo  # Global Z-coordinate of the center of the ring

        axis = self.axis

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
                if axis == 'Y':
                    x = Xo + r1*cos(angle)
                    y = Yo
                    z = Zo + r1*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + r1*sin(angle)
                    z = Zo + r1*cos(angle)
                elif axis == 'Z':
                    x = Xo + r1*sin(angle)
                    y = Yo + r1*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusRingMesh.')
            
            # Generate the outer radius of nodes
            else:
                angle = theta*((i - n) - 1)
                if axis == 'Y':
                    x = Xo + r2*cos(angle)
                    y = Yo 
                    z = Zo + r2*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + r2*sin(angle)
                    z = Zo + r2*cos(angle)
                elif axis == 'Z':
                    x = Xo + r2*sin(angle)
                    y = Yo + r2*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusRingMesh.')
            
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

    def __init__(self, t, E, nu, outer_radius, inner_radius, num_inner_quads, origin=[0, 0, 0],
                 axis='Y', start_node='N1', start_element='Q1'):
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
        self.Xo = origin[0]
        self.Yo = origin[1]
        self.Zo = origin[2]

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

        axis = self.axis

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
                if axis == 'Y':
                    x = Xo + r1*cos(angle)
                    y = Yo
                    z = Zo + r1*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + r1*sin(angle)
                    z = Zo + r1*cos(angle)
                elif axis == 'Z':
                    x = Xo + r1*sin(angle)
                    y = Yo + r1*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusTransRingMesh.')
            
            # Generate the center radius of nodes
            elif i <= 3*n:
                if (i - n) == 1:
                    angle = theta2
                elif (i - n) % 2 == 0:
                    angle += theta2
                else:
                    angle += 2*theta2
                if axis == 'Y':
                    x = Xo + r2*cos(angle)
                    y = Yo 
                    z = Zo + r2*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + r2*sin(angle)
                    z = Zo + r2*cos(angle)
                elif axis == 'Z':
                    x = Xo + r2*sin(angle)
                    y = Yo + r2*cos(angle)
                    z = Zo
            # Generate the outer radius of nodes
            else:
                if (i - 3*n) == 1:
                    angle = 0
                else:
                    angle = theta3*((i - 3*n) - 1)
                if axis == 'Y':
                    x = Xo + r3*cos(angle)
                    y = Yo 
                    z = Zo + r3*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + r3*sin(angle)
                    z = Zo + r3*cos(angle)
                elif axis == 'Z':
                    x = Xo + r3*sin(angle)
                    y = Yo + r3*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for AnnulusTransRingMesh.')
            
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

    def __init__(self, t, E, nu, mesh_size, large_radius, small_radius, height, origin=[0, 0, 0],
                 axis='Y', start_node='N1', start_element='Q1'):
        
        # Create an annulus mesh
        super().__init__(t, E, nu, mesh_size, large_radius, small_radius, origin, start_node, start_element)

        Xo = origin[0]
        Yo = origin[1]
        Zo = origin[2]

        # Adjust the cooridnates of each node to make a frustrum
        for node in self.nodes.values():
            X = node.X
            Y = node.Y
            Z = node.Z
            r = ((X - Xo)**2 + (Z - Zo)**2)**0.5
            if axis == 'Y':
                node.Y += (r - large_radius)/(large_radius - small_radius)*height
            elif axis == 'X':
                node.X += (r - large_radius)/(large_radius - small_radius)*height
            elif axis == 'Z':
                node.Z += (r - large_radius)/(large_radius - small_radius)*height
            else:
                raise Exception('Invalid axis specified for frustrum mesh.')

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

    def __init__(self, t, E, nu, radius, height, num_quads, origin=[0, 0, 0], axis='Y',
                 start_node='N1', start_element='Q1'):

        super().__init__(t, E, nu, start_node=start_node, start_element=start_element)

        self.radius = radius
        self.height = height

        self.num_quads = num_quads

        self.Xo = origin[0]
        self.Yo = origin[1]
        self.Zo = origin[2]

        self.axis = axis

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

        axis = self.axis
        
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
                if axis == 'Y':
                    x = Xo + radius*cos(angle)
                    y = Yo
                    z = Zo + radius*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + radius*sin(angle)
                    z = Zo + radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + radius*sin(angle)
                    y = Yo + radius*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for CylinderRingMesh.')

            # Generate the top nodes of the ring
            else:
                angle = theta*((i - n) - 1)
                if axis == 'Y':
                    x = Xo + radius*cos(angle)
                    y = Yo
                    z = Zo + radius*sin(angle)
                elif axis == 'X':
                    x = Xo
                    y = Yo + radius*sin(angle)
                    z = Zo + radius*cos(angle)
                elif axis == 'Z':
                    x = Xo + radius*sin(angle)
                    y = Yo + radius*cos(angle)
                    z = Zo
                else:
                    raise Exception('Invalid axis specified for CylinderRingMesh.')
            
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
