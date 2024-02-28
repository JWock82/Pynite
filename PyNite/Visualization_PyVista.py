
from json import load
import warnings

from IPython.display import Image
from numpy import array, empty, append, cross, zeros
from numpy.linalg import norm
import pyvista as pv

class Renderer:
    """Used to render finite element models.
    """

    scalar = None

    def __init__(self, model):

        self.model = model

        # Default settings for rendering
        self._annotation_size = 5
        self._deformed_shape = False
        self._deformed_scale = 30
        self._render_nodes = True
        self._render_loads = True
        self._color_map = None
        self._combo_name = 'Combo 1'
        self._case = None
        self._labels = True
        self._scalar_bar = False
        self._scalar_bar_text_size = 24
        self.theme = 'default'

        self.plotter = pv.Plotter()
        self.plotter.set_background('white')  # Setting background color
        # self.plotter.add_logo_widget('./Resources/Full Logo No Buffer - Transparent.png')
        # self.plotter.view_isometric()
        self.plotter.view_xy()
        self.plotter.show_axes()

    @property
    def window_width(self):
        return self.plotter.window_size[0]

    @window_width.setter
    def window_width(self, width):
        height = self.plotter.window_size[1]
        self.plotter.window_size = (width, height)

    @property
    def window_height(self):
        return self.plotter.window_size[1]

    @window_height.setter
    def window_height(self, height):
        width = self.plotter.window_size[0]
        self.plotter.window_size = (width, height)

    @property
    def annotation_size(self):
        return self._annotation_size
    
    @annotation_size.setter
    def annotation_size(self, size):
        self._annotation_size = size
    
    @property
    def deformed_shape(self):
        return self._deformed_shape
    
    @deformed_shape.setter
    def deformed_shape(self, deformed_shape):
        self._deformed_shape = deformed_shape
    
    @property
    def deformed_scale(self):
        return self._deformed_scale
    
    @deformed_scale.setter
    def deformed_scale(self, scale):
        self._deformed_scale = scale
    
    @property
    def render_nodes(self):
        return self._render_nodes
    
    @render_nodes.setter
    def render_nodes(self, render_nodes):
        self._render_nodes = render_nodes

    @property
    def render_loads(self):
        return self._render_loads
    
    @render_loads.setter
    def render_loads(self, render_loads):
        self._render_loads = render_loads
    
    @property
    def color_map(self):
        return self._color_map
    
    @color_map.setter
    def color_map(self, color_map):
        self._color_map = color_map
    
    @property
    def combo_name(self):
        return self._combo_name
    
    @combo_name.setter
    def combo_name(self, combo_name):
        self._combo_name = combo_name
        self._case = None
    
    @property
    def case(self):
        return self._case
    
    @case.setter
    def case(self, case):
        self._case = case
        self._combo_name = None
    
    @property
    def show_labels(self):
        return self._labels
    
    @show_labels.setter
    def show_labels(self, show_labels):
        self._labels = show_labels
    
    @property
    def scalar_bar(self):
        return self._scalar_bar
    
    @scalar_bar.setter
    def scalar_bar(self, scalar_bar):
        self._scalar_bar = scalar_bar
    
    @property
    def scalar_bar_text_size(self):
        return self._scalar_bar_text_size
    
    @scalar_bar_text_size.setter
    def scalar_bar_text_size(self, text_size):
        self._scalar_bar_text_size = text_size

    def render_model(self, interact=True, reset_camera=True):
        """
        Renders the model in a window

        Parameters
        ----------
        interact : bool
            Suppresses interacting with the window if set to `False`. This can be used to capture a
            screenshot without pausing the program for the user to interact. Default is `True`.
        reset_camera : bool
            Resets the camera if set to `True`. Default is `True`.
        """

        # Update the plotter with the latest geometry
        self.update(reset_camera)

        # Render the model (code execution will pause here until the user closes the window)
        self.plotter.show(title='Pynite - Simple Finite Element Analysis for Python', window_size=None, interactive=True, auto_close=None, interactive_update=False, full_screen=None, screenshot=False, return_img=False, cpos=None, jupyter_backend=None, return_viewer=False, return_cpos=None, before_close_callback=None)

    def screenshot(self, filepath='console', interact=True, reset_camera=True):
        """
        Renders the model in a window. When the window is closed a screenshot is captured.

        Parameters
        ----------
        filepath : string
            Sends a screenshot to the specified filepath. The screenshot will be taken when the
            user closes out of the render window.
        interact : bool
            Suppresses interacting with the window if set to `False`. This can be used to capture a
            screenshot without pausing the program for the user to interact. Default is `True`.
        reset_camera : bool
            Resets the camera if set to `True`. Default is `True`.
        """

        # Update the plotter with the latest geometry
        self.update(reset_camera)

        # Show the plotter for interaction
        if interact == True:
            self.plotter.show(auto_close=False)

        # Render the model (code execution will pause here until the user closes the window)
        self.plotter.show(title='Pynite - Simple Finite Element Analysis for Python', window_size=None, interactive=True, auto_close=None, interactive_update=False, full_screen=None, screenshot=filepath, return_img=False, cpos=None, jupyter_backend=None, return_viewer=False, return_cpos=None, before_close_callback=None)

    def update(self, reset_camera=True):
        """
        Builds or rebuilds the pyvista plotter

        Parameters
        ----------
        reset_camera : bool
            Resets the camera if set to `True`. Default is `True`.
        """

        # Input validation
        if self.deformed_shape and self.case != None:
            raise Exception('Deformed shape is only available for load combinations,'
                            ' not load cases.')
        if self.model.LoadCombos == {} and self.render_loads == True and self.case == None:
            self.render_loads = False
            warnings.warn('Unable to render load combination. No load combinations defined.', UserWarning)

        # Clear out the old plot (if any)
        self.plotter.clear()
        
        # Check if nodes are to be rendered
        if self.render_nodes == True:

            if self.theme == 'print':
                color = 'black'
            else:
                color = 'grey'

            # Plot each node in the model
            for node in self.model.Nodes.values():
                self.plot_node(node, self.annotation_size, color)
            
            # Plot each auxiliary node in the model
            for aux_node in self.model.AuxNodes.values():
                self.plot_node(aux_node, self.annotation_size, color)
        
        # Render node labels
        label_points = [[node.X, node.Y, node.Z] for node in self.model.Nodes.values()]
        labels = [node.name for node in self.model.Nodes.values()]
        self.plotter.add_point_labels(label_points, labels, italic=False, bold=False, font_size=24, text_color='grey', font_family=None, shadow=False, show_points=True, point_color=None, point_size=None, name=None, shape_color=None, shape=None, fill_shape=True, margin=3, shape_opacity=1.0, pickable=False, render_points_as_spheres=True, tolerance=0.001, reset_camera=None, always_visible=False, render=True)

        # Create a visual spring for each spring in the model
        # vis_springs = []
        # for spring in self.model.Springs.values():
        #     vis_springs.append(VisSpring(spring, self.model.Nodes, self.annotation_size))    
    
        # Create a visual member for each member in the model
        # vis_members = []
        # for member in self.model.Members.values():
        #     vis_members.append(VisMember(member, self.model.Nodes, self.annotation_size, self.theme))
        
        for member in self.model.Members.values():
            self.plot_member(member)

        # Render the member labels
        label_points = [[(member.i_node.X+member.j_node.X)/2, (member.i_node.Y+member.j_node.Y)/2, (member.i_node.Z+member.j_node.Z)/2] for member in self.model.Members.values()]
        labels = [member.name for member in self.model.Members.values()]
        self.plotter.add_point_labels(label_points, labels, italic=False, bold=False, font_size=24, text_color='grey', font_family=None, shadow=False, show_points=False, point_color=None, point_size=None, name=None, shape_color=None, shape=None, fill_shape=True, margin=3, shape_opacity=1.0, pickable=False, render_points_as_spheres=True, tolerance=0.001, reset_camera=None, always_visible=False, render=True)

        # Add actors for each spring
        # for vis_spring in vis_springs:
        
        #     # Add the spring
        #     self.plotter.add_mesh(vis_spring.mesh)

        #     if self.labels == True:

        #         # Add the actor for the spring label
        #         self.plotter.add_mesh(vis_spring.text_mesh)  

        # Add actors for each member
        # for vis_member in vis_members:
            
        #     # Add the actor for the member
        #     self.plotter.add_mesh(vis_member.mesh)

        #     if self.labels == True:
                
        #         # Add the actor for the member label
        #         self.plotter.add_mesh(vis_member.text_mesh)

        # Render the deformed shape if requested
        if self.deformed_shape == True:

            # Render deformed nodes
            for node in self.model.Nodes.values():
                self.plot_deformed_node(node, self.deformed_scale)
            
            # Render deformed members
            for member in self.model.Members.values():
                self.plot_deformed_member(member, self.deformed_scale)

            # _DeformedShape(self.model, self.deformed_scale, self.annotation_size, self.combo_name, self.render_nodes, self.theme)

        # # Render the loads if requested
        # if (self.combo_name != None or self.case != None) and self.render_loads != False:
        #     _RenderLoads(self.model, renderer, self.annotation_size, self.combo_name, self.case, self.theme)
        
        # Render the plates and quads, if present
        if self.model.Quads or self.model.Plates:
            self.plot_plates(self.deformed_shape, self.deformed_scale, self.color_map, self.scalar_bar, self.scalar_bar_text_size, self.combo_name, self.theme)

        # Set the window's background color
        # if self.theme == 'default':
        #     renderer.SetBackground(0, 0, 128)  # Blue
        # elif self.theme == 'print':
        #     renderer.SetBackground(255, 255, 255)  # White
            
        # Reset the camera if requested by the user
        if reset_camera:
            self.plotter.reset_camera()
    
    def plot_node(self, node, size, color='grey'):
        """Adds a node to the plotter

        :param node: node
        :type node: Node3D
        """
      
        # Get the node's position
        X = node.X # Global X coordinate
        Y = node.Y # Global Y coordinate
        Z = node.Z # Global Z coordinate

        # Generate any supports that occur at the node
        # Check for a fixed suppport
        if node.support_DX and node.support_DY and node.support_DZ and node.support_RX and node.support_RY and node.support_RZ:
            
            # Create a cube using PyVista
            self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y, node.Z),
                                          x_length=size*2,
                                          y_length=size*2,
                                          z_length=size*2))
        
        # Check for a pinned support
        elif node.support_DX and node.support_DY and node.support_DZ and not node.support_RX and not node.support_RY and not node.support_RZ:
            
            # Create a cone using PyVista's Cone function
            self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y-size, node.Z),
                                          direction=(0, 1, 0),
                                          height=size*2,
                                          radius=size*2))

        # Other support conditions
        else:

            # Generate a sphere for the node
            #sphere = pv.Sphere(center=(X, Y, Z), radius=size/2)
            #self.plotter.add_mesh(sphere, name='Node: '+node.name)
            
            # Restrained against X translation
            if node.support_DX:

                # Line showing support direction
                self.plotter.add_mesh((node.X-size, node.Y, node.Z),
                                      (node.X+size, node.Y, node.Z))

                # Cones at both ends
                self.plotter.add_mesh(pv.Cone(center=(node.X-size, node.Y, node.Z),
                                        direction=(1, 0, 0),
                                        height=size*0.6,
                                        radius=size*0.3))
                self.plotter.add_mesh(pv.Cone(center=(node.X+size, node.Y, node.Z),
                                        direction=(-1, 0, 0),
                                        height=size*0.6,
                                        radius=size*0.3))
            
            # Restrained against Y translation
            if node.support_DY:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y-size, node.Z),
                                              (node.X, node.Y+size, node.Z)))

                # Cones at both ends
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y-size, node.Z),
                                              direction=(0, 1, 0),
                                              height=size*0.6,
                                              radius=size*0.3))
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y+size, node.Z),
                                              direction=(0, -1, 0),
                                              height=size*0.6,
                                              radius=size*0.3))
            
            # Restrained against Z translation
            if node.support_DZ:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y, node.Z-size),
                                              (node.X, node.Y, node.Z+size)))

                # Cones at both ends
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y, node.Z-size),
                                              direction=(0, 0, 1),
                                              height=size*0.6,
                                              radius=size*0.3))
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y, node.Z+size),
                                              direction=(0, 0, -1),
                                              height=size*0.6,
                                              radius=size*0.3))
            
            # Restrained against X rotation
            if node.support_RX:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X-1.6*size, node.Y, node.Z),
                                              (node.X+1.6*size, node.Y, node.Z)))

                # Cubes at both ends
                self.plotter.add_mesh(pv.Cube(center=(node.X-1.9*size, node.Y, node.Z),
                                              x_length=size*0.6,
                                              y_length=size*0.6,
                                              z_length=size*0.6))
                self.plotter.add_mesh(pv.Cube(center=(node.X+1.9 *size, node.Y, node.Z),
                                              x_length=size*0.6,
                                              y_length=size*0.6,
                                              z_length=size*0.6))

            # Restrained against rotation about the Y-axis
            if node.support_RY:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y-1.6*size, node.Z),
                                              (node.X, node.Y+1.6*size, node.Z)))

                # Cubes at both ends
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y-1.9*size, node.Z),
                                              x_length=size*0.6,
                                              y_length=size*0.6,
                                              z_length=size*0.6))
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y+1.9*size, node.Z),
                                              x_length=size*0.6,
                                              y_length=size*0.6,
                                              z_length=size*0.6))
            
            # Restrained against rotation about the Z-axis
            if node.support_RZ:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y, node.Z-1.6*size),
                                              (node.X, node.Y, node.Z+1.6*size)))

                # Cubes at both ends
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y, node.Z-1.9*size),
                                              x_length=size*0.6,
                                              y_length=size*0.6,
                                              z_length=size*0.6))
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y, node.Z+1.9*size),
                                              x_length=size*0.6,
                                              y_length=size*0.6,
                                              z_length=size*0.6))

    def plot_member(self, member, theme='default'):
    
        # Generate a line for the member
        line = pv.Line()

        Xi = member.i_node.X
        Yi = member.i_node.Y
        Zi = member.i_node.Z
        line.points[0] = [Xi, Yi, Zi]

        Xj = member.j_node.X
        Yj = member.j_node.Y
        Zj = member.j_node.Z
        line.points[1] = [Xj, Yj, Zj]

        self.plotter.add_mesh(line, color='black', line_width=2)

    def plot_spring(self, spring, size, color='grey'):
        
        # Create points for the line source
        points = zeros((2, 3))  # 2 points in 3D
        line = pv.Line(points)

        # Step through each node in the model and find the position of the i-node and j-node
        i_node = spring.i_node
        j_node = spring.j_node
        Xi, Yi, Zi = i_node.X, i_node.Y, i_node.Z
        line.points[0] = [Xi, Yi, Zi]
        Xj, Yj, Zj = j_node.X, j_node.Y, j_node.Z
        line.points[1] = [Xj, Yj, Zj]

        # Create the PyVista actors
        actor = pv.PolyData(line).tube(radius=0.1)
        label = pv.VectorText(spring.name)

        # Set up an actor for the spring label
        lblActor = pv.Follower(label)
        lblActor.scale(size, size, size)
        lblActor.position = [(Xi+Xj)/2, (Yi+Yj)/2, (Zi+Zj)/2]

        # Add some color
        if color is None:
            actor.point_data["color"] = array([255, 0, 255])  # Magenta
            lblActor.point_data["color"] = array([255, 0, 255])
        elif color == 'black':
            actor.point_data["color"] = array([0, 0, 0])  # Black
            lblActor.point_data["color"] = array([0, 0, 0])

        # Plot the actors
        plotter = pv.Plotter()
        plotter.add_mesh(actor)
        plotter.add_mesh(lblActor)
        plotter.show()
            
    def plot_plates(self, deformed_shape, deformed_scale, color_map, scalar_bar, scalar_bar_text_size, combo_name, theme='default'):
        
        # Start a list of vertices
        plate_vertices = []

        # Start a list of plates (faces) for the mesh.
        plate_faces = []
        
        # `plate_results` will store the results in a list for PyVista
        plate_results = []
        
        # Each element will be assigned a unique element number `i` beginning at 0
        i = 0
        
        # Calculate the smoothed contour results at each node
        _PrepContour(self.model, color_map, combo_name)
        
        # Add each plate and quad in the model to the PyVista dataset
        for item in list(self.model.Plates.values()) + list(self.model.Quads.values()):
            
            # Create a point for each corner (must be in counter clockwise order)
            if deformed_shape:
                p0 = [item.i_node.X + item.i_node.DX[combo_name]*deformed_scale,
                    item.i_node.Y + item.i_node.DY[combo_name]*deformed_scale,
                    item.i_node.Z + item.i_node.DZ[combo_name]*deformed_scale]
                p1 = [item.j_node.X + item.j_node.DX[combo_name]*deformed_scale,
                    item.j_node.Y + item.j_node.DY[combo_name]*deformed_scale,
                    item.j_node.Z + item.j_node.DZ[combo_name]*deformed_scale]
                p2 = [item.m_node.X + item.m_node.DX[combo_name]*deformed_scale,
                    item.m_node.Y + item.m_node.DY[combo_name]*deformed_scale,
                    item.m_node.Z + item.m_node.DZ[combo_name]*deformed_scale]
                p3 = [item.n_node.X + item.n_node.DX[combo_name]*deformed_scale,
                    item.n_node.Y + item.n_node.DY[combo_name]*deformed_scale,
                    item.n_node.Z + item.n_node.DZ[combo_name]*deformed_scale]
            else:
                p0 = [item.i_node.X, item.i_node.Y, item.i_node.Z]
                p1 = [item.j_node.X, item.j_node.Y, item.j_node.Z]
                p2 = [item.m_node.X, item.m_node.Y, item.m_node.Z]
                p3 = [item.n_node.X, item.n_node.Y, item.n_node.Z]
        
            # Add the points to the PyVista dataset
            plate_vertices.append(p0)
            plate_vertices.append(p1)
            plate_vertices.append(p2)
            plate_vertices.append(p3)
            plate_faces.append([4, i*4, i*4 + 1, i*4 + 2, i*4 + 3])

            # Get the contour value for each node
            r0 = item.i_node.contour
            r1 = item.j_node.contour
            r2 = item.m_node.contour
            r3 = item.n_node.contour
            
            # Add plate results to the results list if the user has requested them
            if color_map:

                # Save the results for each corner of the plate - one entry for each corner
                plate_results.append(r0)
                plate_results.append(r1)
                plate_results.append(r2)
                plate_results.append(r3)

            # Move on to the next plate in our lists to repeat the process
            i+=1

        # Add the vertices and the faces to our lists
        plate_vertices = array(plate_vertices)
        plate_faces = array(plate_faces)

        # Create a new PyVista dataset to store plate data
        plate_polydata = pv.PolyData(plate_vertices, plate_faces)

        # Add the results as point data to the PyVista dataset
        if color_map:
            
            plate_polydata = plate_polydata.separate_cells()
            plate_polydata['Contours'] = array(plate_results)
            
            # Add the scalar bar for the contours
            if scalar_bar:
                self.plotter.add_mesh(plate_polydata, scalars='Contours', show_edges=True)
            else:
                self.plotter.add_mesh(plate_polydata)
        else:
            self.plotter.add_mesh(plate_polydata)
      
    def plot_deformed_node(self, node, scale_factor):

        # Calculate the node's deformed position
        newX = node.X + scale_factor * (node.DX[self.combo_name])
        newY = node.Y + scale_factor * (node.DY[self.combo_name])
        newZ = node.Z + scale_factor * (node.DZ[self.combo_name])

        # Generate a sphere source for the node in its deformed position
        sphere = pv.Sphere(radius=0.6 * self.annotation_size, center=[newX, newY, newZ])

        # Add the mesh to the plotter
        self.plotter.add_mesh(sphere)
  
    def plot_deformed_member(self, member, scale_factor):
        
        # Determine if this member is active for each load combination
        if member.active:
        
            L = member.L() # Member length
            T = member.T() # Member local transformation matrix
        
            cos_x = array([T[0, 0:3]]) # Direction cosines of local x-axis
            cos_y = array([T[1, 0:3]]) # Direction cosines of local y-axis
            cos_z = array([T[2, 0:3]]) # Direction cosines of local z-axis

            # Find the initial position of the local i-node
            Xi = member.i_node.X
            Yi = member.i_node.Y
            Zi = member.i_node.Z
        
            # Calculate the local y-axis displacements at 20 points along the member's length
            DY_plot = empty((0, 3))
            for i in range(20):
                    
                # Calculate the local y-direction displacement
                dy_tot = member.deflection('dy', L / 19 * i, self.combo_name)
            
                # Calculate the scaled displacement in global coordinates
                DY_plot = append(DY_plot, dy_tot * cos_y * scale_factor, axis=0)
        
            # Calculate the local z-axis displacements at 20 points along the member's length
            DZ_plot = empty((0, 3)) 
            for i in range(20):
                    
                # Calculate the local z-direction displacement
                dz_tot = member.deflection('dz', L / 19 * i, self.combo_name)
            
                # Calculate the scaled displacement in global coordinates
                DZ_plot = append(DZ_plot, dz_tot * cos_z * scale_factor, axis=0)
        
            # Calculate the local x-axis displacements at 20 points along the member's length
            DX_plot = empty((0, 3)) 
            for i in range(20):
                    
                # Displacements in local coordinates
                dx_tot = [[Xi, Yi, Zi]] + (L / 19 * i + member.deflection('dx', L / 19 * i, self.combo_name) * scale_factor) * cos_x
                    
                # Magnified displacements in global coordinates
                DX_plot = append(DX_plot, dx_tot, axis=0)
            
            # Sum the component displacements to obtain overall displacement
            D_plot = DY_plot + DZ_plot + DX_plot
            
            # Create lines connecting the points
            for i in range(len(D_plot)-1):
                line = pv.Line(D_plot[i], D_plot[i+1])
                self.plotter.add_mesh(line, color='red', line_width=2)
                                  
# class VisDeformedSpring():
    
#     def __init__(self, spring, nodes, scale_factor, combo_name='Combo 1'):

#         # Determine if this spring is active for each load combination
#         self.active = spring.active
        
#         # Generate a line source for the spring
#         self.source = vtk.vtkLineSource()
        
#         # Find the deformed position of the local i-node
#         # Step through each node
#         for node in nodes.values():
      
#             # Check to see if the current node is the i-node
#             if node.name == spring.i_node.name:
#                 Xi = node.X + node.DX[combo_name]*scale_factor
#                 Yi = node.Y + node.DY[combo_name]*scale_factor
#                 Zi = node.Z + node.DZ[combo_name]*scale_factor
#                 self.source.SetPoint1(Xi, Yi, Zi)
        
#             # Check to see if the current node is the i-node
#             if node.name == spring.j_node.name:
#                 Xj = node.X + node.DX[combo_name]*scale_factor
#                 Yj = node.Y + node.DY[combo_name]*scale_factor
#                 Zj = node.Z + node.DZ[combo_name]*scale_factor
#                 self.source.SetPoint2(Xj, Yj, Zj)
        
#         self.source.Update()
    
# class VisPtLoad():
#     '''
#     Creates a point load for the viewer
#     '''
    
#     def __init__(self, position, direction, length, label_text=None, annotation_size=5, color=None):
#         '''
#         Constructor.
      
#         Parameters
#         ----------
#         position : tuple
#           A tuple of X, Y and Z coordinates for the point of the load arrow: (X, Y, Z).
#         direction : tuple
#           A tuple indicating the direction vector for the load arrow: (i, j, k).
#         length : number
#           The length of the load arrow.
#         tip_length : number
#           The height of the arrow head.
#         label_text : string
#           Text that will show up at the tail of the arrow. If set to 'None' no text will be displayed.
#         '''
      
#         # Create a unit vector in the direction of the 'direction' vector
#         unitVector = direction/norm(direction)
      
#         # Create a 'vtkAppendPolyData' filter to append the tip and shaft together into a single dataset
#         self.polydata = vtk.vtkAppendPolyData()
      
#         # Determine if the load is positive or negative
#         if length == 0:
#             sign = 1
#         else:
#             sign = abs(length)/length
      
#         # Generate the tip of the load arrow
#         tip_length = abs(length)/4
#         radius = abs(length)/16
#         tip = vtk.vtkConeSource()
#         tip.SetCenter(position[0] - tip_length*sign*0.5*unitVector[0], \
#                       position[1] - tip_length*sign*0.5*unitVector[1], \
#                       position[2] - tip_length*sign*0.5*unitVector[2])
#         tip.SetDirection([direction[0]*sign, direction[1]*sign, direction[2]*sign])
#         tip.SetHeight(tip_length)
#         tip.SetRadius(radius)
#         tip.Update()
      
#         # Add the arrow tip to the append filter
#         self.polydata.AddInputData(tip.GetOutput())
        
#         # Create the shaft
#         shaft = vtk.vtkLineSource()
#         shaft.SetPoint1(position)
#         shaft.SetPoint2((position[0]-length*unitVector[0], position[1]-length*unitVector[1], position[2]-length*unitVector[2]))
#         shaft.Update()
      
#         # Copy and append the shaft data to the append filter
#         self.polydata.AddInputData(shaft.GetOutput())
#         self.polydata.Update()
      
#         # Create a mapper and actor
#         mapper = vtk.vtkPolyDataMapper()
#         mapper.SetInputConnection(self.polydata.GetOutputPort())
#         self.actor = vtk.vtkActor()
#         if color is None: self.actor.GetProperty().SetColor(0, 255, 0) # Green
#         elif color == 'black': self.actor.GetProperty().SetColor(0, 0, 0)  # Black
#         self.actor.SetMapper(mapper)
      
#         # Create the label if needed
#         if label_text != None:
      
#             # Create the label and set its text
#             self.label = vtk.vtkVectorText()
#             self.label.SetText(label_text)
        
#             # Set up a mapper for the label
#             lblMapper = vtk.vtkPolyDataMapper()
#             lblMapper.SetInputConnection(self.label.GetOutputPort())
        
#             # Set up an actor for the label
#             self.lblActor = vtk.vtkFollower()
#             self.lblActor.SetMapper(lblMapper)
#             self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
#             self.lblActor.SetPosition(position[0] - (length - 0.6*annotation_size)*unitVector[0], \
#                                       position[1] - (length - 0.6*annotation_size)*unitVector[1], \
#                                       position[2] - (length - 0.6*annotation_size)*unitVector[2])
#             if color is None: self.lblActor.GetProperty().SetColor(0, 255, 0)  # Green
#             elif color == 'black': self.lblActor.GetProperty().SetColor(0, 0, 0)  # Black
      
# class VisDistLoad():
#     '''
#     Creates a distributed load for the viewer
#     '''
    
#     def __init__(self, position1, position2, direction, length1, length2, label_text1, label_text2, annotation_size=5, color=None):
#         '''
#         Constructor.
#         '''
      
#         # Calculate the length of the distributed load
#         loadLength = ((position2[0]-position1[0])**2 + (position2[1]-position1[1])**2 + (position2[2]-position1[2])**2)**0.5
      
#         # Find the direction cosines for the line the load acts on
#         lineDirCos = [(position2[0]-position1[0])/loadLength, (position2[1]-position1[1])/loadLength, (position2[2]-position1[2])/loadLength]
      
#         # Find the direction cosines for the direction the load acts in
#         dirDirCos = direction/norm(direction)
      
#         # Create point loads at intervals roughly equal to 75% of the load's largest length (magnitude)
#         # Add text labels to the first and last load arrow
#         if loadLength > 0:
#             num_steps = int(round(0.75*loadLength/max(abs(length1), abs(length2)), 0))
#         else:
#             num_steps = 0

#         num_steps = max(num_steps, 1)
#         step = loadLength/num_steps
#         ptLoads = []
      
#         for i in range(num_steps + 1):
      
#             # Calculate the position (X, Y, Z) of this load arrow's point
#             position = (position1[0] + i*step*lineDirCos[0], position1[1] + i*step*lineDirCos[1], position1[2] + i*step*lineDirCos[2])
        
#             # Determine the length of this load arrow
#             length = length1 + (length2 - length1)/loadLength*i*step
        
#             # Determine the label's text
#             if i == 0:
#                 label_text = label_text1
#             elif i == num_steps:
#                 label_text = label_text2
        
#             # Create the load arrow
#             ptLoads.append(VisPtLoad(position, direction, length, label_text, annotation_size=annotation_size))
          
#         # Draw a line between the first and last load arrow's tails
#         tail_line = vtk.vtkLineSource()
#         tail_line.SetPoint1((position1[0] - length1*dirDirCos[0], position1[1] - length1*dirDirCos[1], position1[2] - length1*dirDirCos[2]))
#         tail_line.SetPoint2((position2[0] - length2*dirDirCos[0], position2[1] - length2*dirDirCos[1], position2[2] - length2*dirDirCos[2]))
      
#         # Combine all the geometry into one 'vtkPolyData' object
#         self.polydata = vtk.vtkAppendPolyData()
#         for arrow in ptLoads:
#             arrow.polydata.Update()
#             self.polydata.AddInputData(arrow.polydata.GetOutput())
        
#         tail_line.Update()
#         self.polydata.AddInputData(tail_line.GetOutput())
#         self.polydata.Update()
      
#         # Create a mapper and actor for the geometry
#         mapper = vtk.vtkPolyDataMapper()
#         mapper.SetInputConnection(self.polydata.GetOutputPort())
#         self.actor = vtk.vtkActor()
#         if color is None: self.actor.GetProperty().SetColor(0, 255, 0)  # Green
#         elif color == 'black': self.actor.GetProperty().SetColor(0, 0, 0)  # Black
#         self.actor.SetMapper(mapper)
      
#         # Get the actors for the labels
#         self.lblActors = [ptLoads[0].lblActor, ptLoads[len(ptLoads) - 1].lblActor]

# class VisMoment():
#     '''
#     Creates a concentrated moment for the viewer
#     '''
    
#     def __init__(self, center, direction, radius, label_text=None, annotation_size=5, color=None):
#         '''
#         Constructor.
      
#         Parameters
#         ----------
#         center : tuple
#           A tuple of X, Y and Z coordinates for center of the moment: (X, Y, Z).
#         direction : tuple
#           A tuple indicating the direction vector for the moment: (i, j, k).
#         radius : number
#           The radius of the moment.
#         tip_length : number
#           The height of the arrow head.
#         label_text : string
#           Text that will show up at the tail of the moment. If set to 'None' no text will be displayed.
#         '''
      
#         # Create an append filter to store load polydata in
#         self.polydata = vtk.vtkAppendPolyData()
        
#         # Find a vector perpendicular to the directional unit vector
#         v1 = direction/norm(direction)  # v1 = The directional unit vector for the moment
#         v2 = _PerpVector(v1)             # v2 = A unit vector perpendicular to v1
#         v3 = cross(v1, v2)
#         v3 = v3/norm(v3)                # v3 = A unit vector perpendicular to v1 and v2
      
#         # Generate an arc for the moment
#         Xc, Yc, Zc = center
#         arc = vtk.vtkArcSource()
#         arc.SetCenter(Xc, Yc, Zc)
#         arc.SetPoint1(Xc + v2[0]*radius, Yc + v2[1]*radius, Zc + v2[2]*radius)
#         arc.SetPoint2(Xc + v3[0]*radius, Yc + v3[1]*radius, Zc + v3[2]*radius)
#         arc.SetNegative(True)
#         arc.SetResolution(20)
#         arc.Update()
#         self.polydata.AddInputData(arc.GetOutput())
      
#         # Generate the arrow tip at the end of the arc
#         tip_length = radius/2
#         cone_radius = radius/8
#         tip = vtk.vtkConeSource()
#         tip.SetCenter(arc.GetPoint1()[0], arc.GetPoint1()[1], arc.GetPoint1()[2])
#         tip.SetDirection(cross(v1, v2))
#         tip.SetHeight(tip_length)
#         tip.SetRadius(cone_radius)
#         tip.Update()
#         self.polydata.AddInputData(tip.GetOutput())
      
#         # Update the polydata one last time now that we're done appending items to it
#         self.polydata.Update()
      
#         # Create the text label
#         label = vtk.vtkVectorText()
#         label.SetText(label_text)
#         lblMapper = vtk.vtkPolyDataMapper()
#         lblMapper.SetInputConnection(label.GetOutputPort())
#         self.lblActor = vtk.vtkFollower()
#         self.lblActor.SetMapper(lblMapper)
#         self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
#         self.lblActor.SetPosition(Xc + v3[0]*(radius + 0.25*annotation_size), \
#                                   Yc + v3[1]*(radius + 0.25*annotation_size), \
#                                   Zc + v3[2]*(radius + 0.25*annotation_size))
#         if color is None: self.lblActor.GetProperty().SetColor(0, 255, 0)  # Green
#         elif color == 'black': self.lblActor.GetProperty().SetColor(0, 0, 0)  # Black

# class VisAreaLoad():
#     '''
#     Creates an area load for the viewer
#     '''
    
#     def __init__(self, position0, position1, position2, position3, direction, length, label_text, annotation_size=5, theme='default'):
#         '''
#         Constructor
#         '''
      
#         # Create a point load for each corner of the area load
#         ptLoads = []
#         ptLoads.append(VisPtLoad(position0, direction, length, label_text, annotation_size=annotation_size))
#         ptLoads.append(VisPtLoad(position1, direction, length, label_text, annotation_size=annotation_size))
#         ptLoads.append(VisPtLoad(position2, direction, length, label_text, annotation_size=annotation_size))
#         ptLoads.append(VisPtLoad(position3, direction, length, label_text, annotation_size=annotation_size))
      
#         # Find the direction cosines for the direction the load acts in
#         dirDirCos = direction/norm(direction)
      
#         # Find the positions of the tails of all the arrows at the corners of the area load. This is
#         # where we will place the polygon.
#         self.p0 = position0 - dirDirCos*length
#         self.p1 = position1 - dirDirCos*length
#         self.p2 = position2 - dirDirCos*length
#         self.p3 = position3 - dirDirCos*length
      
#         # Combine all geometry into one 'vtkPolyData' object
#         self.polydata = vtk.vtkAppendPolyData()
#         for arrow in ptLoads:
#             self.polydata.AddInputData(arrow.polydata.GetOutput())
#         self.polydata.Update()
      
#         # Add a label
#         self.label_actor = ptLoads[0].lblActor

#         # Add color to the area load label
#         if theme == 'print':
#             self.label_actor.GetProperty().SetColor(255/255, 0/255, 0/255)  # red
#         elif theme == 'default':
#             self.label_actor.GetProperty().SetColor(0/255, 255/255, 0/255)  # green

def _PerpVector(v):
    '''
    Returns a unit vector perpendicular to v=[i, j, k]
    '''
    
    i = v[0]
    j = v[1]
    k = v[2]
    
    # Find a vector in a direction perpendicular to <i, j, k>
    if i == 0:
        i2 = 1
        j2 = 0
        k2 = 0
    elif j == 0:
        i2 = 0
        j2 = 1
        k2 = 0
    elif k == 0:
        i2 = 0
        j2 = 0
        k2 = 1
    else:
        i2 = 1
        j2 = 1
        k2 = -(i*i2+j*j2)/k
    
    # Return the unit vector
    return [i2, j2, k2]/norm([i2, j2, k2])

def _PrepContour(model, stress_type='Mx', combo_name='Combo 1'):

    if stress_type != None:
    
        # Erase any previous contours
        for node in model.Nodes.values():
            node.contour = []
    
        # Step through each element in the model
        for element in list(model.Quads.values()) + list(model.Plates.values()):
            
            # Rectangular elements and quadrilateral elements have different local coordinate systems. Rectangles are based on a traditional (x, y) system, while quadrilaterals are based on a 'natural' (r, s) coordinate system. To reduce duplication of code for both these elements we'll define the edges of the plate here for either element using the (r, s) terminology.
            if element.type == 'Rect':
                r_left = 0
                r_right = element.width()
                s_bot = 0
                s_top = element.height()
            else:
                r_left = -1
                r_right = 1
                s_bot = -1
                s_top = 1
      
            # Determine which stress result has been requested by the user
            if stress_type == 'dz':
                # Internally PyNite defines the nodes for a rectangular element in the order (i, j,m, n), while VTK defines the nodes for a quadrilateral element in the order (m, n, i, j)
                if element.type == 'Rect':
                    i, j, m, n = element.d(combo_name)[[2, 8, 14, 20], :]
                else:
                    i, j, m, n = element.d(combo_name)[[14, 20, 2, 8], :]
                element.i_node.contour.append(i)
                element.j_node.contour.append(j)
                element.m_node.contour.append(m)
                element.n_node.contour.append(n)
            elif stress_type == 'Mx':
                element.i_node.contour.append(element.moment(r_left, s_bot, combo_name)[0])
                element.j_node.contour.append(element.moment(r_right, s_bot, combo_name)[0])
                element.m_node.contour.append(element.moment(r_right, s_top, combo_name)[0])
                element.n_node.contour.append(element.moment(r_left, s_top, combo_name)[0])
            elif stress_type == 'My':
                element.i_node.contour.append(element.moment(r_left, s_bot, combo_name)[1])
                element.j_node.contour.append(element.moment(r_right, s_bot, combo_name)[1])
                element.m_node.contour.append(element.moment(r_right, s_top, combo_name)[1])
                element.n_node.contour.append(element.moment(r_left, s_top, combo_name)[1])
            elif stress_type == 'Mxy':
                element.i_node.contour.append(element.moment(r_left, s_bot, combo_name)[2])
                element.j_node.contour.append(element.moment(r_right, s_bot, combo_name)[2])
                element.m_node.contour.append(element.moment(r_right, s_top, combo_name)[2])
                element.n_node.contour.append(element.moment(r_left, s_top, combo_name)[2])
            elif stress_type == 'Qx':
                element.i_node.contour.append(element.shear(r_left, s_bot, combo_name)[0])
                element.j_node.contour.append(element.shear(r_right, s_bot, combo_name)[0])
                element.m_node.contour.append(element.shear(r_right, s_top, combo_name)[0])
                element.n_node.contour.append(element.shear(r_left, s_top, combo_name)[0])
            elif stress_type == 'Qy':
                element.i_node.contour.append(element.shear(r_left, s_bot, combo_name)[1])
                element.j_node.contour.append(element.shear(r_right, s_bot, combo_name)[1])
                element.m_node.contour.append(element.shear(r_right, s_top, combo_name)[1])
                element.n_node.contour.append(element.shear(r_left, s_top, combo_name)[1])
            elif stress_type == 'Sx':
                element.i_node.contour.append(element.membrane(r_left, s_bot, combo_name)[0])
                element.j_node.contour.append(element.membrane(r_right, s_bot, combo_name)[0])
                element.m_node.contour.append(element.membrane(r_right, s_top, combo_name)[0])
                element.n_node.contour.append(element.membrane(r_left, s_top, combo_name)[0])
            elif stress_type == 'Sy':
                element.i_node.contour.append(element.membrane(r_left, s_bot, combo_name)[1])
                element.j_node.contour.append(element.membrane(r_right, s_bot, combo_name)[1])
                element.m_node.contour.append(element.membrane(r_right, s_top, combo_name)[1])
                element.n_node.contour.append(element.membrane(r_left, s_top, combo_name)[1])
            elif stress_type == 'Txy':
                element.i_node.contour.append(element.membrane(r_left, s_bot, combo_name)[2])
                element.j_node.contour.append(element.membrane(r_right, s_bot, combo_name)[2])
                element.m_node.contour.append(element.membrane(r_right, s_top, combo_name)[2])
                element.n_node.contour.append(element.membrane(r_left, s_top, combo_name)[2])                 
      
        # Average the values at each node to obtain a smoothed contour
        for node in model.Nodes.values():
            # Prevent divide by zero errors for nodes with no contour values
            if node.contour != []:
                node.contour = sum(node.contour)/len(node.contour)

# def _DeformedShape(model, scale_factor, annotation_size, combo_name, render_nodes=True, theme='default'):
#     '''
#     Renders the deformed shape of a model.
    
#     Parameters
#     ----------
#     model : FEModel3D
#         Finite element model to be rendered.
#     scale_factor : number
#         The scale factor to apply to the model deformations.
#     annotation_size : number
#         Controls the height of text displayed with the model. The units used for `annotation_size` are
#         the same as those used for lengths in the model. Sizes of other objects (such as nodes) are
#         related to this value.
#     combo_name : string
#         The load case used for rendering the deflected shape.
    
#     Returns
#     -------
#     None.
#     '''
    
#     # Create an empty PolyData
#     polydata = pv.PolyData()
        
#     # Add the springs to the PolyData
#     for spring in model.Springs.values():
#         # Only add the spring if it is active for the given load combination
#         if spring.active[combo_name]:
#             vis_spring = VisDeformedSpring(spring, model.Nodes, scale_factor, combo_name)
#             polydata += vis_spring.source.GetOutput()
            
#     # Adjust the color
#     if theme == 'default':
#         color = [255, 255, 0]  # Yellow
#     elif theme == 'print':
#         color = [26, 26, 26]  # Dark Grey
#     else:
#         color = [255, 255, 255]  # White
    
#     # Create a PyVista actor for the PolyData
#     actor = pv.PolyData(polydata).plot(color=color)


# def _RenderLoads(model, renderer, annotation_size, combo_name, case, theme='default'):

#     # Create an append filter to store all the polydata in. This will allow us to use fewer actors to
#     # display all the loads, which will greatly improve rendering speed as the user interacts. VTK
#     # becomes very slow when a large number of actors are used.
#     polydata = vtk.vtkAppendPolyData()
    
#     # Polygons are treated as cells in VTK. Create a cell array to store all the area load polygons
#     # in. We'll also create a list of points to store the polygon points in. The polydata for these
#     # polygons will be stored separately from the other load data.
#     polygons = vtk.vtkCellArray()
#     polygon_points = vtk.vtkPoints()
#     polygon_polydata = vtk.vtkPolyData()
    
#     # Get the maximum load magnitudes that will be used to normalize the display scale
#     max_pt_load, max_moment, max_dist_load, max_area_load = _MaxLoads(model, combo_name, case)
    
#     # Display the requested load combination, or 'Combo 1' if no load combo or case has been
#     # specified
#     if case == None:
#         # Store model.LoadCombos[combo].factors under a simpler name for use below
#         load_factors = model.LoadCombos[combo_name].factors
#     else:
#         # Set up a load combination dictionary that represents the load case
#         load_factors = {case: 1}
    
#     # Step through each node
#     for node in model.Nodes.values():
    
#         # Step through and display each nodal load
#         for load in node.NodeLoads:
          
#             # Determine if this load is part of the requested LoadCombo or case
#             if load[2] in load_factors:
              
#                 # Calculate the factored value for this load and it's sign (positive or negative)
#                 load_value = load[1]*load_factors[load[2]]
#                 if load_value != 0:
#                     sign = load_value/abs(load_value)
#                 else:
#                     sign = 1
                
#                 # Display the load
#                 if load[0] == 'FX':
#                     ptLoad = VisPtLoad((node.X - 0.6*annotation_size*sign, node.Y, node.Z), [1, 0, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'FY':
#                     ptLoad = VisPtLoad((node.X, node.Y - 0.6*annotation_size*sign, node.Z), [0, 1, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'FZ':
#                     ptLoad = VisPtLoad((node.X, node.Y, node.Z - 0.6*annotation_size*sign), [0, 0, 1], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'MX':
#                     ptLoad = VisMoment((node.X, node.Y, node.Z), (1*sign, 0, 0), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'MY':
#                     ptLoad = VisMoment((node.X, node.Y, node.Z), (0, 1*sign, 0), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'MZ':
#                     ptLoad = VisMoment((node.X, node.Y, node.Z), (0, 0, 1*sign), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
                
#                 polydata.AddInputData(ptLoad.polydata.GetOutput())
#                 renderer.AddActor(ptLoad.lblActor)
#                 ptLoad.lblActor.SetCamera(renderer.GetActiveCamera())
    
#     # Step through each member
#     for member in model.Members.values():
    
#         # Get the direction cosines for the member's local axes
#         dir_cos = member.T()[0:3, 0:3]
      
#         # Get the starting point for the member
#         x_start, y_start, z_start = member.i_node.X, member.i_node.Y, member.i_node.Z
      
#         # Step through each member point load
#         for load in member.PtLoads:
      
#             # Determine if this load is part of the requested load combination
#             if load[3] in load_factors:
        
#                 # Calculate the factored value for this load and it's sign (positive or negative)
#                 load_value = load[1]*load_factors[load[3]]
#                 sign = load_value/abs(load_value)
          
#                 # Calculate the load's location in 3D space
#                 x = load[2]
#                 position = [x_start + dir_cos[0, 0]*x, y_start + dir_cos[0, 1]*x, z_start + dir_cos[0, 2]*x]
          
#                 # Display the load
#                 if load[0] == 'Fx':
#                     ptLoad = VisPtLoad(position, dir_cos[0, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'Fy':
#                     ptLoad = VisPtLoad(position, dir_cos[1, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'Fz':
#                     ptLoad = VisPtLoad(position, dir_cos[2, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'Mx':
#                     ptLoad = VisMoment(position, dir_cos[0, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'My':
#                     ptLoad = VisMoment(position, dir_cos[1, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'Mz':
#                     ptLoad = VisMoment(position, dir_cos[2, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'FX':
#                     ptLoad = VisPtLoad(position, [1, 0, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'FY':
#                     ptLoad = VisPtLoad(position, [0, 1, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'FZ':
#                     ptLoad = VisPtLoad(position, [0, 0, 1], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'MX':
#                     ptLoad = VisMoment(position, [1*sign, 0, 0], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'MY':
#                     ptLoad = VisMoment(position, [0, 1*sign, 0], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
#                 elif load[0] == 'MZ':
#                     ptLoad = VisMoment(position, [0, 0, 1*sign], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size)
            
#                 polydata.AddInputData(ptLoad.polydata.GetOutput())
#                 renderer.AddActor(ptLoad.lblActor)
#                 ptLoad.lblActor.SetCamera(renderer.GetActiveCamera())
    
#         # Step through each member distributed load
#         for load in member.DistLoads:
    
#             # Determine if this load is part of the requested load combination
#             if load[5] in load_factors:
        
#                 # Calculate the factored value for this load and it's sign (positive or negative)
#                 w1 = load[1]*load_factors[load[5]]
#                 w2 = load[2]*load_factors[load[5]]
          
#                 # Calculate the loads location in 3D space
#                 x1 = load[3]
#                 x2 = load[4]
#                 position1 = [x_start + dir_cos[0, 0]*x1, y_start + dir_cos[0, 1]*x1, z_start + dir_cos[0, 2]*x1]
#                 position2 = [x_start + dir_cos[0, 0]*x2, y_start + dir_cos[0, 1]*x2, z_start + dir_cos[0, 2]*x2]
                
#                 # Display the load
#                 if load[0] == 'Fx':
#                     distLoad = VisDistLoad(position1, position2, dir_cos[0, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
#                 elif load[0] == 'Fy':
#                     distLoad = VisDistLoad(position1, position2, dir_cos[1, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
#                 elif load[0] == 'Fz':
#                     distLoad = VisDistLoad(position1, position2, dir_cos[2, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
#                 elif load[0] == 'FX':
#                     distLoad = VisDistLoad(position1, position2, [1, 0, 0], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
#                 elif load[0] == 'FY':
#                     distLoad = VisDistLoad(position1, position2, [0, 1, 0], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
#                 elif load[0] == 'FZ':
#                     distLoad = VisDistLoad(position1, position2, [0, 0, 1], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
               
#                 polydata.AddInputData(distLoad.polydata.GetOutput())
#                 renderer.AddActor(distLoad.lblActors[0])
#                 renderer.AddActor(distLoad.lblActors[1])
#                 distLoad.lblActors[0].SetCamera(renderer.GetActiveCamera())
#                 distLoad.lblActors[1].SetCamera(renderer.GetActiveCamera())   
    
#     # Step through each plate
#     i = 0
#     for plate in list(model.Plates.values()) + list(model.Quads.values()):
    
#         # Get the direction cosines for the plate's local z-axis
#         dir_cos = plate.T()[0:3, 0:3]
#         dir_cos = dir_cos[2]
      
#         # Step through each plate load
#         for load in plate.pressures:
      
#             # Determine if this load is part of the requested load combination
#             if load[1] in load_factors:
        
#                 # Calculate the factored value for this load
#                 load_value = load[0]*load_factors[load[1]]
                
#                 # Find the sign for this load. Intercept any divide by zero errors
#                 if load[0] == 0:
#                     sign = 1
#                 else:
#                     sign = abs(load[0])/load[0]
          
#                 # Find the position of the load's 4 corners
#                 position0 = [plate.i_node.X, plate.i_node.Y, plate.i_node.Z]
#                 position1 = [plate.j_node.X, plate.j_node.Y, plate.j_node.Z]
#                 position2 = [plate.m_node.X, plate.m_node.Y, plate.m_node.Z]
#                 position3 = [plate.n_node.X, plate.n_node.Y, plate.n_node.Z]
          
#                 # Create an area load and get its data
#                 area_load = VisAreaLoad(position0, position1, position2, position3, dir_cos*sign, abs(load_value)/max_area_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
          
#                 # Add the area load's arrows to the overall load polydata
#                 polydata.AddInputData(area_load.polydata.GetOutput())
          
#                 # Add the 4 points at the corners of this area load to the list of points
#                 polygon_points.InsertNextPoint(area_load.p0[0], area_load.p0[1], area_load.p0[2])
#                 polygon_points.InsertNextPoint(area_load.p1[0], area_load.p1[1], area_load.p1[2])
#                 polygon_points.InsertNextPoint(area_load.p2[0], area_load.p2[1], area_load.p2[2])
#                 polygon_points.InsertNextPoint(area_load.p3[0], area_load.p3[1], area_load.p3[2])
          
#                 # Create a polygon based on the four points we just defined.
#                 # The 1st number in `SetId()` is the local point id
#                 # The 2nd number in `SetId()` is the global point id
#                 polygon = vtk.vtkPolygon()
#                 polygon.GetPointIds().SetNumberOfIds(4)
#                 polygon.GetPointIds().SetId(0, i*4)
#                 polygon.GetPointIds().SetId(1, i*4 + 1)
#                 polygon.GetPointIds().SetId(2, i*4 + 2)
#                 polygon.GetPointIds().SetId(3, i*4 + 3)
          
#                 # Add the polygon to the list of polygons
#                 polygons.InsertNextCell(polygon)
                
#                 # Add the load label
#                 renderer.AddActor(area_load.label_actor)
          
#                 # Set the text to follow the camera as the user interacts
#                 area_load.label_actor.SetCamera(renderer.GetActiveCamera())
          
#                 # `i` keeps track of the next polygon's ID. We've just added a polygon, so `i` needs to
#                 # go up 1.
#                 i += 1
        
#                 # Create polygon polydata from all the points and polygons we just defined
#                 polygon_polydata.SetPoints(polygon_points)
#                 polygon_polydata.SetPolys(polygons)
    
#     # Set up an actor and mapper for the loads
#     load_mapper = vtk.vtkPolyDataMapper()
#     load_mapper.SetInputConnection(polydata.GetOutputPort())
#     load_actor = vtk.vtkActor()
#     load_actor.SetMapper(load_mapper)

#     # Colorize the loads
#     if theme == 'default':
#         load_actor.GetProperty().SetColor(0/255, 255/255, 0/255)  # Green
#     elif theme == 'print':
#         load_actor.GetProperty().SetColor(255/255, 0/255, 0/255)  # Red

#     # Add the load actor to the renderer
#     renderer.AddActor(load_actor)
    
#     # Set up an actor and a mapper for the area load polygons
#     polygon_mapper = vtk.vtkPolyDataMapper()
#     polygon_mapper.SetInputData(polygon_polydata)
#     polygon_actor = vtk.vtkActor()

#     # polygon_actor.GetProperty().SetOpacity(0.5)      # 50% opacity
#     polygon_actor.SetMapper(polygon_mapper)
#     renderer.AddActor(polygon_actor)
    
#     # Set the color of the area load polygons
#     if theme == 'default':
#         polygon_actor.GetProperty().SetColor(0/255, 255/255, 0/255)  # Green
#     elif theme == 'print':
#         polygon_actor.GetProperty().SetColor(255/255, 0/255, 0/255)  # Red

def _MaxLoads(model, combo_name=None, case=None):

    max_pt_load = 0
    max_moment = 0
    max_dist_load = 0
    max_area_load = 0
    
    # Find the requested load combination or load case
    if case == None:
    
        # Step through each node
        for node in model.Nodes.values():
      
            # Step through each nodal load to find the largest one
            for load in node.NodeLoads:
              
                # Find the largest loads in the load combination
                if load[2] in model.LoadCombos[combo_name].factors:
                    if load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ':
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[2]]) > max_pt_load:
                            max_pt_load = abs(load[1]*model.LoadCombos[combo_name].factors[load[2]])
                    else:
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[2]]) > max_moment:
                            max_moment = abs(load[1]*model.LoadCombos[combo_name].factors[load[2]])
    
        # Step through each member
        for member in model.Members.values():
      
            # Step through each member point load
            for load in member.PtLoads:
                
                # Find and store the largest point load and moment in the load combination
                if load[3] in model.LoadCombos[combo_name].factors:
          
                    if (load[0] == 'Fx' or load[0] == 'Fy' or load[0] == 'Fz'
                    or  load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ'):
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[3]]) > max_pt_load:
                            max_pt_load = abs(load[1]*model.LoadCombos[combo_name].factors[load[3]])
                    else:
                        if abs(load[1]*model.LoadCombos[combo_name].factors[load[3]]) > max_moment:
                            max_moment = abs(load[1]*model.LoadCombos[combo_name].factors[load[3]])
      
            # Step through each member distributed load
            for load in member.DistLoads:
        
                #Find and store the largest distributed load in the load combination
                if load[5] in model.LoadCombos[combo_name].factors:
          
                    if abs(load[1]*model.LoadCombos[combo_name].factors[load[5]]) > max_dist_load:
                        max_dist_load = abs(load[1]*model.LoadCombos[combo_name].factors[load[5]])
                    if abs(load[2]*model.LoadCombos[combo_name].factors[load[5]]) > max_dist_load:
                        max_dist_load = abs(load[2]*model.LoadCombos[combo_name].factors[load[5]])
      
        # Step through each plate
        for plate in model.Plates.values():
      
            # Step through each plate load
            for load in plate.pressures:
        
                if load[1] in model.LoadCombos[combo_name].factors:
                    if abs(load[0]*model.LoadCombos[combo_name].factors[load[1]]) > max_area_load:
                        max_area_load = abs(load[0]*model.LoadCombos[combo_name].factors[load[1]])
      
        # Step through each quad
        for quad in model.Quads.values():
      
            # Step through each plate load
            for load in quad.pressures:
        
                # Check to see if the load case is in the requested load combination
                if load[1] in model.LoadCombos[combo_name].factors:
                    if abs(load[0]*model.LoadCombos[combo_name].factors[load[1]]) > max_area_load:
                        max_area_load = abs(load[0]*model.LoadCombos[combo_name].factors[load[1]])
      
    # Behavior if case has been specified
    else:
        
        # Step through each node
        for node in model.Nodes.values():
      
            # Step through each nodal load to find the largest one
            for load in node.NodeLoads:
              
                # Find the largest loads in the load case
                if load[2] == case:
                    if load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ':
                        if abs(load[1]) > max_pt_load:
                            max_pt_load = abs(load[1])
                    else:
                        if abs(load[1]) > max_moment:
                            max_moment = abs(load[1])
      
        # Step through each member
        for member in model.Members.values():
      
            # Step through each member point load
            for load in member.PtLoads:
              
                # Find and store the largest point load and moment in the load case
                if load[3] == case:
            
                    if (load[0] == 'Fx' or load[0] == 'Fy' or load[0] == 'Fz'
                    or  load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ'):
                        if abs(load[1]) > max_pt_load:
                            max_pt_load = abs(load[1])
                    else:
                        if abs(load[1]) > max_moment:
                            max_moment = abs(load[1])
        
            # Step through each member distributed load
            for load in member.DistLoads:
          
                # Find and store the largest distributed load in the load case
                if load[5] == case:
            
                    if abs(load[1]) > max_dist_load:
                        max_dist_load = abs(load[1])
                    if abs(load[2]) > max_dist_load:
                        max_dist_load = abs(load[2])
            
            # Step through each plate
            for plate in model.Plates.values():
          
                # Step through each plate load
                for load in plate.pressures:
          
                    if load[1] == case:
            
                        if abs(load[0]) > max_area_load:
                            max_area_load = abs(load[0])
          
        # Step through each quad
        for quad in model.Quads.values():
      
            # Step through each plate load
            for load in quad.pressures:
        
                if load[1] == case:
          
                    if abs(load[0]) > max_area_load:
                        max_area_load = abs(load[0])
        
    # Return the maximum loads in the load combination or load case
    return max_pt_load, max_moment, max_dist_load, max_area_load
