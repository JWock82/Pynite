from __future__ import annotations # Allows more recent type hints features
from json import load
import warnings
from typing import TYPE_CHECKING, Callable, List, Any

from IPython.display import Image
import numpy as np
import pyvista as pv
import math

# Suppress PyVista warnings for clean console output
warnings.filterwarnings('ignore', category=UserWarning, module='pyvista')
warnings.filterwarnings('ignore', message='.*Points is not a float type.*')

# For type checking only - these imports are only used during type checking
if TYPE_CHECKING:
    from typing import List, Union, Tuple, Optional
    from Pynite.Node3D import Node3D
    from Pynite.Member3D import Member3D
    from Pynite.Spring3D import Spring3D
    from Pynite.FEModel3D import FEModel3D

# Allow for 3D interaction within jupyter notebook using trame
try:
    pv.global_theme.trame.jupyter_extension_enabled = True
except:
    # Ignore the exception that is produced if we are not running the code via jupyter
    pass
pv.set_jupyter_backend('trame')

class Renderer:
    """Used to render finite element models.
    """

    scalar: Optional[str] = None

    def __init__(self, model: FEModel3D) -> None:

        self.model: FEModel3D = model

        # Default settings for rendering
        self._annotation_size: float = 5.0
        self._deformed_shape: bool = False
        self._deformed_scale: float = 30.0
        self._render_nodes: bool = True
        self._render_loads: bool = True
        self._color_map: Optional[str] = None
        self._combo_name: Optional[str] = 'Combo 1'
        self._case: Optional[str] = None
        self._labels: bool = True
        self._scalar_bar: bool = False
        self._scalar_bar_text_size: int = 24
        self.theme: str = 'default'

        # Callback list for post-update customization:
        # This is added because `self.update()` clears the plotter, removing user self.plotter configurations.
        # Functions in this list run after Pynite adds actors, allowing further PyVista customizations 
        # (e.g., grid, axes) before render. Each func in this list must accept a `pyvista.Plotter` argument.
        self.post_update_callbacks: List[Callable[[pv.Plotter], None]] = []

        # Create plotter with off_screen mode if pyvista is set globally to OFF_SCREEN
        # This is important for headless CI/testing environments
        self.plotter: pv.Plotter = pv.Plotter(off_screen=pv.OFF_SCREEN)
        self.plotter.set_background('white')  # Setting background color
        # self.plotter.add_logo_widget('./Resources/Full Logo No Buffer - Transparent.png')
        
        # Only set view and axes in interactive mode (renderer must exist for these calls)
        # In off-screen/headless mode, these will be set when update() is called
        if not pv.OFF_SCREEN:
            # self.plotter.view_isometric()
            self.plotter.view_xy()
            self.plotter.show_axes()
            self.plotter.set_viewup((0, 1, 0))  # Set the Y axis to vertical for 3D plots
        
        # Make X button behave like 'q' key - properly exit without destroying plotter
        # Why: By default, PyVista's X button forcefully destroys the render window,
        #      causing warnings and preventing clean shutdown. This observer intercepts
        #      the window close event.
        # How: When ExitEvent fires (X button clicked), we call TerminateApp() on the
        #      interactor, which cleanly exits the event loop just like pressing 'q'.
        self.plotter.iren.add_observer('ExitEvent', lambda obj, event: obj.TerminateApp())

        # Initialize load labels
        self._load_label_points: List[List[float]] = []
        self._load_labels: List[Union[str, float, int]] = []

        # Initialize spring labels
        self._spring_label_points: List[List[float]] = []
        self._spring_labels: List[str] = []

    @property
    def window_width(self) -> int:
        return self.plotter.window_size[0]

    @window_width.setter
    def window_width(self, width: int) -> None:
        height = self.plotter.window_size[1]
        self.plotter.window_size = (width, height)

    @property
    def window_height(self) -> int:
        return self.plotter.window_size[1]

    @window_height.setter
    def window_height(self, height: int) -> None:
        width = self.plotter.window_size[0]
        self.plotter.window_size = (width, height)

    @property
    def annotation_size(self) -> float:
        return self._annotation_size
    
    @annotation_size.setter
    def annotation_size(self, size: float) -> None:
        self._annotation_size = size
    
    @property
    def deformed_shape(self) -> bool:
        return self._deformed_shape
    
    @deformed_shape.setter
    def deformed_shape(self, deformed_shape: bool) -> None:
        self._deformed_shape = deformed_shape
    
    @property
    def deformed_scale(self) -> float:
        return self._deformed_scale
    
    @deformed_scale.setter
    def deformed_scale(self, scale: float) -> None:
        self._deformed_scale = scale
    
    @property
    def render_nodes(self) -> bool:
        return self._render_nodes
    
    @render_nodes.setter
    def render_nodes(self, render_nodes: bool) -> None:
        self._render_nodes = render_nodes

    @property
    def render_loads(self) -> bool:
        return self._render_loads
    
    @render_loads.setter
    def render_loads(self, render_loads: bool) -> None:
        self._render_loads = render_loads
    
    @property
    def color_map(self) -> Optional[str]:
        return self._color_map
    
    @color_map.setter
    def color_map(self, color_map: Optional[str]) -> None:
        self._color_map = color_map
    
    @property
    def combo_name(self) -> Optional[str]:
        return self._combo_name
    
    @combo_name.setter
    def combo_name(self, combo_name: Optional[str]) -> None:
        self._combo_name = combo_name
        self._case = None

    @property
    def case(self) -> Optional[str]:
        return self._case

    @case.setter
    def case(self, case: Optional[str]) -> None:
        self._case = case
        self._combo_name = None

    @property
    def show_labels(self) -> bool:
        return self._labels

    @show_labels.setter
    def show_labels(self, show_labels: bool) -> None:
        self._labels = show_labels

    @property
    def scalar_bar(self) -> bool:
        return self._scalar_bar

    @scalar_bar.setter
    def scalar_bar(self, scalar_bar: bool) -> None:
        self._scalar_bar = scalar_bar

    @property
    def scalar_bar_text_size(self) -> int:
        return self._scalar_bar_text_size

    @scalar_bar_text_size.setter
    def scalar_bar_text_size(self, text_size: int) -> None:
        self._scalar_bar_text_size = text_size

    def render_model(self, reset_camera: bool = True, off_screen: bool = False) -> None:
        """
        Renders the model in a window

        Parameters
        ----------
        reset_camera : bool
            Resets the camera if set to `True`. Default is `True`.
        off_screen : bool
            Renders off-screen without displaying a window. Useful for testing or
            generating images in headless environments. Default is `False`.
        """

        # Set off-screen mode if requested (for testing/headless environments)
        # This must be done BEFORE calling update() so the plotter is in off-screen mode
        # when adding meshes, preventing PyVista camera initialization issues
        if off_screen:
            self.plotter.off_screen = True

        # Update the plotter with the latest geometry
        self.update(reset_camera)

        # Render the model (code execution will pause here until the user closes the window)
        try:
            self.plotter.show(title='Pynite - Simple Finite Element Analysis for Python', auto_close=False)
        finally:
            # Explicitly deep clean the plotter to prevent garbage collection issues
            # Why: PyVista/VTK objects can have circular references that cause AttributeErrors
            #      during Python's garbage collection phase when the program exits.
            # How: deep_clean() forces immediate cleanup of all VTK objects, meshes, and actors
            #      before Python's garbage collector runs, preventing the "'NoneType' object has
            #      no attribute 'check_attribute'" exception that would otherwise appear.
            self.plotter.deep_clean()

    def screenshot(self, filepath: str = './Pynite_Image.png', interact: bool = True, reset_camera: bool = False) -> None:
        """
        Saves a screenshot of the rendered model.

        In non-Jupyter notebook environments, if `interact` is set to `True`, the plotter will show a window that allows you to set the scene prior to the screenshot. Press `q` or click the `X` button when you are ready to capture the screenshot and close the window.

        For Jupyter notebooks, the scene must be set using `render_model` prior to using `screenshot`.

        :param filepath: The filepath to write the image to.
        :type filepath: str, optional
        :param interact: When set to `True` the user can set the scene before the screenshot is taken. Once the scene is set, press 'q' to take the screenshot. Defaults to `True`.
        :type interact: bool, optional
        :param reset_camera: Resets the plotter's camera, forcing it to return to the original rendered view for the screenshot. Note that resetting the camera will ignore any changes you've made to the scene via interaction. Defaults to `False`.
        :type reset_camera: bool, optional
        """

        # Update the plotter with the latest geometry
        self.update(reset_camera)

        # In a Jupyter notebook pyvista will always take the screenshot prior to allowing user interaction. In order to get interaction before taking the screenshot, `render_model` must be called in a cell before `screenshot`. Since `render_model` will show the plotter, there is no need to show it again when using `screenshot`. This next line prevents showing the plotter twice.
        self.plotter.notebook = False

        # For non-Jupyter environments, determine if the user should interact with the window before capturing the screenshot
        if interact == False: self.plotter.off_screen = True
        
        # Save the screenshot to the specified filepath. Note that `auto_close` shuts down the entire plotter after the screenshot is taken, rather than just closing the window. We'll set `auto_close=False` to allow the plotter to remain active. Note that the window must be closed by pressing `q`. Closing it with the 'X' button in the window's corner will close the whole plotter down.
        self.plotter.show(
            title='Pynite - Simple Finite Element Anlaysis for Python',
            screenshot=filepath,
            auto_close=False
            )

    def update(self, reset_camera: bool = True) -> None:
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
        if self.model.load_combos == {} and self.render_loads == True and self.case == None:
            self.render_loads = False
            warnings.warn('Unable to render load combination. No load combinations defined.', UserWarning)

        # Clear out the old plot (if any)
        self.plotter.clear()
        
        # Set up view and axes (works for both interactive and off-screen modes)
        try:
            self.plotter.view_xy()
            self.plotter.show_axes()
            self.plotter.set_viewup((0, 1, 0))
        except:
            # Silently fail if not supported in this context
            pass

        # Clear out internally stored labels (if any)
        self._load_label_points = []
        self._load_labels = []

        self._spring_label_points = []
        self._spring_labels = []
        
        # Check if nodes are to be rendered
        if self.render_nodes == True:

            if self.theme == 'print':
                color = 'black'
            else:
                color = 'grey'

            # Plot each node in the model
            for node in self.model.nodes.values():
                self.plot_node(node, color)
        
        # Render node labels
        label_points = [[node.X, node.Y, node.Z] for node in self.model.nodes.values()]
        labels = [node.name for node in self.model.nodes.values()]
        
        self.plotter.add_point_labels(label_points, labels, bold=False, text_color='black', show_points=True, point_color='grey', point_size=5, shape=None, render_points_as_spheres=True)

        # Check if there are springs in the model
        if self.model.springs:

            # Render the springs
            for spring in self.model.springs.values():
                self.plot_spring(spring, 'grey')
            
            # Render the spring labels
            self.plotter.add_point_labels(self._spring_label_points, self._spring_labels, text_color='black', bold=False, shape=None, render_points_as_spheres=False)
        
        # Render the members
        for member in self.model.members.values():
            self.plot_member(member)

        # Render the member labels
        label_points = [[(member.i_node.X+member.j_node.X)/2, (member.i_node.Y+member.j_node.Y)/2, (member.i_node.Z+member.j_node.Z)/2] for member in self.model.members.values()]
        labels = [member.name for member in self.model.members.values()]
        self.plotter.add_point_labels(label_points, labels, bold=False, text_color='black', show_points=False, shape=None, render_points_as_spheres=False)

        # Render the deformed shape if requested
        if self.deformed_shape == True:

            # Render deformed nodes
            # for node in self.model.nodes.values():
            #     self.plot_deformed_node(node, self.deformed_scale)
            
            # Render deformed members
            for member in self.model.members.values():
                self.plot_deformed_member(member, self.deformed_scale)
            
            # Render deformed springs
            for spring in self.model.springs.values():
                self.plot_spring(spring, 'red', deformed=True)

            # _DeformedShape(self.model, self.deformed_scale, self.annotation_size, self.combo_name, self.render_nodes, self.theme)

        # Render the loads if requested
        if (self.combo_name != None or self.case != None) and self.render_loads != False:

            # Plot the loads
            self.plot_loads()

            # Plot the load labels
            self.plotter.add_point_labels(self._load_label_points, self._load_labels, bold=False, text_color='green', show_points=False, shape=None, render_points_as_spheres=False)
        
        # Render the plates and quads, if present
        if self.model.quads or self.model.plates:
            self.plot_plates(self.deformed_shape, self.deformed_scale, self.color_map, self.combo_name)
        
        # Determine whether to show or hide the scalar bar
        # if self._scalar_bar == False:
        #     self.plotter.scalar_bar.VisibilityOff()

        # Execute user-defined post-update callbacks.
        # Allows plotter customization after internal Pynite updates. See __init__.
        if hasattr(self, 'post_update_callbacks') and self.post_update_callbacks:
            for func in self.post_update_callbacks:
                if callable(func):
                    try:
                        # Pass the plotter instance to the user's function
                        func(self.plotter)
                    except Exception as e:
                        warnings.warn(f"Error executing post-update callback {func.__name__}: {e}")
                else:
                    warnings.warn(f"Item in post_update_callbacks is not callable: {func}")

        # Reset the camera if requested by the user
        if reset_camera:
            self.plotter.reset_camera()


    def plot_node(self, node: Node3D, color: str = 'grey') -> None:
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
                                          x_length=self.annotation_size*2,
                                          y_length=self.annotation_size*2,
                                          z_length=self.annotation_size*2),
                                  color=color)
        
        # Check for a pinned support
        elif node.support_DX and node.support_DY and node.support_DZ and not node.support_RX and not node.support_RY and not node.support_RZ:
            
            # Create a cone using PyVista's Cone function
            self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y - self.annotation_size, node.Z),
                                          direction=(0, 1, 0),
                                          height=self.annotation_size*2,
                                          radius=self.annotation_size*2),
                                  color=color)

        # Other support conditions
        else:

            # Generate a sphere for the node
            # sphere = pv.Sphere(center=(X, Y, Z), radius=0.4*self.annotation_size)
            # self.plotter.add_mesh(sphere, name='Node: '+ node.name, color=color)
            
            # Restrained against X translation
            if node.support_DX:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X - self.annotation_size, node.Y, node.Z),
                                              (node.X + self.annotation_size, node.Y, node.Z)),
                                      color=color)

                # Cones at both ends
                self.plotter.add_mesh(pv.Cone(center=(node.X - self.annotation_size, node.Y,
                                                      node.Z),
                                              direction=(1, 0, 0), height=self.annotation_size*0.6,
                                              radius=self.annotation_size*0.3),
                                      color=color)
                self.plotter.add_mesh(pv.Cone(center=(node.X + self.annotation_size, node.Y,
                                                      node.Z),
                                              direction=(-1, 0, 0),
                                              height=self.annotation_size*0.6,
                                              radius=self.annotation_size*0.3),
                                      color=color)
            
            # Restrained against Y translation
            if node.support_DY:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y - self.annotation_size, node.Z),
                                              (node.X, node.Y + self.annotation_size, node.Z)),
                                      color=color)

                # Cones at both ends
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y - self.annotation_size,
                                                      node.Z), direction=(0, 1, 0),
                                                      height=self.annotation_size*0.6,
                                                      radius=self.annotation_size*0.3),
                                      color=color)
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y + self.annotation_size,
                                                      node.Z),
                                                      direction=(0, -1, 0),
                                                      height=self.annotation_size*0.6,
                                                      radius=self.annotation_size*0.3),
                                      color=color)
            
            # Restrained against Z translation
            if node.support_DZ:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y, node.Z-self.annotation_size),
                                              (node.X, node.Y, node.Z+self.annotation_size)),
                                      color=color)

                # Cones at both ends
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y, node.Z-self.annotation_size),
                                              direction=(0, 0, 1),
                                              height=self.annotation_size*0.6,
                                              radius=self.annotation_size*0.3),
                                      color=color)
                self.plotter.add_mesh(pv.Cone(center=(node.X, node.Y, node.Z+self.annotation_size),
                                              direction=(0, 0, -1),
                                              height=self.annotation_size*0.6,
                                              radius=self.annotation_size*0.3),
                                      color=color)
            
            # Restrained against X rotation
            if node.support_RX:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X-1.6*self.annotation_size, node.Y, node.Z),
                                              (node.X+1.6*self.annotation_size, node.Y, node.Z)),
                                      color=color)

                # Cubes at both ends
                self.plotter.add_mesh(pv.Cube(center=(node.X-1.9*self.annotation_size, node.Y,
                                                      node.Z),
                                              x_length=self.annotation_size*0.6,
                                              y_length=self.annotation_size*0.6,
                                              z_length=self.annotation_size*0.6),
                                      color=color)
                self.plotter.add_mesh(pv.Cube(center=(node.X+1.9 *self.annotation_size, node.Y,
                                                      node.Z),
                                              x_length=self.annotation_size*0.6,
                                              y_length=self.annotation_size*0.6,
                                              z_length=self.annotation_size*0.6),
                                      color=color)

            # Restrained against rotation about the Y-axis
            if node.support_RY:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y-1.6*self.annotation_size, node.Z),
                                              (node.X, node.Y+1.6*self.annotation_size, node.Z)),
                                      color=color)

                # Cubes at both ends
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y-1.9*self.annotation_size,
                                                      node.Z),
                                              x_length=self.annotation_size*0.6,
                                              y_length=self.annotation_size*0.6,
                                              z_length=self.annotation_size*0.6),
                                      color=color)
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y+1.9*self.annotation_size,
                                                      node.Z),
                                              x_length=self.annotation_size*0.6,
                                              y_length=self.annotation_size*0.6,
                                              z_length=self.annotation_size*0.6),
                                      color=color)
            
            # Restrained against rotation about the Z-axis
            if node.support_RZ:

                # Line showing support direction
                self.plotter.add_mesh(pv.Line((node.X, node.Y, node.Z-1.6*self.annotation_size),
                                              (node.X, node.Y, node.Z+1.6*self.annotation_size)),
                                      color=color)

                # Cubes at both ends
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y,
                                                      node.Z-1.9*self.annotation_size),
                                              x_length=self.annotation_size*0.6,
                                              y_length=self.annotation_size*0.6,
                                              z_length=self.annotation_size*0.6),
                                      color=color)
                self.plotter.add_mesh(pv.Cube(center=(node.X, node.Y,
                                                      node.Z+1.9*self.annotation_size),
                                              x_length=self.annotation_size*0.6,
                                              y_length=self.annotation_size*0.6,
                                              z_length=self.annotation_size*0.6),
                                      color=color)

    def plot_member(self, member: Member3D, theme: str = 'default') -> None:
        """
        Adds a member to the plotter. This method generates a line representing a structural member between two nodes, and adds it to the plotter with specified theme settings.
        
        Parameters
        ==========
        :param member: The structural member to be plotted, containing information about its end nodes.
        :type member: Member
        :param theme: The theme for plotting the member. Default is 'default'.
        :type theme: str
        """
    
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

    def plot_spring(self, spring: Spring3D, color: str = 'grey', deformed: bool = False) -> None:
        """
        Adds a spring to the plotter. This method generates a zig-zag line representing a spring between two nodes, and adds it to the plotter with specified theme settings.
        """
        
        # Scale the spring's zigzags
        size = self.annotation_size

        # Find the spring's i-node and j-node
        i_node = spring.i_node
        j_node = spring.j_node

        # Find the spring's node coordinates
        Xi, Yi, Zi = i_node.X, i_node.Y, i_node.Z
        Xj, Yj, Zj = j_node.X, j_node.Y, j_node.Z

        # Determine if the spring should be plotted in its deformed shape
        if deformed:
            Xi = Xi + i_node.DX[self.combo_name]*self.deformed_scale
            Yi = Yi + i_node.DY[self.combo_name]*self.deformed_scale
            Zi = Zi + i_node.DZ[self.combo_name]*self.deformed_scale
            Xj = Xj + j_node.DX[self.combo_name]*self.deformed_scale
            Yj = Yj + j_node.DY[self.combo_name]*self.deformed_scale
            Zj = Zj + j_node.DZ[self.combo_name]*self.deformed_scale

        # Calculate the spring direction vector and length
        direction = np.array([Xj, Yj, Zj]) - np.array([Xi, Yi, Zi])
        length = np.linalg.norm(direction)

        # Normalize the direction vector
        direction = direction / length

        # Calculate perpendicular vectors for zig-zag plane
        arbitrary_vector = np.array([1, 0, 0])
        if np.allclose(direction, arbitrary_vector) or np.allclose(direction, -arbitrary_vector):
            arbitrary_vector = np.array([0, 1, 0])
        perp_vector1 = np.cross(direction, arbitrary_vector)
        perp_vector1 /= np.linalg.norm(perp_vector1)
        perp_vector2 = np.cross(direction, perp_vector1)
        perp_vector2 /= np.linalg.norm(perp_vector2)

        # Define the length of the straight segments
        straight_segment_length = length / 10
        zigzag_length = length - 2 * straight_segment_length

        # Generate points for the zig-zag line
        num_zigs = 4
        num_points = num_zigs * 2
        amplitude = size
        t = np.linspace(0, zigzag_length, num_points)
        zigzag_pattern = amplitude * np.tile([1, -1], num_zigs)
        zigzag_points = np.outer(t, direction) + np.outer(zigzag_pattern, perp_vector1)

        # Add the straight segments to the start and end
        start_point = np.array([Xi, Yi, Zi])
        end_point = np.array([Xj, Yj, Zj])
        start_segment = start_point + direction * straight_segment_length
        end_segment = end_point - direction * straight_segment_length

        # Adjust the zigzag points to the correct position
        zigzag_points += start_segment

        # Combine the points
        points = np.vstack([start_point, start_segment, zigzag_points, end_segment, end_point])

        # Add lines connecting the points
        num_points = len(points)
        lines = np.zeros((num_points - 1, 3), dtype=int)
        lines[:, 0] = 2
        lines[:, 1] = np.arange(num_points - 1, dtype=int)
        lines[:, 2] = np.arange(1, num_points, dtype=int)

        # Create a PolyData object for the zig-zag line
        zigzag_line = pv.PolyData(points, lines=lines)

        # Create a plotter and add the zig-zag line
        self.plotter.add_mesh(zigzag_line, color=color, line_width=2)

        # Add the spring label to the list of labels
        self._spring_labels.append(spring.name)
        self._spring_label_points.append([(Xi + Xj) / 2, (Yi + Yj) / 2, (Zi + Zj) / 2])

            
    def plot_plates(self, deformed_shape: bool, deformed_scale: float, color_map: Optional[str], combo_name: Optional[str]) -> None:
        
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
        for item in list(self.model.plates.values()) + list(self.model.quads.values()):
            
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
        plate_vertices = np.array(plate_vertices)
        plate_faces = np.array(plate_faces)

        # Create a new PyVista dataset to store plate data
        plate_polydata = pv.PolyData(plate_vertices, plate_faces)

        # Add the results as point data to the PyVista dataset
        if color_map:
            
            plate_polydata = plate_polydata.separate_cells()
            plate_polydata['Contours'] = np.array(plate_results)
            
            # Add the scalar bar for the contours
            if self._scalar_bar == True:
                self.plotter.add_mesh(plate_polydata, scalars='Contours', show_edges=True)
            else:
                self.plotter.add_mesh(plate_polydata)

        else:
            self.plotter.add_mesh(plate_polydata)
      
    def plot_deformed_node(self, node: Node3D, scale_factor: float, color: str = 'grey') -> None:

        # Calculate the node's deformed position
        newX = node.X + scale_factor * (node.DX[self.combo_name])
        newY = node.Y + scale_factor * (node.DY[self.combo_name])
        newZ = node.Z + scale_factor * (node.DZ[self.combo_name])

        # Generate a sphere source for the node in its deformed position
        sphere = pv.Sphere(radius=0.4*self.annotation_size, center=[newX, newY, newZ])

        # Add the mesh to the plotter
        self.plotter.add_mesh(sphere, color=color)
  
    def plot_deformed_member(self, member: Member3D, scale_factor: float) -> None:
        
        # Determine if this member is active for each load combination
        if member.active:
        
            L = member.L() # Member length
            T = member.T() # Member local transformation matrix
        
            cos_x = np.array([T[0, 0:3]]) # Direction cosines of local x-axis
            cos_y = np.array([T[1, 0:3]]) # Direction cosines of local y-axis
            cos_z = np.array([T[2, 0:3]]) # Direction cosines of local z-axis

            # Find the initial position of the local i-node
            Xi = member.i_node.X
            Yi = member.i_node.Y
            Zi = member.i_node.Z
        
            # Calculate the local y-axis displacements at 20 points along the member's length
            DY_plot = np.empty((0, 3))
            for i in range(20):
                    
                # Calculate the local y-direction displacement
                dy_tot = member.deflection('dy', L / 19 * i, self.combo_name)
            
                # Calculate the scaled displacement in global coordinates
                DY_plot = np.append(DY_plot, dy_tot * cos_y * scale_factor, axis=0)
        
            # Calculate the local z-axis displacements at 20 points along the member's length
            DZ_plot = np.empty((0, 3)) 
            for i in range(20):
                    
                # Calculate the local z-direction displacement
                dz_tot = member.deflection('dz', L / 19 * i, self.combo_name)
            
                # Calculate the scaled displacement in global coordinates
                DZ_plot = np.append(DZ_plot, dz_tot * cos_z * scale_factor, axis=0)
        
            # Calculate the local x-axis displacements at 20 points along the member's length
            DX_plot = np.empty((0, 3)) 
            for i in range(20):
                    
                # Displacements in local coordinates
                dx_tot = [[Xi, Yi, Zi]] + (L / 19 * i + member.deflection('dx', L / 19 * i, self.combo_name) * scale_factor) * cos_x
                    
                # Magnified displacements in global coordinates
                DX_plot = np.append(DX_plot, dx_tot, axis=0)
            
            # Sum the component displacements to obtain overall displacement
            D_plot = DY_plot + DZ_plot + DX_plot
            
            # Create lines connecting the points
            for i in range(len(D_plot)-1):
                line = pv.Line(D_plot[i], D_plot[i+1])
                self.plotter.add_mesh(line, color='red', line_width=2)

    def plot_pt_load(self, position: Tuple[float, float, float], direction: Union[Tuple[float, float, float], np.ndarray], 
                    length: float, label_text: Optional[Union[str, float, int]] = None, color: str = 'green') -> None:

        # Create a unit vector in the direction of the 'direction' vector
        unitVector = direction/np.linalg.norm(direction)
        
        # Determine if the load is positive or negative
        if length == 0:
            sign = 1
        else:
            sign = abs(length)/length

        # Generate the tip of the load arrow
        tip_length = abs(length) / 4
        radius = abs(length) / 16
        tip = pv.Cone(center=(position[0] - tip_length*sign*unitVector[0]/2,
                              position[1] - tip_length*sign*unitVector[1]/2,
                              position[2] - tip_length*sign*unitVector[2]/2),
                              direction=(direction[0]*sign, direction[1]*sign, direction[2]*sign),
                              height=tip_length, radius=radius)

        # Plot the tip
        self.plotter.add_mesh(tip, color=color)

        # Create the shaft (you'll need to specify the second point)
        X_tail = position[0] - unitVector[0]*length
        Y_tail = position[1] - unitVector[1]*length
        Z_tail = position[2] - unitVector[2]*length
        shaft = pv.Line(pointa=position, pointb=(X_tail, Y_tail, Z_tail))
        
        # Save the data necessary to create the load's label
        if label_text is not None:
            self._load_labels.append(sig_fig_round(label_text, 3))
            self._load_label_points.append([X_tail, Y_tail, Z_tail])
        
        # Plot the shaft
        self.plotter.add_mesh(shaft, line_width=2, color=color)                         

    def plot_dist_load(self, position1: Tuple[float, float, float], position2: Tuple[float, float, float], 
                      direction: Union[np.ndarray, Tuple[float, float, float]], length1: float, length2: float,
                      label_text1: Optional[Union[str, float, int]], label_text2: Optional[Union[str, float, int]], 
                      color: str = 'green') -> None:

        # Calculate the length of the distributed load
        load_length = ((position2[0] - position1[0])**2 + (position2[1] - position1[1])**2 + (position2[2] - position1[2])**2)**0.5

        # Find the direction cosines for the line the load acts on
        line_dir_cos = [(position2[0] - position1[0])/load_length,
                        (position2[1] - position1[1])/load_length,
                        (position2[2] - position1[2])/load_length]

        # Find the direction cosines for the direction the load acts in
        dir_dir_cos = direction/np.linalg.norm(direction)

        # Create point loads at intervals roughly equal to 75% of the load's largest length
        # Add text labels to the first and last load arrow
        if load_length > 0:
            num_steps = int(round(0.75 * load_length/max(abs(length1), abs(length2)), 0))
        else:
            num_steps = 0

        num_steps = max(num_steps, 1)
        step = load_length/num_steps

        for i in range(num_steps + 1):

            # Calculate the position (X, Y, Z) of this load arrow's point
            position = (position1[0] + i*step*line_dir_cos[0],
                        position1[1] + i*step*line_dir_cos[1],
                        position1[2] + i*step*line_dir_cos[2])

            # Determine the length of this load arrow
            length = length1 + (length2 - length1)/load_length*i*step

            # Determine the label's text
            if i == 0:
                label_text = label_text1
            elif i == num_steps:
                label_text = label_text2
            else:
                label_text = None

            # Plot the load arrow
            self.plot_pt_load(position, dir_dir_cos, length, label_text, color)

        # Draw a line between the first and last load arrow's tails (using cylinder here for better visualization)
        tail_line = pv.Line(position1 - dir_dir_cos*length1, position2 - dir_dir_cos*length2)

        # Combine all geometry into a single PolyData object
        self.plotter.add_mesh(tail_line, color=color)

    def plot_moment(self, center: Tuple[float, float, float], direction: Union[Tuple[float, float, float], np.ndarray], 
                    radius: float, label_text: Optional[Union[str, float, int]] = None, color: str = 'green') -> None:

        # Convert the direction vector into a unit vector
        v1 = direction/np.linalg.norm(direction)

        # Find any vector perpendicular to the moment direction vector. This will serve as a
        # vector from the center of the arc pointing to the tail of the moment arc.
        v2 = _PerpVector(v1)
        
        # Generate the arc for the moment
        arc = pv.CircularArcFromNormal(center, resolution=20, normal=v1, angle=215, polar=v2*radius)
        
        # Add the arc to the plot
        self.plotter.add_mesh(arc, line_width=2, color=color)

        # Generate the arrow tip at the end of the arc
        tip_length = radius/4
        cone_radius = radius/16
        cone_direction = -np.cross(v1, arc.center - arc.points[-1])
        tip = pv.Cone(center=arc.points[-1], direction=cone_direction, height=tip_length,
                      radius=cone_radius)

        # Add the tip to the plot
        self.plotter.add_mesh(tip, color=color)

        # Create the text label
        if label_text:
            text_pos = center + (radius + 0.25*self.annotation_size)*v2
            self._load_label_points.append(text_pos)
            self._load_labels.append(label_text)

    def plot_area_load(self, position0, position1, position2, position3, direction, length, label_text, color='green'):

        # Find the direction cosines for the direction the load acts in
        dir_dir_cos = direction / np.linalg.norm(direction)

        # Find the positions of the tails of all the arrows at the corners
        self.p0 = position0 - dir_dir_cos * length
        self.p1 = position1 - dir_dir_cos * length
        self.p2 = position2 - dir_dir_cos * length
        self.p3 = position3 - dir_dir_cos * length

        # Plot the area load arrows
        self.plot_pt_load(position0, dir_dir_cos, length, label_text, color)
        self.plot_pt_load(position1, dir_dir_cos, length, color=color)
        self.plot_pt_load(position2, dir_dir_cos, length, color=color)
        self.plot_pt_load(position3, dir_dir_cos, length, color=color)

        # Create the area load polygon (quad)
        quad = pv.Quadrilateral([self.p0, self.p1, self.p2, self.p3])

        self.plotter.add_mesh(quad, color=color)

    def _calc_max_loads(self):

        max_pt_load = 0
        max_moment = 0
        max_dist_load = 0
        max_area_load = 0
        
        # Find the requested load combination or load case
        if self.case == None:
        
            # Step through each node
            for node in self.model.nodes.values():
        
                # Step through each nodal load to find the largest one
                for load in node.NodeLoads:
                
                    # Find the largest loads in the load combination
                    if load[2] in self.model.load_combos[self.combo_name].factors:
                        if load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ':
                            if abs(load[1]*self.model.load_combos[self.combo_name].factors[load[2]]) > max_pt_load:
                                max_pt_load = abs(load[1]*self.model.load_combos[self.combo_name].factors[load[2]])
                        else:
                            if abs(load[1]*self.model.load_combos[self.combo_name].factors[load[2]]) > max_moment:
                                max_moment = abs(load[1]*self.model.load_combos[self.combo_name].factors[load[2]])
        
            # Step through each member
            for member in self.model.members.values():
        
                # Step through each member point load
                for load in member.PtLoads:
                    
                    # Find and store the largest point load and moment in the load combination
                    if load[3] in self.model.load_combos[self.combo_name].factors:
            
                        if (load[0] == 'Fx' or load[0] == 'Fy' or load[0] == 'Fz'
                        or  load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ'):
                            if abs(load[1]*self.model.load_combos[self.combo_name].factors[load[3]]) > max_pt_load:
                                max_pt_load = abs(load[1]*self.model.load_combos[self.combo_name].factors[load[3]])
                        else:
                            if abs(load[1]*self.model.load_combos[self.combo_name].factors[load[3]]) > max_moment:
                                max_moment = abs(load[1]*self.model.load_combos[self.combo_name].factors[load[3]])
        
                # Step through each member distributed load
                for load in member.DistLoads:
            
                    #Find and store the largest distributed load in the load combination
                    if load[5] in self.model.load_combos[self.combo_name].factors:
            
                        if abs(load[1]*self.model.load_combos[self.combo_name].factors[load[5]]) > max_dist_load:
                            max_dist_load = abs(load[1]*self.model.load_combos[self.combo_name].factors[load[5]])
                        if abs(load[2]*self.model.load_combos[self.combo_name].factors[load[5]]) > max_dist_load:
                            max_dist_load = abs(load[2]*self.model.load_combos[self.combo_name].factors[load[5]])
        
            # Step through each plate
            for plate in self.model.plates.values():
        
                # Step through each plate load
                for load in plate.pressures:
            
                    if load[1] in self.model.load_combos[self.combo_name].factors:
                        if abs(load[0]*self.model.load_combos[self.combo_name].factors[load[1]]) > max_area_load:
                            max_area_load = abs(load[0]*self.model.load_combos[self.combo_name].factors[load[1]])
        
            # Step through each quad
            for quad in self.model.quads.values():
        
                # Step through each plate load
                for load in quad.pressures:
            
                    # Check to see if the load case is in the requested load combination
                    if load[1] in self.model.load_combos[self.combo_name].factors:
                        if abs(load[0]*self.model.load_combos[self.combo_name].factors[load[1]]) > max_area_load:
                            max_area_load = abs(load[0]*self.model.load_combos[self.combo_name].factors[load[1]])
        
        # Behavior if case has been specified
        else:
            
            # Step through each node
            for node in self.model.nodes.values():
        
                # Step through each nodal load to find the largest one
                for load in node.NodeLoads:
                
                    # Find the largest loads in the load case
                    if load[2] == self.case:
                        if load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ':
                            if abs(load[1]) > max_pt_load:
                                max_pt_load = abs(load[1])
                        else:
                            if abs(load[1]) > max_moment:
                                max_moment = abs(load[1])
        
            # Step through each member
            for member in self.model.members.values():
        
                # Step through each member point load
                for load in member.PtLoads:
                
                    # Find and store the largest point load and moment in the load case
                    if load[3] == self.case:
                
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
                    if load[5] == self.case:
                
                        if abs(load[1]) > max_dist_load:
                            max_dist_load = abs(load[1])
                        if abs(load[2]) > max_dist_load:
                            max_dist_load = abs(load[2])
                
                # Step through each plate
                for plate in self.model.plates.values():
            
                    # Step through each plate load
                    for load in plate.pressures:
            
                        if load[1] == self.case:
                
                            if abs(load[0]) > max_area_load:
                                max_area_load = abs(load[0])
            
            # Step through each quad
            for quad in self.model.quads.values():
        
                # Step through each plate load
                for load in quad.pressures:
            
                    if load[1] == self.case:
            
                        if abs(load[0]) > max_area_load:
                            max_area_load = abs(load[0])
            
        # Prevent division by zero errors by ensuring max values are never zero
        # If a load type has no loads, set it to 1 to avoid crashes during normalization
        if max_pt_load == 0:
            max_pt_load = 1
        if max_moment == 0:
            max_moment = 1
        if max_dist_load == 0:
            max_dist_load = 1
        if max_area_load == 0:
            max_area_load = 1
        
        # Return the maximum loads for the load combo or load case
        return max_pt_load, max_moment, max_dist_load, max_area_load
            
    def plot_loads(self):
      
        # Get the maximum load magnitudes that will be used to normalize the display scale
        max_pt_load, max_moment, max_dist_load, max_area_load = self._calc_max_loads()
        
        # Display the requested load combination, or 'Combo 1' if no load combo or case has been
        # specified
        if self.case is None:
            # Store model.load_combos[combo].factors under a simpler name for use below
            load_factors = self.model.load_combos[self.combo_name].factors
        else:
            # Set up a load combination dictionary that represents the load case
            load_factors = {self.case: 1}
        
        # Step through each node
        for node in self.model.nodes.values():
        
            # Step through and display each nodal load
            for load in node.NodeLoads:
            
                # Determine if this load is part of the requested LoadCombo or case
                if load[2] in load_factors:
                
                    # Calculate the factored value for this load and it's sign (positive or
                    # negative)
                    load_value = load[1]*load_factors[load[2]]
                    if load_value != 0:
                        sign = load_value/abs(load_value)
                    else:
                        sign = 1
                    
                    # Determine the direction of this load
                    if load[0] == 'FX' or load[0] == 'MX': direction = (sign, 0, 0)
                    elif load[0] == 'FY' or load[0] == 'MY': direction = (0, sign, 0)
                    elif load[0] == 'FZ' or load[0] == 'MZ': direction = (0, 0, sign)

                    # Display the load
                    if load[0] in {'FX', 'FY', 'FZ'}:
                        self.plot_pt_load((node.X, node.Y, node.Z), direction,
                                          abs(load_value/max_pt_load)*5*self.annotation_size,
                                          load_value, 'green')
                    elif load[0] in {'MX', 'MY', 'MZ'}:
                        self.plot_moment((node.X, node.Y, node.Z), direction, abs(load_value/max_moment)*2.5*self.annotation_size, str(load_value), 'green')
        
        # Step through each member
        for member in self.model.members.values():
        
            # Get the direction cosines for the member's local axes
            dir_cos = member.T()[0:3, 0:3]
        
            # Get the starting point for the member
            x_start, y_start, z_start = member.i_node.X, member.i_node.Y, member.i_node.Z
        
            # Step through each member point load
            for load in member.PtLoads:
        
                # Determine if this load is part of the requested load combination
                if load[3] in load_factors:
            
                    # Calculate the factored value for this load and it's sign (positive or negative)
                    load_value = load[1]*load_factors[load[3]]
                    sign = load_value/abs(load_value)
            
                    # Calculate the load's location in 3D space
                    x = load[2]
                    position = [x_start + dir_cos[0, 0]*x, y_start + dir_cos[0, 1]*x, z_start + dir_cos[0, 2]*x]
            
                    # Display the load
                    if load[0] == 'Fx':
                        self.plot_pt_load(position, dir_cos[0, :], load_value/max_pt_load*5*self.annotation_size, load_value)
                    elif load[0] == 'Fy':
                        self.plot_pt_load(position, dir_cos[1, :], load_value/max_pt_load*5*self.annotation_size, load_value)
                    elif load[0] == 'Fz':
                        self.plot_pt_load(position, dir_cos[2, :], load_value/max_pt_load*5*self.annotation_size, load_value)
                    elif load[0] == 'Mx':
                        self.plot_moment(position, dir_cos[0, :]*sign, abs(load_value)/max_moment*2.5*self.annotation_size, str(load_value))
                    elif load[0] == 'My':
                        self.plot_moment(position, dir_cos[1, :]*sign, abs(load_value)/max_moment*2.5*self.annotation_size, str(load_value))
                    elif load[0] == 'Mz':
                        self.plot_moment(position, dir_cos[2, :]*sign, abs(load_value)/max_moment*2.5*self.annotation_size, str(load_value))
                    elif load[0] == 'FX':
                        self.plot_pt_load(position, [1, 0, 0], load_value/max_pt_load*5*self.annotation_size, load_value)
                    elif load[0] == 'FY':
                        self.plot_pt_load(position, [0, 1, 0], load_value/max_pt_load*5*self.annotation_size, load_value)
                    elif load[0] == 'FZ':
                        self.plot_pt_load(position, [0, 0, 1], load_value/max_pt_load*5*self.annotation_size, load_value)
                    elif load[0] == 'MX':
                        self.plot_moment(position, [1*sign, 0, 0], abs(load_value)/max_moment*2.5*self.annotation_size, str(load_value))
                    elif load[0] == 'MY':
                        self.plot_moment(position, [0, 1*sign, 0], abs(load_value)/max_moment*2.5*self.annotation_size, str(load_value))
                    elif load[0] == 'MZ':
                        self.plot_moment(position, [0, 0, 1*sign], abs(load_value)/max_moment*2.5*self.annotation_size, str(load_value))
        
            # Step through each member distributed load
            for load in member.DistLoads:
        
                # Determine if this load is part of the requested load combination
                if load[5] in load_factors:
            
                    # Calculate the factored value for this load and it's sign (positive or negative)
                    w1 = load[1]*load_factors[load[5]]
                    w2 = load[2]*load_factors[load[5]]
            
                    # Calculate the loads location in 3D space
                    x1 = load[3]
                    x2 = load[4]
                    position1 = [x_start + dir_cos[0, 0]*x1, y_start + dir_cos[0, 1]*x1, z_start + dir_cos[0, 2]*x1]
                    position2 = [x_start + dir_cos[0, 0]*x2, y_start + dir_cos[0, 1]*x2, z_start + dir_cos[0, 2]*x2]
                    
                    # Display the load
                    if load[0] in {'Fx', 'Fy', 'Fz', 'FX', 'FY', 'FZ'}:

                        # Determine the load direction
                        if load[0] == 'Fx': direction = dir_cos[0, :]
                        elif load[0] == 'Fy': direction = dir_cos[1, :]
                        elif load[0] == 'Fz': direction = dir_cos[2, :]
                        elif load[0] == 'FX': direction = [1, 0, 0]
                        elif load[0] == 'FY': direction = [0, 1, 0]
                        elif load[0] == 'FZ': direction = [0, 0, 1]

                        # Plot the distributed load
                        self.plot_dist_load(position1, position2, direction, w1/max_dist_load*5*self.annotation_size, w2/max_dist_load*5*self.annotation_size, str(sig_fig_round(w1, 3)), str(sig_fig_round(w2, 3)), 'green')
        
        # Step through each plate
        for plate in list(self.model.plates.values()) + list(self.model.quads.values()):
        
            # Get the direction cosines for the plate's local z-axis
            dir_cos = plate.T()[0:3, 0:3]
            dir_cos = dir_cos[2]
        
            # Step through each plate load
            for load in plate.pressures:
        
                # Determine if this load is part of the requested load combination
                if load[1] in load_factors:
            
                    # Calculate the factored value for this load
                    load_value = load[0]*load_factors[load[1]]
                    
                    # Find the sign for this load. Intercept any divide by zero errors
                    if load[0] == 0:
                        sign = 1
                    else:
                        sign = abs(load[0])/load[0]
            
                    # Find the position of the load's 4 corners
                    position0 = [plate.i_node.X, plate.i_node.Y, plate.i_node.Z]
                    position1 = [plate.j_node.X, plate.j_node.Y, plate.j_node.Z]
                    position2 = [plate.m_node.X, plate.m_node.Y, plate.m_node.Z]
                    position3 = [plate.n_node.X, plate.n_node.Y, plate.n_node.Z]
            
                    # Create an area load and get its data
                    self.plot_area_load(position0, position1, position2, position3, dir_cos*sign, load_value/max_area_load*5*self.annotation_size, str(sig_fig_round(load_value, 3)), color='green')

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
    return [i2, j2, k2]/np.linalg.norm([i2, j2, k2])

def _PrepContour(model, stress_type='Mx', combo_name='Combo 1'):

    if stress_type != None:
    
        # Erase any previous contours
        for node in model.nodes.values():
            node.contour = []
        
        # Check for global stresses:
        if stress_type in ['MX', 'MY', 'MZ', 'QX', 'QY', 'QZ', 'SX', 'SY']:
            local = False
        else:
            local = True
    
        # Step through each element in the model
        for element in list(model.quads.values()) + list(model.plates.values()):
            
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
                i, j, m, n = element.d(combo_name)[[2, 8, 14, 20], :]
                element.i_node.contour.append(i)
                element.j_node.contour.append(j)
                element.m_node.contour.append(m)
                element.n_node.contour.append(n)
            elif stress_type.upper() == 'MX':
                element.i_node.contour.append(element.moment(r_left, s_bot, local, combo_name)[0])
                element.j_node.contour.append(element.moment(r_right, s_bot, local, combo_name)[0])
                element.m_node.contour.append(element.moment(r_right, s_top, local, combo_name)[0])
                element.n_node.contour.append(element.moment(r_left, s_top, local, combo_name)[0])
            elif stress_type.upper() == 'MY':
                element.i_node.contour.append(element.moment(r_left, s_bot, local, combo_name)[1])
                element.j_node.contour.append(element.moment(r_right, s_bot, local, combo_name)[1])
                element.m_node.contour.append(element.moment(r_right, s_top, local, combo_name)[1])
                element.n_node.contour.append(element.moment(r_left, s_top, local, combo_name)[1])
            elif stress_type.upper() == 'MXY':
                element.i_node.contour.append(element.moment(r_left, s_bot, local, combo_name)[2])
                element.j_node.contour.append(element.moment(r_right, s_bot, local, combo_name)[2])
                element.m_node.contour.append(element.moment(r_right, s_top, local, combo_name)[2])
                element.n_node.contour.append(element.moment(r_left, s_top, local, combo_name)[2])
            elif stress_type.upper() == 'QX':
                element.i_node.contour.append(element.shear(r_left, s_bot, local, combo_name)[0])
                element.j_node.contour.append(element.shear(r_right, s_bot, local, combo_name)[0])
                element.m_node.contour.append(element.shear(r_right, s_top, local, combo_name)[0])
                element.n_node.contour.append(element.shear(r_left, s_top, local, combo_name)[0])
            elif stress_type.upper() == 'QY':
                element.i_node.contour.append(element.shear(r_left, s_bot, local, combo_name)[1])
                element.j_node.contour.append(element.shear(r_right, s_bot, local, combo_name)[1])
                element.m_node.contour.append(element.shear(r_right, s_top, local, combo_name)[1])
                element.n_node.contour.append(element.shear(r_left, s_top, local, combo_name)[1])
            elif stress_type.upper() == 'SX':
                element.i_node.contour.append(element.membrane(r_left, s_bot, local, combo_name)[0])
                element.j_node.contour.append(element.membrane(r_right, s_bot, local, combo_name)[0])
                element.m_node.contour.append(element.membrane(r_right, s_top, local, combo_name)[0])
                element.n_node.contour.append(element.membrane(r_left, s_top, local, combo_name)[0])
            elif stress_type.upper() == 'SY':
                element.i_node.contour.append(element.membrane(r_left, s_bot, local, combo_name)[1])
                element.j_node.contour.append(element.membrane(r_right, s_bot, local, combo_name)[1])
                element.m_node.contour.append(element.membrane(r_right, s_top, local, combo_name)[1])
                element.n_node.contour.append(element.membrane(r_left, s_top, local, combo_name)[1])
            elif stress_type.upper() == 'TXY':
                element.i_node.contour.append(element.membrane(r_left, s_bot, local, combo_name)[2])
                element.j_node.contour.append(element.membrane(r_right, s_bot, local, combo_name)[2])
                element.m_node.contour.append(element.membrane(r_right, s_top, local, combo_name)[2])
                element.n_node.contour.append(element.membrane(r_left, s_top, local, combo_name)[2])                 
      
        # Average the values at each node to obtain a smoothed contour
        for node in model.nodes.values():
            # Prevent divide by zero errors for nodes with no contour values
            if node.contour != []:
                node.contour = sum(node.contour)/len(node.contour)

def sig_fig_round(number, sig_figs):

    # Check for strings or other convertible data types
    if not isinstance(number, (float, int)):
        try:
            number = float(number)
        except:
            raise ValueError(f"{number} is not a number. Ensure that `number` is numeric.")

    if number == 0:
        return 0

    # Calculate the magnitude of the number
    magnitude = math.floor(math.log10(abs(number)))

    # Calculate the number of decimal places to round to
    decimal_places = sig_figs - 1 - magnitude

    # Round the number to the specified number of decimal places
    rounded_number = round(number, decimal_places)

    return rounded_number
