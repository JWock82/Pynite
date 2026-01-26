"""PyVista-based visualization for PyNite finite element models.

This module provides classes and methods for rendering 3D models with support
for undeformed/deformed geometry, load glyphs, contours, and internal diagrams.
Docstrings follow reStructuredText field-list format for Sphinx autodoc.
"""

from __future__ import annotations # Allows more recent type hints features
from json import load
import warnings
from typing import TYPE_CHECKING, Callable, List, Any, Optional, Union, Tuple

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
    """Render finite element models using PyVista.

    This class provides methods for rendering 3D models (nodes, members, springs,
    plates, quads) with deformed shapes, loads, contours, and diagrams.
    """

    scalar: Optional[str] = None

    def __init__(self, model: FEModel3D) -> None:
        """Initialize the renderer with a finite element model.

        :param FEModel3D model: Finite element model to render.
        """
        self.model: FEModel3D = model

        # Default settings for rendering
        self._annotation_size: Optional[float] = None  # None means auto-calculate
        self._annotation_size_manual: bool = False  # Track if user manually set the size
        self._annotation_size_cached: Optional[float] = None  # Cache to avoid recalculating 2600+ times per render
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
        self._member_diagrams: Optional[str] = None  # Options: None, 'Fy', 'Fz', 'My', 'Mz', 'Fx', 'Tx'
        self._diagram_scale: float = 30.0
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
        """Size of text annotations and visual elements in model units.

        Auto-calculated as 5% of shortest distance between nodes if not
        manually set.
        """
        if self._annotation_size is None or not self._annotation_size_manual:
            # Return cached value if available; calculate once on first access
            if self._annotation_size_cached is None:
                self._annotation_size_cached = self._calculate_auto_annotation_size()
            return self._annotation_size_cached
        return self._annotation_size

    @annotation_size.setter
    def annotation_size(self, size: float) -> None:
        self._annotation_size = size
        self._annotation_size_manual = True  # Mark as manually set

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

    @property
    def member_diagrams(self) -> Optional[str]:
        """Member diagram type to display.

        Options: ``None``, ``'Fy'``, ``'Fz'``, ``'My'``, ``'Mz'``, ``'Fx'``, ``'Tx'``.
        """
        return self._member_diagrams

    @member_diagrams.setter
    def member_diagrams(self, diagram_type: Optional[str]) -> None:
        valid_options = [None, 'Fy', 'Fz', 'My', 'Mz', 'Fx', 'Tx']
        if diagram_type not in valid_options:
            raise ValueError(f"member_diagrams must be one of {valid_options}, got '{diagram_type}'")
        self._member_diagrams = diagram_type

    @property
    def diagram_scale(self) -> float:
        """Scale factor for member diagram visualization."""
        return self._diagram_scale

    @diagram_scale.setter
    def diagram_scale(self, scale: float) -> None:
        self._diagram_scale = scale

    def _calculate_auto_annotation_size(self) -> float:
        """Calculate automatic annotation size as 5% of shortest node distance.

        Uses vectorized NumPy operations for fast computation on large meshes.
        Result is cached by the annotation_size property to avoid recalculation.

        :returns: Annotation size in model units (``5.0`` fallback if <2 nodes).
        :rtype: float
        """
        from numpy import asarray
        
        nodes = list(self.model.nodes.values())

        # Need at least 2 nodes to calculate distance
        if len(nodes) < 2:
            return 5.0  # Default fallback

        # Extract node coordinates as numpy array (much faster than nested loops)
        coords = asarray([[node.X, node.Y, node.Z] for node in nodes])
        
        # Calculate pairwise distances using vectorized operations
        # For each node, compute distance to all other nodes, then find minimum
        min_distance = float('inf')
        for i in range(len(coords)):
            # Vectorized distance calculation for node i to all others
            diffs = coords[i+1:] - coords[i]
            distances = (diffs**2).sum(axis=1)**0.5
            
            # Find minimum distance from this node to subsequent nodes
            valid_distances = distances[distances > 0]
            if len(valid_distances) > 0:
                local_min = valid_distances.min()
                if local_min < min_distance:
                    min_distance = local_min

        # If all nodes are at same location, use default
        if min_distance == float('inf') or min_distance == 0:
            return 5.0

        # Return 5% of shortest distance
        return min_distance * 0.05

    def render_model(self, reset_camera: bool = True, off_screen: bool = False) -> None:
        """Render the model in a window.

        :param bool reset_camera: Reset the camera before rendering (default ``True``).
        :param bool off_screen: Render off-screen without displaying a window
            (useful for headless environments, default ``False``).
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
        """Save a screenshot of the rendered model.

        In non-Jupyter environments, if ``interact=True`` the user can adjust
        the view before taking the screenshot (press ``'q'`` to proceed).
        For Jupyter, use ``render_model`` to set the scene before ``screenshot``.

        :param str filepath: File path to write PNG image (default ``'./Pynite_Image.png'``).
        :param bool interact: Allow user interaction before screenshot (default ``True``).
        :param bool reset_camera: Reset camera to original view (default ``False``).
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
        """Rebuild the PyVista plotter with current settings.

        Updates all visual elements (geometry, loads, labels, contours, diagrams)
        and calls post-update callbacks.

        :param bool reset_camera: Reset camera to fit model (default ``True``).
        """

        # Clear annotation size cache to recalculate if model has changed
        self._annotation_size_cached = None

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

        # Build visual helper objects. These classes encapsulate geometry and label bookkeeping
        # so that the renderer logic stays readable while still leveraging PyVista efficiently.
        vis_nodes: List[VisNode] = []
        vis_springs: List[VisSpring] = []
        vis_members: List[VisMember] = []

        if self.render_nodes:
            node_color = 'black' if self.theme == 'print' else 'grey'
            vis_nodes = [VisNode(node, self.annotation_size, node_color) for node in self.model.nodes.values()]
            for vis_node in vis_nodes:
                vis_node.add_to_plotter(self.plotter)

        if self.model.springs:
            # Undeformed springs are added here; deformed ones come later if requested.
            vis_springs = [VisSpring(spring, self.annotation_size, 'grey', False, self.deformed_scale, self.combo_name) for spring in self.model.springs.values()]
            for vis_spring in vis_springs:
                vis_spring.add_to_plotter(self.plotter)

        if self.model.members:
            vis_members = [VisMember(member, self.theme) for member in self.model.members.values()]
            for vis_member in vis_members:
                vis_member.add_to_plotter(self.plotter)

        # Render the deformed shape if requested. Deformed visuals are built separately so the
        # undeformed geometry remains visible alongside deformed overlays when desired.
        if self.deformed_shape:
            for member in self.model.members.values():
                vis_def_member = VisDeformedMember(member, self.deformed_scale, self.combo_name)
                vis_def_member.add_to_plotter(self.plotter)

            for spring in self.model.springs.values():
                vis_def_spring = VisSpring(spring, self.annotation_size, 'red', True, self.deformed_scale, self.combo_name)
                vis_def_spring.add_to_plotter(self.plotter)

        # Render labels if requested. We gather label locations from the visual helper classes to
        # avoid duplicating coordinate math here.
        if self.show_labels and vis_nodes:
            label_points = [vis_node.label_point for vis_node in vis_nodes]
            labels = [vis_node.label for vis_node in vis_nodes]
            self.plotter.add_point_labels(label_points, labels, bold=False, text_color='black', show_points=True, point_color='grey', point_size=5, shape=None, render_points_as_spheres=True)

        if self.show_labels and vis_springs:
            self._spring_label_points = [vis_spring.label_point for vis_spring in vis_springs]
            self._spring_labels = [vis_spring.label for vis_spring in vis_springs]
            self.plotter.add_point_labels(self._spring_label_points, self._spring_labels, text_color='black', bold=False, shape=None, render_points_as_spheres=False)

        if self.show_labels and vis_members:
            label_points = [vis_member.label_point for vis_member in vis_members]
            labels = [vis_member.label for vis_member in vis_members]
            self.plotter.add_point_labels(label_points, labels, bold=False, text_color='black', show_points=False, shape=None, render_points_as_spheres=False)

        # Render the loads if requested
        if (self.combo_name != None or self.case != None) and self.render_loads != False:

            # Plot the loads
            self.plot_loads()

            # Plot the load labels
            self.plotter.add_point_labels(self._load_label_points, self._load_labels, bold=False, text_color='green', show_points=False, shape=None, render_points_as_spheres=False)

        # Render the plates and quads, if present
        if self.model.quads or self.model.plates:
            self.plot_plates(self.deformed_shape, self.deformed_scale, self.color_map, self.combo_name)

        # Render member diagrams if requested
        if self.member_diagrams and (self.combo_name is not None or self.case is not None):
            self.plot_member_diagrams()

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
        """Add a node and its support geometry to the plotter.

        :param Node3D node: Node to visualize (support conditions displayed).
        :param str color: Color for nodes/supports (default ``'grey'``).
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
        """Add a member to the plotter.

        Generates a line representing the structural member between its end nodes.

        :param Member3D member: Structural member to plot.
        :param str theme: Rendering theme (default ``'default'``).
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
        """Add a spring to the plotter.

        Generates a zig-zag line representing the spring between its end nodes.

        :param Spring3D spring: Spring to visualize.
        :param str color: Line color (default ``'grey'``).
        :param bool deformed: Plot deformed shape if ``True`` (default ``False``).
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
        """Add plate/quad elements to the plotter with optional contours.

        :param bool deformed_shape: Plot deformed geometry if ``True``.
        :param float deformed_scale: Scale factor for deformations.
        :param str | None color_map: Contour result type (``None`` disables contours).
        :param str | None combo_name: Load combination name for results.
        """

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
        """Add a node in its deformed position to the plotter.

        :param Node3D node: Node to visualize in deformed state.
        :param float scale_factor: Scale factor for displacements.
        :param str color: Sphere color (default ``'grey'``).
        """

        # Calculate the node's deformed position
        newX = node.X + scale_factor * (node.DX[self.combo_name])
        newY = node.Y + scale_factor * (node.DY[self.combo_name])
        newZ = node.Z + scale_factor * (node.DZ[self.combo_name])

        # Generate a sphere source for the node in its deformed position
        sphere = pv.Sphere(radius=0.4*self.annotation_size, center=[newX, newY, newZ])

        # Add the mesh to the plotter
        self.plotter.add_mesh(sphere, color=color)

    def plot_deformed_member(self, member: Member3D, scale_factor: float) -> None:
        """Add a member in its deformed configuration to the plotter.

        :param Member3D member: Member to visualize in deformed state.
        :param float scale_factor: Scale factor for displacements.
        """

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
        """Add a point-load arrow to the plotter.

        :param tuple position: Arrow tip coordinates ``(X, Y, Z)``.
        :param tuple|ndarray direction: Direction vector (normalized).
        :param float length: Arrow length; sign determines direction.
        :param str|float|int|None label_text: Label text (optional).
        :param str color: Arrow color (default ``'green'``).
        """

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
        """Add a linearly varying distributed load to the plotter.

        :param tuple position1: Start point of the load.
        :param tuple position2: End point of the load.
        :param ndarray|tuple direction: Load direction vector (normalized).
        :param float length1: Arrow length at start (sign indicates direction).
        :param float length2: Arrow length at end (sign indicates direction).
        :param str|float|int|None label_text1: Label at start.
        :param str|float|int|None label_text2: Label at end.
        :param str color: Arrow color (default ``'green'``).
        """

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
        """Add a concentrated moment to the plotter.

        :param tuple center: Center point of the moment arc.
        :param tuple|ndarray direction: Direction vector for the moment axis.
        :param float radius: Arc radius.
        :param str|float|int|None label_text: Label text (optional).
        :param str color: Arc and arrow color (default ``'green'``).
        """

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
        """Add an area load (quad with arrows) to the plotter.

        :param tuple position0: First corner of the loaded area.
        :param tuple position1: Second corner of the loaded area.
        :param tuple position2: Third corner of the loaded area.
        :param tuple position3: Fourth corner of the loaded area.
        :param tuple|ndarray direction: Load direction vector (normalized).
        :param float length: Arrow length; sign determines orientation.
        :param str label_text: Label text shown at one corner.
        :param str color: Polygon and arrow color (default ``'green'``).
        """

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
        """Calculate maximum load magnitudes for normalization.

        :returns: Tuple of (max_pt_load, max_moment, max_dist_load, max_area_load)
            with zero values replaced by 1 to avoid division errors.
        :rtype: tuple[float, float, float, float]
        """

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
        """Add all loads (nodal, member, plate) to the plotter.

        Renders point loads, moments, distributed loads, and area loads
        with appropriate glyphs and labels.
        """

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

    def _get_max_internal_forces(self, diagram_type: str, combo_name: str) -> float:
        """Calculate maximum internal force/moment magnitudes across all members for consistent scaling.

        :param str diagram_type: Diagram type (``'Fy'``, ``'Fz'``, ``'My'``, ``'Mz'``, ``'Fx'``, ``'Tx'``).
        :param str combo_name: Load combination name.
        :returns: Maximum absolute value of the internal force/moment.
        :rtype: float
        """
        
        max_value = 0
        
        # Iterate through all members to find the global maximum
        for member in self.model.members.values():
            
            # Skip inactive members for this combo
            if combo_name not in member.active or not member.active[combo_name]:
                continue
            
            try:
                # Get the maximum value for this member based on diagram type
                if diagram_type == 'Fy':
                    member_max = member.max_shear('Fy', combo_name)
                    member_min = member.min_shear('Fy', combo_name)
                elif diagram_type == 'Fz':
                    member_max = member.max_shear('Fz', combo_name)
                    member_min = member.min_shear('Fz', combo_name)
                elif diagram_type == 'My':
                    member_max = member.max_moment('My', combo_name)
                    member_min = member.min_moment('My', combo_name)
                elif diagram_type == 'Mz':
                    member_max = member.max_moment('Mz', combo_name)
                    member_min = member.min_moment('Mz', combo_name)
                elif diagram_type == 'Fx':
                    member_max = member.max_axial(combo_name)
                    member_min = member.min_axial(combo_name)
                elif diagram_type == 'Tx':
                    member_max = member.max_torque(combo_name)
                    member_min = member.min_torque(combo_name)
                else:
                    continue
                
                # Track the global maximum absolute value
                max_value = max(max_value, abs(member_max), abs(member_min))
            
            except Exception:
                # Skip members that don't have the requested results
                continue
        
        # Prevent division by zero
        if max_value == 0:
            max_value = 1
        
        return max_value

    def plot_member_diagrams(self) -> None:
        """Add internal force/moment diagrams to members.

        Displays member diagrams (Fy, Fz, My, Mz, Fx, Tx) based on the
        current ``member_diagrams`` setting.
        """
        
        # Determine which combo/case to use
        combo_name = self.combo_name if self.combo_name is not None else self.case
        
        # Calculate global maximum internal force/moment for consistent scaling across all members
        global_max = self._get_max_internal_forces(self.member_diagrams, combo_name)
        
        # Create diagrams for each active member
        for member in self.model.members.values():
            
            # Check if member is active for the specified combo
            if combo_name not in member.active or not member.active[combo_name]:
                continue
            
            try:
                # Get member information
                i_node = member.i_node
                j_node = member.j_node
                Xi, Yi, Zi = i_node.X, i_node.Y, i_node.Z
                Xj, Yj, Zj = j_node.X, j_node.Y, j_node.Z
                L = member.L()
                
                # Get transformation matrix for local coordinates
                T = member.T()
                cos_x = np.array([T[0, 0:3]])  # Local x-axis (along member)
                cos_y = np.array([T[1, 0:3]])  # Local y-axis
                cos_z = np.array([T[2, 0:3]])  # Local z-axis
                
                # Member base line
                member_start = np.array([Xi, Yi, Zi])
                member_dir = np.array([Xj - Xi, Yj - Yi, Zj - Zi])
                
                # Determine perpendicular direction for diagram offset
                if self.member_diagrams in ['Fy', 'Mz']:
                    perp_dir = cos_y[0]  # Use y direction for offset
                elif self.member_diagrams in ['Fz', 'My']:
                    perp_dir = cos_z[0]  # Use z direction for offset
                else:
                    perp_dir = cos_y[0]  # Default to y direction
                
                # Get result values at points along member
                n_points = 20
                x_array = np.linspace(0, L, n_points)
                
                if self.member_diagrams == 'Fy':
                    results = member.shear_array('Fy', n_points, combo_name, x_array)
                    y_values = results[1]
                    max_value = member.max_shear('Fy', combo_name)
                    min_value = member.min_shear('Fy', combo_name)
                elif self.member_diagrams == 'Fz':
                    results = member.shear_array('Fz', n_points, combo_name, x_array)
                    y_values = results[1]
                    max_value = member.max_shear('Fz', combo_name)
                    min_value = member.min_shear('Fz', combo_name)
                elif self.member_diagrams == 'My':
                    results = member.moment_array('My', n_points, combo_name, x_array)
                    y_values = results[1]
                    max_value = member.max_moment('My', combo_name)
                    min_value = member.min_moment('My', combo_name)
                elif self.member_diagrams == 'Mz':
                    results = member.moment_array('Mz', n_points, combo_name, x_array)
                    y_values = results[1]
                    max_value = member.max_moment('Mz', combo_name)
                    min_value = member.min_moment('Mz', combo_name)
                elif self.member_diagrams == 'Fx':
                    results = member.axial_array(n_points, combo_name, x_array)
                    y_values = results[1]
                    max_value = member.max_axial(combo_name)
                    min_value = member.min_axial(combo_name)
                elif self.member_diagrams == 'Tx':
                    results = member.torque_array(n_points, combo_name, x_array)
                    y_values = results[1]
                    max_value = member.max_torque(combo_name)
                    min_value = member.min_torque(combo_name)
                else:
                    continue
                
                # Normalize values for better visualization
                # Use global_max for consistent scaling across all members
                if global_max > 0:
                    normalized_values = y_values / global_max
                else:
                    normalized_values = y_values
                
                # Create baseline points
                baseline_points = []
                for i, x in enumerate(x_array):
                    pos_along_member = member_start + (x / L) * member_dir
                    baseline_points.append(pos_along_member)
                
                # Create diagram points (displaced from member axis)
                diagram_points = []
                for i, x in enumerate(x_array):
                    pos_along_member = member_start + (x / L) * member_dir
                    diag_displacement = (normalized_values[i] * self.diagram_scale * 0.5) * perp_dir
                    diagram_pt = pos_along_member + diag_displacement
                    diagram_points.append(diagram_pt)
                
                # Create baseline line
                baseline_pts = np.array(baseline_points)
                baseline_lines = np.zeros((len(baseline_pts)-1, 3), dtype=int)
                baseline_lines[:, 0] = 2
                baseline_lines[:, 1] = np.arange(len(baseline_pts)-1, dtype=int)
                baseline_lines[:, 2] = np.arange(1, len(baseline_pts), dtype=int)
                baseline_polydata = pv.PolyData(baseline_pts, lines=baseline_lines)
                
                # Create diagram line
                diagram_pts = np.array(diagram_points)
                diagram_lines = np.zeros((len(diagram_pts)-1, 3), dtype=int)
                diagram_lines[:, 0] = 2
                diagram_lines[:, 1] = np.arange(len(diagram_pts)-1, dtype=int)
                diagram_lines[:, 2] = np.arange(1, len(diagram_pts), dtype=int)
                diagram_polydata = pv.PolyData(diagram_pts, lines=diagram_lines)
                
                # Create connector lines (vertical lines from baseline to diagram)
                connector_pts_list = []
                connector_lines_list = []
                pt_idx = 0
                for i in range(len(x_array)):
                    connector_pts_list.append(baseline_points[i])
                    connector_pts_list.append(diagram_points[i])
                    connector_lines_list.append([2, pt_idx, pt_idx + 1])
                    pt_idx += 2
                
                if connector_pts_list:
                    connector_pts = np.array(connector_pts_list)
                    connector_lines = np.array(connector_lines_list)
                    connector_polydata = pv.PolyData(connector_pts, lines=connector_lines)
                    
                    # Set color based on theme
                    if self.theme == 'default':
                        color = (0, 1, 1)  # Cyan
                    else:
                        color = (0, 0, 0)  # Black for print
                    
                    # Add all diagram components to plotter
                    self.plotter.add_mesh(baseline_polydata, color=color, line_width=1)
                    self.plotter.add_mesh(diagram_polydata, color=color, line_width=1)
                    self.plotter.add_mesh(connector_polydata, color=color, line_width=1)
                    
                    # Add value labels at max and min points
                    if self.show_labels:
                        max_idx = int(y_values.argmax())
                        min_idx = int(y_values.argmin())
                        
                        # Place labels with offset proportional to annotation_size
                        label_offset = 0.1 * self.annotation_size
                        max_pos = np.array(diagram_points[max_idx]) + label_offset * perp_dir
                        min_pos = np.array(diagram_points[min_idx]) + label_offset * perp_dir
                        
                        max_str = f"{max_value:.3g}"
                        min_str = f"{min_value:.3g}"
                        
                        # Add text labels (show_points defaults to True)
                        self.plotter.add_point_labels([max_pos], [max_str], text_color=color, point_color=color, bold=False, shape=None, render_points_as_spheres=False)
                        self.plotter.add_point_labels([min_pos], [min_str], text_color=color, point_color=color, bold=False, shape=None, render_points_as_spheres=False)
            
            except Exception:
                # Silently skip members that fail to render diagrams
                pass


# === Visualization helper classes ===
# These PyVista-focused classes mirror the VTK helpers used in Visualization.py. They keep
# geometry and label construction encapsulated, which keeps Renderer.update readable and makes
# it easy to swap draw strategies without touching the orchestrating logic.


class VisNode:
    """Visual wrapper for a Node3D using PyVista primitives.

    Encapsulates construction of node and support-condition glyphs for the plotter.
    """

    def __init__(self, node: 'Node3D', annotation_size: float, color: str) -> None:
        """Build visual elements for a node.

        :param Node3D node: Node to visualize.
        :param float annotation_size: Base size for support glyphs.
        :param str color: Color for all glyphs.
        """
        self.node = node
        self.annotation_size = annotation_size
        self.color = color
        self.label = node.name
        self.label_point = [node.X, node.Y, node.Z]
        self.meshes: List[pv.PolyData] = []
        self._build_geometry()

    def _build_geometry(self) -> None:
        """Assemble PyVista meshes for the node and support conditions.

        Creates cubes, cones, or lines depending on the node's restraint flags.
        """
        n = self.node
        s = self.annotation_size

        # Fixed support: represent with a cube for a compact glyph.
        if n.support_DX and n.support_DY and n.support_DZ and n.support_RX and n.support_RY and n.support_RZ:
            self.meshes.append(pv.Cube(center=(n.X, n.Y, n.Z), x_length=s*2, y_length=s*2, z_length=s*2))
            return

        # Pinned support: use a single cone pointed up.
        if n.support_DX and n.support_DY and n.support_DZ and not n.support_RX and not n.support_RY and not n.support_RZ:
            self.meshes.append(pv.Cone(center=(n.X, n.Y - s, n.Z), direction=(0, 1, 0), height=s*2, radius=s*2))
            return

        # For partial restraints, build individual glyphs for each restraint direction/rotation.
        if n.support_DX:
            self.meshes.append(pv.Line((n.X - s, n.Y, n.Z), (n.X + s, n.Y, n.Z)))
            self.meshes.append(pv.Cone(center=(n.X - s, n.Y, n.Z), direction=(1, 0, 0), height=s*0.6, radius=s*0.3))
            self.meshes.append(pv.Cone(center=(n.X + s, n.Y, n.Z), direction=(-1, 0, 0), height=s*0.6, radius=s*0.3))

        if n.support_DY:
            self.meshes.append(pv.Line((n.X, n.Y - s, n.Z), (n.X, n.Y + s, n.Z)))
            self.meshes.append(pv.Cone(center=(n.X, n.Y - s, n.Z), direction=(0, 1, 0), height=s*0.6, radius=s*0.3))
            self.meshes.append(pv.Cone(center=(n.X, n.Y + s, n.Z), direction=(0, -1, 0), height=s*0.6, radius=s*0.3))

        if n.support_DZ:
            self.meshes.append(pv.Line((n.X, n.Y, n.Z - s), (n.X, n.Y, n.Z + s)))
            self.meshes.append(pv.Cone(center=(n.X, n.Y, n.Z - s), direction=(0, 0, 1), height=s*0.6, radius=s*0.3))
            self.meshes.append(pv.Cone(center=(n.X, n.Y, n.Z + s), direction=(0, 0, -1), height=s*0.6, radius=s*0.3))

        if n.support_RX:
            self.meshes.append(pv.Line((n.X - 1.6*s, n.Y, n.Z), (n.X + 1.6*s, n.Y, n.Z)))
            self.meshes.append(pv.Cube(center=(n.X - 1.9*s, n.Y, n.Z), x_length=s*0.6, y_length=s*0.6, z_length=s*0.6))
            self.meshes.append(pv.Cube(center=(n.X + 1.9*s, n.Y, n.Z), x_length=s*0.6, y_length=s*0.6, z_length=s*0.6))

        if n.support_RY:
            self.meshes.append(pv.Line((n.X, n.Y - 1.6*s, n.Z), (n.X, n.Y + 1.6*s, n.Z)))
            self.meshes.append(pv.Cube(center=(n.X, n.Y - 1.9*s, n.Z), x_length=s*0.6, y_length=s*0.6, z_length=s*0.6))
            self.meshes.append(pv.Cube(center=(n.X, n.Y + 1.9*s, n.Z), x_length=s*0.6, y_length=s*0.6, z_length=s*0.6))

        if n.support_RZ:
            self.meshes.append(pv.Line((n.X, n.Y, n.Z - 1.6*s), (n.X, n.Y, n.Z + 1.6*s)))
            self.meshes.append(pv.Cube(center=(n.X, n.Y, n.Z - 1.9*s), x_length=s*0.6, y_length=s*0.6, z_length=s*0.6))
            self.meshes.append(pv.Cube(center=(n.X, n.Y, n.Z + 1.9*s), x_length=s*0.6, y_length=s*0.6, z_length=s*0.6))

        # If no supports are present, add a small sphere so the node is still visible.
        if not self.meshes:
            self.meshes.append(pv.Sphere(center=(n.X, n.Y, n.Z), radius=0.4*s))

    def add_to_plotter(self, plotter: pv.Plotter) -> None:
        """Add all prepared meshes to the plotter.

        :param pv.Plotter plotter: Plotter to receive the meshes.
        """
        for mesh in self.meshes:
            plotter.add_mesh(mesh, color=self.color)


class VisSpring:
    """Visual wrapper for a Spring3D with optional deformed rendering."""

    def __init__(self, spring: 'Spring3D', annotation_size: float, color: str, deformed: bool, scale: float, combo_name: Optional[str]) -> None:
        """Build visual elements for a spring.

        :param Spring3D spring: Spring to visualize.
        :param float annotation_size: Scale for zig-zag amplitude.
        :param str color: Line color.
        :param bool deformed: If ``True`` use deformed coordinates.
        :param float scale: Deformation scale factor.
        :param str | None combo_name: Load combination for displacements.
        """
        self.spring = spring
        self.annotation_size = annotation_size
        self.color = color
        self.deformed = deformed
        self.scale = scale
        self.combo_name = combo_name
        self.mesh: Optional[pv.PolyData] = None
        self.label = spring.name
        self.label_point: Optional[List[float]] = None
        self._build_geometry()

    def _build_geometry(self) -> None:
        i_node = self.spring.i_node
        j_node = self.spring.j_node

        # Pull coordinates, optionally in deformed configuration.
        Xi, Yi, Zi = i_node.X, i_node.Y, i_node.Z
        Xj, Yj, Zj = j_node.X, j_node.Y, j_node.Z
        if self.deformed and self.combo_name is not None:
            Xi += i_node.DX[self.combo_name]*self.scale
            Yi += i_node.DY[self.combo_name]*self.scale
            Zi += i_node.DZ[self.combo_name]*self.scale
            Xj += j_node.DX[self.combo_name]*self.scale
            Yj += j_node.DY[self.combo_name]*self.scale
            Zj += j_node.DZ[self.combo_name]*self.scale

        direction = np.array([Xj, Yj, Zj]) - np.array([Xi, Yi, Zi])
        length = np.linalg.norm(direction)
        if length == 0:
            # Avoid degenerate geometry; skip drawing this spring.
            return

        direction = direction / length

        # Build a zig-zag polyline to suggest a spring.
        arbitrary_vector = np.array([1, 0, 0])
        if np.allclose(direction, arbitrary_vector) or np.allclose(direction, -arbitrary_vector):
            arbitrary_vector = np.array([0, 1, 0])
        perp_vector1 = np.cross(direction, arbitrary_vector)
        perp_vector1 /= np.linalg.norm(perp_vector1)
        perp_vector2 = np.cross(direction, perp_vector1)
        perp_vector2 /= np.linalg.norm(perp_vector2)

        straight_segment_length = length / 10
        zigzag_length = length - 2 * straight_segment_length

        num_zigs = 4
        num_points = num_zigs * 2
        amplitude = self.annotation_size
        t = np.linspace(0, zigzag_length, num_points)
        zigzag_pattern = amplitude * np.tile([1, -1], num_zigs)
        zigzag_points = np.outer(t, direction) + np.outer(zigzag_pattern, perp_vector1)

        start_point = np.array([Xi, Yi, Zi])
        end_point = np.array([Xj, Yj, Zj])
        start_segment = start_point + direction * straight_segment_length
        end_segment = end_point - direction * straight_segment_length
        zigzag_points += start_segment

        points = np.vstack([start_point, start_segment, zigzag_points, end_segment, end_point])
        lines = np.zeros((len(points) - 1, 3), dtype=int)
        lines[:, 0] = 2
        lines[:, 1] = np.arange(len(points) - 1, dtype=int)
        lines[:, 2] = np.arange(1, len(points), dtype=int)

        self.mesh = pv.PolyData(points, lines=lines)
        self.label_point = [(Xi + Xj) / 2, (Yi + Yj) / 2, (Zi + Zj) / 2]

    def add_to_plotter(self, plotter: pv.Plotter) -> None:
        """Add the spring mesh to the plotter.

        :param pv.Plotter plotter: Plotter to receive the mesh.
        """
        if self.mesh is not None:
            plotter.add_mesh(self.mesh, color=self.color, line_width=2)


class VisMember:
    """Visual wrapper for a Member3D as a simple line."""

    def __init__(self, member: 'Member3D', theme: str) -> None:
        """Build visual elements for a member.

        :param Member3D member: Member to visualize.
        :param str theme: Rendering theme (``'default'`` or ``'print'``).
        """
        self.member = member
        self.theme = theme
        self.mesh: pv.PolyData = self._build_geometry()
        self.label = member.name
        self.label_point = [(member.i_node.X + member.j_node.X) / 2,
                            (member.i_node.Y + member.j_node.Y) / 2,
                            (member.i_node.Z + member.j_node.Z) / 2]

    def _build_geometry(self) -> pv.PolyData:
        line = pv.Line()
        line.points[0] = [self.member.i_node.X, self.member.i_node.Y, self.member.i_node.Z]
        line.points[1] = [self.member.j_node.X, self.member.j_node.Y, self.member.j_node.Z]
        return line

    def add_to_plotter(self, plotter: pv.Plotter) -> None:
        """Add the member line to the plotter.

        :param pv.Plotter plotter: Plotter to receive the mesh.
        """
        color = 'black' if self.theme == 'print' else 'black'
        plotter.add_mesh(self.mesh, color=color, line_width=2)


class VisDeformedMember:
    """Visual wrapper for a deformed Member3D polyline."""

    def __init__(self, member: 'Member3D', scale_factor: float, combo_name: Optional[str]) -> None:
        """Build a deformed member polyline.

        :param Member3D member: Member to visualize in deformed state.
        :param float scale_factor: Scale factor for displacements.
        :param str | None combo_name: Load combination name for results.
        """
        self.member = member
        self.scale_factor = scale_factor
        self.combo_name = combo_name
        self.mesh: Optional[pv.PolyData] = None
        self._build_geometry()

    def _build_geometry(self) -> None:
        if self.combo_name is None:
            return
        if not self.member.active:
            return

        L = self.member.L()
        T = self.member.T()
        cos_x = np.array([T[0, 0:3]])
        cos_y = np.array([T[1, 0:3]])
        cos_z = np.array([T[2, 0:3]])

        Xi = self.member.i_node.X
        Yi = self.member.i_node.Y
        Zi = self.member.i_node.Z

        DY_plot = np.empty((0, 3))
        for i in range(20):
            dy_tot = self.member.deflection('dy', L / 19 * i, self.combo_name)
            DY_plot = np.append(DY_plot, dy_tot * cos_y * self.scale_factor, axis=0)

        DZ_plot = np.empty((0, 3))
        for i in range(20):
            dz_tot = self.member.deflection('dz', L / 19 * i, self.combo_name)
            DZ_plot = np.append(DZ_plot, dz_tot * cos_z * self.scale_factor, axis=0)

        DX_plot = np.empty((0, 3))
        for i in range(20):
            dx_tot = [[Xi, Yi, Zi]] + (L / 19 * i + self.member.deflection('dx', L / 19 * i, self.combo_name) * self.scale_factor) * cos_x
            DX_plot = np.append(DX_plot, dx_tot, axis=0)

        D_plot = DY_plot + DZ_plot + DX_plot
        if len(D_plot) == 0:
            return

        self.mesh = pv.lines_from_points(D_plot, close=False)

    def add_to_plotter(self, plotter: pv.Plotter) -> None:
        """Add the deformed member line to the plotter.

        :param pv.Plotter plotter: Plotter to receive the mesh.
        """
        if self.mesh is not None:
            plotter.add_mesh(self.mesh, color='red', line_width=2)

def _PerpVector(v):
    """Return a unit vector perpendicular to ``v``.

    :param v: Input vector ``[i, j, k]``.
    :type v: array-like
    :returns: Unit vector perpendicular to ``v``.
    :rtype: list
    """

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
    """Populate nodal ``contour`` values for plate/quad results.

    :param FEModel3D model: Model whose plates/quads are processed.
    :param str stress_type: Contour type (``'Mx'``, ``'My'``, ``'Mxy'``,
        ``'Qx'``, ``'Qy'``, ``'Sx'``, ``'Sy'``, ``'Txy'``, or ``'dz'``).
    :param str combo_name: Load combination name (default ``'Combo 1'``).
    """
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
