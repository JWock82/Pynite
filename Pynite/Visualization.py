"""VTK-based visualization helpers for PyNite models.

This module provides a collection of helper classes and functions used by
``Renderer`` to display undeformed/deformed geometry, loads, contours, and
member diagrams. The docstrings follow reStructuredText/Numpy style so they can
be rendered cleanly by Sphinx/ReadTheDocs.
"""

from __future__ import annotations # Allows more recent type hints features
import warnings

from IPython.display import Image
from numpy import array, empty, append, cross, asarray, linspace
from numpy.linalg import norm
import vtk

class Renderer():
    """Renderer for visualizing finite element models using VTK.

    This class provides a flexible interface for rendering 3D finite element
    models (nodes, members, springs, plates, quads). It can display undeformed
    and deformed geometry, load glyphs, contours, and internal force diagrams.

    :param FEModel3D model: Finite element model to be rendered.
    :ivar float annotation_size: Text/marker scale (auto-computed when ``None``).
    :ivar bool deformed_shape: Toggle deformed-shape display.
    :ivar float deformed_scale: Scale factor for deformations.
    :ivar bool render_nodes: Toggle node rendering.
    :ivar bool render_loads: Toggle load glyph rendering.
    :ivar str | None color_map: Plate/quad contour result key (``'Qx'``, ``'Mx'``, etc.).
    :ivar str | None combo_name: Load combination name (exclusive with ``case``).
    :ivar str | None case: Load case name (exclusive with ``combo_name``).
    :ivar bool labels: Toggle text labels for elements.
    :ivar bool scalar_bar: Toggle scalar bar for contours.
    :ivar int scalar_bar_text_size: Font size for scalar bar labels.
    :ivar str theme: ``'default'`` or ``'print'`` color theme.
    :ivar tuple[int, int] window_size: Render window dimensions.
    """

    scalar = None

    def __init__(self, model):
        """Initialize the renderer with a finite element model.

        :param FEModel3D model: The finite element model to render.
        """
        self.model = model

        # Default settings for rendering
        self._annotation_size = None  # None means auto-calculate
        self._annotation_size_manual = False  # Track if user manually set the size
        self._annotation_size_cached = None  # Cache to avoid recalculating
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
        self._theme = 'default'
        self._show_load_info = True
        self._member_diagrams = None  # Options: None, 'Fy', 'Fz', 'My', 'Mz', 'Fx', 'Tx'
        self._diagram_scale = 30  # Scale factor for diagram visualization

        # Initialize VTK objects
        self.renderer = vtk.vtkRenderer()
        self.window = vtk.vtkRenderWindow()
        self.window.SetWindowName('Pynite - Simple Finite Element Analysis in Python')
        self.window.AddRenderer(self.renderer)

    @property
    def window_size(self):
        """Window size as (width, height) tuple."""
        return self.window.GetSize()

    @window_size.setter
    def window_size(self, size):
        """Set window size from tuple or list (width, height)."""
        width, height = size
        self.window.SetSize(width, height)

    @property
    def combo_name(self):
        """Load combination name. When set, case is automatically set to None."""
        return self._combo_name

    @combo_name.setter
    def combo_name(self, value):
        self._combo_name = value
        if value is not None:
            self._case = None

    @property
    def case(self):
        """Load case name. When set, combo_name is automatically set to None."""
        return self._case

    @case.setter
    def case(self, value):
        self._case = value
        if value is not None:
            self._combo_name = None

    @property
    def annotation_size(self):
        """Size of text annotations and visual elements in model units.
        
        If not manually set, automatically calculates as 5% of the shortest
        distance between nodes in the model.
        """
        if self._annotation_size is None or not self._annotation_size_manual:
            # Return cached value if available; calculate once on first access
            if self._annotation_size_cached is None:
                self._annotation_size_cached = self._calculate_auto_annotation_size()
            return self._annotation_size_cached
        return self._annotation_size

    @annotation_size.setter
    def annotation_size(self, value):
        self._annotation_size = value
        self._annotation_size_manual = True  # Mark as manually set

    @property
    def deformed_shape(self):
        """Whether to render the deformed shape of the model."""
        return self._deformed_shape

    @deformed_shape.setter
    def deformed_shape(self, value):
        self._deformed_shape = value

    @property
    def deformed_scale(self):
        """Scale factor for deformation visualization."""
        return self._deformed_scale

    @deformed_scale.setter
    def deformed_scale(self, value):
        self._deformed_scale = value

    @property
    def render_nodes(self):
        """Whether to render nodes in the visualization."""
        return self._render_nodes

    @render_nodes.setter
    def render_nodes(self, value):
        self._render_nodes = value

    @property
    def render_loads(self):
        """Whether to render applied loads."""
        return self._render_loads

    @render_loads.setter
    def render_loads(self, value):
        self._render_loads = value

    @property
    def color_map(self):
        """Type of stress/force contour to display on plates/quads.

        Valid options: 'Qx', 'Qy', 'Mx', 'My', 'Mxy', 'Sx', 'Sy', 'Txy'
        - Qx, Qy: Out-of-plane shear forces
        - Mx, My, Mxy: Local out-of-plane bending moments
        - Sx, Sy: Membrane forces
        - Txy: In-plane shear force
        """
        return self._color_map

    @color_map.setter
    def color_map(self, value):
        self._color_map = value

    @property
    def labels(self):
        """Whether to display text labels for elements."""
        return self._labels

    @labels.setter
    def labels(self, value):
        self._labels = value

    @property
    def scalar_bar(self):
        """Whether to display a scalar bar legend for contour plots."""
        return self._scalar_bar

    @scalar_bar.setter
    def scalar_bar(self, value):
        self._scalar_bar = value

    @property
    def scalar_bar_text_size(self):
        """Font size for scalar bar text."""
        return self._scalar_bar_text_size

    @scalar_bar_text_size.setter
    def scalar_bar_text_size(self, value):
        self._scalar_bar_text_size = value

    @property
    def theme(self):
        """Visual theme: 'default' (dark background) or 'print' (white background)."""
        return self._theme

    @theme.setter
    def theme(self, value):
        self._theme = value

    @property
    def show_load_info(self):
        """Whether to display load case/combo information in top left corner."""
        return self._show_load_info

    @show_load_info.setter
    def show_load_info(self, value):
        self._show_load_info = value

    @property
    def member_diagrams(self):
        """Type of internal force/moment diagram to display on members.

        Valid options:
        - None: No diagrams (default)
        - 'Fy': Shear force in local y-direction
        - 'Fz': Shear force in local z-direction
        - 'My': Bending moment about local y-axis
        - 'Mz': Bending moment about local z-axis
        - 'Fx': Axial force along member
        - 'Tx': Torsional moment about member axis
        """
        return self._member_diagrams

    @member_diagrams.setter
    def member_diagrams(self, value):
        valid_options = [None, 'Fy', 'Fz', 'My', 'Mz', 'Fx', 'Tx']
        if value not in valid_options:
            raise ValueError(f"member_diagrams must be one of {valid_options}, got {value}")
        self._member_diagrams = value

    @property
    def diagram_scale(self):
        """Scale factor for internal force/moment diagram visualization."""
        return self._diagram_scale

    @diagram_scale.setter
    def diagram_scale(self, value):
        self._diagram_scale = value

    def _calculate_auto_annotation_size(self):
        """Calculate automatic annotation size as 5% of shortest node distance.

        Uses vectorized NumPy operations for fast computation on large meshes.
        Result is cached by the annotation_size property to avoid recalculation.

        :returns: Annotation size in model units (``5.0`` fallback if <2 nodes
            or co-located nodes).
        :rtype: float
        """
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

    def render_model(self, interact=True, reset_camera=True):
        """Render the model in a window.

        :param bool interact: If ``False`` suppresses the VTK interactor (useful
            for scripted screenshots). Defaults to ``True``.
        :param bool reset_camera: If ``True`` resets the camera before rendering
            (default ``True``).
        :returns: The VTK render window.
        :rtype: vtk.vtkRenderWindow
        """

        # Get the render window
        window = self.window

        # Update the renderer
        self.update(reset_camera)

        # Render the window
        window.Render()

        # Handle user interaction if requested by the user
        if interact:

            # Set up an interactor. The interactor style determines how user interactions affect the
            # view. The trackball camera style behaves much like popular commercial CAD programs.
            interactor = vtk.vtkRenderWindowInteractor()
            style = vtk.vtkInteractorStyleTrackballCamera()
            interactor.SetInteractorStyle(style)
            interactor.SetRenderWindow(self.window)

            # Add coordinate axes in the bottom left corner
            axes = vtk.vtkAxesActor()  # Create a 3D axes actor showing X, Y, and Z axes
            axes.SetTotalLength(1.5, 1.5, 1.5)  # Set the length of each axis arrow to 1.5 units
            axes.SetShaftTypeToLine()  # Use simple lines instead of cylinders for the axis shafts (cleaner appearance)
            axes.SetAxisLabels(1)  # Enable the display of axis labels (X, Y, Z)
            axes.SetCylinderRadius(0.02)  # Set the radius of the arrow tips to 0.02 units

            # Create the orientation marker widget to display the axes in the corner
            axes_widget = vtk.vtkOrientationMarkerWidget()  # Create a widget to display the axes actor
            axes_widget.SetOrientationMarker(axes)  # Attach the axes actor to the widget
            axes_widget.SetInteractor(interactor)  # Connect the widget to the render window interactor
            axes_widget.SetViewport(0.0, 0.0, 0.2, 0.2)  # Position in bottom left corner (x_min, y_min, x_max, y_max as fractions of window)
            axes_widget.SetEnabled(1)  # Enable the widget so it displays in the render window
            axes_widget.InteractiveOff()  # Disable interaction with the axes widget (it will only rotate with the camera)

            # Start the interactor. Code execution will pause here until the user closes the window.
            interactor.Start()

            # Finalize the render window once the user closes out of it. I don't understand everything
            # this does, but I've found screenshots will cause the program to crash if this line is
            # omitted. I have noticed it will shut down the interactor.
            window.Finalize()

        return window

    def screenshot(self, filepath='console', interact=True, reset_camera=True):
        """Render the model and capture a screenshot.

        :param str filepath: Destination for the PNG. ``'console'`` returns an
            IPython ``Image``; ``'BytesIO'`` returns a ``BytesIO`` object; any
            other value is treated as a filesystem path. Defaults to ``'console'``.
        :param bool interact: When ``True`` (default) start the VTK interactor
            so the user can orbit/zoom before closing the window.
        :param bool reset_camera: When ``True`` (default) resets the camera
            before rendering.
        :returns: ``Image`` when ``filepath='console'``, ``BytesIO`` when
            ``filepath='BytesIO'``, otherwise ``None`` after writing to disk.
        :rtype: IPython.display.Image | io.BytesIO | None
        """

        # Render the model in a window and save the window
        window = self.render_model(interact, reset_camera)

        # Screenshot code
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(window)
        w2if.SetInputBufferTypeToRGB()
        w2if.ReadFrontBufferOff()

        # These next two lines are in the examples and documentation for VTK, but don't seem to do
        # anything. I've left them here in case I find a bug somewhere down the line that needs
        # fixing.
        # w2if.Update()
        # w2if.Modified()

        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(w2if.GetOutputPort())

        if filepath == 'console' or filepath == 'BytesIO':

            writer.SetWriteToMemory(1)
            writer.Write()
            fig_file = memoryview(writer.GetResult()).tobytes()

            # Now that we're done with the render window, finalize it
            window.Finalize()

            if filepath == 'console':
                return Image(fig_file)
            elif filepath == 'BytesIO':
                from io import BytesIO
                return BytesIO(fig_file)
        else:

            writer.SetFileName(filepath)
            writer.Write()

            # Now that we're done with the render window, finalize it
            window.Finalize()

            return

    def update(self, reset_camera=True):
        """Rebuild the VTK renderer with current settings.

        Updates deformed shapes, loads, contours, labels, and diagrams. Called
        automatically by ``render_model`` and ``screenshot``.

        :param bool reset_camera: If ``True`` resets the camera to fit the model
            (default ``True``).
        :raises Exception: If invalid configuration is detected.
        """

        # Clear annotation size cache to recalculate if model has changed
        self._annotation_size_cached = None

        # Input validation
        if self.deformed_shape and self.case != None:
            self.deformed_shape = False
            warnings.warn('Deformed shape is only available for load combinations. Deformed shape will not be rendered for load cases.', UserWarning)
        if self.model.load_combos == {} and self.render_loads == True and self.case == None:
            self.render_loads = False
            warnings.warn('Unable to render load combination. No load combinations defined.', UserWarning)

        # Check if nodes are to be rendered
        if self.render_nodes == True:

            # Create a visual node for each node in the model
            vis_nodes = []
            for node in self.model.nodes.values():
                vis_nodes.append(VisNode(node, self.annotation_size))

        # Create a visual spring for each spring in the model
        vis_springs = []
        for spring in self.model.springs.values():
            vis_springs.append(VisSpring(spring, self.model.nodes, self.annotation_size))

        # Create a visual member for each member in the model
        vis_members = []
        for member in self.model.members.values():
            vis_members.append(VisMember(member, self.model.nodes, self.annotation_size, self.theme))

        # Get the renderer
        renderer = self.renderer

        # Clear out all the old actors from any previous renderings
        for actor in renderer.GetActors():
            renderer.RemoveActor(actor)

        # Clear out all the old view props (2D annotations, etc.)
        for prop in renderer.GetViewProps():
            renderer.RemoveViewProp(prop)

        # Add actors for each spring
        for vis_spring in vis_springs:

            # Add the actor for the spring
            renderer.AddActor(vis_spring.actor)

            if self.labels == True:
                # Add the actor for the spring label
                renderer.AddActor(vis_spring.lblActor)

                # Set the text to follow the camera as the user interacts. This will
                # require a reset of the camera (see below)
                vis_spring.lblActor.SetCamera(renderer.GetActiveCamera())

        # Add actors for each member
        for vis_member in vis_members:

            # Add the actor for the member
            renderer.AddActor(vis_member.actor)

            if self.labels == True:

                # Add the actor for the member label
                renderer.AddActor(vis_member.lblActor)

                # Set the text to follow the camera as the user interacts. This will
                # require a reset of the camera (see below)
                vis_member.lblActor.SetCamera(renderer.GetActiveCamera())

        # Check if nodes are to be rendered
        if self.render_nodes == True:

            # Combine the polydata from each node

            # Create an append filter for combining node polydata
            node_polydata = vtk.vtkAppendPolyData()

            for vis_node in vis_nodes:

                # Add the node's polydata
                node_polydata.AddInputData(vis_node.polydata.GetOutput())

                if self.labels == True:

                    if self.theme == 'print':

                        # Adjust the node label's color
                        vis_node.lblActor.GetProperty().SetColor(0, 0, 1)  # Blue

                    # Add the actor for the node label
                    renderer.AddActor(vis_node.lblActor)

                    # Set the text to follow the camera as the user interacts. This will
                    # require a reset of the camera (see below)
                    vis_node.lblActor.SetCamera(renderer.GetActiveCamera())

            # Update the node polydata in the append filter
            node_polydata.Update()

            # Create a mapper and actor for the nodes
            node_mapper = vtk.vtkPolyDataMapper()
            node_mapper.SetInputConnection(node_polydata.GetOutputPort())
            node_actor = vtk.vtkActor()
            node_actor.SetMapper(node_mapper)

            # Adjust the color of all the nodes.
            if self.theme == 'print':
                node_actor.GetProperty().SetColor(0, 0, 1)  # Blue

            # Add the node actor to the renderer
            renderer.AddActor(node_actor)

        # Render the deformed shape if requested
        if self.deformed_shape == True:
            _DeformedShape(self.model, renderer, self.deformed_scale, self.annotation_size, self.combo_name, self.render_nodes, self.theme)

        # Render the loads if requested
        if (self.combo_name != None or self.case != None) and self.render_loads != False:
            _RenderLoads(self.model, renderer, self.annotation_size, self.combo_name, self.case, self.theme)

        # Render the plates and quads, if present
        if self.model.quads or self.model.plates:
            _RenderContours(self.model, renderer, self.deformed_shape, self.deformed_scale,
                            self.color_map, self.scalar_bar, self.scalar_bar_text_size,
                            self.combo_name, self.theme)

        # Set the window's background color
        if self.theme == 'default':
            renderer.SetBackground(0, 0, 0.5)  # Blue
        elif self.theme == 'print':
            renderer.SetBackground(1, 1, 1)  # White

        # Add text overlay in the top left corner showing load case/combo information
        if self.show_load_info:

            # Determine which load case or combo is being displayed
            text_str = None
            if self.case is not None:
                text_str = f"Load Case: {self.case}"
            elif self.combo_name is not None:
                # Check if this is a modal combination
                combo = self.model.load_combos.get(self.combo_name)
                if combo and combo.combo_tags and 'modal' in combo.combo_tags:
                    # Extract mode number and get frequency
                    mode_num = int(self.combo_name.split()[1]) - 1  # "Mode 3" -> index 2
                    if hasattr(self.model, 'frequencies') and mode_num < len(self.model.frequencies):
                        freq = self.model.frequencies[mode_num]
                        text_str = f"Load Combo: {self.combo_name} - {freq:.3f} Hz"
                    else:
                        text_str = f"Load Combo: {self.combo_name}"
                else:
                    text_str = f"Load Combo: {self.combo_name}"

            # Create the corner annotation if there's text to display
            if text_str:
                corner_annotation = vtk.vtkCornerAnnotation()
                corner_annotation.SetText(2, text_str)  # 2 = upper left corner
                corner_annotation.GetTextProperty().SetFontSize(12)
                corner_annotation.GetTextProperty().SetFontFamilyToArial()

                # Set color based on theme
                if self.theme == 'print':
                    corner_annotation.GetTextProperty().SetColor(0, 0, 0)  # Black text for print theme
                else:
                    corner_annotation.GetTextProperty().SetColor(1, 1, 1)  # White text for default theme

                # Add the corner annotation to the renderer
                renderer.AddViewProp(corner_annotation)

        # Render member internal force/moment diagrams if requested
        if self.member_diagrams is not None:
            _RenderMemberDiagrams(self.model, renderer, self.member_diagrams,
                                 self.diagram_scale, self.combo_name, self.case, self.theme, self.annotation_size)

        # Reset the camera
        if reset_camera: renderer.ResetCamera()


# Converts a node object into a node for the viewer
class VisNode():
    """VTK representation of a node and its supports."""

    # Constructor
    def __init__(self, node, annotation_size=5):
        """Create the VTK actors for a node.

        :param Node3D node: Node to visualize (support flags are displayed).
        :param float annotation_size: Base scale for spheres, cones, cubes, and
            label size (default ``5``).
        """

        # Create an append filter to append all the sources related to the node into a single 'PolyData' object
        self.polydata = vtk.vtkAppendPolyData()

        # Get the node's position
        X = node.X  # Global X coordinate
        Y = node.Y  # Global Y coordinate
        Z = node.Z  # Global Z coordinate

        # Generate a sphere source for the node
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(X, Y, Z)
        sphere.SetRadius(0.6*annotation_size)
        sphere.Update()
        self.polydata.AddInputData(sphere.GetOutput())

        # Create the text for the node label
        label = vtk.vtkVectorText()
        label.SetText(node.name)

        # Set up a mapper for the node label
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())

        # Set up an actor for the node label
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition(X + 0.6*annotation_size, Y + 0.6*annotation_size, Z)

        # Generate any supports that occur at the node
        # Check for a fixed suppport
        if (node.support_DX == True and node.support_DY == True and node.support_DZ == True and node.support_RX == True and node.support_RY == True and node.support_RZ == True):

            # Create the fixed support
            support = vtk.vtkCubeSource()
            support.SetCenter(node.X, node.Y, node.Z)
            support.SetXLength(annotation_size*1.2)
            support.SetYLength(annotation_size*1.2)
            support.SetZLength(annotation_size*1.2)

            # Copy and append the support data to the append filter
            support.Update()
            self.polydata.AddInputData(support.GetOutput())

        # Check for a pinned support
        elif node.support_DX == True and node.support_DY == True and node.support_DZ == True \
        and node.support_RX == False and node.support_RY == False and node.support_RZ == False:

            # Create the pinned support
            support = vtk.vtkConeSource()
            support.SetCenter(node.X, node.Y-0.6*annotation_size, node.Z)
            support.SetDirection((0, 1, 0))
            support.SetHeight(annotation_size*1.2)
            support.SetRadius(annotation_size*1.2)

            # Copy and append the support data to the append filter
            support.Update()
            self.polydata.AddInputData(support.GetOutput())

        # Other support conditions
        else:

            # Restrained against X translation
            if node.support_DX == True:

                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X-annotation_size, node.Y, node.Z)
                support1.SetPoint2(node.X+annotation_size, node.Y, node.Z)

                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())

                support2 = vtk.vtkConeSource()
                support2.SetCenter(node.X-annotation_size, node.Y, node.Z)
                support2.SetDirection((1, 0, 0))
                support2.SetHeight(annotation_size*0.6)
                support2.SetRadius(annotation_size*0.3)

                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())

                support3 = vtk.vtkConeSource()
                support3.SetCenter(node.X+annotation_size, node.Y, node.Z)
                support3.SetDirection((-1, 0, 0))
                support3.SetHeight(annotation_size*0.6)
                support3.SetRadius(annotation_size*0.3)

                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())

            # Restrained against Y translation
            if node.support_DY == True:

                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y-annotation_size, node.Z)
                support1.SetPoint2(node.X, node.Y+annotation_size, node.Z)

                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())

                support2 = vtk.vtkConeSource()
                support2.SetCenter(node.X, node.Y-annotation_size, node.Z)
                support2.SetDirection((0, 1, 0))
                support2.SetHeight(annotation_size*0.6)
                support2.SetRadius(annotation_size*0.3)

                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())

                support3 = vtk.vtkConeSource()
                support3.SetCenter(node.X, node.Y+annotation_size, node.Z)
                support3.SetDirection((0, -1, 0))
                support3.SetHeight(annotation_size*0.6)
                support3.SetRadius(annotation_size*0.3)

                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())

            # Restrained against Z translation
            if node.support_DZ == True:

                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y, node.Z-annotation_size)
                support1.SetPoint2(node.X, node.Y, node.Z+annotation_size)

                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())

                support2 = vtk.vtkConeSource()
                support2.SetCenter(node.X, node.Y, node.Z-annotation_size)
                support2.SetDirection((0, 0, 1))
                support2.SetHeight(annotation_size*0.6)
                support2.SetRadius(annotation_size*0.3)

                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())

                support3 = vtk.vtkConeSource()
                support3.SetCenter(node.X, node.Y, node.Z+annotation_size)
                support3.SetDirection((0, 0, -1))
                support3.SetHeight(annotation_size*0.6)
                support3.SetRadius(annotation_size*0.3)

                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())

            # Restrained against rotation about the X-axis
            if node.support_RX == True:

                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X-1.6*annotation_size, node.Y, node.Z)
                support1.SetPoint2(node.X+1.6*annotation_size, node.Y, node.Z)

                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())

                support2 = vtk.vtkCubeSource()
                support2.SetCenter(node.X-1.9*annotation_size, node.Y, node.Z)
                support2.SetXLength(annotation_size*0.6)
                support2.SetYLength(annotation_size*0.6)
                support2.SetZLength(annotation_size*0.6)

                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())

                support3 = vtk.vtkCubeSource()
                support3.SetCenter(node.X+1.9*annotation_size, node.Y, node.Z)
                support3.SetXLength(annotation_size*0.6)
                support3.SetYLength(annotation_size*0.6)
                support3.SetZLength(annotation_size*0.6)

                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())

            # Restrained against rotation about the Y-axis
            if node.support_RY == True:

                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y-1.6*annotation_size, node.Z)
                support1.SetPoint2(node.X, node.Y+1.6*annotation_size, node.Z)

                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())

                support2 = vtk.vtkCubeSource()
                support2.SetCenter(node.X, node.Y-1.9*annotation_size, node.Z)
                support2.SetXLength(annotation_size*0.6)
                support2.SetYLength(annotation_size*0.6)
                support2.SetZLength(annotation_size*0.6)

                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())

                support3 = vtk.vtkCubeSource()
                support3.SetCenter(node.X, node.Y+1.9*annotation_size, node.Z)
                support3.SetXLength(annotation_size*0.6)
                support3.SetYLength(annotation_size*0.6)
                support3.SetZLength(annotation_size*0.6)

                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())

            # Restrained against rotation about the Z-axis
            if node.support_RZ == True:

                # Create the support
                support1 = vtk.vtkLineSource()  # The line showing the support direction
                support1.SetPoint1(node.X, node.Y, node.Z-1.6*annotation_size)
                support1.SetPoint2(node.X, node.Y, node.Z+1.6*annotation_size)

                # Copy and append the support data to the append filter
                support1.Update()
                self.polydata.AddInputData(support1.GetOutput())

                support2 = vtk.vtkCubeSource()
                support2.SetCenter(node.X, node.Y, node.Z-1.9*annotation_size)
                support2.SetXLength(annotation_size*0.6)
                support2.SetYLength(annotation_size*0.6)
                support2.SetZLength(annotation_size*0.6)

                # Copy and append the support data to the append filter
                support2.Update()
                self.polydata.AddInputData(support2.GetOutput())

                support3 = vtk.vtkCubeSource()
                support3.SetCenter(node.X, node.Y, node.Z+1.9*annotation_size)
                support3.SetXLength(annotation_size*0.6)
                support3.SetYLength(annotation_size*0.6)
                support3.SetZLength(annotation_size*0.6)

                # Copy and append the support data to the append filter
                support3.Update()
                self.polydata.AddInputData(support3.GetOutput())

        # Update the append filter
        self.polydata.Update()

        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.polydata.GetOutputPort())
        self.actor = vtk.vtkActor()

        # Set the mapper for the node's actor
        self.actor.SetMapper(mapper)


class VisSpring():
    """VTK representation of a spring element."""

    def __init__(self, spring, nodes, annotation_size=5, color=None):
        """Create the VTK actors for a spring.

        :param Spring3D spring: Spring instance to visualize.
        :param dict[str, Node3D] nodes: Model nodes used to locate the i/j ends.
        :param float annotation_size: Base size for labels and arrow heads
            (default ``5``).
        :param str | None color: ``'black'`` for black lines; ``None`` uses magenta.
        """

        # Generate a line source for the spring
        line = vtk.vtkLineSource()

        # Step through each node in the model and find the position of the
        # i-node and j-node
        for node in nodes.values():

            # Check to see if the current node is the i-node
            if node.name == spring.i_node.name:
                Xi = node.X
                Yi = node.Y
                Zi = node.Z
                line.SetPoint1(Xi, Yi, Zi)

            # Check to see if the current node is the j-node
            elif node.name == spring.j_node.name:
                Xj = node.X
                Yj = node.Y
                Zj = node.Z
                line.SetPoint2(Xj, Yj, Zj)

        # Set up a mapper for the spring
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(line.GetOutputPort())

        # Set up an actor and a mapper for the spring
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)

        # Create the text for the spring label
        label = vtk.vtkVectorText()
        label.SetText(spring.name)

        # Set up a mapper for the spring label
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())

        # Set up an actor for the spring label
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition((Xi+Xj)/2, (Yi+Yj)/2, (Zi+Zj)/2)

        # Add some color
        if color is None:
            self.actor.GetProperty().SetColor(255, 0, 255)     # Magenta
            self.lblActor.GetProperty().SetColor(255, 0, 255)
        elif color == 'black':
            self.actor.GetProperty().SetColor(0, 0, 0)         # Black
            self.lblActor.GetProperty().SetColor(0, 0, 0)

# Converts a member object into a member for the viewer
class VisMember():
    """VTK representation of a frame/beam member."""

    # Constructor
    def __init__(self, member, nodes, annotation_size=5, theme='default'):
        """Create the VTK actors for a member.

        :param Member3D member: Member to visualize.
        :param dict[str, Node3D] nodes: Model nodes used to find end coordinates.
        :param float annotation_size: Base size for labels (default ``5``).
        :param str theme: ``'default'`` or ``'print'`` to control colors.
        """

        # Generate a line for the member
        line = vtk.vtkLineSource()

        # Step through each node in the model and find the position of the i-node and j-node
        for node in nodes.values():

            # Check to see if the current node is the i-node
            if node.name == member.i_node.name:
                Xi = node.X
                Yi = node.Y
                Zi = node.Z
                line.SetPoint1(Xi, Yi, Zi)

            # Check to see if the current node is the j-node
            elif node.name == member.j_node.name:
                Xj = node.X
                Yj = node.Y
                Zj = node.Z
                line.SetPoint2(Xj, Yj, Zj)

        # Set up a mapper for the member
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(line.GetOutputPort())

        # Set up an actor for the member
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)

        # Create the text for the member label
        label = vtk.vtkVectorText()
        label.SetText(member.name)

        # Set up a mapper for the member label
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())

        # Set up an actor for the member label
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition((Xi+Xj)/2, (Yi+Yj)/2, (Zi+Zj)/2)

        # Adjust the color of the member if the theme is 'print'
        if theme == 'print':
            self.actor.GetProperty().SetColor(0/255, 0/255, 0/255)  # Black
            self.lblActor.GetProperty().SetColor(0/255, 0/255, 0/255) # Black

# Converts a node object into a node in its deformed position for the viewer
class VisDeformedNode():
    """Sphere representing a node in its deformed position."""

    def __init__(self, node, scale_factor, annotation_size=5, combo_name='Combo 1'):
        """Build the VTK geometry for a deformed node.

        :param Node3D node: Node whose displacement results are used.
        :param float scale_factor: Multiplier applied to translations for display.
        :param float annotation_size: Radius scale for the sphere (default ``5``).
        :param str combo_name: Load combination name used to fetch displacements.
        """

        # Calculate the node's deformed position
        newX = node.X + scale_factor*(node.DX[combo_name])
        newY = node.Y + scale_factor*(node.DY[combo_name])
        newZ = node.Z + scale_factor*(node.DZ[combo_name])

        # Generate a sphere source for the node in its deformed position
        self.source = vtk.vtkSphereSource()
        self.source.SetCenter(newX, newY, newZ)
        self.source.SetRadius(0.6*annotation_size)
        self.source.Update()

class VisDeformedMember():
    """Polyline representing a member's deformed shape."""

    def __init__(self, member, nodes, scale_factor, combo_name='Combo 1'):
        """Create the deformed polyline for a member.

        :param Member3D member: Member whose deformation is plotted.
        :param dict[str, Node3D] nodes: Model nodes to locate the member start.
        :param float scale_factor: Multiplier applied to displacement values.
        :param str combo_name: Load combination name used for results.
        """

        # Determine if this member is active for each load combination
        self.active = member.active

        L = member.L() # Member length
        T = member.T() # Member local transformation matrix

        cos_x = array([T[0,0:3]]) # Direction cosines of local x-axis
        cos_y = array([T[1,0:3]]) # Direction cosines of local y-axis
        cos_z = array([T[2,0:3]]) # Direction cosines of local z-axis

        # Find the initial position of the local i-node
        # Step through each node
        for node in nodes.values():

            # Check to see if the current node is the i-node
            if node.name == member.i_node.name:
                Xi = node.X
                Yi = node.Y
                Zi = node.Z

        # Calculate the local y-axis displacements at 20 points along the member's
        # length
        DY_plot = empty((0, 3))
        for i in range(20):

            # Calculate the local y-direction displacement
            dy_tot = member.deflection('dy', L/19*i, combo_name)

            # Calculate the scaled displacement in global coordinates
            DY_plot = append(DY_plot, dy_tot*cos_y*scale_factor, axis=0)

        # Calculate the local z-axis displacements at 20 points along the member's
        # length
        DZ_plot = empty((0, 3))
        for i in range(20):

            # Calculate the local z-direction displacement
            dz_tot = member.deflection('dz', L/19*i, combo_name)

            # Calculate the scaled displacement in global coordinates
            DZ_plot = append(DZ_plot, dz_tot*cos_z*scale_factor, axis=0)

        # Calculate the local x-axis displacements at 20 points along the member's
        # length
        DX_plot = empty((0, 3))
        for i in range(20):

            # Displacements in local coordinates
            dx_tot = [[Xi, Yi, Zi]] + (L/19*i + member.deflection('dx', L/19*i, combo_name)*scale_factor)*cos_x

            # Magnified displacements in global coordinates
            DX_plot = append(DX_plot, dx_tot, axis=0)

        # Sum the component displacements to obtain overall displacement
        D_plot = DY_plot + DZ_plot + DX_plot

        # Generate vtk points
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(len(D_plot))

        for i in range(len(D_plot)):
            points.SetPoint(i, D_plot[i, 0], D_plot[i, 1], D_plot[i, 2])

        # Generate vtk lines
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(len(D_plot))

        for i in range(len(D_plot)):
            lines.InsertCellPoint(i)

        # Create a polyline source from the defined points and lines
        self.source = vtk.vtkPolyData()
        self.source.SetPoints(points)
        self.source.SetLines(lines)

class VisDeformedSpring():
    """Line representation of a spring in its deformed position."""

    def __init__(self, spring, nodes, scale_factor, combo_name='Combo 1'):
        """Create the deformed line for a spring.

        :param Spring3D spring: Spring whose deformation is plotted.
        :param dict[str, Node3D] nodes: Model nodes to locate end coordinates.
        :param float scale_factor: Multiplier applied to displacement values.
        :param str combo_name: Load combination name used for results.
        """

        # Determine if this spring is active for each load combination
        self.active = spring.active

        # Generate a line source for the spring
        self.source = vtk.vtkLineSource()

        # Find the deformed position of the local i-node
        # Step through each node
        for node in nodes.values():

            # Check to see if the current node is the i-node
            if node.name == spring.i_node.name:
                Xi = node.X + node.DX[combo_name]*scale_factor
                Yi = node.Y + node.DY[combo_name]*scale_factor
                Zi = node.Z + node.DZ[combo_name]*scale_factor
                self.source.SetPoint1(Xi, Yi, Zi)

            # Check to see if the current node is the i-node
            if node.name == spring.j_node.name:
                Xj = node.X + node.DX[combo_name]*scale_factor
                Yj = node.Y + node.DY[combo_name]*scale_factor
                Zj = node.Z + node.DZ[combo_name]*scale_factor
                self.source.SetPoint2(Xj, Yj, Zj)

        self.source.Update()

class VisPtLoad():
    """Arrow representation of a concentrated load."""

    def __init__(self, position, direction, length, label_text: str | None = None, annotation_size=5, theme: str = 'default'):
        """Create a point-load arrow and optional label.

        :param tuple position: Coordinates of the arrow tip ``(X, Y, Z)``.
        :param tuple direction: Direction vector for the arrow (normalized).
        :param float length: Arrow length; sign controls pointing direction.
        :param str | None label_text: Text near the tail; ``None`` hides it.
        :param float annotation_size: Base scale for text/geometry (default ``5``).
        :param str theme: ``'default'`` or ``'print'`` color scheme.
        """

        # Create a unit vector in the direction of the 'direction' vector
        unitVector = direction/norm(direction)

        # Create a 'vtkAppendPolyData' filter to append the tip and shaft together into a single dataset
        self.polydata = vtk.vtkAppendPolyData()

        # Determine if the load is positive or negative
        if length == 0:
            sign = 1
        else:
            sign = abs(length)/length

        # Generate the tip of the load arrow
        tip_length = abs(length)/4
        radius = abs(length)/16
        tip = vtk.vtkConeSource()
        tip.SetCenter(position[0] - tip_length*sign*0.5*unitVector[0], \
                      position[1] - tip_length*sign*0.5*unitVector[1], \
                      position[2] - tip_length*sign*0.5*unitVector[2])
        tip.SetDirection([direction[0]*sign, direction[1]*sign, direction[2]*sign])
        tip.SetHeight(tip_length)
        tip.SetRadius(radius)
        tip.Update()

        # Add the arrow tip to the append filter
        self.polydata.AddInputData(tip.GetOutput())

        # Create the shaft
        shaft = vtk.vtkLineSource()
        shaft.SetPoint1(position)
        shaft.SetPoint2((position[0]-length*unitVector[0], position[1]-length*unitVector[1], position[2]-length*unitVector[2]))
        shaft.Update()

        # Copy and append the shaft data to the append filter
        self.polydata.AddInputData(shaft.GetOutput())
        self.polydata.Update()

        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.polydata.GetOutputPort())
        self.actor = vtk.vtkActor()

        # Set the colors to match the theme
        if theme == 'default':
            self.actor.GetProperty().SetColor(0, 1, 0) # Green
        elif theme == 'print':
            self.actor.GetProperty().SetColor(0, 0.75, 0)  # Dark Green

        self.actor.SetMapper(mapper)

        # Create the label if needed
        if label_text != None:

            # Create the label and set its text
            self.label = vtk.vtkVectorText()
            self.label.SetText(label_text)

            # Set up a mapper for the label
            lblMapper = vtk.vtkPolyDataMapper()
            lblMapper.SetInputConnection(self.label.GetOutputPort())

            # Set up an actor for the label
            self.lblActor = vtk.vtkFollower()
            self.lblActor.SetMapper(lblMapper)
            self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
            self.lblActor.SetPosition(position[0] - (length - 0.6*annotation_size)*unitVector[0], \
                                      position[1] - (length - 0.6*annotation_size)*unitVector[1], \
                                      position[2] - (length - 0.6*annotation_size)*unitVector[2])

            if theme == 'default':
                self.lblActor.GetProperty().SetColor(0, 1, 0)  # Green
            elif theme == 'print':
                self.lblActor.GetProperty().SetColor(0, 0.75, 0)  # Black


class VisDistLoad():
    """Series of arrows representing a distributed load."""

    def __init__(self, position1, position2, direction, length1, length2, label_text1, label_text2, annotation_size=5, theme = 'default'):
        """Create the VTK actors for a linearly varying load.

        :param tuple position1: Coordinates of the start point of the load.
        :param tuple position2: Coordinates of the end point of the load.
        :param tuple direction: Direction vector for each arrow (normalized).
        :param float length1: Arrow length at ``position1`` (sign indicates direction).
        :param float length2: Arrow length at ``position2`` (sign indicates direction).
        :param str label_text1: Label text shown on the first arrow.
        :param str label_text2: Label text shown on the last arrow.
        :param float annotation_size: Base scale for geometry (default ``5``).
        :param str theme: Color theme, ``'default'`` or ``'print'``.
        """

        # Calculate the length of the distributed load
        loadLength = ((position2[0]-position1[0])**2 + (position2[1]-position1[1])**2 + (position2[2]-position1[2])**2)**0.5

        # Find the direction cosines for the line the load acts on
        lineDirCos = [(position2[0]-position1[0])/loadLength, (position2[1]-position1[1])/loadLength, (position2[2]-position1[2])/loadLength]

        # Find the direction cosines for the direction the load acts in
        dirDirCos = direction/norm(direction)

        # Create point loads at intervals roughly equal to 75% of the load's largest length (magnitude)
        # Add text labels to the first and last load arrow
        if loadLength > 0:
            num_steps = int(round(0.75*loadLength/max(abs(length1), abs(length2)), 0))
        else:
            num_steps = 0

        num_steps = max(num_steps, 1)
        step = loadLength/num_steps
        ptLoads = []

        for i in range(num_steps + 1):

            # Calculate the position (X, Y, Z) of this load arrow's point
            position = (position1[0] + i*step*lineDirCos[0], position1[1] + i*step*lineDirCos[1], position1[2] + i*step*lineDirCos[2])

            # Determine the length of this load arrow
            length = length1 + (length2 - length1)/loadLength*i*step

            # Determine the label's text
            if i == 0:
                label_text = label_text1
            elif i == num_steps:
                label_text = label_text2

            # Create the load arrow
            ptLoads.append(VisPtLoad(position, direction, length, label_text, annotation_size, theme))

        # Draw a line between the first and last load arrow's tails
        tail_line = vtk.vtkLineSource()
        tail_line.SetPoint1((position1[0] - length1*dirDirCos[0], position1[1] - length1*dirDirCos[1], position1[2] - length1*dirDirCos[2]))
        tail_line.SetPoint2((position2[0] - length2*dirDirCos[0], position2[1] - length2*dirDirCos[1], position2[2] - length2*dirDirCos[2]))

        # Combine all the geometry into one 'vtkPolyData' object
        self.polydata = vtk.vtkAppendPolyData()
        for arrow in ptLoads:
            arrow.polydata.Update()
            self.polydata.AddInputData(arrow.polydata.GetOutput())

        tail_line.Update()
        self.polydata.AddInputData(tail_line.GetOutput())
        self.polydata.Update()

        # Create a mapper and actor for the geometry
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.polydata.GetOutputPort())
        self.actor = vtk.vtkActor()

        # Set the color
        if theme == 'default':
            self.actor.GetProperty().SetColor(0, 1, 0)  # Green
        elif theme == 'print':
            self.actor.GetProperty().SetColor(0, 0.75, 0)  # Dark Green

        self.actor.SetMapper(mapper)

        # Get the actors for the labels
        self.lblActors = [ptLoads[0].lblActor, ptLoads[len(ptLoads) - 1].lblActor]

class VisMoment():
    """Arc-and-arrow representation of a concentrated moment."""

    def __init__(self, center, direction, radius, label_text=None, annotation_size=5, theme='default'):
        """Create the VTK actors for a moment glyph.

        :param tuple center: Center point of the moment arc ``(X, Y, Z)``.
        :param tuple direction: Direction vector indicating the moment axis.
        :param float radius: Radius of the circular arc.
        :param str | None label_text: Text shown near the arrow head; ``None`` hides it.
        :param float annotation_size: Base scale for text/geometry (default ``5``).
        :param str theme: ``'default'`` or ``'print'`` to control colors.
        """

        # Create an append filter to store load polydata in
        self.polydata = vtk.vtkAppendPolyData()

        # Find a vector perpendicular to the directional unit vector
        v1 = direction/norm(direction)  # v1 = The directional unit vector for the moment
        v2 = _PerpVector(v1)             # v2 = A unit vector perpendicular to v1
        v3 = cross(v1, v2)
        v3 = v3/norm(v3)                # v3 = A unit vector perpendicular to v1 and v2

        # Generate an arc for the moment
        Xc, Yc, Zc = center
        arc = vtk.vtkArcSource()
        arc.SetCenter(Xc, Yc, Zc)
        arc.SetPoint1(Xc + v2[0]*radius, Yc + v2[1]*radius, Zc + v2[2]*radius)
        arc.SetPoint2(Xc + v3[0]*radius, Yc + v3[1]*radius, Zc + v3[2]*radius)
        arc.SetNegative(True)
        arc.SetResolution(20)
        arc.Update()
        self.polydata.AddInputData(arc.GetOutput())

        # Generate the arrow tip at the end of the arc
        tip_length = radius/2
        cone_radius = radius/8
        tip = vtk.vtkConeSource()
        tip.SetCenter(arc.GetPoint1()[0], arc.GetPoint1()[1], arc.GetPoint1()[2])
        tip.SetDirection(cross(v1, v2))
        tip.SetHeight(tip_length)
        tip.SetRadius(cone_radius)
        tip.Update()
        self.polydata.AddInputData(tip.GetOutput())

        # Update the polydata one last time now that we're done appending items to it
        self.polydata.Update()

        # Create the text label
        label = vtk.vtkVectorText()
        label.SetText(label_text)
        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())
        self.lblActor = vtk.vtkFollower()
        self.lblActor.SetMapper(lblMapper)
        self.lblActor.SetScale(annotation_size, annotation_size, annotation_size)
        self.lblActor.SetPosition(Xc + v3[0]*(radius + 0.25*annotation_size), \
                                  Yc + v3[1]*(radius + 0.25*annotation_size), \
                                  Zc + v3[2]*(radius + 0.25*annotation_size))

        if theme == 'default':
            self.lblActor.GetProperty().SetColor(0, 1, 0)  # Green
        elif theme == 'print':
            self.lblActor.GetProperty().SetColor(0, 0.75, 0)  # Dark Green

class VisAreaLoad():
    """Polygon and arrows used to visualize a uniform area load."""

    def __init__(self, position0, position1, position2, position3, direction, length, label_text, annotation_size=5, theme='default'):
        """Create the VTK actors for an area load.

        :param tuple position0: First corner of the loaded area.
        :param tuple position1: Second corner of the loaded area.
        :param tuple position2: Third corner of the loaded area.
        :param tuple position3: Fourth corner of the loaded area.
        :param tuple direction: Direction vector for the load arrows (normalized).
        :param float length: Arrow length; sign determines arrow orientation.
        :param str label_text: Text displayed at a corner arrow.
        :param float annotation_size: Base scale for text/geometry (default ``5``).
        :param str theme: Color theme, ``'default'`` or ``'print'``.
        """

        # Create a point load for each corner of the area load
        ptLoads = []
        ptLoads.append(VisPtLoad(position0, direction, length, label_text, annotation_size, theme))
        ptLoads.append(VisPtLoad(position1, direction, length, label_text, annotation_size, theme))
        ptLoads.append(VisPtLoad(position2, direction, length, label_text, annotation_size, theme))
        ptLoads.append(VisPtLoad(position3, direction, length, label_text, annotation_size, theme))

        # Find the direction cosines for the direction the load acts in
        dirDirCos = direction/norm(direction)

        # Find the positions of the tails of all the arrows at the corners of the area load. This is
        # where we will place the polygon.
        self.p0 = position0 - dirDirCos*length
        self.p1 = position1 - dirDirCos*length
        self.p2 = position2 - dirDirCos*length
        self.p3 = position3 - dirDirCos*length

        # Combine all geometry into one 'vtkPolyData' object
        self.polydata = vtk.vtkAppendPolyData()
        for arrow in ptLoads:
            self.polydata.AddInputData(arrow.polydata.GetOutput())
        self.polydata.Update()

        # Add a label
        self.label_actor = ptLoads[0].lblActor

        # Add color to the area load label
        if theme == 'default':
            self.label_actor.GetProperty().SetColor(0, 1, 0)  # Green
        elif theme == 'print':
            self.label_actor.GetProperty().SetColor(0, 0.75, 0)  # Dark Green

def _PerpVector(v):
    """Return a unit vector perpendicular to ``v``.

    :param array-like v: Input vector ``[i, j, k]``.
    :returns: Unit vector perpendicular to ``v``.
    :rtype: numpy.ndarray
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
    return [i2, j2, k2]/norm([i2, j2, k2])

def _PrepContour(model, stress_type='Mx', combo_name='Combo 1'):
    """Populate nodal ``contour`` values for plate/quad results.

    :param FEModel3D model: Model whose plate and quad results are used.
    :param str stress_type: Contour to compute (``'Mx'``, ``'My'``, ``'Mxy'``,
        ``'Qx'``, ``'Qy'``, ``'Sx'``, ``'Sy'``, ``'Txy'``, or ``'dz'``).
    :param str combo_name: Load combination name (default ``'Combo 1'``).
    """

    if stress_type != None:

        # Erase any previous contours
        for node in model.nodes.values():
            node.contour = []

        # Step through each element in the model
        for element in list(model.quads.values()) + list(model.plates.values()):

            # Rectangular elements and quadrilateral elements have different local coordinate systems.
            # Rectangles are based on a traditional (x, y) system, while quadrilaterals are based on a
            # 'natural' (r, s) coordinate system. To reduce duplication of code for both these elements
            # we'll define the edges of the plate here for either element using the (r, s) terminology.
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
                # Internally Pynite defines the nodes for a rectangular element in the order (i, j, m, n),
                # while it defines the nodes for a quadrilateral element in the order (m, n, i, j)
                if element.type == 'Rect':
                    i, j, m, n = element.d(combo_name)[[2, 8, 14, 20], :]
                else:
                    i, j, m, n = element.d(combo_name)[[14, 20, 2, 8], :]
                element.i_node.contour.append(i)
                element.j_node.contour.append(j)
                element.m_node.contour.append(m)
                element.n_node.contour.append(n)
            elif stress_type == 'Mx':
                element.i_node.contour.append(element.moment(r_left, s_bot, True, combo_name)[0])
                element.j_node.contour.append(element.moment(r_right, s_bot, True, combo_name)[0])
                element.m_node.contour.append(element.moment(r_right, s_top, True, combo_name)[0])
                element.n_node.contour.append(element.moment(r_left, s_top, True, combo_name)[0])
            elif stress_type == 'My':
                element.i_node.contour.append(element.moment(r_left, s_bot, True, combo_name)[1])
                element.j_node.contour.append(element.moment(r_right, s_bot, True, combo_name)[1])
                element.m_node.contour.append(element.moment(r_right, s_top, True, combo_name)[1])
                element.n_node.contour.append(element.moment(r_left, s_top, True, combo_name)[1])
            elif stress_type == 'Mxy':
                element.i_node.contour.append(element.moment(r_left, s_bot, True, combo_name)[2])
                element.j_node.contour.append(element.moment(r_right, s_bot, True, combo_name)[2])
                element.m_node.contour.append(element.moment(r_right, s_top, True, combo_name)[2])
                element.n_node.contour.append(element.moment(r_left, s_top, True, combo_name)[2])
            elif stress_type == 'Qx':
                element.i_node.contour.append(element.shear(r_left, s_bot, True, combo_name)[0])
                element.j_node.contour.append(element.shear(r_right, s_bot, True, combo_name)[0])
                element.m_node.contour.append(element.shear(r_right, s_top, True, combo_name)[0])
                element.n_node.contour.append(element.shear(r_left, s_top, True, combo_name)[0])
            elif stress_type == 'Qy':
                element.i_node.contour.append(element.shear(r_left, s_bot, True, combo_name)[1])
                element.j_node.contour.append(element.shear(r_right, s_bot, True, combo_name)[1])
                element.m_node.contour.append(element.shear(r_right, s_top, True, combo_name)[1])
                element.n_node.contour.append(element.shear(r_left, s_top, True, combo_name)[1])
            elif stress_type == 'Sx':
                element.i_node.contour.append(element.membrane(r_left, s_bot, True, combo_name)[0])
                element.j_node.contour.append(element.membrane(r_right, s_bot, True, combo_name)[0])
                element.m_node.contour.append(element.membrane(r_right, s_top, True, combo_name)[0])
                element.n_node.contour.append(element.membrane(r_left, s_top, True, combo_name)[0])
            elif stress_type == 'Sy':
                element.i_node.contour.append(element.membrane(r_left, s_bot, True, combo_name)[1])
                element.j_node.contour.append(element.membrane(r_right, s_bot, True, combo_name)[1])
                element.m_node.contour.append(element.membrane(r_right, s_top, True, combo_name)[1])
                element.n_node.contour.append(element.membrane(r_left, s_top, True, combo_name)[1])
            elif stress_type == 'Txy':
                element.i_node.contour.append(element.membrane(r_left, s_bot, True, combo_name)[2])
                element.j_node.contour.append(element.membrane(r_right, s_bot, True, combo_name)[2])
                element.m_node.contour.append(element.membrane(r_right, s_top, True, combo_name)[2])
                element.n_node.contour.append(element.membrane(r_left, s_top, True, combo_name)[2])

        # Average the values at each node to obtain a smoothed contour
        for node in model.nodes.values():
            # Prevent divide by zero errors for nodes with no contour values
            if node.contour != []:
                node.contour = (sum(node.contour)/len(node.contour))[0]  # The [0] converts it from an array to a float

def _DeformedShape(model, vtk_renderer, scale_factor, annotation_size, combo_name, render_nodes=True, theme='default'):
    """Render the deformed shape of a model.

    :param FEModel3D model: Finite element model to be rendered.
    :param vtk.vtkRenderer vtk_renderer: Renderer receiving the actor.
    :param float scale_factor: Scale factor applied to deformations.
    :param float annotation_size: Base size for spheres/text.
    :param str combo_name: Load combination name used for displacements.
    :param bool render_nodes: If ``True`` include deformed nodes (default).
    :param str theme: ``'default'`` or ``'print'`` color scheme.
    """

    # Create an append filter to add all the shape polydata to
    append_filter = vtk.vtkAppendPolyData()

    # Check if nodes are to be rendered
    if render_nodes == True:

        # Add the deformed nodes to the append filter
        for node in model.nodes.values():

            vis_node = VisDeformedNode(node, scale_factor, annotation_size, combo_name)
            append_filter.AddInputData(vis_node.source.GetOutput())

    # Add the springs to the append filter
    for spring in model.springs.values():

        # Only add the spring if it is active for the given load combination
        if spring.active[combo_name] == True:

            vis_spring = VisDeformedSpring(spring, model.nodes, scale_factor, combo_name)
            append_filter.AddInputData(vis_spring.source.GetOutput())

    # Add the members to the append filter
    for member in model.members.values():

        # Only add the member if it is active for the given load combination.
        if member.active[combo_name] == True:

            vis_member = VisDeformedMember(member, model.nodes, scale_factor, combo_name)
            append_filter.AddInputData(vis_member.source)

    # Create a mapper and actor for the append filter
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(append_filter.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Adjust the color
    if theme == 'default':
        actor.GetProperty().SetColor(255/255, 255/255, 0/255)  # Yellow
    elif theme == 'print':
        actor.GetProperty().SetColor(26/255, 26/255, 26/255)  # Dark Grey

    # Add the actor to the renderer
    vtk_renderer.AddActor(actor)

def _RenderLoads(model, renderer, annotation_size, combo_name, case, theme='default'):
    """Add load glyphs to a VTK renderer.

    :param FEModel3D model: Model containing loads to draw.
    :param vtk.vtkRenderer renderer: Renderer receiving the load actors.
    :param float annotation_size: Base scale for text and arrow geometry.
    :param str | None combo_name: Load combination name to display (ignored
        when ``case`` is provided).
    :param str | None case: Load case name to display (mutually exclusive with
        ``combo_name``).
    :param str theme: ``'default'`` or ``'print'`` color scheme.
    """

    # Create an append filter to store all the polydata in. This will allow us to use fewer actors to
    # display all the loads, which will greatly improve rendering speed as the user interacts. VTK
    # becomes very slow when a large number of actors are used.
    polydata = vtk.vtkAppendPolyData()

    # Polygons are treated as cells in VTK. Create a cell array to store all the area load polygons
    # in. We'll also create a list of points to store the polygon points in. The polydata for these
    # polygons will be stored separately from the other load data.
    polygons = vtk.vtkCellArray()
    polygon_points = vtk.vtkPoints()
    polygon_polydata = vtk.vtkPolyData()

    # Track whether any loads have been added to avoid VTK errors with empty append filters
    has_loads = False

    # Get the maximum load magnitudes that will be used to normalize the display scale
    max_pt_load, max_moment, max_dist_load, max_area_load = _MaxLoads(model, combo_name, case)

    # Display the requested load combination, or 'Combo 1' if no load combo or case has been
    # specified
    if case == None:
        # Store model.load_combos[combo].factors under a simpler name for use below
        load_factors = model.load_combos[combo_name].factors
    else:
        # Set up a load combination dictionary that represents the load case
        load_factors = {case: 1}

    # Step through each node
    for node in model.nodes.values():

        # Step through and display each nodal load
        for load in node.NodeLoads:

            # Determine if this load is part of the requested LoadCombo or case
            if load[2] in load_factors:

                # Calculate the factored value for this load and it's sign (positive or negative)
                load_value = load[1]*load_factors[load[2]]
                if load_value != 0:
                    sign = load_value/abs(load_value)
                else:
                    sign = 1

                # Display the load
                if load[0] == 'FX':
                    ptLoad = VisPtLoad((node.X - 0.6*annotation_size*sign, node.Y, node.Z), [1, 0, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'FY':
                    ptLoad = VisPtLoad((node.X, node.Y - 0.6*annotation_size*sign, node.Z), [0, 1, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'FZ':
                    ptLoad = VisPtLoad((node.X, node.Y, node.Z - 0.6*annotation_size*sign), [0, 0, 1], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'MX':
                    ptLoad = VisMoment((node.X, node.Y, node.Z), (1*sign, 0, 0), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'MY':
                    ptLoad = VisMoment((node.X, node.Y, node.Z), (0, 1*sign, 0), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'MZ':
                    ptLoad = VisMoment((node.X, node.Y, node.Z), (0, 0, 1*sign), abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)

                polydata.AddInputData(ptLoad.polydata.GetOutput())
                renderer.AddActor(ptLoad.lblActor)
                ptLoad.lblActor.SetCamera(renderer.GetActiveCamera())
                has_loads = True

    # Step through each member
    for member in model.members.values():

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
                    ptLoad = VisPtLoad(position, dir_cos[0, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'Fy':
                    ptLoad = VisPtLoad(position, dir_cos[1, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'Fz':
                    ptLoad = VisPtLoad(position, dir_cos[2, :], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'Mx':
                    ptLoad = VisMoment(position, dir_cos[0, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'My':
                    ptLoad = VisMoment(position, dir_cos[1, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'Mz':
                    ptLoad = VisMoment(position, dir_cos[2, :]*sign, abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'FX':
                    ptLoad = VisPtLoad(position, [1, 0, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'FY':
                    ptLoad = VisPtLoad(position, [0, 1, 0], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'FZ':
                    ptLoad = VisPtLoad(position, [0, 0, 1], load_value/max_pt_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'MX':
                    ptLoad = VisMoment(position, [1*sign, 0, 0], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'MY':
                    ptLoad = VisMoment(position, [0, 1*sign, 0], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)
                elif load[0] == 'MZ':
                    ptLoad = VisMoment(position, [0, 0, 1*sign], abs(load_value)/max_moment*2.5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)

                polydata.AddInputData(ptLoad.polydata.GetOutput())
                renderer.AddActor(ptLoad.lblActor)
                ptLoad.lblActor.SetCamera(renderer.GetActiveCamera())
                has_loads = True

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
                if load[0] == 'Fx':
                    distLoad = VisDistLoad(position1, position2, dir_cos[0, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'Fy':
                    distLoad = VisDistLoad(position1, position2, dir_cos[1, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'Fz':
                    distLoad = VisDistLoad(position1, position2, dir_cos[2, :], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'FX':
                    distLoad = VisDistLoad(position1, position2, [1, 0, 0], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'FY':
                    distLoad = VisDistLoad(position1, position2, [0, 1, 0], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)
                elif load[0] == 'FZ':
                    distLoad = VisDistLoad(position1, position2, [0, 0, 1], w1/max_dist_load*5*annotation_size, w2/max_dist_load*5*annotation_size, '{:.3g}'.format(w1), '{:.3g}'.format(w2), annotation_size)

                polydata.AddInputData(distLoad.polydata.GetOutput())
                renderer.AddActor(distLoad.lblActors[0])
                renderer.AddActor(distLoad.lblActors[1])
                distLoad.lblActors[0].SetCamera(renderer.GetActiveCamera())
                distLoad.lblActors[1].SetCamera(renderer.GetActiveCamera())
                has_loads = True

    # Step through each plate
    i = 0
    for plate in list(model.plates.values()) + list(model.quads.values()):

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
                area_load = VisAreaLoad(position0, position1, position2, position3, dir_cos*sign, abs(load_value)/max_area_load*5*annotation_size, '{:.3g}'.format(load_value), annotation_size, theme)

                # Add the area load's arrows to the overall load polydata
                polydata.AddInputData(area_load.polydata.GetOutput())
                has_loads = True

                # Add the 4 points at the corners of this area load to the list of points
                polygon_points.InsertNextPoint(area_load.p0[0], area_load.p0[1], area_load.p0[2])
                polygon_points.InsertNextPoint(area_load.p1[0], area_load.p1[1], area_load.p1[2])
                polygon_points.InsertNextPoint(area_load.p2[0], area_load.p2[1], area_load.p2[2])
                polygon_points.InsertNextPoint(area_load.p3[0], area_load.p3[1], area_load.p3[2])

                # Create a polygon based on the four points we just defined.
                # The 1st number in `SetId()` is the local point id
                # The 2nd number in `SetId()` is the global point id
                polygon = vtk.vtkPolygon()
                polygon.GetPointIds().SetNumberOfIds(4)
                polygon.GetPointIds().SetId(0, i*4)
                polygon.GetPointIds().SetId(1, i*4 + 1)
                polygon.GetPointIds().SetId(2, i*4 + 2)
                polygon.GetPointIds().SetId(3, i*4 + 3)

                # Add the polygon to the list of polygons
                polygons.InsertNextCell(polygon)

                # Add the load label
                renderer.AddActor(area_load.label_actor)

                # Set the text to follow the camera as the user interacts
                area_load.label_actor.SetCamera(renderer.GetActiveCamera())

                # `i` keeps track of the next polygon's ID. We've just added a polygon, so `i` needs to
                # go up 1.
                i += 1

                # Create polygon polydata from all the points and polygons we just defined
                polygon_polydata.SetPoints(polygon_points)
                polygon_polydata.SetPolys(polygons)

    # Only create and add actors if there are loads to render (avoids VTK errors with empty append filters)
    if has_loads:
        # Set up an actor and mapper for the loads
        load_mapper = vtk.vtkPolyDataMapper()
        load_mapper.SetInputConnection(polydata.GetOutputPort())
        load_actor = vtk.vtkActor()
        load_actor.SetMapper(load_mapper)

        # Colorize the loads
        if theme == 'default':
            load_actor.GetProperty().SetColor(0, 1, 0)  # Green
        elif theme == 'print':
            load_actor.GetProperty().SetColor(0, 0.75, 0)  # Dark Green

        # Add the load actor to the renderer
        renderer.AddActor(load_actor)

    # Only create polygon actors if there are area loads (i > 0 means we added polygons)
    if i > 0:
        # Set up an actor and a mapper for the area load polygons
        polygon_mapper = vtk.vtkPolyDataMapper()
        polygon_mapper.SetInputData(polygon_polydata)
        polygon_actor = vtk.vtkActor()

        # polygon_actor.GetProperty().SetOpacity(0.5)      # 50% opacity
        polygon_actor.SetMapper(polygon_mapper)
        renderer.AddActor(polygon_actor)

        # Set the color of the area load polygons
        if theme == 'default':
            polygon_actor.GetProperty().SetColor(0, 1, 0)  # Green
        elif theme == 'print':
            polygon_actor.GetProperty().SetColor(0, 0.75, 0)  # Dark Green

def _RenderContours(model, renderer, deformed_shape, deformed_scale, color_map, scalar_bar, scalar_bar_text_size, combo_name, theme='default'):
    """Render plate/quad contours and optional scalar bar.

    :param FEModel3D model: Model containing plates/quads.
    :param vtk.vtkRenderer renderer: Renderer receiving the contour actors.
    :param bool deformed_shape: If ``True`` use deformed node positions.
    :param float deformed_scale: Scale factor for deformations when drawn.
    :param str | None color_map: Result component to map to colors; ``None``
        disables contours.
    :param bool scalar_bar: If ``True`` add a scalar bar legend.
    :param int scalar_bar_text_size: Font size for scalar bar labels.
    :param str combo_name: Load combination name used for contour values.
    :param str theme: ``'default'`` or ``'print'`` color choices.
    """

    # Create a new `vtkCellArray` object to store the elements
    plates = vtk.vtkCellArray()

    # Create a `vtkPoints` object to store the coordinates of the corners of the elements
    plate_points = vtk.vtkPoints()

    # Create 2 lists to store plate result
    # `results` will store the results in a Python iterable list
    # `plate_results` will store the results in a `vtkDoubleArray` for VTK
    results = []
    plate_results = vtk.vtkDoubleArray()
    plate_results.SetNumberOfComponents(1)

    # Each element will be assigned a unique element number `i` beginning at 0
    i = 0

    # Calculate the smoothed contour results at each node
    _PrepContour(model, color_map, combo_name)

    # Add each plate and quad in the model to the cell array we just created
    for item in list(model.plates.values()) + list(model.quads.values()):

        # Create a point for each corner (must be in counter clockwise order)
        if deformed_shape == True:
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

        # Add the points to the `vtkPoints` object we created earlier
        plate_points.InsertNextPoint(p0)
        plate_points.InsertNextPoint(p1)
        plate_points.InsertNextPoint(p2)
        plate_points.InsertNextPoint(p3)

        # Create a `vtkQuad` based on the four points we just defined
        # The 1st number in `SetId()` is the local point id
        # The 2nd number in `SetId()` is the global point id
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0, i*4)
        quad.GetPointIds().SetId(1, i*4 + 1)
        quad.GetPointIds().SetId(2, i*4 + 2)
        quad.GetPointIds().SetId(3, i*4 + 3)

        # Get the contour value for each node
        # Convert to scalar for NumPy >= 2.4.0 compatibility
        # Handle cases where contour might be empty list, list with values, or scalar
        r0 = item.i_node.contour
        r1 = item.j_node.contour
        r2 = item.m_node.contour
        r3 = item.n_node.contour

        if color_map != None:

            # Save the results to the Python list of results we created earlier
            results.append(r0)
            results.append(r1)
            results.append(r2)
            results.append(r3)

            # Save the results to the `vtkDoubleArray` list of results for VTK
            plate_results.InsertNextTuple([r0])
            plate_results.InsertNextTuple([r1])
            plate_results.InsertNextTuple([r2])
            plate_results.InsertNextTuple([r3])

        # Insert the quad into the cell array
        plates.InsertNextCell(quad)

        # Increment `i` for the next plate
        i += 1

    # Create a `vtkPolyData` object to store plate data in
    plate_polydata = vtk.vtkPolyData()

    # Add the points and plates to the dataset
    plate_polydata.SetPoints(plate_points)
    plate_polydata.SetPolys(plates)

    # Setup actor and mapper for the plates
    plate_mapper = vtk.vtkPolyDataMapper()
    plate_mapper.SetInputData(plate_polydata)
    plate_actor = vtk.vtkActor()
    plate_actor.SetMapper(plate_mapper)

    # Map the results to the plates
    if color_map != None:

        plate_polydata.GetPointData().SetScalars(plate_results)

        # Create a `vtkLookupTable` for the colors used to map results
        lut = vtk.vtkLookupTable()
        lut.SetTableRange(min(results), max(results))
        lut.SetNumberOfColors(256)
        # The commented code below can be uncommented and modified to change the color scheme
        # ctf = vtk.vtkColorTransferFunction()
        # ctf.SetColorSpaceToDiverging()
        # ctf.AddRGBPoint(min(results), 255, 0, 255)  # Purple
        # ctf.AddRGBPoint(max(results), 255, 0, 0)    # Red
        # for i in range(256):
        #   rgb = list(ctf.GetColor(float(i)/256))
        #   rgb.append(1.0)
        #   lut.SetTableValue(i, *rgb)
        plate_mapper.SetLookupTable(lut)
        plate_mapper.SetUseLookupTableScalarRange(True)
        plate_mapper.SetScalarModeToUsePointData()
        lut.Build()

        # Add the scalar bar for the contours.
        if scalar_bar:

            if Renderer.scalar == None:
                Renderer.scalar = vtk.vtkScalarBarActor()

            scalar = Renderer.scalar

            # This next group of lines controls the font on the scalar bar
            scalar.SetUnconstrainedFontSize(1)
            scalar_text = vtk.vtkTextProperty()
            scalar_text.SetFontSize(max(int(scalar_bar_text_size), 1))
            scalar_text.SetBold(1)

            # The `vtkTextProperty` object is white by default
            if theme == 'print':
                scalar_text.SetColor(0, 0, 0)  # Black

            scalar.SetLabelTextProperty(scalar_text)

            scalar.SetMaximumWidthInPixels(100)

            scalar.SetTextPositionToPrecedeScalarBar()

            scalar.SetLookupTable(lut)

            renderer.AddActor(scalar)

    # Add the actor for the plates
    renderer.AddActor(plate_actor)

def _MaxLoads(model, combo_name=None, case=None):
    """Return maximum load magnitudes for scaling glyphs.

    :param FEModel3D model: Model containing loads.
    :param str | None combo_name: Load combination name to evaluate.
    :param str | None case: Load case name to evaluate.
    :returns: Maximum point load, moment, distributed load, and area load
        magnitudes (zeros replaced by ones to avoid division by zero).
    :rtype: tuple[float, float, float, float]
    """

    max_pt_load = 0
    max_moment = 0
    max_dist_load = 0
    max_area_load = 0

    # Find the requested load combination or load case
    if case == None:

        # Step through each node
        for node in model.nodes.values():

            # Step through each nodal load to find the largest one
            for load in node.NodeLoads:

                # Find the largest loads in the load combination
                if load[2] in model.load_combos[combo_name].factors:
                    if load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ':
                        if abs(load[1]*model.load_combos[combo_name].factors[load[2]]) > max_pt_load:
                            max_pt_load = abs(load[1]*model.load_combos[combo_name].factors[load[2]])
                    else:
                        if abs(load[1]*model.load_combos[combo_name].factors[load[2]]) > max_moment:
                            max_moment = abs(load[1]*model.load_combos[combo_name].factors[load[2]])

        # Step through each member
        for member in model.members.values():

            # Step through each member point load
            for load in member.PtLoads:

                # Find and store the largest point load and moment in the load combination
                if load[3] in model.load_combos[combo_name].factors:

                    if (load[0] == 'Fx' or load[0] == 'Fy' or load[0] == 'Fz'
                    or  load[0] == 'FX' or load[0] == 'FY' or load[0] == 'FZ'):
                        if abs(load[1]*model.load_combos[combo_name].factors[load[3]]) > max_pt_load:
                            max_pt_load = abs(load[1]*model.load_combos[combo_name].factors[load[3]])
                    else:
                        if abs(load[1]*model.load_combos[combo_name].factors[load[3]]) > max_moment:
                            max_moment = abs(load[1]*model.load_combos[combo_name].factors[load[3]])

            # Step through each member distributed load
            for load in member.DistLoads:

                #Find and store the largest distributed load in the load combination
                if load[5] in model.load_combos[combo_name].factors:

                    if abs(load[1]*model.load_combos[combo_name].factors[load[5]]) > max_dist_load:
                        max_dist_load = abs(load[1]*model.load_combos[combo_name].factors[load[5]])
                    if abs(load[2]*model.load_combos[combo_name].factors[load[5]]) > max_dist_load:
                        max_dist_load = abs(load[2]*model.load_combos[combo_name].factors[load[5]])

        # Step through each plate
        for plate in model.plates.values():

            # Step through each plate load
            for load in plate.pressures:

                if load[1] in model.load_combos[combo_name].factors:
                    if abs(load[0]*model.load_combos[combo_name].factors[load[1]]) > max_area_load:
                        max_area_load = abs(load[0]*model.load_combos[combo_name].factors[load[1]])

        # Step through each quad
        for quad in model.quads.values():

            # Step through each plate load
            for load in quad.pressures:

                # Check to see if the load case is in the requested load combination
                if load[1] in model.load_combos[combo_name].factors:
                    if abs(load[0]*model.load_combos[combo_name].factors[load[1]]) > max_area_load:
                        max_area_load = abs(load[0]*model.load_combos[combo_name].factors[load[1]])

    # Behavior if case has been specified
    else:

        # Step through each node
        for node in model.nodes.values():

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
        for member in model.members.values():

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
            for plate in model.plates.values():

                # Step through each plate load
                for load in plate.pressures:

                    if load[1] == case:

                        if abs(load[0]) > max_area_load:
                            max_area_load = abs(load[0])

        # Step through each quad
        for quad in model.quads.values():

            # Step through each plate load
            for load in quad.pressures:

                if load[1] == case:

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

    # Return the maximum loads in the load combination or load case
    return max_pt_load, max_moment, max_dist_load, max_area_load

def _MaxInternalForces(model, diagram_type, combo_name):
    """Calculate maximum internal force/moment magnitudes across all members for consistent scaling.

    :param FEModel3D model: The finite element model.
    :param str diagram_type: Diagram type (``'Fy'``, ``'Fz'``, ``'My'``, ``'Mz'``, ``'Fx'``, ``'Tx'``).
    :param str combo_name: Load combination name.
    :returns: Maximum absolute value of the internal force/moment.
    :rtype: float
    """

    max_value = 0

    # Iterate through all members to find the global maximum
    for member in model.members.values():

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

class VisMemberDiagram():
    """Creates a visual representation of internal forces/moments along a member."""

    def __init__(self, member, nodes, diagram_type, scale_factor, combo_name='Combo 1', theme='default', n_points=20, annotation_size=5, global_max=None):
        """Build a VTK polyline diagram for a member.

        :param Member3D member: Member to diagram.
        :param dict[str, Node3D] nodes: Model nodes.
        :param str diagram_type: One of ``'Fy'``, ``'Fz'``, ``'My'``, ``'Mz'``,
            ``'Fx'``, ``'Tx'``.
        :param float scale_factor: Scale factor for diagram magnitude.
        :param str combo_name: Load combination name (default ``'Combo 1'``).
        :param str theme: Visual theme (``'default'`` or ``'print'``).
        :param int n_points: Number of sample points along member.
        :param float annotation_size: Text/label scale.
        :param float | None global_max: Global maximum value for consistent scaling across all members.
            If ``None``, scales based on individual member maximum.
        """

        self.polydata = vtk.vtkAppendPolyData()
        self.label_actors = []  # Store text label actors
        
        # Get member information
        i_node = member.i_node
        j_node = member.j_node
        Xi, Yi, Zi = i_node.X, i_node.Y, i_node.Z
        Xj, Yj, Zj = j_node.X, j_node.Y, j_node.Z
        L = member.L()
        
        # Get transformation matrix for local coordinates
        T = member.T()
        cos_x = array([T[0, 0:3]])  # Local x-axis (along member)
        cos_y = array([T[1, 0:3]])  # Local y-axis
        cos_z = array([T[2, 0:3]])  # Local z-axis
        
        # Member base line
        member_start = array([Xi, Yi, Zi])
        member_dir = array([Xj - Xi, Yj - Yi, Zj - Zi])
        member_unit = member_dir / norm(member_dir)
        
        # Determine perpendicular direction for diagram offset
        if diagram_type in ['Fy', 'Mz']:
            perp_dir = cos_y[0]  # Use y direction for offset
        elif diagram_type in ['Fz', 'My']:
            perp_dir = cos_z[0]  # Use z direction for offset
        else:
            perp_dir = cos_y[0]  # Default to y direction
        
        # Get result values at points along member
        from numpy import linspace
        x_array = linspace(0, L, n_points)
        
        if diagram_type == 'Fy':
            results = member.shear_array('Fy', n_points, combo_name, x_array)
            y_values = results[1]
            label = 'Fy'
            max_value = member.max_shear('Fy', combo_name)
            min_value = member.min_shear('Fy', combo_name)
        elif diagram_type == 'Fz':
            results = member.shear_array('Fz', n_points, combo_name, x_array)
            y_values = results[1]
            label = 'Fz'
            max_value = member.max_shear('Fz', combo_name)
            min_value = member.min_shear('Fz', combo_name)
        elif diagram_type == 'My':
            results = member.moment_array('My', n_points, combo_name, x_array)
            y_values = results[1]
            label = 'My'
            max_value = member.max_moment('My', combo_name)
            min_value = member.min_moment('My', combo_name)
        elif diagram_type == 'Mz':
            results = member.moment_array('Mz', n_points, combo_name, x_array)
            y_values = results[1]
            label = 'Mz'
            max_value = member.max_moment('Mz', combo_name)
            min_value = member.min_moment('Mz', combo_name)
        elif diagram_type == 'Fx':
            results = member.axial_array(n_points, combo_name, x_array)
            y_values = results[1]
            label = 'Fx'
            max_value = member.max_axial(combo_name)
            min_value = member.min_axial(combo_name)
        elif diagram_type == 'Tx':
            results = member.torque_array(n_points, combo_name, x_array)
            y_values = results[1]
            label = 'Tx'
            max_value = member.max_torque(combo_name)
            min_value = member.min_torque(combo_name)
        else:
            return
        
        # Find indices of max and min values in the results array
        max_idx = int(y_values.argmax())
        min_idx = int(y_values.argmin())
        
        # Normalize values for better visualization
        # Use global_max if provided for consistent scaling across all members
        normalization_max = global_max if global_max is not None else max(abs(y_values))
        if normalization_max > 0:
            normalized_values = y_values / normalization_max
        else:
            normalized_values = y_values
        
        # Create diagram line (positive and negative sides)
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        
        point_idx = 0
        
        # Create baseline and diagram lines
        for i, x in enumerate(x_array):
            # Position along member
            pos_along_member = member_start + (x / L) * member_dir
            # Baseline point (on member axis)
            pts_idx = points.InsertNextPoint(pos_along_member)
            point_idx += 1
        
        # Add diagram points (displaced from member axis)
        for i, x in enumerate(x_array):
            pos_along_member = member_start + (x / L) * member_dir
            diag_displacement = (normalized_values[i] * scale_factor * 0.5) * perp_dir
            diagram_pt = pos_along_member + diag_displacement
            points.InsertNextPoint(diagram_pt)
            point_idx += 1
        
        # Create baseline line
        baseline_line = vtk.vtkCellArray()
        baseline_line.InsertNextCell(len(x_array))
        for i in range(len(x_array)):
            baseline_line.InsertCellPoint(i)
        
        # Create diagram line
        diagram_line = vtk.vtkCellArray()
        diagram_line.InsertNextCell(len(x_array))
        for i in range(len(x_array)):
            diagram_line.InsertCellPoint(len(x_array) + i)
        
        # Create connection lines (vertical lines from baseline to diagram)
        connections = vtk.vtkCellArray()
        for i in range(len(x_array)):
            connections.InsertNextCell(2)
            connections.InsertCellPoint(i)
            connections.InsertCellPoint(len(x_array) + i)
        
        # Create polydata
        self.source = vtk.vtkPolyData()
        self.source.SetPoints(points)
        self.source.SetLines(baseline_line)
        
        # Add diagram line
        diagram_lines = vtk.vtkPolyData()
        diagram_lines.SetPoints(points)
        diagram_lines.SetLines(diagram_line)
        self.polydata.AddInputData(diagram_lines)
        
        # Add connection lines
        connections_polydata = vtk.vtkPolyData()
        connections_polydata.SetPoints(points)
        connections_polydata.SetLines(connections)
        self.polydata.AddInputData(connections_polydata)
        
        # Add baseline
        baseline_polydata = vtk.vtkPolyData()
        baseline_polydata.SetPoints(points)
        baseline_polydata.SetLines(baseline_line)
        self.polydata.AddInputData(baseline_polydata)
        
        self.polydata.Update()
        
        # Create mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.polydata.GetOutputPort())
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)
        
        # Set color based on theme
        if theme == 'default':
            self.actor.GetProperty().SetColor(0, 1, 1)  # Cyan for diagrams
        else:
            self.actor.GetProperty().SetColor(0, 0, 0)  # Black for print theme
        
        # Set line width
        self.actor.GetProperty().SetLineWidth(1)
        
        # Create text labels for max and min values
        self._create_value_labels(diagram_type, max_idx, min_idx, max_value, min_value, x_array, points, 
                                 normalized_values, scale_factor, perp_dir, theme, member_start, member_dir, L, annotation_size, member, combo_name)
    
    def _create_value_labels(self, diagram_type, max_idx, min_idx, max_value, min_value, x_array, points, 
                            normalized_values, scale_factor, perp_dir, theme, member_start, member_dir, L, annotation_size, member, combo_name):
        """Create text labels showing max and min values on the diagram."""
        
        # Format values for display (remove trailing zeros)
        max_str = f"{max_value:.3g}"
        min_str = f"{min_value:.3g}"
        
        # Get the actual diagram point positions (which are stored in the VTK points object)
        # The diagram points start at index n_points in the points array
        n_points = len(x_array)
        max_pt_idx = n_points + max_idx
        min_pt_idx = n_points + min_idx
        
        # Get the world coordinates of the diagram points
        max_pos = points.GetPoint(max_pt_idx)
        min_pos = points.GetPoint(min_pt_idx)
        
        # Add displacement for visibility, proportional to annotation_size
        label_offset = 0.1 * annotation_size
        max_pos = array(max_pos) + label_offset * perp_dir
        min_pos = array(min_pos) + label_offset * perp_dir
        
        # Text size should scale with annotation_size
        text_scale = annotation_size
        
        # Create max value label using vtkFollower (follows camera)
        max_text = vtk.vtkFollower()
        max_label = vtk.vtkVectorText()
        max_label.SetText(max_str)
        
        max_mapper = vtk.vtkPolyDataMapper()
        max_mapper.SetInputConnection(max_label.GetOutputPort())
        max_text.SetMapper(max_mapper)
        
        # Position the label at the max diagram point
        max_text.SetPosition(max_pos[0], max_pos[1], max_pos[2])
        max_text.SetScale(text_scale, text_scale, text_scale)  # Scale text with annotation_size
        
        # Configure properties for max
        max_prop = max_text.GetProperty()
        if theme == 'default':
            max_prop.SetColor(0, 1, 1)  # Cyan text
        else:
            max_prop.SetColor(0, 0, 0)  # Black text
        max_prop.SetLineWidth(1)
        
        # Create min value label
        min_text = vtk.vtkFollower()
        min_label = vtk.vtkVectorText()
        min_label.SetText(min_str)
        
        min_mapper = vtk.vtkPolyDataMapper()
        min_mapper.SetInputConnection(min_label.GetOutputPort())
        min_text.SetMapper(min_mapper)
        
        # Position the label at the min point
        min_text.SetPosition(min_pos[0], min_pos[1], min_pos[2])
        min_text.SetScale(text_scale, text_scale, text_scale)  # Scale text with annotation_size
        
        # Configure properties for min
        min_prop = min_text.GetProperty()
        if theme == 'default':
            min_prop.SetColor(0, 1, 1)  # Cyan text
        else:
            min_prop.SetColor(0, 0, 0)  # Black text
        min_prop.SetLineWidth(1)
        
        # Store the text actors for renderer integration
        self.label_actors = [
            (max_text, max_pos),
            (min_text, min_pos)
        ]


def _RenderMemberDiagrams(model, renderer, diagram_type, scale_factor, combo_name=None, case=None, theme='default', annotation_size=5):
    """Render internal force/moment diagrams on members.

    :param FEModel3D model: The finite element model.
    :param vtk.vtkRenderer renderer: The VTK renderer.
    :param str diagram_type: Diagram to render (``'Fy'``, ``'Fz'``, ``'My'``,
        ``'Mz'``, ``'Fx'``, ``'Tx'``).
    :param float scale_factor: Scale factor for visualization.
    :param str | None combo_name: Load combination name.
    :param str | None case: Load case name.
    :param str theme: Visual theme.
    :param float annotation_size: Size for annotations/labels.
    """

    # Determine which combo/case to use
    if case is None and combo_name is None:
        combo_name = 'Combo 1'
    elif case is not None:
        # For load cases, use the case name as combo_name in the member methods
        combo_name = case
    
    # Calculate global maximum internal force/moment for consistent scaling across all members
    global_max = _MaxInternalForces(model, diagram_type, combo_name)
    
    # Create diagrams for each active member
    for member in model.members.values():
        
        # Check if member is active for the specified combo
        if combo_name not in member.active or not member.active[combo_name]:
            continue
        
        try:
            vis_diagram = VisMemberDiagram(member, model.nodes, diagram_type, 
                                          scale_factor, combo_name, theme, annotation_size=annotation_size, 
                                          global_max=global_max)
            renderer.AddActor(vis_diagram.actor)
            
            # Add text label actors
            for label_actor, world_pos in vis_diagram.label_actors:
                # For vtkFollower actors, set the camera so they always face the viewer
                label_actor.SetCamera(renderer.GetActiveCamera())
                
                # Add the label actor to the renderer
                renderer.AddActor(label_actor)

        except Exception:
            # Skip members that can't create diagrams (e.g., inactive members)
            pass