=========
Rendering
=========

`Pynite` provides two complementary visualization libraries to help you visualize and analyze your finite element models. Both renderers support interactive visualization, load glyphs, deformed shapes, contour plots, and internal force/moment diagrams. You can choose the renderer that best fits your workflow.

Choosing a Renderer
===================

**PyVista Renderer (Recommended for most users)**

The ``Rendering`` module uses `PyVista <https://docs.pyvista.org/>`_, a high-level Python interface to `VTK <https://www.vtk.org/>`_ by Kitware. PyVista is modern, actively maintained, and provides an intuitive API. Use this renderer for:

- Interactive 3D visualization and manipulation
- Jupyter notebook integration with Trame backend
- Publication-quality renderings
- Modern Python workflows

.. code-block:: python

    from Pynite.Rendering import Renderer
    my_renderer = Renderer(my_model)
    my_renderer.render_model()

**VTK Renderer (Alternative)**

The ``Visualization`` module provides direct VTK rendering via ``Renderer`` class. This renderer offers similar features with a slightly different API and may be preferred for:

- Lightweight environments with minimal dependencies
- Specific VTK customization needs
- Legacy code compatibility

.. code-block:: python

    from Pynite.Visualization import Renderer
    my_renderer = Renderer(my_model)
    my_renderer.render_model()

The Renderer Object
===================

Both renderers expose a ``Renderer`` class that controls visualization:

.. code-block:: python

    # PyVista renderer (recommended)
    from Pynite.Rendering import Renderer
    my_rndr = Renderer(my_model)

    # VTK renderer (alternative)
    from Pynite.Visualization import Renderer
    my_rndr = Renderer(my_model)

Each renderer provides the same core visualization capabilities through a consistent API.


`Renderer` Properties and Attributes
====================================

The ``Renderer`` object has many properties you can adjust to customize your visualization. Below is a comprehensive guide:

Geometry Visibility
-------------------

Control what parts of the model are displayed:

.. code-block:: python

    my_rndr.render_nodes = True       # Show/hide nodes (default: True)
    my_rndr.render_loads = True       # Show/hide load glyphs (default: True)
    my_rndr.labels = True             # Show/hide text labels (default: True, PyVista uses show_labels)

Deformation and Scale
---------------------

Display and control the deformed shape of the model:

.. code-block:: python

    my_rndr.deformed_shape = True     # Toggle deformed shape display (default: False)
    my_rndr.deformed_scale = 30       # Scale factor for deformations (default: 30)

Larger scale factors exaggerate deformations for visibility in stiff models.

Annotation and Label Sizing
---------------------------

Control the size of visual elements relative to your model:

.. code-block:: python

    my_rndr.annotation_size = 5       # Size of text, nodes, and support glyphs (default: auto)

When set to ``None`` (or not set), annotation size is automatically computed as 5% of the shortest distance between nodes. You can manually set it to override:

.. code-block:: python

    my_rndr.annotation_size = 2       # Force smaller annotations

Load Combination and Case Selection
-----------------------------------

Choose which load results to display:

.. code-block:: python

    my_rndr.combo_name = 'Load Combo 1'  # Render results for a load combination
    my_rndr.case = None                  # (Mutually exclusive with combo_name)
    
    # To switch to a load case:
    my_rndr.case = 'Dead Load'           # Automatically clears combo_name
    my_rndr.combo_name = None

**Note:** Setting ``combo_name`` clears ``case`` and vice versa. Deformed shapes are only available for load combinations, not load cases.

Plate and Quad Contours
-----------------------

Display stress and result contours on plate and quad elements:

.. code-block:: python

    my_rndr.color_map = None          # No contours (default)
    my_rndr.color_map = 'Mx'          # Moment for bending along (not about) x-axis
    my_rndr.color_map = 'My'          # Moment for bending along (not about) y-axis
    my_rndr.color_map = 'Mxy'         # Twisting moment
    my_rndr.color_map = 'Qx'          # Shear for bending along (not about) the x-axis
    my_rndr.color_map = 'Qy'          # Shear for bending along (not about) the y-axis
    my_rndr.color_map = 'Sx'          # Membrane stress in x-direction
    my_rndr.color_map = 'Sy'          # Membrane stress in y-direction
    my_rndr.color_map = 'Txy'         # In-plane shear stress

Scalar Bar for Contours
-----------------------

Display a legend showing contour value ranges:

.. code-block:: python

    my_rndr.scalar_bar = False             # Hide scalar bar (default)
    my_rndr.scalar_bar = True              # Show scalar bar
    my_rndr.scalar_bar_text_size = 24      # Font size for scalar bar labels (default: 24)

Visual Theme
------------

Switch between rendering themes:

.. code-block:: python

    my_rndr.theme = 'default'         # Dark background (default)
    my_rndr.theme = 'print'           # White background (better for publications)

The ``'print'`` theme uses black lines and text on white, ideal for generating publication-quality images.

Window Size
-----------

Control the render window dimensions for both renderers:

.. code-block:: python

    my_rndr.window_width = 1920
    my_rndr.window_height = 1080

Both VTK and PyVista renderers support ``window_width`` and ``window_height`` properties to set the initial window size.

Internal Force and Moment Diagrams
==================================

You can visualize internal forces and moments directly on members using member diagrams. This is useful for quick visual inspection of shear, moment, axial, and torque distributions along members. Each diagram automatically displays labels showing the maximum and minimum values.

.. code-block:: python

    # Display a moment diagram in the z-direction
    my_rndr.member_diagrams = 'Mz'
    my_rndr.diagram_scale = 30      # Adjust diagram size (default: 30)
    my_rndr.render_model()

Available Diagram Types
-----------------------

The ``member_diagrams`` property accepts the following values:

- ``'Fy'`` - Shear force in local y-direction
- ``'Fz'`` - Shear force in local z-direction
- ``'My'`` - Bending moment about local y-axis
- ``'Mz'`` - Bending moment about local z-axis
- ``'Fx'`` - Axial force along member
- ``'Tx'`` - Torsional moment about member axis
- ``None`` - No diagrams (default)

Diagram Display Features
------------------------

- **Diagram Lines**: Visual representation of internal force/moment distribution along members
- **Value Labels**: Maximum and minimum values are automatically labeled on the diagram
- **Color Coding**: Diagrams use cyan in default theme and black in print theme for visibility

Diagram Customization
---------------------

The diagram visualization can be customized with the following property:

- ``diagram_scale`` (float, default=30) - Controls the size of the diagrams. Larger values make diagrams more visually prominent and easier to read labels.

.. code-block:: python

    # Example: Display axial force diagrams with custom scale
    my_rndr.member_diagrams = 'Fx'
    my_rndr.diagram_scale = 50       # Larger diagrams with larger labels
    my_rndr.render_model()

Example: Viewing Different Diagrams
------------------------------------

You can easily switch between different diagram types:

.. code-block:: python

    # View moment diagram with max/min labels
    my_rndr.member_diagrams = 'Mz'
    my_rndr.render_model()

    # Switch to shear diagram
    my_rndr.member_diagrams = 'Fz'
    my_rndr.render_model()

    # Display axial forces
    my_rndr.member_diagrams = 'Fx'
    my_rndr.render_model()

    # Display torque
    my_rndr.member_diagrams = 'Tx'
    my_rndr.render_model()

    # Turn off diagrams
    my_rndr.member_diagrams = None
    my_rndr.render_model()

Rendering the Model and Screenshots
===================================

Core Rendering Operations
---------------------------

**render_model()**

Renders the finite element model with the current visualization settings in an interactive window or Jupyter notebook.

.. code-block:: python

    my_rndr.render_model()

**screenshot(file_path)**

Captures the current render window as an image and saves it to disk.

.. code-block:: python

    my_rndr.screenshot('path/to/output.png')
    my_rndr.screenshot('C:\\Users\\username\\Desktop\\model.png')

**update()**

Refreshes the visualization after changing renderer properties without re-rendering the entire model.

.. code-block:: python

    # Change visualization properties
    my_rndr.deformed_shape = True
    my_rndr.deformed_scale = 50
    
    # Quickly refresh the display
    my_rndr.update()

Use ``update()`` when you've changed properties like ``deformed_scale``, ``annotation_size``, ``color_map``, or load combination and want to see changes immediately without rebuilding the entire scene.

Interactive Workflow Example
------------------------------

.. code-block:: python

    from Pynite.Rendering import Renderer  # or Pynite.Visualization
    
    # Create renderer
    my_rndr = Renderer(my_model)
    
    # View undeformed shape
    my_rndr.render_model()
    
    # Add deformation and quickly refresh
    my_rndr.deformed_shape = True
    my_rndr.deformed_scale = 20
    my_rndr.update()
    
    # Switch diagram type
    my_rndr.member_diagrams = 'Mz'
    my_rndr.update()
    
    # Save result
    my_rndr.screenshot('mz_diagram.png')

Renderer Class Reference
========================

.. autoclass:: Pynite.Rendering.Renderer
   :members:
   :undoc-members:

.. autoclass:: Pynite.Visualization.Renderer
   :members:
   :undoc-members:

Common Visualization Tasks
=========================

**View Deformed Shape with Bending Moments**

.. code-block:: python

    my_rndr.deformed_shape = True
    my_rndr.deformed_scale = 25
    my_rndr.member_diagrams = 'Mz'
    my_rndr.diagram_scale = 20
    my_rndr.render_model()

**Generate Publication-Quality Image**

.. code-block:: python

    my_rndr.theme = 'print'           # White background
    my_rndr.annotation_size = 3       # Smaller annotations for clarity
    my_rndr.scalar_bar = True         # Show contour legend
    my_rndr.member_diagrams = None    # Remove diagrams for clean look
    my_rndr.render_model()
    my_rndr.screenshot('figure_01.png')

**Compare Multiple Load Cases**

.. code-block:: python

    load_cases = ['Dead Load', 'Live Load', 'Wind Load']
    
    for case_name in load_cases:
        my_rndr.case = case_name
        my_rndr.member_diagrams = 'Fx'
        my_rndr.render_model()
        my_rndr.screenshot(f'{case_name}.png')

**Visualize Plate Stresses**

.. code-block:: python

    my_rndr.color_map = 'Mxy'         # Twisting moment
    my_rndr.scalar_bar = True
    my_rndr.deformed_shape = False
    my_rndr.render_nodes = False
    my_rndr.render_model()