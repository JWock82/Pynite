=========
Rendering
=========

`Pynite` has a built-in ``Rendering`` library to help you visualize your model. This library can be used to produce basic rendererings of the model, loads, the deformed shape, and can generate interactive plots in `Jupyter Notebooks`.

The ``Rendering`` library relies on `Pyvista` which in turn relies on `The Visualization Toolkit (VTK)` by `Kitware`. This library was chosen because it is cross-platform, widely used, regularly maintained, and reliable.

The `Renderer` object
=====================

`Pynite` uses a ``Renderer`` object in the ``Rendering`` library to do the heavy lifting for you. To create a ``Renderer`` for your model you just create a new instance of the ``Renderer`` class and pass your model to it:

.. code-block:: python

    # Import the `Renderer` object
    from Pynite.Rendering import Renderer

    # Create a `Renderer` for this model
    my_rndr = Renderer(my_model)


`Renderer` Options
==================

The ``Renderer`` object has some default properties you may wish to change. Below are examples of these options:

.. code-block:: python
    
    my_rndr.deformed_shape = True      # Turns on the deformed shape during rendering
    my_rndr.annotation_size = 2        # Scales model annotations to 2 model length units
    my_rndr.deformed_scale = 30        # Sets the scale of deformations
    my_rndr.render_nodes = False       # Turns off the nodes during rendering
    my_rndr.render_loads = False       # Turns off load rendering
    my_rndr.color_map = 'Txy'          # Changes plate contours to show shear stresses
    my_rndr.combo_name = '1.4D'        # The name of the load combination to render
    my_rndr.case = None                # The load case to render - ignored if a load combination is specified
    my_rndr.labels = False             # Turns off label renderings
    my_rndr.scalar_bar = True          # Turns on the scalar bar for plate contours
    my_rndr.scalar_bar_text_size = 24  # Adjusts the text height for the scalar bar

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

Once you've got your settings the way you want them, you can render your model using ``render_model``:

.. code-block:: python

    my_rndr.render_model()

You can also create a screenshot using ``screenshot``. Pass the ``filepath`` for the screenshot as a string. Use ``interact=True`` if you want to position the screenshot before capturing it, then it ``q`` to capture it. If you close out of the window instead of hitting ``q`` it will capture the original positioning instead (this is a ``Pyvista`` nuance).

.. code-block:: python

    my_rndr.screeenshot(filepath, interact=True)

Renderer Class Reference
========================

.. autoclass:: Pynite.Rendering.Renderer
   :members:
   :undoc-members: