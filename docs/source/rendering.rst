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

.. autoclass:: PyNite.Rendering.Renderer
   :members:
   :undoc-members: