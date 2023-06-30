=====
Nodes
=====

Each node has 6 degrees of freedom (3 translations and 3 rotations). Each of these degrees of freedom needs to be restrained for the model to be stable. For example, applying full end releases to all the members framing into a node will allow the node to "spin". The node would either have to be attached to at least one of the members rotationally, or supported rotationally by a support.

Creating a New Node
===================

Nodes can be added to an existing model using the ``FEModel3D.add_node`` method:

.. code-block:: python

    # Add a node named 'N1' at location (0, 0, 0)
    my_model.add_node('N1', 0, 0, 0)

Once added to the model, the node will be stored in the ``Nodes`` dictionary of the model for easy access later on.

.. code-block:: python

    # Get the node we just created from the model
    node_1 = my_model.Nodes['N1']

Adding Supports to a Node
=========================

You can define supports for any node you've created using the ``def_support`` method. You'll need to specify the node name, and a ``True`` or ``False`` value for each support at the node. The supports are specified in the following order: X-translation, Y-translation, Z-translation, rotation about X, rotation about Y, rotation about Z.

You can also use keyword arguments if you only want to specify specific supports. Any supports you do not specify will be set to ``False``.

If supports have already been defined for a node they will be overwritten when ``def_support`` is called again for that node.

.. code-block:: python

    # Add a pinned support at node `N1`
    my_model.def_support('N1', True, True, True, False, False, False)

    # Another way to add a pinned support at node `N1`
    my_model.def_support('N1', 1, 1, 1, 0, 0, 0)

    # Yet another way to add a pinned support at node 'N1'
    my_model.def_support('N1', support_DX=True, support_DY=True, support_DZ=True)

    # Turn node 'N1' into a rotational support about the Y-axis.
    my_model.def_support('N1', support_RY=True)

Adding Nodal Spring Supports
============================

Nodal spring supports can be defined using the ``def_support_spring`` method. Nodal spring supports can be tension-only, compression-only, or two-way. Compression-only springs can be very useful for modeling foundation problems. When using a tension-only or compression-only spring the solution is iterative, so be sure not to use the ``FEModel3D.analyze_linear()`` solver.

You'll need to specify which node to apply the spring to, which degree of freedome to apply it to (``'DX'``, ``'DY'``, ``'DZ'``, ``'RX'``, ``'RY'``, or ``'RZ'``), the stiffness value, and a direction (``'-'`` = tension-only, ``'+'`` = compression-only, ``None`` = two-way).

Note that support springs are defined individually for each degree of freedom instead of for all degrees of freedom at the same time. Unlike with regular supports, defining subsequent support springs for the same node will not override prior assignments unless the assignment is to the same degree of freedom.

.. code-block:: python

    # Add a compression only spring to node `N1` with a stiffness of 5 in the Y-direction
    my_model.def_support_spring('N1', dof='DY', stiffness=5, direction='+')

    # Add a rotational two-way spring to node `N1`
    my_model.def_support_spring('N1', dof='RX', stiffness=3, direction=None)

Adding Nodal Loads
==================

Use the ``FEModel3D.add_node_load`` method to add nodal loads to a model.

.. code-block:: python

    # Add a moment about the global Z axis to node `N1` for load case 'D'
    my_model.add_node_load('N1', 'MZ', 20, 'D')

    # Add a force in the global X direction to node 'N1' for load case 'E'
    my_model.add_node_load('N1', 'FX', 30, 'E')

## Adding Node Displacements (e.g. Support Settlements)
Use the `AddNodeDisplacement` method to model a known nodal displacement, such as a support settlement.

Syntax:

    AddNodeDisplacement (self, Node, Direction, Magnitude): 
    
    Node : string
        The name of the node where the nodal displacement is being applied.
    Direction : string
        'DX' = Displacement in the global X-direction
        'DY' = Displacement in the global Y-direction
        'DZ' = Displacement in the global Z-direction
        'RX' = Rotation about the global X-axis
        'RY' = Rotation about the global Y-axis
        'RZ' = Rotation about the global Z-axis
    Magnitude : number
        The magnitude of the displacement.

Example:
    
    # Add a nodal displacement of -1.5 at node N4 in the global Y-direction
    myModel.AddNodeDisplacement('N4', 'DY', -1.5)

## Getting Node Results
The `FEModel3D` class stores nodes in a Python dictionary. Nodes can be accessed using the sytax `FEModel3D.Nodes['node_name']`.

Once you've retrieved a node you can access its reactions and displacements as node class attributes. Reactions and displacements are also stored in dictionary format, with the keys being the load combination names.

Examples:

    # Printing the Y-reaction and the reaction moment about the Z-axis at nodes "N2" and "N3" respectively
    print(myModel.Nodes['N2'].RxnFY['1.2D+1.0W'])
    print(myModel.Nodes['N3'].RxnMZ['1.2D+1.0W'])