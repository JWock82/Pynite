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

Once added to the model, the node will be stored in the ``nodes`` dictionary of the model for easy access later on.

.. code-block:: python

    # Get the node we just created from the model
    node_1 = my_model.nodes['N1']

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

Enforced Nodal Displacements (e.g., Support Settlements)
=======================================================================
Use the ``FEModel3D.def_node_disp`` method to model a known nodal displacement or rotation such as a support settlement.

.. code-block:: python

    # Add an enforced displacement of -1.5 at node 'N4' in global Y
    my_model.def_node_disp('N4', direction='DY', magnitude=-1.5)

Getting Node Results
====================

After running an analysis, each node stores its results in dictionaries keyed by load
combination name. Access a node via ``model.nodes['Name']`` and then query the result
you need by combination.

Available nodal results
-----------------------

- Displacements (global)
  - ``DX['Combo']``, ``DY['Combo']``, ``DZ['Combo']``
- Rotations (global)
  - ``RX['Combo']``, ``RY['Combo']``, ``RZ['Combo']``
- Reactions (global forces)
  - ``RxnFX['Combo']``, ``RxnFY['Combo']``, ``RxnFZ['Combo']``
- Reactions (global moments)
  - ``RxnMX['Combo']``, ``RxnMY['Combo']``, ``RxnMZ['Combo']``

Notes
-----

- Results are populated only after a successful call to `analyze_linear()`, `analyze()`,
  or `analyze_PDelta()`. Each combination analyzed will appear as a key in the
  per-result dictionaries.
- Reaction values are reported in the global axes and represent the support reactions
  required to enforce supports and enforced displacements/rotations at the node.
- Enforced displacements/rotations are inputs (`EnforcedDX/DY/DZ/RX/RY/RZ` on the node); the
  corresponding reaction response will appear in `RxnF*`/`RxnM*` after analysis.
- Spring supports are defined per-DOF and stored as `spring_D* = [stiffness, direction, active]`
  where `direction` is `'+'` (compression-only), `'-'` (tension-only), or `None` (two-way).

Examples
--------

.. code-block:: python

    n = my_model.nodes['N2']

    # Global displacements/rotations for combo 'D+L'
    ux = n.DX['D+L']
    uy = n.DY['D+L']
    rz = n.RZ['D+L']

    # Support reactions for combo '1.2D+1.0W'
    vy = n.RxnFY['1.2D+1.0W']
    mz = n.RxnMZ['1.2D+1.0W']

    # Inspect loads assigned to the node (list of tuples: (Direction, value, case))
    for (dirn, val, case) in n.NodeLoads:
        print(dirn, val, case)
The ``FEModel3D`` class stores nodes in a Python dictionary. Nodes can be accessed using the syntax ``FEModel3D.nodes['node_name']``.

Once you've retrieved a node you can access its reactions and displacements as node class attributes. Reactions and displacements are also stored in dictionary format, with the keys being the load combination names.

.. code-block:: python

    # Print the Y-reaction and the reaction moment about the Z-axis at nodes 'N2' and 'N3'
    print(my_model.nodes['N2'].RxnFY['1.2D+1.0W'])
    print(my_model.nodes['N3'].RxnMZ['1.2D+1.0W'])

Node API Quick Reference
========================

- Create: ``FEModel3D.add_node(name, X, Y, Z) -> str``
- Supports: ``FEModel3D.def_support(node_name, support_DX=False, support_DY=False, support_DZ=False, support_RX=False, support_RY=False, support_RZ=False)``
- Spring supports: ``FEModel3D.def_support_spring(node_name, dof, stiffness, direction=None)``
- Enforced disp/rot: ``FEModel3D.def_node_disp(node_name, direction, magnitude)``
- Loads: ``FEModel3D.add_node_load(node_name, direction, P, case='Case 1')``
- Results on a node object: ``DX/DY/DZ/RX/RY/RZ['Combo']``, ``RxnFX/RxnFY/RxnFZ/RxnMX/RxnMY/RxnMZ['Combo']``

Worked Example
==============

.. code-block:: python

    from Pynite import FEModel3D

    model = FEModel3D()
    n1 = model.add_node('N1', 0, 0, 0)
    n2 = model.add_node('N2', 0, 10, 0)

    # Pinned base and vertical settlement at the top
    model.def_support('N1', support_DX=True, support_DY=True, support_DZ=True)
    model.def_node_disp('N2', 'DY', -0.25)

    # Add a small horizontal force
    model.add_node_load('N2', 'FX', 2.0, case='S')
    model.add_load_combo('S', {'S': 1.0})
    model.analyze_linear()

    print(model.nodes['N2'].DX['S'])     # horizontal drift
    print(model.nodes['N1'].RxnFY['S'])  # base vertical reaction

API Cross-Links and Reference
=============================

- Class: :class:`~Pynite.Node3D.Node3D`

.. autoclass:: Pynite.Node3D.Node3D
   :members:
   :undoc-members:
   :show-inheritance:


