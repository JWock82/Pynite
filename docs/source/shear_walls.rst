===========
Shear Walls
===========

.. _shear_walls:

The ``ShearWall`` helper creates and manages an in-plane shell-element model of
a structural shear wall. It generates the wall mesh, optional wall returns,
supports, and story loads; identifies wall piers and coupling beams around
openings; and reports component forces after analysis.

Create a shear wall through :meth:`Pynite.FEModel3D.FEModel3D.add_shear_wall`.
The returned helper is stored in ``model.shear_walls`` under the provided name.

Coordinate System
=================

Shear-wall dimensions and all wall-specific coordinates use a local x-y system:

* Local ``x`` runs along the wall length, from its left edge.
* Local ``y`` runs upward, from the wall base.
* Local ``z`` is normal to the wall plane.

Set ``plane='XY'``, ``'XZ'``, or ``'YZ'`` and use ``origin`` to position the
local system in the global model. In-plane shear acts in the local x direction.

Workflow
========

Define the wall geometry and features before analysis. Calling ``generate()``
explicitly is optional: the model regenerates a shear wall automatically when
it has changed before an analysis.

.. code-block:: python

    from Pynite import FEModel3D

    model = FEModel3D()
    model.add_material('Concrete', E=3600 * 144, G=1500 * 144, nu=0.2, rho=0.150)

    wall = model.add_shear_wall(
        'Wall 1',
        mesh_size=1.0,
        length=20.0,
        height=15.0,
        thickness=1.0,
        material_name='Concrete',
        ky_mod=0.35,
        plane='XY',
        origin=[0.0, 0.0, 0.0],
    )

    # Add a 14 ft by 12 ft opening, leaving 3 ft boundary piers and a 3 ft beam.
    wall.add_opening('Opening 1', x_start=3.0, y_start=0.0, width=14.0, height=12.0)
    wall.add_support(elevation=0.0)
    wall.add_story('Roof', elevation=15.0)
    wall.add_shear('Roof', force=100.0, case='E')

    model.add_load_combo('1.0E', {'E': 1.0})
    model.analyze_linear()

Openings, Flanges, and Materials
================================

Use :meth:`~Pynite.ShearWall.ShearWall.add_opening` to place rectangular
openings in local coordinates. The wall automatically detects piers and
coupling beams after it is generated. Identifiers are assigned based on the
resulting wall regions, so use the ``piers`` and ``coupling_beams`` dictionaries
or the layout plots to confirm the generated IDs.

Use :meth:`~Pynite.ShearWall.ShearWall.add_flange` to create a wall return at a
local x coordinate. Use
:meth:`~Pynite.ShearWall.ShearWall.asign_material` to apply a different
material and thickness to a rectangular region of the wall. The method name is
spelled ``asign_material`` in the public API.

Stories, Loads, and Stiffness
=============================

A story defines a horizontal diaphragm level and the portion of the wall it
loads. ``add_shear`` distributes an in-plane force across that story; positive
force acts along positive local x. ``add_axial`` distributes a downward axial
force across the story.

Each call to ``add_story`` also creates a 100-unit in-plane test load used by
:meth:`~Pynite.ShearWall.ShearWall.stiffness` after analysis:

.. code-block:: python

    roof_stiffness = wall.stiffness('Roof')

Component Results
=================

After analysis, each automatically identified pier is available in
``wall.piers`` and each coupling beam in ``wall.coupling_beams``. Both report
``(P, M, V, shear_span_ratio)``. The force sign convention matches
:class:`~Pynite.PhysMember.PhysMember` internal member results.

Piers
-----

Piers report at their ``'bottom'`` and ``'top'`` edges. The final value is
``M/(V*L)``, where ``L`` is the pier length.

.. code-block:: python

    pier = wall.piers['P1']
    P_bottom, M_bottom, V_bottom, ratio_bottom = pier.sum_forces('1.0E', 'bottom')
    P_top, M_top, V_top, ratio_top = pier.sum_forces('1.0E', 'top')

    wall.print_piers('1.0E')  # prints one Bottom and one Top row per pier

Coupling Beams
--------------

Coupling beams report at their ``'left'`` and ``'right'`` edges. The final
value is ``M/(V*H)``, where ``H`` is the coupling-beam height.

.. code-block:: python

    beam = wall.coupling_beams['B1']
    P_left, M_left, V_left, ratio_left = beam.sum_forces('1.0E', 'left')
    P_right, M_right, V_right, ratio_right = beam.sum_forces('1.0E', 'right')

    wall.print_coupling_beams('1.0E')  # prints one Left and one Right row per beam

Layout Plots and Screenshots
============================

Use ``draw_piers`` and ``draw_coupling_beams`` to inspect the regions identified
by the wall helper. By default they return a pyplot object that can be saved;
pass ``show=True`` to display the layout directly.

.. code-block:: python

    wall.draw_piers(show=True)
    wall.draw_coupling_beams(show=True)
    wall.screenshots('1.0E', dir_path='results', renderer_backend='vtk')

``screenshots`` writes two shear-contour images plus the pier and coupling-beam
layout images. It supports the ``'vtk'`` renderer by default or ``'pyvista'``.

API Reference
=============

.. autoclass:: Pynite.ShearWall.ShearWall
   :members:
   :undoc-members:

.. autoclass:: Pynite.ShearWall.Pier
   :members:

.. autoclass:: Pynite.ShearWall.CouplingBeam
   :members:
