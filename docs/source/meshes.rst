
======
Meshes
======

.. _meshes:

PyNite provides several utilities for generating 2D shell/plate **meshes** and
curved surface meshes (annuli, cylinders, frustra) that automatically create
:class:`~Pynite.Node3D.Node3D` nodes and :class:`~Pynite.Quad3D.Quad3D` (or
:class:`~Pynite.Plate3D.Plate3D`) elements, add them to your
:class:`~Pynite.FEModel3D.FEModel3D` model, and keep track of the mesh's
membership for post-processing.

This page explains the mesh classes, common options, openings, numbering,
and convenience result helpers (:meth:`max_moment`, :meth:`min_moment`,
:meth:`max_shear`, :meth:`min_shear`, :meth:`max_membrane`, :meth:`min_membrane`).

Quick Start
===========

Create a rectangular quad mesh on the global XY plane and add it to a model::

    from Pynite.FEModel3D import FEModel3D
    from Pynite.Mesh import RectangleMesh

    model = FEModel3D()
    model.add_material('Conc', E=3600*144, nu=0.2, rho=0.0)  # example units

    # Add a rectangular mesh to the model
    model.add_rectangle_mesh(
      name='Slab',
      mesh_size = 2.0,
      width = 20.0,
      height = 25.0,
      thickness = 0.5,
      material_name = 'Conc',
      kx_mod = 1.0, ky_mod = 1.0,
      origin = [0, 0, 0],
      plane = 'XZ'
    )

    # Rectangular meshes are unique, in that you can add openings to them
    model.meshes['Slab'].add_rect_opening(
      name = 'Opening 1',
      x_left = 5.0,  # The mesh's local x-y coordinate system
      y_bott = 6.0,
      width = 5.0,
      height = 7.0
    )

    # Once all openings are defined the mesh can be generated
    model.meshes['Slab'].generate()  # builds nodes/elements and adds them to the model

    # Now that the nodes and elements in the mesh have been generated, we can add supports.
    # Let's add fixed supports to the perimeter
    for node in model.meshes['Slab'].nodes.values():
      # Check if the node lies on the perimeter
      if node.X == 0.0 or node.X == 20.0 or node.Z == 0.0 or node.Z == 25.0:  # The mesh was generated in the XZ plane
        model.def_support(node.name, True, True, True, True, True, True)  # Fix all degrees of freedom at the node

    # Solve the model
    model.analyze()

    # Solve the model, then query peak values by combo/tag
    slab = model.meshes['Slab']; mx_max = slab.max_moment('Mx', combo_tags=['service', 'strength'])  # Check tagged combinations
    qy_min = slab.min_shear('Qy', combo_tags='Combo 1')  # Check just one combination

Why Use Mesh Classes?
=====================

* They **parametrically** create nodes and elements with consistent numbering.
* They **add** created nodes/elements directly to ``model.nodes``, ``model.quads``
  and/or ``model.plates``. If a proposed name already exists in the model,
  duplicates are automatically **renamed** to unique IDs to avoid collisions.
* They retain their own ``nodes`` and ``elements`` dicts, so you can post‑process
  a region without sifting through the entire model.
* They include **result helpers** that scan corner and center points for max/min
  membrane stresses, bending moments, and shears per combination/tag filters.

Common Concepts
===============

Coordinate Planes / Axes
~~~~~~~~~~~~~~~~~~~~~~~~

Most mesh generators create geometry on a named plane or around a principal axis:

* Planar meshes: ``plane='XY'`` | ``'YZ'`` | ``'XZ'``
* Axisymmetric/curved families (annulus, cylinder, frustrum): ``axis='X'|'Y'|'Z'``

Element Type
~~~~~~~~~~~~

``element_type='Quad'`` produces :class:`~Pynite.Quad3D.Quad3D` (DKMQ + membrane)
elements (recommended). ``'Rect'`` produces :class:`~Pynite.Plate3D.Plate3D`
(rectangular plate bending + membrane).

Orthotropy Modifiers
~~~~~~~~~~~~~~~~~~~~

``kx_mod`` and ``ky_mod`` scale in‑plane stiffness along the element's **local**
x and y axes (membrane only). ``1.0`` means isotropic behavior.

Numbering & Names
~~~~~~~~~~~~~~~~~

* Provide ``start_node='N1'`` and ``start_element='Q1'`` (or ``'R1'``) to seed
  names. The mesh will increment from those numbers.
* If a generated name already exists in the model, it is automatically renamed
  to a **unique** name and both the mesh and model dicts are updated to match.

Openings (RectangleMesh)
~~~~~~~~~~~~~~~~~~~~~~~~

Rectangular openings can be declared **before** ``generate()`` and are respected
during meshing. Nodes and elements fully inside an opening are removed, and
orphans are cleaned up::

    model.meshes['Slab'].add_rect_opening('Stair 1', x_left=6.0, y_bott=3.0, width=3.0, height=4.0)
    model.meshes['Slab'].generate()

Load Combinations / Tag Filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The result helpers accept either:

* a **single** combination name: ``combo_tags='Combo 1'``
* a list of **tags**: ``combo_tags=['service', 'wind']`` (any combo having at
  least one of the tags is evaluated)

After solution, results are sampled at the **four corners** plus the **center**
of each element for the requested component and compared across all eligible
combos.

API Overview
============

Base Class: :class:`Pynite.Mesh.Mesh`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: Pynite.Mesh.Mesh
   :members:
   :undoc-members:
   :show-inheritance:

Concrete Classes
~~~~~~~~~~~~~~~~

RectangleMesh
^^^^^^^^^^^^^

A structured grid on a rectangular region, optionally with internal control
lines and rectangular openings.

.. autoclass:: Pynite.Mesh.RectangleMesh
   :members:
   :undoc-members:
   :show-inheritance:

**Parameters (key):**

* ``mesh_size`` – target element edge size (controls i/j divisions between control points)
* ``width, height`` – dimensions in the mesh's **local** x/y
* ``origin`` – local (0,0,0) in global coords
* ``plane`` – ``'XY'|'YZ'|'XZ'``
* ``x_control, y_control`` – optional internal control coordinates to be
  honored exactly (deduped with tolerance)
* ``element_type`` – ``'Quad'`` or ``'Rect'``

**Openings:** See :meth:`add_rect_opening` and :class:`RectOpening`.

RectOpening
^^^^^^^^^^^

.. autoclass:: Pynite.Mesh.RectOpening
   :members:
   :undoc-members:
   :show-inheritance:

AnnulusMesh
^^^^^^^^^^^

Builds concentric **rings** of quads from an inner radius to an outer radius.
The algorithm increases mesh density as needed; transitions are handled with
:class:`AnnulusTransRingMesh` rings.

.. autoclass:: Pynite.Mesh.AnnulusMesh
   :members:
   :undoc-members:
   :show-inheritance:

AnnulusRingMesh
^^^^^^^^^^^^^^^

Single annular ring of quads between two radii and about an axis.

.. autoclass:: Pynite.Mesh.AnnulusRingMesh
   :members:
   :undoc-members:
   :show-inheritance:

AnnulusTransRingMesh
^^^^^^^^^^^^^^^^^^^^

Annular ring that **triples** element count between inner and mid/outer radii
to smoothly transition to a finer mesh.

.. autoclass:: Pynite.Mesh.AnnulusTransRingMesh
   :members:
   :undoc-members:
   :show-inheritance:

FrustrumMesh
^^^^^^^^^^^^

A conical frustum generated by first creating an :class:`AnnulusMesh` then
adjusting node heights linearly with radius to achieve the specified frustum
height along the chosen ``axis``.

.. autoclass:: Pynite.Mesh.FrustrumMesh
   :members:
   :undoc-members:
   :show-inheritance:

CylinderMesh
^^^^^^^^^^^^

A vertical (or X/Z‐axis) cylindrical surface meshed into circumferential
courses and vertical divisions. Choose ``element_type='Quad'`` or ``'Rect'``.

.. autoclass:: Pynite.Mesh.CylinderMesh
   :members:
   :undoc-members:
   :show-inheritance:

Result Helpers (per Mesh)
=========================

All helpers scan **corner** and **center** points per element for eligible
combinations.

* :meth:`max_moment(direction='Mx', combo_tags=...)` – returns **max** of
  ``Mx``, ``My``, ``Mxy`` or global ``MX, MY, MZ``.
* :meth:`min_moment(direction='Mx', combo_tags=...)` – returns **min**.
* :meth:`max_shear(direction='Qx', combo_tags=...)` – returns **max** of
  ``Qx, Qy`` (local) or ``QX, QY`` (global).
* :meth:`min_shear(direction='Qx', combo_tags=...)` – returns **min**.
* :meth:`max_membrane(direction='Sx', combo_tags=...)` – returns **max** of
  membrane ``Sx, Sy, Sxy`` (local) or ``SX, SY`` (global mem. components).
* :meth:`min_membrane(direction='Sx', combo_tags=...)` – returns **min**.

**Notes**

* Directions ending in **uppercase** (``'MX'``, ``'QY'``, ``'SX'``/``'SY'``)
  indicate **global** components where implemented; lowercase are **local**.
* If no eligible combo is found, helpers return ``0.0``.

Examples
========

Rectangular slab with openings and control lines::

    model = FEModel3D()

    model.add_rectangle_mesh(
      name = 'Slab',
      mesh_size = 1.0,
      width = 36.0, height = 24.0,
      thickness = 0.667,
      material_name = 'Conc',
      kx_mod = 1.0, ky_mod = 1.0,
      origin = [0, 12, 0],
      plane = 'XZ',
      x_control=[12.0, 24.0], y_control=[12.0],
      element_type = 'Quad'
    )

    model.meshes['Slab'].add_rect_opening('Elevator', x_left=8.0, y_bott=6.0, width=4.0, height=5.0)
    model.meshes['Slab'].add_rect_opening('Stair',    x_left=22.0, y_bott=6.0, width=5.0, height=9.0)
    model.meshes['Slab'].generate()

    # After solving combinations with tags 'service' and 'strength'
    mxy_max = model.meshes['Slab'].max_moment('Mxy', combo_tags=['service', 'strength'])

Annular tank roof shell (about global Y‑axis)::

    model.add_annulus_mesh(name='AM1',
                           mesh_size = 1.5,
                           outer_radius = 25.0,
                           inner_radius = 5.0,
                           thickness = 0.375,
                           material_name = 'Steel',
                           kx_mod = 1.0, ky_mod = 1.0,
                           origin = [0, 0, 0],
                           axis = 'Y'
                           )
    
    model.meshes['AM1'].generate()

Cylindrical wall (about global Z-axis)::

    model.add_cylinder_mesh(name = 'CM1',
                            mesh_size = 2.0,
                            radius = 20.0,
                            height = 18.0,
                            thickness = 0.4,
                            material_name = 'Conc',
                            kx_mod = 1.0, ky_mod = 1.0,
                            origin = [0, 0, 0],
                            axis = 'Z',
                            num_elements = 30
                            element_type = 'Quad'
                            )

    model.meshes['CM1'].generate()

Implementation Details & Behavior
=================================

* **Duplicate Names:** During ``generate()``, if any node/element name is found
  to already exist in the model, it is **renamed** (e.g., by calling
  ``model.unique_name(...)``) and the mesh’s own ``nodes/elements`` dicts are
  updated so both remain consistent.

* **Model Integration:** After generation, all mesh nodes/elements are inserted
  into the owning model’s dictionaries (``model.nodes``, ``model.quads`` and/or
  ``model.plates``). The mesh sets ``is_generated=True`` when finished.

* **Element Local Coordinates:** Quads compute local axes from i→j and the
  element plane; results sampling uses natural coordinates ``(xi, eta)=(-1,1)``
  at corners and their average for the center.

* **Orthotropy:** ``kx_mod`` / ``ky_mod`` affect only the **membrane** (plane
  stress) stiffness matrix ``Cm``; shear modulus is unaffected.

Troubleshooting & Tips
======================

* For **very skewed** quads, monitor Jacobian signs; extremely distorted shapes
  can trigger warnings in element routines.
* Use ``num_elements`` on curved meshes (e.g., :class:`CylinderMesh`) to make
  sure seams align with adjacent meshes.
* When combining meshes, prefer consistent ``start_node`` / ``start_element``
  ranges to keep naming readable; duplicates will auto‑rename if needed.
* For openings close to edges or control lines, consider mesh sizes that avoid
  sliver elements.

Reference
=========

This page documents the mesh API and typical usage patterns. See also:

* :mod:`Pynite.Quad3D` — DKMQ + membrane quadrilateral formulation
* :mod:`Pynite.Plate3D` — rectangular plate + membrane element
* :mod:`Pynite.Node3D` — nodal DOFs and result storage


