============
FEModel3D
============

The ``FEModel3D`` class is the primary entry point for building and analyzing structural models in PyNite. This page provides an overview, common workflows, and links to the full API reference.

Quick start
-----------

.. code-block:: python

    # Import the modeling class
    from Pynite import FEModel3D

    # Create a new model
    model = FEModel3D()

    # Geometry
    model.add_node('N1', 0.0, 0.0, 0.0)
    model.add_node('N2', 3.0, 0.0, 0.0)

    model.add_material('A36', E=29_000_000.0, G=11_200_000.0, nu=0.3, rho=0.283)
    model.add_section('Wsect', A=10.0, Iy=100.0, Iz=200.0, J=5.0)

    m1 = model.add_member('M1', i_node='N1', j_node='N2', material_name='A36', section_name='Wsect')

    # Supports and loads
    model.def_support('N1', support_DX=True, support_DY=True, support_DZ=True, support_RY=True, support_RZ=True)
    model.add_node_load('N2', direction='FZ', P=-5.0, case='D')

    # Load combinations and analysis
    model.add_load_combo('1.0D', {'D': 1.0})
    model.analyze_linear(log=False)

    # Results
    uz = model.nodes['N2'].DZ['1.0D']
    rxn = model.nodes['N1'].RxnFZ['1.0D']

Core containers
---------------

Internally, models are organized into Python dictionaries for fast lookup by name. Keys are the user-provided names returned by the corresponding ``add_...`` methods. The most commonly used containers are:

- ``model.nodes``: Node objects keyed by node name.
- ``model.members``: Physical members keyed by member name. Each physical member may contain one or more sub-members for analysis.
- ``model.springs``: Discrete spring elements between nodes.
- ``model.plates``: Rectangular plate elements.
- ``model.quads``: General quadrilateral shell elements (MITC4 bending, isoparametric membrane).
- ``model.materials``: Material definitions.
- ``model.sections``: Section property definitions, including steel sections.
- ``model.meshes``: Parametric mesh helpers you can generate and then query for created nodes/elements.
- ``model.shear_walls``: Shear wall helpers that generate wall meshes and keep track of membership.
- ``model.load_combos``: Load combinations and factors.

Useful properties and utilities
-------------------------------

- ``model.load_cases``: Property returning a sorted list of all unique load-case names found in the model.
- ``model.solution``: String describing the last analysis type that ran (e.g., ``'Linear'``, ``'P-Delta'``).
- ``model.merge_duplicate_nodes(tolerance=0.001)``: Merges coincident nodes and rewires connected objects; returns removed node names.
- ``model.delete_node|delete_member|delete_spring(name)``: Remove items and invalidate results.

Building models
---------------

- Nodes: ``add_node(name, X, Y, Z)``
- Materials: ``add_material(name, E, G, nu, rho, fy=None)``
- Sections: ``add_section(name, A, Iy, Iz, J)`` or ``add_steel_section(name, A, Iy, Iz, J, Zy, Zz, material_name)``
- Members: ``add_member(name, i_node, j_node, material_name, section_name, rotation=0.0, tension_only=False, comp_only=False)``
- Plates/Quads:
  - ``add_plate(name, i_node, j_node, m_node, n_node, t, material_name, kx_mod=1.0, ky_mod=1.0)``
  - ``add_quad(name, i_node, j_node, m_node, n_node, t, material_name, kx_mod=1.0, ky_mod=1.0)``
- Springs: ``add_spring(name, i_node, j_node, ks, tension_only=False, comp_only=False)``

Meshing helpers
---------------

- ``add_rectangle_mesh(name, mesh_size, width, height, thickness, material_name, ..., element_type='Quad')``
- ``add_annulus_mesh(name, mesh_size, outer_radius, inner_radius, thickness, material_name, ...)``
- ``add_frustrum_mesh(name, mesh_size, large_radius, small_radius, height, thickness, material_name, ...)``
- ``add_cylinder_mesh(name, mesh_size, radius, height, thickness, material_name, ..., element_type='Quad')``
- ``add_shear_wall(name, mesh_size, length, height, thickness, material_name, ky_mod=0.35, plane='XY', origin=[0,0,0])``

Supports, releases, and imposed displacements
------------------------------------------------------------

- Supports: ``def_support(node_name, support_DX=False, support_DY=False, support_DZ=False, support_RX=False, support_RY=False, support_RZ=False)``
- Support springs: ``def_support_spring(node_name, dof, stiffness, direction=None)`` where ``dof`` ∈ {``DX``, ``DY``, ``DZ``, ``RX``, ``RY``, ``RZ``} and ``direction`` ∈ {``+``, ``-``, ``None``}.
- Member end releases: ``def_releases(member_name, Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=False, Rzi=False, Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=False, Rzj=False)``
- Enforced nodal displacements/rotations: ``def_node_disp(node_name, direction, magnitude)``

Loading
-------

- Nodal loads: ``add_node_load(node_name, direction, P, case='Case 1')``
- Member loads:
  - Point: ``add_member_pt_load(member_name, direction, P, x, case='Case 1')``
  - Distributed: ``add_member_dist_load(member_name, direction, w1, w2, x1=None, x2=None, case='Case 1')``
  - Self-weight (members only): ``add_member_self_weight(global_direction, factor, case='Case 1')``
- Surface pressures:
  - Plates: ``add_plate_surface_pressure(plate_name, pressure, case='Case 1')``
  - Quads: ``add_quad_surface_pressure(quad_name, pressure, case='Case 1')``
- Clear all loads: ``delete_loads()``

Load combinations
-----------------

- ``add_load_combo(name, factors, combo_tags=None)`` creates a named combination from load-case factors. Tags are useful for filtering during analysis (e.g., design groups).

Analysis options
----------------

- Linear elastic: ``analyze_linear(log=False, check_stability=True, check_statics=False, sparse=True, combo_tags=None)``
- General elastic with tension/compression-only iteration and optional load stepping: ``analyze(log=False, check_stability=True, check_statics=False, max_iter=30, sparse=True, combo_tags=None, spring_tolerance=0, member_tolerance=0, num_steps=1)``
- P-Delta geometric nonlinearity: ``analyze_PDelta(log=False, check_stability=True, max_iter=30, sparse=True, combo_tags=None)``

Results access
--------------

Most results are stored on the individual objects and keyed by load combination name. Examples:

- Nodal displacements/rotations: ``model.nodes['N2'].DX['Combo']``, ``.DY``, ``.DZ``, ``.RX``, ``.RY``, ``.RZ``
- Nodal reactions: ``model.nodes['N1'].RxnFX['Combo']``, ``.RxnFY``, ``.RxnFZ``, ``.RxnMX``, ``.RxnMY``, ``.RxnMZ``
- Member forces and diagrams are available on member objects; see member-specific documentation.

Naming and uniqueness
---------------------

If you pass ``None`` or an empty string for a name in most ``add_...`` methods, the model will automatically assign a unique name (e.g., ``N#``, ``M#``, ``Q#``, etc.). The helper ``unique_name(dictionary, prefix)`` underlies this behavior for internal generators.

Tips and common patterns
------------------------

- Prefer global load directions (``FX``, ``FY``, ``FZ``) for self-weight and simple vertical/gravity loads; use local directions (``Fx``, ``Fy``, ``Fz``) for member-centric loading.
- Use ``combo_tags`` to select subsets of combinations when calling analyzers.
- For tension-only/comp-only models, use ``analyze(...)`` with multiple ``num_steps`` for better convergence.
- ``sparse=True`` generally provides faster solutions for larger models.

API reference
-------------

.. autoclass:: Pynite.FEModel3D
   :members:
   :undoc-members:
   :show-inheritance:
