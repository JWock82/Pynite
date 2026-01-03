======
Modal
======

Modal analysis computes a structure's natural frequencies and mode shapes. Pynite provides
the ``FEModel3D.analyze_modal()`` method which assembles the global stiffness and mass
matrices, solves the generalized eigenvalue problem, and stores the resulting mode shapes
as displacement results in modal load combinations.

How it works
============

- Assembles the global stiffness matrix ``[K]`` and mass matrix ``[M]`` using the
  specified ``mass_combo_name`` and ``mass_direction``.
- Partitions out supported degrees of freedom before solving to avoid singularities.
- Solves the generalized eigenproblem ``[K]{φ} = λ[M]{φ}`` where ``λ = ω²``. Frequencies
  are returned in Hz as ``f = ω / (2π)``.
- Mode shapes (eigenvectors) are expanded back into the model's full DOF set and
  stored in load combinations named ``Mode 1``, ``Mode 2``, etc.

Usage
=====

``analyze_modal`` signature
---------------------------

``FEModel3D.analyze_modal(num_modes: int = 12, mass_combo_name: str = 'Combo 1', mass_direction: str = 'Y', gravity: float = 1.0, log=False, check_stability=True)``

Parameters
----------

- ``num_modes``: Number of modes to calculate (default: ``12``).
- ``mass_combo_name``: Name of the load combination used to convert loads to masses.
  Defaults to ``'Combo 1'``. If no load combos exist, a default is created automatically.
- ``mass_direction``: Direction to use when converting loads to equivalent mass
  (``'X'``, ``'Y'``, or ``'Z'``). Loads applied in this direction will be treated as
  masses (default: ``'Y'``).
- ``gravity``: Acceleration used when converting loads to mass (default: ``1.0``).
- ``log``: If ``True``, prints progress messages to the console.
- ``check_stability``: If ``True``, checks the stiffness matrix for unsupported DOFs
  and will raise an exception if instabilities are found.

Requirements and notes
======================

- SciPy is required for modal analysis. The method uses SciPy's sparse eigen-solver
  internally; if SciPy is not installed the method will raise an exception.
- The solver filters modes by solving for a limited number of eigenpairs. Increase
  ``num_modes`` if additional modes are required.
- Ensure materials have density (``rho``) and/or a suitable ``mass_combo_name`` with
  applied loads so that a nonzero mass matrix is assembled. An exception will be
  raised if no mass terms are found.

Example
-------

Create a simple model, add mass-producing loads (or define material densities), and
run modal analysis::

   from Pynite.FEModel3D import FEModel3D

   model = FEModel3D()
   # add nodes, members, materials, etc.
   # ensure material densities or a load combo exists to create mass terms
   freqs = model.analyze_modal(num_modes=6, mass_direction='Y', gravity=9.81, log=True)

The returned value is a list (NumPy array) of frequencies in Hz and the mode
displacement results are stored in load combinations named ``Mode 1``, ``Mode 2``, etc.

``` 
