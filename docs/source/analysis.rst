========
Analysis
========

`Pynite` offers several analysis options you can choose from. This section will help you identify the anlysis options you want to run.

Sparse vs. Dense Solvers
========================

Each of the analysis options in Pynite allows you to use either a sparse or dense matrix solver. The dense matrix solver stores all the values in the stiffness matrix, whereas the sparse solver only stores non-zero values. Both solvers provide correct solutions, but each will perform differently depending on your model.

The sparse solver is the default solver used by Pynite. It is well suited to large models. It uses less memory and solves large models faster. It solves small models slower, but usually the difference is not noticed because the models are small anyway. In order to use it you'll need to have ``Scipy`` installed.

You can switch to the dense solver by passing ``sparse=False`` to the analysis method you are using. ``Scipy`` does not need to be installed to use the dense solver. The dense solver is well suited to small models. If you need to repeatedly solve a small model, the dense solver can offer some performance advantages over the sparse solver.

General Analysis
================

Use the ``FEModel3D.analyze()`` method to run a general analysis. This analysis is iterative if there are tension-only or compression-only elements or supports in the model.

Linear Analysis
===============

Linear analysis can be performed by using the ``FEModel3D.analyze_linear()`` method. This method of analysis is very fast, but is limited to models without nonlinear features such as tension-only or compression-only elements and supports, or P-:math:`\Delta` effects. This method only needs to assemble the stiffness matrix once because it uses analytical superposition of forces to generate load combinations. Superposition requires a linear model.

P-:math:`\Delta` Analysis
=========================

P-:math:`\Delta` analysis is required by many building codes for frame structures. It covered in great detail here: :doc:`PDelta`

Other Useful Options
====================

You can pass a few other parameters to each analysis. To check statics use `check_statics=True` in the analysis command. To check stability use the `check_stability=True` option. Note that this option does slow down solution speed.
