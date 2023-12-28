====================
P-:math:`\Delta` and P-:math:`\delta` Analysis
====================

When a structure is loaded, it displaces. Once a structure has displaced, the forces act on a geometry that is different than the geometry that was initially analyzed. Secondary forces develop due to eccentricities created by these displacements. These secondary forces are known collectively as P-:math:`\Delta` and P-:math:`\delta` effects.

P-:math:`\Delta` effects are secondary forces that develop due to the displacements of nodes in the model, while P-:math:`\delta` effects are secondary forces that develop due to member displacements.

There are traditionally 2 types of procedures used to run P-:math:`\Delta` analysis:
1. The iterative procedure - A basic analysis is rerun considering the effects of the displaced goemetry, leading to further displacements. These further displacements require further iterations until the displacements either converge or diverge.
2. Use of the geometric stiffness matrix `Kg`- A stiffness matrix that adjusts member stiffnesses for the effects of member axial loads is added to the overall stiffness matrix. This method allows for a direct linear solution if the axial load on the members is known in advance.

P-:math:`\delta` effects can be captured by modeling additional nodes along the members length. This forces the analysis to track member internal displacements at the nodes during analysis. In most cases only 2 or 3 intermediate nodes are required. AISC has published an good engineering journal article discussing this method.

Procedure
=========

Pynite can perform perform P-:math:`\Delta` analysis, and with the use of internal nodes added along members P-:math:`\delta` analysis. Please note that in Pynite P-:math:`\Delta` and P-:math:`\delta` effects are not considered for plates.

To run a P-:math:`\Delta` analysis use the ``FEModel3D.anapyze_PDelta()`` method as follows:

.. code-block:: python

    # Run a P-Delta analysis on `my_model`
    my_model.analyze_PDelta()

For a more detailed description of the options available for P-:math:`\Delta` analysis see the :doc:`FEModel3D </FEModel3D>` documentation.

The procedure Pynite uses is as follows:

1. Run a load combination using a simple linear-elastic analysis in order to calculate member axial loads.
2. Perform tension-only and compression-only iterations. If tension/compression-only analysis did not converge, undo the analysis, deactivate/reactivate members as needed, and go back to step 1.
3. Run P-:math:`\Delta` analysis using the axial loads determined from step 1.
4. Perform tension-only and compression-only iterations. If tension/compression-only analysis did not converge, undo the P-:math:`\Delta` analysis, deactivate/reactivate members as needed, and go back to step 3.
5. Go back to step 1 for the next load combination. Repeat until all applicable load combinations have been evaluated.

Validation
==========

AISC has provided benchmark tests found in AISC 360-16 commentary C2.1 that can be used as a guage to test a program's P-:math:`\Delta` and P-:math:`\delta` analysis capabilities. Pynite's P-:math:`\Delta` and P-:math:`\delta` analysis pass the AISC benchmark tests when two intermediate nodes are introduced at the third points of the members to pick up P-:math:`\delta` effects. The number of nodes necessary to model P-:math:`\delta` effects may vary from one model to the next.

Limitations
===========
*Note that P-:math:`\Delta` and P-:math:`\delta` analysis is just one part of an overall second-order analysis. See the building code for additional requirements that may be applicable, such as stiffness reductions and notional loads.

*Loads are applied in a single load step. This is sufficient for most common cases.

*P-:math:`\Delta` effects are not considered for plate elements.
