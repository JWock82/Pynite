=================
Pushover Analysis
=================

Pushover analysis is a powerful way to consider nonlinear behavior in a structure. A pushover load pattern is incrementally added to each load combination analyzed. P-Delta effects and tension/compression-only anlysis is included.

Here's how it works, step by step, in Pynite:

1. Run a load combination using a simple linear-elastic analysis.
2. Perform P-Delta analysis and tension/compression-only analysis with the load we've just applied to capture geometric nonlinear effects in the load combination.
3. Apply a user-defined fraction of the pushover load and adjust the stiffness matrix for inelastic effects at member ends. Adjust the total calculated displacements accordingly.
4. Repeat steps 2 and 3 until the full pushover load has been applied.
5. Repeat steps 2 and 3, but without any additional pushover load, until the displacements converge (or diverge).
5. Go back to step 1 for the next load combination. Repeat until all applicable load combinations have been evaluated for the pushover load.

Limitations
===========
*Nonlinear member behavior is limited to locations where there are nodes in your model. This is usually sufficient for most cases.

*Nonlinear behavior of plates and quads is not considered in the analysis.

*Only the pushover load pattern is applied incrementally. Other loads defined by your load combinations are first run assuming elastic behavior. The pushover load is then added incrementally in subsequent analysis iterations. What this means is that if you have non-linear behavior without applying the pushover loads, your analysis may be incorrect. This is done to ensure that only the pushover load is being incrementally applied. For example, gravity loads should usually not vary as a seismic load is incrementally applied.

Defining Members for Use with Pushover Analysis
===============================================
Capturing nonlinear behavior requires information regarding the section properties of members. Use of Pynite's `Sections` feature allows Pynite to calculate plastic section properties. Pynite will assume that every member can develop full cross-sectional strength at every node. Be cautioned that forces developing at member ends may require connections that are as strong as the member.

Determining an Appropriate Load Step
====================================
The load step used for nonlinear analysis is critical. Too large of a load step may lead to multiple plastic mechanisms forming at once, and/or innacurate estimations of plastic stiffness reductions. Too small of a load step will lead to excessive solve time. Generally applying the load 10-20% at a time will be sufficient, but the user is advised to verify that a smaller load step does not lead to significant changes in results.

Creating a Pushover Load Combination & Viewing Results
======================================================

Pushover analysis will run through your load combinations like any of Pynite's other solvers, except that it will incrementally add a static pushover load. In addition to the load combinations you would normally set up, you'll need to add a pushover load combination that defines the final pushover load, as well as the load step used to get there.

Since the pushover load combination is added to each of the other combinations, results will not be available for the pushover combination itself. Rather, pushover effects will show up in the results for the other load combinations.

.. code-block:: python

    # Add a pushover load combination
    my_model.add_load_combo('Push Combo' {'E': 0.1})

In the example above, any loads in the earthquake load case 'E' will be included in the pushover combination 'Push Combo'. The use of the 0.1 factor on 'E' indicates that the load will be increased 10% for each iteration until a final load factor of 1.0 is applied. The use of 'E' and 'Push Combo' are arbitrary. You can use any names you like for load cases and load combinations defining a pushover analysis.

With a pushover load combination defined, we are ready to run a pushover analysis. Below is an example using the 'Push Combo' we defined above.

.. code-block:: python

    # Run the pushover anlaysis
    my_model.analyze_pushover(log=False, check_stability=True, push_load='Push Combo', max_iter=30, tol=0.01, sparse=True, combo_tags=None)
