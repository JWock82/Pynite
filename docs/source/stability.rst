=========
Stability
=========

Perhaps the trickiest part of finite element modeling is making sure a model is stable. That may
sound trivial, but in finite element analysis it is nuanced. If you are getting "singularity"
errors in your model, it may not be as stable as you think it is. Remember that just because a
model is static does not mean it is stable.

Degrees of Freedom
==================
In 3 dimensions there are 6 "degrees of freedom" for motion to occur in. An object can translate
in either the X, Y or Z direction, and it can rotate about the X, Y or Z axis. These 3 translations
and 3 rotations together sum to 6 degrees of freedom. For an object to be stable, all 6 degrees
of freedom need to be stabilized

Instability Type 1: Rigid Body Motion
=====================================
Rigid body motion occurs when the entire structure can move as a unit in one of the 6 degrees of
freedom. This occurs when there is insufficient support for the structure. When rigid body motion
occurs, the stiffness matrix is "singular", and unsolvable. Adding supports can eliminate this type
of singularity.

Instability Type 2: Nodal Instability
=====================================
Once the global structure is stabilized, there may still be localized instabilities within it. Just
as all 6 degrees of freedom need to be constrained at the global level, so do instabilities at the
nodal (local) level need to be constrained. These types of instabilities also show up as
"singularities" in the stiffness matrix.

Consider a truss with multiple hinged-ended members framing into a node. What's to stop the node
from spinning if all the members it's attached to are hinged? The correct way to model this joint
would be to leave one of the member ends at the node unhinged. Rotations at the node would be then
be restrained by that member. It may be that there are no rotations at the node, but the degree of freedom
still needs to be tethered to a member (or a support) to make the node stable.

When selecting which member to tether the rotation to, consider what the connection looks like in
real life. If you were to apply a moment to the joint, where would it go? Unless you're dealing with
true pins, it will go somewhere. Even true pins are not frictionless (I've never seen them
perpetually rotating in the real world), and you may as well pick which member you want to assign
stray forces at that degree of freedom to.

Pynite can help you isolate these types of instabilities by passing ``check_stability=True`` to
any of the analysis methods (see below for an example). Note that this does slow down analysis, so
you usually only want to do this when you are debugging a model.

.. code-block:: python

    my_model.analyze(check_stability=True)

Instability Type 3: Second Order Effects
========================================
Once the Type 1 and Type 2 instabilities are rooted out of your model, there may be a third type of instability lurking quietly inside it. These types of instabilities can go unnoticed without careful attention to the analysis procedure.

Most material codes recognize that a first-order analysis is often not enough. As a structure deforms under load, the geometry changes, making the analysis results invalid. This requires us to reanalyze the structure under the new deformed geometry.

The P-:math:`\Delta` analysis feature is one tool PyNite provides to help check for this type of instability. P-:math:`\Delta` analysis is but one piece of the overall picture, which is usually defined by the material's building code. Additional considerations are usually necessary to capture this behavior, such as P-:math:`\delta` analysis, material stiffness reductions and notional loads. Note also that Pynite only performs P-:math:`\Delta` analysis for members. P-:math:`\Delta` effects are not considered for plates. See :doc:`PDelta` for more information.

Second order effects are usually only critical for models with slender members, or members with high axial loads in the presence of bending moments. If you don't have these types of features in your model, second order effects may be negligible. Consult your building code for specific criteria in making this determination. Either way, running a P-Delta analysis can help identify erorrs in the model you may have been unaware of. A model that cannot pass a P-Delta analysis is often something you don't really want to build.
