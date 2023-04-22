=======
Members
=======
All members in Pynite are beam-column elements, meaning they can handle axial and transverse loads.
By default all members are also "physical members", meaning they automatically segment themselves
at any internal nodes.

Local Coordinate System
=======================
Each member starts at its i-node and ends at its j-node. The local x-axis for the member is defined
by a vector going from the i-node to the j-node.

By default the local z-axis is always horizontal in Pynite, and is on the right-hand side of the
member.

The local y-axis is defined as the cross-product of the local z-axis with the local x-axis. In
other words, the local y-axis is always perpendicular to the member and to the local z-axis.

End Releases
============
End releases can be applied to each end of a member to simulate pinned connections. End releases
can be applied using the ```FEmodel3D.def_release()``` method. See below for an example. By
applying rotational end releases to both ends of a member you can simulate two-way truss members.

.. code-block:: python

    # The following line turns member M1 into a pin-ended member
    my_model.def_release('M1', Dxi=False, Dyi=False, Dzi=False, Rxi=False, Ryi=True, Rzi=True, Dxj=False, Dyj=False, Dzj=False, Rxj=False, Ryj=True, Rzj=True)

    # This next line does the same thing as the previous line - just simplified
    my_model.def_release('M1', False, False, False, False, True, True, False, False, False, False, True, True)

    # This next line is yet another simple way to do the same thing
    my_model.def_release('M1', 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1)

Note that in the code above, ```Dxi``` stands for displacement in the local x direction at the
i-node, ```Rjz```` stands for rotation about the local z axis at the j-node, and so forth.
