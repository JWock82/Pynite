============
Fundamentals
============

The ``FEModel3D`` class does most of the work for you in Pynite. Below is a detailed reference of
all the members available within the class. To begin modeling create a new instance of the
``FEModel3D`` class as follows:

.. code-block:: python

    from PyNite import FEModel3D
    my_model = FEModel3D()

That's it! now you have a finite element model named ``my_model``. Usually you'll want to give it
a name that makes sense to you rather than ``my_model``. For the purposes of this documentation
anytime you see ``my_model`` remember that we're just referring to a ``FEModel3D`` object that we
have created.

The ``FEModel3D`` object has many "methods" built into it to help you create and work with your
model. For example, to add a node to the model we can simply use the ``add_node`` method:

.. code-block:: python

    my_model.add_node('N1', 0, 0, 0)

The methods available through the ``FEModel3D`` object are documented below.

As you create your model, the ``FEModel3D`` object organizes and stores the information you give it.
This information is stored primarily in Python `dictionaries`. You can access these dictionaries and
retrieve information from them if you recall the name of the object stored in the dictionary. For
example, to retrieve the X-coordinate for the node we just created we could access it in the
``Nodes`` dictionary as follows:

.. code-block:: python

    my_model.Nodes['N1'].X

This returns the X coordinate for the node. Notice the "dot" operator above. As you use the
dictionaries a good code environment can prompt you to see what information is available when you
use the "dot" operator.

Avaliable dictionaries are:

* ``my_model.Nodes``       A dictionary containing all the nodes in the model
* ``my_model.AuxNodes``    A dictionary containing all the auxiliary nodes in the model
* ``my_model.Members``     A dictionary containing all the members in the model
* ``my_model.Springs``     A dictionary containing all the springs in the model
* ``my_model.Plates``      A dictionary containing all the rectangular plates in the model
* ``my_model.Quads``       A dictionary containing all the quadrilaterals in the model
* ``my_model.Materials``   A dictionary containing all the materials in the model
* ``my_model.LoadCombos``  A dictionary containing all the load combinations in the model
* ``my_model.Meshes``      A dictionary containing all the meshes in the model

FEModel3D Class Reference
=========================

.. autoclass:: PyNite.FEModel3D
   :members:
