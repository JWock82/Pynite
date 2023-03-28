.. Pynite documentation master file, created by
   sphinx-quickstart on Tue Mar 21 19:42:11 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pynite's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. toctree::
   quickstart

Pynite provides simple finite element analysis for Python.

This documentation is just getting started. Keep checking back for updates.

Installation
============

The easiest way to install pynite is via `pip`:
::
   
   $pip install PyNiteFEA

Be sure to install `PyNiteFEA` rather than `PyNite`.

Modeling Basics
===============

You can start a new model using the `FEModel3D`` class:
::

   from PyNite import FEModel3D
   my_model = FEModel3D()

Once you have a model started you can add elements to it using the dot operator:
::

   my_model.add_node('N1', 0, 0, 0)

To use PyNite effectively you should become familiar with the methods available within the
`FEModel3D` class:

.. automodule:: PyNite.FEModel3D
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
