============
Installation
============

`Pynite` can be installed via `pip` (Python's default package installer). Run the following command from the command line interface:

.. code-block:: console

    $ pip install PyniteFEA[all]

Be sure to install ``PyniteFEA`` rather than ``Pynite``. The second one is a different package that has nothing to do with finite element analysis. The ``[all]`` in the line above is only needed if you want to use `Pynite's` built in visualization and reporting features. For a lighter-weight installation you can omit ``[all]``.

For engineers who may be new to `python`, the ``$`` symbol in the line above represents any directory in your system that ``pip`` is accessible from. If you checked the box to place `python` on your system ``PATH`` during installation (recommended) you can run ``pip`` from any directory (see image below). Otherwise it's sitting in the ``Scripts`` folder in the directory where Python was installed, and you'll need to either add it to your system ``PATH`` or navigate to that directory via the console and run the ``pip`` command above.

.. figure:: ../img/pip_installation_example.png

Adding the [all] command is equivalent to running the following additional installation commands:

.. code-block:: console

    $ pip install matplotlib
    $ pip install scipy
    $ pip install pdfkit
    $ pip install jinja2
    $ pip install vtk
    $ pip install pyvista[all,trame]
    $ pip install ipywidgets
    $ pip install trame_jupyter_extension

Note that if you use [all] without `Jupyter Lab` installed, it will install `Jupyter Lab` for you. If you don't want `Jupyter Lab` installed it's best to install `Pynite` without the [all] command, and individually install the optional dependences listed above that you do want.

With `Jupyter Lab` installed you'll be able to view some of the documents used to help derive `Pynite` located in the `GitHub` repository, if you also install `Sympy`. It can be installed as follows:

.. code-block:: console
    
    $ pip install sympy

Dependencies
============

Below is a detailed description of how `Pynite` uses each of its dependencies to help you decide which (if any) installation customizations you may want to make.

Required Dependencies Installed No Matter What
----------------------------------------------

* `numpy`: Used for matrix algebra and dense matrix solver
* `PrettyTable`: Used to format tabular output
* `scipy`: Used for sparse matrix solver to improve solution speed and memory management. In most cases you'll want to install scipy. It is highly recommended for large models.

Additional Optional Dependencies Included with ``[all]``
--------------------------------------------------------

* `matplotlib`: used for plotting member diagrams.
* `PDFKit`: Used for generating pdf reports. In order to generate pdf reports, PDFKit requires you to have wkhtmltopdf installed on your computer. This is a free program available for download at https://wkhtmltopdf.org/downloads.html. Once installed, you'll need to help Pynite find it. On Windows, this can be done by setting your PATH environment variable to include the path to "wkhtmltopdf.exe" after installation. For example, mine is installed at "C:\Program Files\wkhtmltopdf\bin"
* `jinja2`: Used by `Pynite` for templating reports into HTML prior to HTML-to-pdf conversion. Only needed if you plan to use the PDF reporting features.
* `pyvista[all, trame]`: Used for interactive visualization in `Jupyter` notebooks.
* `ipywidgets`: Used for interactive visualization in `Jupyter` notebooks.
* `trame_jupyter_extension`: Used for interactive visualization in `Jupyter` notebooks.
* `VTK`: (legacy) Used for visualization - Note that VTK is a little picky about which version of Python you are running. You must run a 64 bit installation of Python, rather than a 32 bit version. VTK is published by Kitware. I've noticed Kitware takes a little time updating VTK to be compatible anytime a new version of Python is released. If you're having trouble installing VTK, you can see which versions of Python are supported by visiting https://pypi.org/project/vtk/#files. VTK does not need to be installed if you don't plan to use the visualization tools built into Pynite. Pynite is moving away from `VTK` and toward `pyvista` in an effort to simplify.
* `jupyterlab`: Installed by default with some of the other libraries if not already installed.

Other Dependencies You May Want to consider
-------------------------------------------
* `sympy`: Only needed if you want to view the derivations used to build Pynite.
