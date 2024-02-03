====================
Creating PDF Reports
====================

Requirements
============

In order to use Pynite's reporting features you'll need to install the following packages. These requirements only need to be installed if you plan on using the reporting features:

* PDFKit - ``PDFKit`` converts HTML files to PDF files via ``wkhtmltopdf`` (see discussion on ``wkhtmltopdf`` below).
* Jinja2 - ``Jinja2`` is a templating engine used for creating HTML files using Python.

You'll also need to install ``wkhtmltopdf``, which is a program, rather than a Python library. Installers for various operating systems can be freely downloaded here: `download wkhtmltopdf <https://wkhtmltopdf.org/downloads.html>`_.

Configuring ``wkhtmltopdf`` on Windows
======================================

Once ``wkhtmltopdf`` is installed, you'll need to set your environment PATH variable to include the path to ``wkhtmltopdf.exe``. This step ensures PDFKit can find and use ``wkhtmltopdf``. On Windows this can generally be done by going to ``Control Panel -> System & Security -> System -> Advanced System Settings -> Environment Variables``. On Windows 10 you can also just type "env" in the search bar to bring it up.

From there you can edit the "Path" variable to include the file path to the folder containing ``wkhtmltopdf.exe``. The path you need to add will most likely be "C:\\Program Files\\wkhtmltopdf\\bin" or something similar, depending on where you installed it.

Configuring ``wkhtmltopdf`` on Other Operating Systems
======================================================

I've only ever done this for Windows, so for other operating systems I recommend going to the PDFKit documentation on PyPI here: `PDFKit on PyPI <https://pypi.org/project/pdfkit/>`_. It explains PKFKit configuration options that "should" work, thought I've never done it myself. Good luck!

Generating & Customizing Reports
================================

Generate reports by importing the ``Reporting`` module and then calling the ``create_report`` method. By default, the report prints everything, but it can be customized by setting the various keyword arguments in the ``create_report`` method to ``False``. Also by default, the PDF report is saved in the ``Pynite`` folder unless another filepath is specified using the ``output_filepath`` argument.

Syntax:

.. code-block:: python

    create_report(model, output_filepath='.//PyNite Report.pdf', nodes=True,
                  members=True, plates=True, member_releases=True,
                  node_reactions=True, node_displacements=True,
                  member_end_forces=True, member_internal_forces=True,
                  plate_corner_forces=True, plate_center_forces=True,
                  plate_corner_membrane=True, plate_center_membrane=True)

Example:

.. code-block:: python

    # Import the reporting module
    from PyNite import Reporting
    
    # Create the report
    Reporting.create_report(my_model, output_filepath='./My Report.pdf')

Reporting Class Reference
=========================
.. automodule:: PyNite.Reporting
    :members: create_report