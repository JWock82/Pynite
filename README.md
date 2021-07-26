# PyNite
[![Build Status](https://travis-ci.com/JWock82/PyNite.svg?branch=master)](https://travis-ci.com/JWock82/PyNite)
![PyPI - Downloads](https://img.shields.io/pypi/dm/PyNiteFEA)
<img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/JWock82/PyNite">
![GitHub last commit](https://img.shields.io/github/last-commit/JWock82/PyNite)
![GitHub](https://img.shields.io/github/license/JWock82/PyNite)

An easy to use elastic 3D structural engineering finite element analysis library for Python.

<div>
<img src="https://github.com/JWock82/PyNite/blob/master/Resources/PyNite%20Bin.PNG" width="300"/>
</div>

# Current Capabilities
* 3D static analysis of elastic structures.
* P-&Delta; analysis of frame type structures.
* Member point loads, linearly varying distributed loads, and nodal loads are supported.
* Classify loads by load case and create load combinations from load cases.
* Produces shear, moment, and deflection results and diagrams for each member.
* Tension-only and compression-only elements.
* Springs: two-way, tension-only, and compression-only.
* Quadrilateral plate elements (based on an isoparametric formulation).
* Rectangular plate elements (based on a 12-term polynomial formulation).
* Basic meshing algorithms for some common shapes.
* Reports support reactions.
* Rendering of model geometry, supports, load cases, load combinations, and deformed shapes.
* Generates PDF reports for models and model results.

# Project Objectives
As I've gotten into the structural engineering profession, I've found there's a need for an easy to use open-source finite element package. I hope to help fill that need by prioritizing the following:

1. Accuracy: There are no guarantees PyNite is error free, but accuracy and correctness are a priority. When bugs or errors are identified, top priority will be given to eliminate them. PyNite's code is frequently reviewed, and its output is regularly being tested against problems with known solutions to isolate errors. If you find an error, please report it as an issue.

2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in outdated languages, or require extensive knowledge of finite element theory and/or element formulations to use. PyNite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package.

4. Improvement: I plan to continue supporting PyNite for many years to come. There are a lot of pieces I'd like to add to PyNite going forward. There's a lot of potential to create extensions as well to solve all kinds of engineering problems. There are more problems to solve than I have time for, so some priorities will have to be made. The plan is to keep PyNite mainstream, adding core functionality first. Occasionally however I may just add what interests me at the time.

5. Collaboration: The intent is to keep PyNite free and open source. This will encourage future development and contributions. Keeping it open source will allow anyone to inspect and improve the code it runs on. If you see an area you think you can help PyNite improve in you are encouraged to contribute. I'd like to get PyNite doing a lot more. Don't be offended if I'm a little slow to accept your contributions. FEA is a very technical subject and accuracy is extremely important to me. Sometimes I'm a little slow understanding both FEA and Python and it takes some time for me to comprehend what's being proposed. I also have a young family to take care of that takes first priority.

# Support
Whether you just need help getting started with PyNite, or are looking to build something more complex, there are a few resources available:
* The examples in the "Examples" folder in this repository cover a variety of simple problems. The comments in the examples provide additional guidance on how PyNite works.
* The wiki in this repository provides basic guidance on how to use PyNite.
* If you're looking for more direct guidance on using PyNite, or for help coding a project, I am available on a private consulting basis. You can reach out to me directly at Building.Code@outlook.com to discuss options.

# Dependencies
PyNite depends on the following packages:
## Required Dependencies
* numpy: used for matrix algebra and dense matrix solver
* scipy: used for sparse matrix solver to improve solution speed
* matplotlib: used for plotting member diagrams
* PrettyTable : used to format tabular output

## Optional Dependencies
* VTK: Used for visualization - Note that VTK is a little picky about which version of Python you are running. You must run a 64 bit installation of Python, rather than a 32 bit version. VTK is published by Kitware. I've noticed Kitware takes a little time updating VTK to be compatible anytime a new version of Python is released. If you're having trouble installing VTK, you can see which versions of Python are supported by visiting https://pypi.org/project/vtk/#files. VTK does not need to be installed if you don't plan to use the visualization tools built into PyNite.
* PDFKit: Used for generating pdf reports. In order to generate pdf reports, PDFKit requires you to have wkhtmltopdf installed on your computer. This is a free program available for download at https://wkhtmltopdf.org/downloads.html. Once installed, you'll need to help PyNite find it. On Windows, this can be done by setting your PATH environment variable to include the path to "wkhtmltopdf.exe" after installation. For example, mine is installed at "C:\Program Files\wkhtmltopdf\bin"
* jinja2: Used for templating reports into HTML prior to HTML-to-pdf conversion.
* jupyterlab: Only needed if you want to view the derivations used to build PyNite.
* sympy: Only needed if you want to view the derivations used to build PyNite.

# Example Projects
Here's a list of projects that run on PyNite:

* Building Code (https://building-code.herokuapp.com/) - This one is my personal side project.
* Standard Solver (https://www.standardsolver.com/)

# What's New?
Version 0.0.40
* Added support springs.
* Reorganized visualization code.

Version 0.0.39
* Made console output optional. By default it is now turned off.
* Fixed a warning that occured during P-Delta analysis. Sparse matrix solution efficiency was poor.
* Fixed deformed shape rendering for plates and quads. The scale factor was not being applied.
* Added a basic screenshot feature for rendering. If `screenshot` is set to a filepath in the `RenderModel` method a .PNG screenshot will be saved to the filepath and the render window will be closed automatically.

Version 0.0.38
* Bug fix for load vector calculation. When multiple load cases were present on a plate or quad, only the last entered load case was being considered.
* Additional fix for an exception that occured when unused load cases were in a model that was being rendered.

Version 0.0.37
* Fixed rectangular mesh control point feature. It was not working.
* Added functions to simplify getting localized max/min moment and shear results from meshes.
* Fixed an exception that occured when unused load cases were in a model that was being rendered.

Version 0.0.36
* Correction to sign convention to rectangular plate loads. They were being applied in the opposite direction than was specified by the user. Quadrilaterals were not affected.
* Changed sign convention on rectangular plate to match the sign convention for quadrilaterals. These two elements are derived using different bending sign conventions and it seemed to be appropriate to make them behave the same way for the end user.

Version 0.0.35
* Issues with rectangular plate elements have been fixed. Membrane stiffness terms were being placed in the wrong location in the element's global stifness matrix. This bug was identical to the one that had been occuring in quadrilateral elements prior to v0.0.28.
* Bug fix for `dz` countours for rectangular plate elements. The individual plate corner displacements were being mapped to the wrong corners. Only the `dz` contours were affected by this issue.

Version 0.0.34
* Bug fix for contour smoothing on rectangular plates. The plate results were correct, but the contours were slightly off after adding the contour smoothing feature.
* Added rectangular meshes. More to come on this feature.

Versions 0.0.32 and 0.0.33
* Started work in integrating Travis-CI into GitHub for testing PyNite code.

