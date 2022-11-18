<div align="center">
  <img src="https://github.com/JWock82/PyNite/raw/master/Resources/website_logo_solid_background.png" width=40% align="center"/>
  <br>
  <h1>Simple Finite Element Analysis in Python</h1>
</div>

![Build Status](https://github.com/JWock82/PyNite/actions/workflows/build-and-test.yml/badge.svg)
![PyPI - Downloads](https://img.shields.io/pypi/dm/PyNiteFEA)
<img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/JWock82/PyNite">
![GitHub last commit](https://img.shields.io/github/last-commit/JWock82/PyNite)
![GitHub](https://img.shields.io/github/license/JWock82/PyNite)

An easy to use elastic 3D structural engineering finite element analysis library for Python.

# Installation
The easiest way to install Pynite is with pip: `pip install PyNiteFEA`

# Current Capabilities
* 3D static analysis of elastic structures.
* P-&Delta; analysis of frame type structures.
* Member point loads, linearly varying distributed loads, and nodal loads are supported.
* Classify loads by load case and create load combinations from load cases.
* Produces shear, moment, and deflection results and diagrams for each member.
* Automatic handling of internal nodes on frame members (physical members).
* Tension-only and compression-only elements.
* Spring elements: two-way, tension-only, and compression-only.
* Spring supports: two-way and one-way.
* Quadrilateral plate elements (MITC4 formulation).
* Rectangular plate elements (12-term polynomial formulation).
* Basic meshing algorithms for some common shapes and for openings in rectangular walls.
* Reports support reactions.
* Rendering of model geometry, supports, load cases, load combinations, and deformed shapes.
* Generates PDF reports for models and model results.

# Project Objectives
As I've gotten into the structural engineering profession, I've found there's a need for an easy to use open-source finite element package. I hope to help fill that need by prioritizing the following:

1. Accuracy: There are no guarantees PyNite is error free, but accuracy and correctness are a priority. When bugs or errors are identified, top priority will be given to eliminate them. PyNite's code is frequently reviewed, and its output is tested against a suite of textbook problems with known solutions using continuous integration (CI) anytime a change to the code base is made. If you do happen to find an error, please report it as an issue.

2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in difficult languages, or require extensive knowledge of finite element theory and/or element formulations to use. PyNite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package.

4. Improvement: PyNite is getting better at what it does. Each new feature provides leverage to build upon previous features in more elaborate ways. Key to improvement is (a) maintaining the core features that other features rely on, and (b) adding new features that are solid stepping stones to other features. Improvement most often happens by getting the small and simple things right incrementally, rather than with sweeping overhauls all at once.

5. Collaboration: The intent is to keep PyNite free and open source. This will encourage future development and contributions. Keeping it open source will allow anyone to inspect and improve the code it runs on. If you see an area you can help PyNite improve in you are encouraged to contribute. Please follow the contributing guidelines given in this repository.

# Support
Whether you just need help getting started with PyNite, or are looking to build something more complex, there are a few resources available:
* The examples in the "Examples" folder in this repository cover a variety of simple problems. The comments in the examples provide additional guidance on how PyNite works.
* The wiki in this repository provides basic guidance on how to use PyNite, although it has fallen behind the development of the program.
* If you're looking for more direct guidance on using PyNite, or for help coding a project, I am available on a private consulting basis. You can reach out to me directly at Building.Code@outlook.com to discuss options.

# Dependencies
PyNite depends on the following packages:
## Required Dependencies
* numpy: used for matrix algebra and dense matrix solver
* matplotlib: used for plotting member diagrams
* PrettyTable : used to format tabular output

## Optional Dependencies
* scipy: Used for sparse matrix solver to improve solution speed and memory management. In most cases you'll want to install scipy. 
* VTK: Used for visualization - Note that VTK is a little picky about which version of Python you are running. You must run a 64 bit installation of Python, rather than a 32 bit version. VTK is published by Kitware. I've noticed Kitware takes a little time updating VTK to be compatible anytime a new version of Python is released. If you're having trouble installing VTK, you can see which versions of Python are supported by visiting https://pypi.org/project/vtk/#files. VTK does not need to be installed if you don't plan to use the visualization tools built into PyNite.
* PDFKit: Used for generating pdf reports. In order to generate pdf reports, PDFKit requires you to have wkhtmltopdf installed on your computer. This is a free program available for download at https://wkhtmltopdf.org/downloads.html. Once installed, you'll need to help PyNite find it. On Windows, this can be done by setting your PATH environment variable to include the path to "wkhtmltopdf.exe" after installation. For example, mine is installed at "C:\Program Files\wkhtmltopdf\bin"
* IPython: Used for displaying screenshots from VTK.
* jinja2: Used for templating reports into HTML prior to HTML-to-pdf conversion.
* jupyterlab: Only needed if you want to view the derivations used to build PyNite.
* sympy: Only needed if you want to view the derivations used to build PyNite.

# What's New?

v0.0.69
* Bug fix for rendering screenshots. The ability to interact with the render window was being disabled after the first screenshot had been taken, forcing subsequent screenshots to use the same view as the first one.

v0.0.68
* Bug fix for member distributed loads on physical members. Added a unit test to check for this error going forward. This bug only affected physical members (new as of v0.0.67) that had distributed loads and internal nodes.

v0.0.67
* Added physical members. Members now automatically detect internal nodes and subdivide themselves and their loads.
* Refactoring: deprecated old method names for member results. You may now have some errors show up if you still try to get member results using the old method names.
* Bug fix for P-Delta analysis. Global displacements were correct, but member internal forces were neglecting the geometric stiffness matrix. The impact of this bug was minimal, since the strain induced by correct global displacements was still being considered prior to this update. You should see a slight change to member P-Delta results.
* Code simplification for P-Delta analysis.

v0.0.66
* Code simplification and bug fix for merging duplicate nodes.
* When nodes are merged, support conditions for the deleted node are now assigned to the remaining node.
* Added a linear solver for faster analysis of simple models. If you don't need P-Delta analysis or tension/compression-only analysis this solver saves time by only assembling the global stiffness matrix once.

v0.0.65
* Improved the `merge_duplicate_nodes` method. It seemed to be working, but it was hard to follow, and there may have been cases where it didn't work as expected. Simplified the code for this method to make it clear what it was doing, and to make it more efficient. Added comments explaining each step.
* Screenshot size is now adjustable when rendering.
* Fixed a bug for `RectangleMesh` where it could not be used repeatedly.
* Refactoring: changed `Name` to `name` throughout code. For example, `Node3D.Name` is now `Node3D.name`.
* Fixed obsolete method names that had not been updated.
* Scalar bar text size can now be controlled. It had strange behavior before. It would change with the window size (until the window size was too small).
* More work on the new `Renderer` class. This class is being built to give the user more control over the appearance and behavior of renderings.
* Bug fix for nodal springs applied in the 'RY' and 'RZ' direction. Exceptions were being thrown in some cases.

v0.0.63 thru v0.0.64
* Fixed the `add_mesh` method. It was not working properly after version 0.0.62.
* Made stability checks optional. Stability checks add significant solve time. If you are confident your model is stable, you can skip the stability check by toggling `check_stability` to `False` in your call to your analysis command.

v0.0.62
* PyNite now checks for nodal instabilities when analyzing a model. If nodal instabilities are found, PyNite will output the unstable nodes and directions to the console, and will throw an exception.
* Added a method called `rename` to the `FEModel3D` class for quickly renaming all the nodes and elements in the model in sequential order.
* Added a `last_node` and `last_element` attribute to the `mesh` class. These methods can be used to get the name of the last node or element in a mesh.
* Improved the reliability of the `add_mesh` method. It now can handle adding meshes containing node and element names already defined in the model. It automatically resolves the duplicate names.

# Example Projects
Here's a list of projects that run on PyNite:

* Building Code (https://building-code.herokuapp.com/) - This one is my personal side project.
* Standard Solver (https://www.standardsolver.com/)
