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

# Current Capabilities
* 3D static analysis of elastic structures.
* P-&Delta; analysis of frame type structures.
* Member point loads, linearly varying distributed loads, and nodal loads are supported.
* Classify loads by load case and create load combinations from load cases.
* Produces shear, moment, and deflection results and diagrams for each member.
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

2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in outdated languages, or require extensive knowledge of finite element theory and/or element formulations to use. PyNite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package.

4. Improvement: I plan to continue supporting PyNite for many years to come. There are a lot of pieces I'd like to add to PyNite going forward. There's a lot of potential to create extensions as well to solve all kinds of engineering problems. There are more problems to solve than I have time for, so some priorities will have to be made. The plan is to keep PyNite mainstream, adding core functionality first. Occasionally however I may just add what interests me at the time.

5. Collaboration: The intent is to keep PyNite free and open source. This will encourage future development and contributions. Keeping it open source will allow anyone to inspect and improve the code it runs on. If you see an area you think you can help PyNite improve in you are encouraged to contribute. I'd like to get PyNite doing a lot more. Don't be offended if I'm a little slow to accept your contributions. FEA is a very technical subject and accuracy is extremely important to me. Sometimes I'm a little slow understanding both FEA and Python and it takes some time for me to comprehend what's being proposed. I also have a young family to take care of that takes first priority.

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
v0.0.63 thru v0.0.64
* Fixed the `add_mesh` method. It was not working properly after version 0.0.62.
* Made stability checks optional. Stability checks add significant solve time. If you are confident your model is stable, you can skip the stability check by toggling `check_stability` to `False` in your call to your analysis command.

v0.0.62
* PyNite now checks for nodal instabilities when analyzing a model. If nodal instabilities are found, PyNite will output the unstable nodes and directions to the console, and will throw an exception.
* Added a method called `rename` to the `FEModel3D` class for quickly renaming all the nodes and elements in the model in sequential order.
* Added a `last_node` and `last_element` attribute to the `mesh` class. These methods can be used to get the name of the last node or element in a mesh.
* Improved the reliability of the `add_mesh` method. It now can handle adding meshes containing node and element names already defined in the model. It automatically resolves the duplicate names.

v0.0.57 thru 0.0.61
* Fixed a stubborn bug that wouldn't create openings if they didn't have a node inside them. This prevented openings from showing up in some meshes.
* Added unit testing to help identify similar problems with mesh openings in the future.

v0.0.56
* Bug fix for orthotropic rectangular plates. The stiffness was slightly off on rectangular plates in version 0.0.55. Prior versions were not affected.
* Added unit testing for in-plane (membrane) behavior, and stiffness modification factors. PyNite solutions matched theoretical solutions within 0.1% for uncracked sections, and within 2% for cracked sections.

v0.0.55
* Added stiffness modification factors for rectangular plate and quarilateral elements. Orthotropic in-plane behavior can now be modeled. This can be used to model the cracked stiffness of concrete and masonry for in-plane loads. Exercise caution when using this feature. These factors only apply to in-plane stiffnesses in the element's local x and y directions. Out-of-plane stiffnesses are not modified. Ensure element local axes are aligned to the directions you want to apply the stiffness modifications to.
* The arguments for many methods relating to rect plates and quads have been reorganized with the addition of the stiffness modification factors. With this update any plate/quad models you have will likely need to be updated.
* Reformulated rectangular plate element in-plane (membrane) forces. These elements are now isoparametric just like the quad elements. This was done to make them orthotropic.
* Improved docstrings for plates and quads.

v0.0.54
* Bug fix for member point loads applied in a global direction. They were overiding the member end forces due to an ambiguous variable name.
* Bug fix for rendering. Several variable names had underscores that didn't belong there. It was causing rendering errors. See version 0.0.53 notes for a full list of what's been changing with regard to rendering.

v0.0.53
* `scipy` is now optional for the dense matrix solver. By default `PyNite` uses the sparse matrix solver. If you wish to use the dense matrix solver, you'll need to set the `sparse` parameter equal to `False` when you analyze a model. Generally, the sparse matrix solver is faster and uses less memory. This feature was added to improve compatibility with other programs that may not work well with `scipy`.
* Started work on a `Renderer` class for visualization to better organize the code.
* Renamed several variables in the `render_model` function to make them more descriptive of what they control. If you use keyword arguments in your calls to `render_model` you may need to change some of the variable names.

v0.0.52
* Major bug fix for quadrilaterals. In-plane stiffnesses were 1/4 of what they should have been.
* Added more unit testing for plates and quadrilaterals.
* Fixed a bug where cylindrical quad meshes didn't calculate the number of quads that fit in the circumference correctly. This bug only affected models where the number was not explicitly specified.
* Added plate/quad surface load validation to ensure loads are not added to ficticious plates/quads.
* Improved the `remove_duplicate_nodes` method. It wasn't working properly for springs. It also now returns a list of the names of the nodes that were removed. Also added a docstring for this method.
* Added the `orphaned_nodes` method that checks the model for orphaned nodes and returns a list of the names of orphaned nodes.

v0.0.51
* Internal changes to some matrix operations. This was in response to Issue #102 where statics were not checking out for some 3D plate models.
* Reduced stiffness of plate/quad element drilling degree of freedom. Again in response to Issue #102.

v0.0.50
* Bug fix for cylindrical meshes.
* Added tension/compression only support springs. This required some reworking of how support springs are implemented. See the "Beam on Elastic Foundation" example for an example of how to properly use support springs.
* More PEP8 style changes. This shouldn't change the way you use the program. Some of the internal calls to function names still referenced the old function names. It was leading to deprecation warnings during runtime.
* Fixed a divide by zero error for plotting contours on models with nodes not attached to plates/quads.

v0.0.49
* Major speed boost. A lot of time was being wasted reading/writing to a SciPy `lil_matrix` format. The stiffness matrix is now assembled as a `coo_format` first, and then switched to lil format later on to facilitate slicing.
* Refactored node names to match PEP8 style guide. Instead of `iNode`, `jNode`, etc. PyNite now uses `i_node`, `j_node` etc. This should only affect you if you were using the keyword arguments for node names.
* Removed test modules from package distribution. These files were unnecessary for users to have installed on their computers.

v0.0.48 - Global Load Directions
* Global load directions can now be used for member loads now. Use capital notation to apply the member load in the global X, Y, or Z direction (e.g. 'FX', 'MY', 'FZ'...). Use lower case notation to apply the member load in the member's local x, y, or z direction (e.g. 'Fx', 'My', 'Fz'...).
* Bug fix for rendering models without plates/quads (special thanks to tamalone1 for this pull request!)

# Example Projects
Here's a list of projects that run on PyNite:

* Building Code (https://building-code.herokuapp.com/) - This one is my personal side project.
* Standard Solver (https://www.standardsolver.com/)
