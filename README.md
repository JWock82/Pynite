# PyNite
[![Build Status](https://travis-ci.com/JWock82/PyNite.svg?branch=master)](https://travis-ci.com/JWock82/PyNite)
![PyPI - Downloads](https://img.shields.io/pypi/dm/PyNiteFEA)
<img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/JWock82/PyNite">
![GitHub](https://img.shields.io/github/license/JWock82/PyNite)

An easy to use linear elastic 3D structural engineering finite element analysis library for Python.

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

4. Improvement: I plan to continue supporting PyNite for many years to come. There are a lot of pieces I'd like to add to PyNite going forward: triangular plates, plate meshing and loading algorithms, dynamics, pushover anlysis, etc. There's a lot of potential to create extensions as well to solve all kinds of engineering problems. There are more problems to solve than I have time for, so some priorities will have to be made. The plan is to keep PyNite mainstream, adding core functionality first. Occasionally however I may just add what interests me at the time.

5. Collaboration: The intent is to keep PyNite free and open source. This will encourage future development and contributions. Keeping it open source will allow anyone to inspect and improve the code it runs on. If you see an area you think you can help PyNite improve in you are encouraged to contribute. I'd like to get PyNite doing a lot more. Don't be offended if I'm a little slow to accept your contributions. FEA is a very technical subject and accuracy is extremely important to me. Sometimes I'm a little slow understanding both FEA and Python and it takes some time for me to comprehend what's being proposed. I also have a young family to take care of that takes first priority.

# Dependencies
PyNite depends on the following packages:
## Required Dependencies
* numpy: used for matrix algebra and dense matrix solver
* scipy: used for sparse matrix solver to improve solution speed
* matplotlib: used for plotting member diagrams
* PrettyTable : used to format tabular output

## Optional Dependencies
* VTK: used for visualization - Note that VTK is a little picky about which version of Python you are running. You must run a 64 bit installation of Python, rather than a 32 bit version. VTK does not need to be installed if you don't plan to use the visualization tools built into PyNite.
* PDFKit: Used for generating pdf reports. In order to generate pdf reports, PDFKit requires you to have wkhtmltopdf installed on your computer. This is a free program available for download at https://wkhtmltopdf.org/downloads.html. Once installed, you'll need to help PyNite find it. On Windows, this can be done by setting your PATH environment variable to include the path to "wkhtmltopdf.exe" after installation. For example, mine is installed at "C:\Program Files\wkhtmltopdf\bin"
* jinja2: Used for templating reports into HTML prior to HTML-to-pdf conversion.
* jupyterlab: Only needed if you want to view the derivations used to build PyNite.
* sympy: Only needed if you want to view the derivations used to build PyNite.

# Example Projects
Here's a list of projects that run on PyNite:

* Building Code (https://building-code.herokuapp.com/) - This one is my personal side project.
* Standard Solver (https://www.standardsolver.com/)

# What's New?
Versions 0.0.32 and 0.0.33
* Started work in integrating Travis-CI into GitHub for testing PyNite code.

Version 0.0.31
* Added contour smoothing for better plate and quad contour plots.
* Greatly reduced memory usage. The stiffness matrix is now stored as a sparse matrix by default.
  If the user opts not to use the sparse solver it is converted to a dense matrix later on.
* Revised circular meshes so that the Y-axis is upward instead of the Z-axis.
* Bug fix - Rendering springs was broken after the recent change to storing springs in dictionaries.

Version 0.0.30
* Fixed a bug where quad and plate displacement contours weren't being plotted correctly when user-defined load combinations were being used.
* Added a few basic meshing features for quads.
* Added the ability to merge duplicate nodes.

Version 0.0.29
* Fixed a bug where load combinations were being ignored on plates and quadrilaterals. The default load combination was being used for these items instead of user defined load combinations.
* Fixed a bug where auxiliary nodes were causing renderings to crash in some cases.
* Fixed a bug where area loads were always pointing in the positive direction during rendering. The loads were being applied properly in the model even though they were showing up incorrectly in the rendering.
* Quad membrane stress contours can now be rendered.
* A few meshing algorithms have been added to simplify building models from quadrilaterals. Rings, annulusses and frustrums can now be meshed automatically. This feature will continue to be expanded and examples will follow.
* Added an option to turn off rendering of labels in the `RenderModel` method. This can greatly speed up interaction on large plate/quad models.

Version 0.0.28
* Issues with quadrilateral elements have been fixed. Membrane stiffness terms were being placed in the wrong location in the element's global stifness matrix.
* Nodes, members, plates, quads, springs, and auxiliary nodes are now stored in dictionaries for faster computing and easier user access. For example, instead of using the syntax `FEModel3D.GetNode('node_name')` to retrieve a node, you can now alternatively access a node directly from the `Nodes` dictionary using the syntax `FEModel3D.Nodes['node_name']`.

