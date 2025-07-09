<div align="center">
  <img src="https://github.com/JWock82/Pynite/raw/main/Resources/Full Logo No Buffer.png" width=40% align="center"/>
  <br>
  <h1>Simple Finite Element Analysis in Python</h1>
</div>

![Build Status](https://github.com/JWock82/Pynite/actions/workflows/build-and-test.yml/badge.svg)
[![codecov](https://codecov.io/gh/JWock82/Pynite/branch/main/graph/badge.svg?token=ZH18US3A7P)](https://codecov.io/gh/JWock82/Pynite)
![PyPI - Downloads](https://img.shields.io/pypi/dm/PyNiteFEA)
<img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/JWock82/Pynite">
![GitHub last commit](https://img.shields.io/github/last-commit/JWock82/Pynite)
![GitHub](https://img.shields.io/github/license/JWock82/Pynite)
[![Documentation Status](https://readthedocs.org/projects/pynite/badge/?version=latest)](https://pynite.readthedocs.io/en/latest/?badge=latest)

An easy to use elastic 3D structural engineering finite element analysis library for Python.

# Installation
The easiest way to install Pynite is with pip: `pip install PyniteFEA[all]`.

For a more detailed discussion on installation options and dependencies see https://pynite.readthedocs.io/en/latest/installation.html

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
* Quadrilateral plate elements (DKMQ formulation).
* Rectangular plate elements (12-term polynomial formulation).
* Basic meshing algorithms for some common shapes and for openings in rectangular meshes.
* Reports support reactions.
* Rendering of model geometry, supports, load cases, load combinations, and deformed shapes.
* Advanced tools for modeling and analyzing complex shear walls.
* Generates PDF reports for models and model results.

# Project Objectives
As I've gotten into the structural engineering profession, I've found there's a need for an easy to use open-source finite element package. I hope to help fill that need by prioritizing the following:

1. Accuracy: There are no guarantees Pynite is error free, but accuracy and correctness are a priority. When bugs or errors are identified, top priority will be given to eliminate them. Pynite's code is frequently reviewed, and its output is tested against a suite of textbook problems with known solutions using continuous integration (CI) anytime a change to the code base is made. If you do happen to find an error, please report it as an issue.

2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in difficult languages, or require extensive knowledge of finite element theory and/or element formulations to use. Pynite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package.

4. Improvement: Pynite is getting better at what it does. Each new feature provides leverage to build upon previous features in more elaborate ways. Key to improvement is (a) maintaining the core features that other features rely on, and (b) adding new features that are solid stepping stones to other features. Improvement most often happens by getting the small and simple things right incrementally, rather than with sweeping overhauls all at once.

5. Collaboration: If you see a way to improve Pynite, you are encouraged to contribute. There are many simple ways to contribute that don't take much effort. Issue reports and pull requests can be very helpful. One easy way to contribute is to add to or improve the library of simple example problems in the `Examples` folder. Please keep them relatively simple. Most users learn Pynite from following these simple examples. Another way to help Pynite without having to know too much about its internal workings is to help with the documentation files in the `docs\source` folder of this repository. These files help new users learn Pynite. If you are able to make a bigger commitment, and would like to become a regular contributor to the project, please reach out about becoming a collaborator.

# Support
Whether you just need help getting started with Pynite, or are looking to build something more complex, there are a few resources available:
* The examples in the "Examples" folder in this repository cover a variety of simple problems. The comments in the examples provide additional guidance on how Pynite works.
* Documentation is a work in progress and can be found on readthedocs here: https://pynite.readthedocs.io/en/latest/index.html.
* If you're looking for more direct guidance on using Pynite, or for help coding a project, I am available on a private consulting basis. You can reach out to me directly at Building.Code@outlook.com to discuss options.

# Sponsors
* A special thanks to @edson-gusmoes for sponsoring `Pynite`!

# Example Projects
Here's a list of projects that use Pynite:

* Building Code (https://building-code.herokuapp.com/) - This one is my personal side project.
* Standard Solver (https://www.standardsolver.com/)
* Civils.ai (https://civils.ai/1/free-3D-finite-element-structural-analysis)
* Phaenotyp (https://github.com/bewegende-Architektur/Phaenotyp) (https://youtu.be/shloSw9HjVI)

# What's New?
v1.2.0
* Added the ability to apply loads in steps in the `FEModel3D.analyze` method via the new `num_steps` keyword argument. `num_steps` defaults to 1, but if a higher number is used the load will be split into the specified number of steps and be applied incrementally. This can be very helpful when dealing with complex tension/compression-only scenarios that otherwise may have had trouble converging. Consider using 10 to 20 load steps when dealing with complex T/C-only models. The ability to perform step-wise analysis better approximates nonlinear load-displacement curves and sets `Pynite` up for future nonlinear analysis features.
* Added P-&Delta; effects to vectorized results. Prior to this fix, member end forces and internal forces were being calculated correctly, but the P-&Delta; flag was not triggered in the methods that used vectorized results. Moment plots and moment arrays were the only items that were affected. Spot checking moments, and checking max/min moments were not affected.

v1.1.2
* Corrected a long-standing issue transforming quad local bending and membrane stresses to global coordinates. This did not affect quad local stress results. Only if you were converting the results to global coordinates would the issue arise.
* Array results no longer include extra points at discontinuities. This was affecting some users who were using array results. Note that without the extra points, you may need to use a larger number of points to identify max/min values at discontinuities.
* Added axial stiffness adjustments to P-&Delta; analysis. P-&Delta; effects primarily affect bending stiffness, but axial stiffness is affected too. While these axial stiffness adjustments are often negligible, they can be important in some cases. `Pynite` now considers axial stiffness adjustments when running P-&Delta; analysis.
* Greatly simplified code for P-&Delta; analysis. P-&Delta; analysis should converge quicker with fewer unnecessary iterations. `Pynite` was not taking full advantage of the fact that the geometric stiffness matrix eliminates the need to iterate when solving for P-&Delta; effects. Extra solutions were being run for no reason.
* Cleaned up many "code smells".

v1.1.1
* Bug fix for member array results. Exceptions were being thrown in some cases.
* Simplified code for array functions by removing redundant algorithms.

v1.1.0
* Added a `VTKWriter` class to allow for easy exporting to `Paraview`.
* Improved type hints for a simpler user experience.
* Fixed a bug in the pdf reports. Load combos were not being sent to the report template, which was preventing any load results from being displayed in the report.
* Simplified code for vector extraction of member segment results, and improved results reporting at load discontinuities along beams. Two results values are possible at mathematical discontinues, and only one was being reported. Now both are reported. This issue was only noticable if array results along a member were requested at large intervals rather than small intervals.
* Bug fix for array results format for physical members. This only affected users who used the array results functions.

v1.0.1
* Changes to testing code coverage less than 2% no longer trigger build failure.
* Code cleanup: removed `DKMQ.py` from the repository that was no longer in use. The code for the DKMQ element now lives in `Quad3D.py` instead.
* Frustrum meshes generated about the global X and Z axes are now being generated correctly.
* Fixed a bug that was not letting plate contours render for `VTK` users using load combinations other than 'Combo 1'. This bug was introduced recently when the option for global stress results was added. Global stress results for `VTK` rendering are still not supported yet. Users are urged to switch to using `pyvista` rendering instead as `VTK` rendering is on its way out of `Pynite` and may be only minimally maintained going forward.

v1.0.0
* v1.0 is here! I feel the program is stable enough and has been around long enough to be battle tested and to call it v1.0.
* Important!!! - Changed all calls from `PyNite` to `Pynite` this matches the logo, and made more sense. I'm not sure why I ever capitalized that N to begin with, but going forward from v1.0, `Pynite` has a lowercase n. I've been wanting to make this change as part of the v1.0 release.
* Added a new `ShearWall` class that assists you in constructing and analyzing shear walls. This tool automatically detects piers and coupling beams, and finds the forces inside them and calculates their ascpect ratios, which can be handy for seismic design. It reports stiffness of multi-story shear walls at each story to help with rigid diaphragm analysis. It allows for modeling walls with openings, steps, and partial depth diaphragm loading.
* `vtk` and `pyvista` are now optional dependencies. This change streamlines installation for users who don't rely on `Pynite's` built-in visualization tools. From now on, `Pynite` should be installed using `$ pip install PyNiteFEA[all]` for most users.
* `Pynite` no longer uses auxiliary nodes to define member cross-section rotation. You can now directly specify the rotation (in degrees) when you define a member using the `rotation` argument.

v0.0.98-100
* Bug fix for `FEModel3D.add_section`. It was throwing exceptions and had not been updated to match the examples.
* Improvements to spring rendering in `pyvista`. Up until this point spring elements were being rendered as lines. They now render as zigzag lines in `pyvista`. There is still more work for improvement on spring rendering, but this is a good start.

v0.0.97
* Fixed physical member load and deflection diagrams. Physical members are a newer feature. Member internal results were being reported correctly, but the diagrams for these members had not been revised to plot correctly. The old method for plain members was still being used. Physical members were not considering that a physical member was made from multiple submembers, and results for each span needed to be combined to get the whole plot.
* Switched some commonly used python libraries to be installed by default with `Pynite`. Most `Pynite` users will want these libraries installed for full-featured use of `Pynite`. These libraries help with `Pynite` visualizations, plotting, the sparse solver, and `Jupyter Lab` functionality. This is just easier for new python users. I was getting a lot of questions about how to set up libraries, and this takes the guesswork away. This is part of `Pynite's` objective to stay easy to use.

v0.0.96
* Changed quad elements from MITC4 formulation to DKMQ formulation. This greatly improves plate results at corners and increases the speed with which the plate's stiffness matrix is assembled. MITC4 element code has been retained as legacy code, but is no longer used by the program.
* ***Breaking Changes***: Implemented snake-case for dictionary names (e.g. `FEModel3D.Nodes` is now `FEModel3D.nodes`). These changes were made to prepare `Pynite` for a v1.0 release that is consistent with the `PEP8` style guide for `python`.
* Bug fix for tension/conpression-only member internal results. While global results were correct, member internal results were showing results from the first tension/compression only iteration.
* Member results arrays can now be customized to pick up user defined points. Member results arrays generate results much faster now too.

v0.0.95
* Bug fix for rendering negative point loads via `Pyvista`. They were being rendered as positive loads. The analysis was not impacted by this bug.

v0.0.94
* Added rendering via `Pyvista`. This greatly simplified the rendering code and provided a fresh look to the rendereings. Renderings in jupyter are now interactive. Global axes are also now shown in rendereings. To use `Pyvista` instead of `VTK`, use the new `Rendering` library rather than the old `Visualization` library. Rendering via `VTK` directly is still available.
* Bug fix for member self-weight. The program was throwing exceptions instead of calculating member self-weight. Added a unit test to help prevent this issue from occuring again as code changes.
* Refactored `material` to be `material_name` in the code. The prior naming convention caused confusion which led to the self-weight bug.

v0.0.93
* Fixed phantom reactions showing up at unsupported nodes. If there was a support defined at a node, the program was summing reactions for all directions at the node, rather than just the supported directions. This caused the program to report "extra" reaction directions at any supported node (if the user queried them). Element forces/stresses were not affected as this was a post-processing reaction summing issue. Reactions for supported directions were summed correctly, except in the case of nodes with both spring supports and other supports. Only unsupported directions, and nodes with both spring supports and other supports, were showing phantom reactions. This bug also caused statics checks to fail from time to time.
* Reorganized physical member code to match member code more consistently.
