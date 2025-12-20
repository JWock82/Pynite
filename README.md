<div align="center">
  <img src="https://github.com/JWock82/Pynite/raw/main/Resources/Full Logo No Buffer.png" width=40% align="center"/>
  <br>
  <h1>Simple Finite Element Analysis in Python</h1>
</div>

![Build Status](https://github.com/JWock82/Pynite/actions/workflows/build-and-test.yml/badge.svg)
[![codecov](https://codecov.io/gh/JWock82/Pynite/branch/main/graph/badge.svg?token=ZH18US3A7P)](https://codecov.io/gh/JWock82/Pynite)
![PyPI - Downloads](https://img.shields.io/pypi/dm/PyNiteFEA?cacheSeconds=86400)
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
* Modal analysis of frame type structures.
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
* Advanced tools for modeling and analyzing complex shear walls and mat foundations.
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
v2.0.2
* Added docstrings to the VTK `Renderer` class to help the user.
* Enforced use of properties instead of attributes in the VTK `Renderer` and added docstrings to properties help the user make decisions.
* Removed the VTK `Renderer`'s properties that began with "set_". These were redundant and caused confusion. The prefix "set_" is no longer used to access or set any of the `Renderer`'s properties.

v2.0.1
* Pynite no longer struggles rendering load cases and combinations when no loads are present. This is especially helpful for rendering modal load combinations, which don't have loads (only masses).
* The pyvista plotter's "X" button now works just like pressing "q". This was an annoying pyvista-ism that led to error messages when closing the window with the "X" button.

v2.0.0
* Added modal analysis. You can now find the natural frequencies and mode shapes of frame-type structures (plates not yet supported). A special thanks to @boustrephon for heading up this feature.
* Dropping support for python < 3.10.

v1.6.2
* Added helper method to `FEModel3D` for mat foundations.

v1.6.1
* Added mesh control for mat foundations.
* Added soil pressure calculatiosn for mat foundations.
* Added docstrings for mat foundations.
* Mat foundations should still be considered a "beta" feature until further testing and development is complete.

v1.6.0
* Added mat foundations to `Pynite`. This is still a "beta" feature. Mat foundations are functional, but they are still fairly basic and are largely untested. They don't yet report soil pressures, or allow for you to set effective strip widths. They do support point loads applied to the mat, and they automatically assign soil springs when you supply a subgrade modulus of reaction. Use at your own risk.
* Improved speed for descritizing models with many physical members.
* Adjusted a few display colors for the VTK Renderer's 'print' theme in an effort to make more printer friendly visualizations.
* Added more tests to expand code coverage of visualization.
* Major improvements to documentation program wide.

v1.5.0
* Extended max/min enveloping to plate/quad meshes. Now you can quickly check the max/min mesh forces accross multiple load combinations in one simple command by passing a list of `combo_tags` to the functions.
* Added max/min enveloping to torsion.
* Refactored shear walls to be more `Pynite`-esque. They work similar to other items in `Pynite` now within the `FEModel3D` class, rather than being a stand-alone type of model. Multiple shear walls can be in one model now, and they can be defined parallel to any major global plane and originating at any position. They are much more user-friendly than they were before.
* Fixed an error in member rotations that has persisted since v1.0.0. Rotations about the member's local x-axis were not being applied correctly. They are now applied using the Rodrigues formula.
* Changed the way $P-\delta$ member internal deflections are calculated. Before the program would iterate, sometimes without a solution. Now a closed-form solution for $P-\delta$ member internal deflections is used. This should speed up deflection calculations for members, and be more robust.

v1.4.0
* Added the ability to check across multiple combinations for max/min member force and deflection results using combo tags. Simply substitute a list of combo tags to the functions instead of a load combination name and the functions will envelope results across any combinations having the given tags. No more searching all combinations manually!
* Improved the efficiency of member max/min result functions by removing redundant segmentation routine calls.
* Set the Y-axis to vertical for 3D rendering in Pyvista. This allows the "isometric" button to display the model in the correct orientation in Jupyter.
* Continued working toward pushover analysis. Development of this feature has been slower than anticipated. Some major hurdles have been overcome, but the code for this feature is still largely experimental, and may be so for sometime.

v1.3.1
* Reverted to an older Pynite convention for comparing boolean operators. Users who used `1` and `0` for boolean inputs in methods instead of `True` and `False` were experiencing issues. This change allows for multiple ways for users to input booleans.

v1.3.0
* Fixed shear wall thickness issues. The new `ShearWall` feature was not applying the wall thickness correctly. This has been fixed.
* Fixed a bug in reaction calculations for positive nodal springs. Only reactions for positive springs (e.g. tension-only) were affected.
* Refactored `AxialDeflection` to `axial_deflection` to be consistent with the rest of the code.
* Fixed a bug in reporting axial deflection arrays. The program was throwing exceptions when axial deflection array results were requested.
* Fixed a bug preventing reports from generating.
* Added support for html reporting, and better error messages for getting `pdfkit` and `wkhtmltopdf` setup for reporting.

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
