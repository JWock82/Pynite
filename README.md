# PyNite
A linear elastic 3D structural engineering finite element library for Python.

# Current Capabilities
* 3D static analysis of elastic structures.
* P-&Delta; analysis of frame type structures.
* Member point loads, linearly varying distributed loads, and nodal loads are supported.
* Classify loads by load case and create load combinations from load cases.
* Produces shear, moment, and deflection results and diagrams for each member.
* Rectangular plate elements.
* Reports support reactions.
* Rendering of model geometry, supports, load cases, load combinations, and deformed shapes.

# Project Objectives
As I've gotten into the structural engineering profession, I've found there's a need for an easy to use open-source finite element package. I hope to help fill that need by prioritizing the following:

1. Accuracy: There are no guarantees PyNite is error free, but accuracy and correctness are a priority. When bugs or errors are identified, top priority will be given to eliminate them. PyNite's code is frequently reviewed, and its output is regularly being tested against problems with known solutions to isolate errors. If you find an error, please report it as an issue.

2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in outdated languages, or require extensive knowledge of finite element theory and/or element formulations to use. PyNite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package.

4. Improvement: I'm not inexperienced, but I'm still relatively young in my career with many working years ahead of me. I plan to continue developing/improving PyNite for many years to come. There are a lot of pieces I'd like to add to PyNite going forward: improvements to plates, dynamics, pushover anlysis, etc. There's a lot of potential to create extensions as well to solve all kinds of engineering problems. There are more problems to solve than I have time for, so some priorities will have to be made. The plan is to keep PyNite mainstream, adding core functionality first. Occasionally however I may just add what interests me at the time.

5. Collaboration: The intent is to keep PyNite free and open source. This will encourage future development and contributions. Keeping it open source will allow anyone to inspect and improve the code it runs on. If you see an area you think you can help PyNite improve in you are encouraged to contribute. I'd like to get PyNite doing a lot more. Don't be offended if I'm a little slow to accept your contributions. FEA is a very technical subject and accuracy is extremely important to me. Sometimes I'm a little slow understanding both FEA and Python and it takes some time for me to comprehend what's being proposed. I also have a young family to take care of that takes first priority.

# Dependencies
PyNite depends on the following packages:
* numpy: used for matrix algebra
* matplotlib: used for plotting member diagrams
* tabulate : used to format tabular output
* vtk: used for visualization - Note that VTK is a little picky about which version of Python you are running. You must run a 64 bit installation of Python, rather than a 32 bit version. The latest version of Python 3.8+ is not supported by VTK yet. I recommend running a 64 bit version of Python 3.7+. VTK does not need to be installed if you don't plan to use the visualization tools built into PyNite.

# What's New?
Version 0.0.12
Load rendering has been added for load cases and combinations! It should be noted that the `Visualization.DeformedShape` function has been dapricated. All deformed shapes are now accessed through the `RenderModel` function by setting the `deformed_shape` variable to `True`. This was done to keep PyNite simple, as it allowed for a lot of redundant code to be removed.

Version 0.0.11:
Load cases and load combinations have been added! Each load can now be assigned a load case (e.g. 'D', 'L', 'S', or 'any string you want'), and you can specify multiple load combinations (e.g. '1.2D + 1.6L + 0.5S') with different load factors for each load case within the combination. You may notice slightly different behavior from PyNite now that load combinations have been added. Most results are now presented by load combination in a Python dictionary format. A few of the examples have been updated to demonstrate this feature. More examples will be updated soon. A full description of how to use this feature will also be posted shortly.
