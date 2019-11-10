# PyNite
A linear elastic 3D frame, truss, and beam finite element library for Python.

# Current Capabilities
* 3D static analysis of linear elastic frame structures.
* Member point loads, linearly varying distributed loads, and nodal loads are supported.
* Produces shear, moment, and deflection results and diagrams for each member.
* Reports support reactions.
* Basic rendering of model geometry.

# Project Objectives
After 12 years as a structural engineer, I've found there's a need for an easy to use open-source finite element package. I hope to help fill that need by prioritizing the following:

1. Accuracy: There are no guarantees PyNite is error free, but accuracy and correctness are a priority. When bugs or errors are identified, top priority will be given to eliminate them. PyNite's code is frequently reviewed, and its output is regularly being tested against problems with known solutions to isolate errors. If you find an error, please report it as an issue.

2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in outdated languages, or require extensive knowledge of finite element theory and/or element formulations to use. PyNite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package.

4. Improvement: I plan to continue developing/improving PyNite for many years to come. There are a lot of pieces I'd like to add to PyNite going forward: plates, P-Delta analysis, dynamics, pushover anlysis, etc. There's a lot of potential to create extensions as well to solve all kinds of engineering problems. There are more problems to sove than I have time for, so some priorities will have to be made. The plan is to keep PyNite mainstream, adding core functionality first. Occasionally however I may just add what interests me at the time.

5. Collaboration: This project is ambitious for me to tackle as a lone engineer, but it could really take off with the help of the structural/mechanical engineering community.

# Dependencies
PyNite depends on the following packages:
* numpy: used for matrix algebra
* matplotlib: used for plotting member diagrams
* vtk: used for visualization - note that vtk requires a 64 bit installation of python. vtk does not need to be installed if you don't plan to used the visualization tools in PyNite.

# Contributions Welcome
If you would like to contribute to this project to expand, debug or improve it, I welcome your contributions. I'd like to get the library doing a lot more than just solving frames. Don't be offended if I'm a little slow to accept your contributions. FEA is a very technical subject and accuracy is extremely important to me. Sometimes I'm a little slow understanding both FEA and Python and it takes some time for me to comprehend what's being proposed. Also I have 4 little kids to take care of that divert a lot of my attention :)
