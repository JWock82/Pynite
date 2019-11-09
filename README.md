# PyNite
A linear elastic 3D frame, truss, and beam finite element library for Python.

# Project Objectives
After 12 years as a structural engineer, I've found there's a need for an easy to use open-source finite element package. I hope to help fill that need by prioritizing the following:

1. Accuracy: Innacurate finite element solutions can have serious consequences. There are no guarantees PyNite is error free, but accuracy and correctness are a priority. When bugs or errors are identified, top priority will be given to eliminate them. PyNite's code is frequently reviewed, and its output is regularly being tested against problems with known solutions to isolate errors. If you find an error, please report it.
2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in outdated languages, or require extensive knowledge of finite element theory and element formulations to use. PyNite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package. While languages like C++ could have provided more horsepower, Python was selected as the language because of its simplicity and its excellent scientific packages like numpy and matplotlib.
3. Improvement: I'm relatively young in my career as a structural engineer. I plan on being an engineer for many more years to come. PyNite is a tool I use to improve my understanding as an engineer, and I plan to continue learning new things, bringin PyNite along with me. There are a lot of pieces I'd like to add to PyNite going forward: plates, P-Delta analysis, dynamics, pushover anlysis, etc. There's a lot of potential to create extensions as well to solve all kinds of engineering problems. I may not have time to add all these things myself, but I'm commited to improvement.
4. Community Involvement: This project could really take off with the help of the structural/mechanical engineering community. I'd love to see people who know what they're doing better than I do help provide guidance to this project. I'm a design engineer, not a theorist.

# Dependencies
PyNite depends on the following packages:
* numpy: used for matrix algebra
* matplotlib: used for plotting member diagrams
* vtk: used for visualization - note that vtk requires a 64 bit installation of python. vtk does not need to be installed if you don't plan to used the visualization tools in PyNite.

# Contributions Welcome
If you would like to contribute to this project to expand, debug or improve it, I welcome your contributions. I'd like to get the library doing a lot more than just solving frames. Don't be offended if I'm a little slow to accept your contributions. FEA is a very technical subject and accuracy is extremely important to me. Sometimes I'm a little slow understanding both FEA and Python and it takes some time for me to comprehend what's being proposed. Also I have 4 little kids to take care of that divert a lot of my attention :)
