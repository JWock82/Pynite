# PyNite
A linear elastic 3D structural engineering finite element analysis library for Python.

# Current Capabilities
* 3D static analysis of elastic structures.
* P-&Delta; analysis of frame type structures.
* Member point loads, linearly varying distributed loads, and nodal loads are supported.
* Classify loads by load case and create load combinations from load cases.
* Produces shear, moment, and deflection results and diagrams for each member.
* Tension-only and compression-only elements.
* Springs: two-way, tension-only, and compression-only.
* Rectangular and quadrilateral plate elements.
* Reports support reactions.
* Rendering of model geometry, supports, load cases, load combinations, and deformed shapes.
* Generates PDF reports for models and model results.

# Project Objectives
As I've gotten into the structural engineering profession, I've found there's a need for an easy to use open-source finite element package. I hope to help fill that need by prioritizing the following:

1. Accuracy: There are no guarantees PyNite is error free, but accuracy and correctness are a priority. When bugs or errors are identified, top priority will be given to eliminate them. PyNite's code is frequently reviewed, and its output is regularly being tested against problems with known solutions to isolate errors. If you find an error, please report it as an issue.

2. Simplicity: There are other finite element alternatives out there with many more capabilities, but they are often lacking in documentation, written in outdated languages, or require extensive knowledge of finite element theory and/or element formulations to use. PyNite is not intended to be the most technically advanced solver out there. Rather, the goal is to provide a robust yet simple general purpose package.

4. Improvement: I plan to continue supporting PyNite for many years to come. There are a lot of pieces I'd like to add to PyNite going forward: improvements to plates, dynamics, pushover anlysis, etc. There's a lot of potential to create extensions as well to solve all kinds of engineering problems. There are more problems to solve than I have time for, so some priorities will have to be made. The plan is to keep PyNite mainstream, adding core functionality first. Occasionally however I may just add what interests me at the time.

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

# What's New?
Version 0.0.24: Corrections to reactions for plate and quad elements. Fixed member concentrated loads that were not rendering with the rest of the model.

Version 0.0.23: Bug fix for concentrated moments. Fixed end reactions on one side of the Member3D element were being calculated incorrectly for concentrated moments.

Version 0.0.22: This version makes some significant changes. Major improvements include:
* Solution speed has been greatly improved for large models. The code has been profiled and many slow spots in the solution process have been isolated and removed.
* A sparse matrix solver has been added as the default solver to improve solution speed. The sparse solver may be slower on smaller models, but those models usually solve faster anyway. The dense matrix solver can be used if this is a problem.
* Quadrilateral elements have been implemented.
* Rectangular plate nodes are now defined clockwise instead of counter-clockwise. This change was made in order to be consitent with the new quad element formulation.
* Rendering of plate and quad stresses is now supported.
* Plate surface pressure loads are now supported.
