# PyNite
A linear elastic 3D structural engineering finite element library for Python.

# Current Capabilities
* 3D static analysis of elastic structures.
* P-&Delta; analysis of frame type structures.
* Member point loads, linearly varying distributed loads, and nodal loads are supported.
* Classify loads by load case and create load combinations from load cases.
* Produces shear, moment, and deflection results and diagrams for each member.
* Tension-only and compression-only elements.
* Springs: two-way, tension-only, and compression-only.
* Rectangular plate elements.
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
* numpy: used for matrix algebra
* matplotlib: used for plotting member diagrams
* PrettyTable : used to format tabular output
## Optional Dependencies
* vtk: used for visualization - Note that VTK is a little picky about which version of Python you are running. You must run a 64 bit installation of Python, rather than a 32 bit version. The latest version of Python 3.8+ is not supported by VTK yet. I recommend running a 64 bit version of Python 3.7+. VTK does not need to be installed if you don't plan to use the visualization tools built into PyNite.
* pdfkit: Used for generating pdf reports. In order to generate pdf reports, PDFKit requires you to have wkhtmltopdf installed on your computer. This is a free program available for download at https://wkhtmltopdf.org/downloads.html. Once installed, you'll need to help PyNite find it. On Windows, this can be done by setting your PATH environment variable to include the path to "wkhtmltopdf.exe" after installation. For example, mine is installed at "C:\Program Files\wkhtmltopdf\bin"
* jinja2: Used for templating reports into HTML prior to HTML-to-pdf conversion.

# Support and Feature Requests
I am a structural engineer first and foremost. If you're looking for help understanding, applying, interpreting, or extending the engineering behind PyNite, I'd be happy to help you on a private consulting basis. I have almost 14 years experience designing steel, concrete, aluminum, and masonry structures for many large scale multi-disciplinary projects. I have a strong background in U.S. building code requirements. You can contact me directly at <Building.Code@outlook.com> to discuss options.

I get more feature requests than I have time to handle. I chip away at them as I find time, but it's a slow process. However, if there is a feature you really want expedited, especially a complex one that takes considerable effort, I'm generally willing to devolop it on a consulting basis. I'd prefer to keep any new features I develop free for all users under the MIT License however. You are still welcome to use it for your proprietary software. That promotes expansion of PyNite's capabilities and more opportunities for me to improve PyNite for everyone. Again, you can contact me at <Building.Code@outlook.com> to discuss options.

# What's New?
Version 0.0.15:
Added springs to PyNite. Springs can be two-way, tension-only, or compression-only. Also, fixed numerous bugs relating to tension-only/compression-only analysis features in version 0.0.14.

Version 0.0.14:
Clean looking pdf reports can now be generated by PyNite. See "Simple Beam - Point Load.py" and "Out-of-Plane Wall Panel.py" for examples on how to use this feature. The wiki also provides a detailed explanation of how to implement PDF reporting. Input validation has also been improved thanks to the help of tamalone1.

Members can now be defined as tension-only or compression-only. In the process of adding this feature, a bug was removed from the P-Delta analysis that used to cause P-Delta analyis to crash if any of the members didn't have axial load on them. A new example demonstrating how to perform tension-only analysis has been added to the `Examples` folder.

Version 0.0.13:
Plate in-plane (membrane) stress results are now available. Before only out-of-plane stresses were reported. Documentation will be added to the wiki soon.

PyNite now uses `PrettyTable` instead of `tabulate`. The switch was made because `tabulate` has stricter licensing (GPLv3). `PrettyTable` is BSD (3 clause) licensed. What this means is that `tabulate` could not be used in proprietary software, whereas `PrettyTable` can.

Version 0.0.12:
Load rendering has been added for load cases and combinations! It should be noted that the `Visualization.DeformedShape` function has been dapricated. All deformed shapes are now accessed through the `RenderModel` function by setting the `deformed_shape` variable to `True`. This was done to keep PyNite simple, as it allowed for a lot of redundant code to be removed. Check out the wiki for more clarification.

Version 0.0.11:
Load cases and load combinations have been added! Each load can now be assigned a load case (e.g. 'D', 'L', 'S', or 'any string you want'), and you can specify multiple load combinations (e.g. '1.2D + 1.6L + 0.5S') with different load factors for each load case within the combination. You may notice slightly different behavior from PyNite now that load combinations have been added. Most results are now presented by load combination in a Python dictionary format. A few of the examples have been updated to demonstrate this feature. More examples will be updated soon. A full description of how to use this feature will also be posted shortly.
