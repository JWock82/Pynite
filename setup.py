import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyNiteFEA",
    version="2.0.2",
    author="D. Craig Brinck, PE, SE",
    author_email="Building.Code@outlook.com",
    description="A simple elastic 3D structural finite element library for Python.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JWock82/Pynite.git",
    packages=setuptools.find_packages(include=['Pynite', 'Pynite.*']),
    package_data = {'Pynite': ['*html', '*.css', '*Full Logo No Buffer - Transparent.png']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'numpy',
        'PrettyTable',
        'scipy',
        'matplotlib',
    ],
    extras_require = {
        'all': ['Ipython', 'vtk>=9.3.0', 'pyvista[all,trame]>=0.43.4', 'trame_jupyter_extension', 'ipywidgets', 'pdfkit', 'Jinja2'],
        'vtk':  ['IPython', 'vtk>=9.3.0'],
        'pyvista': ['pyvista[all,trame]>=0.43.4', 'trame_jupyter_extension', 'ipywidgets'],
        'reporting': ['pdfkit', 'Jinja2'],
        'derivations': ['jupyterlab', 'sympy']
    },
    include_package_data = True,
    python_requires = '>=3.10',
)
