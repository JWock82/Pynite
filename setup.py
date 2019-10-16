import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyNiteFEA",
    version="0.0.3",
    author="D. Craig Brinck, SE",
    author_email="Building.Code@outlook.com",
    description="A 3D frame, beam, and truss finite element package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JWock82/PyNite.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)