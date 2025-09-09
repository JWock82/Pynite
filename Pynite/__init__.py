# Select libraries that will be imported into Pynite for the user
from Pynite.FEModel3D import FEModel3D
from Pynite.ShearWall import ShearWall
import Pynite

from pip._vendor import pkg_resources

def get_version(package):
    package = package.lower()
    return next((p.version for p in pkg_resources.working_set if p.project_name.lower() == package), "No match")

__version__= get_version("PyniteFEA")
