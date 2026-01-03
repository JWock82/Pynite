# Select libraries that will be imported into Pynite for the user
from Pynite.FEModel3D import FEModel3D
from Pynite.ShearWall import ShearWall
import Pynite

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("PyniteFEA")
except PackageNotFoundError:
    # Package not installed, use development version
    __version__ = "dev"
