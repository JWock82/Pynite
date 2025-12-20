# Select libraries that will be imported into Pynite for the user
from Pynite.FEModel3D import FEModel3D
from Pynite.ShearWall import ShearWall
import Pynite

from importlib.metadata import version

__version__ = version("PyniteFEA")
