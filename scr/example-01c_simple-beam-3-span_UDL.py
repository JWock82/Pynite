# PyNite Example
# Example: Simply supported beam with a point load
# Units: inches and kips, use consistent units

from operator import ge
import numpy as np
from PyNite import FEModel3D
from PyNite import Reporting
from PyNite import Visualization
from pathlib import Path
import matplotlib.pyplot as plt
from itertools import accumulate

# Updated parameters for this plot
plt.style.use("../plt/figstylefile.mplstyle")


class Section:
    def __init__(
        self,
        section_name,
        section_area,
        moment_inertia_zaxis,
        moment_inertia_yaxis,
        moment_inertia_polar=None,
        young_mod=29000,
        shear_mod=11400,
    ):
        """Section Properites

        Parameters
        ----------

        section_name : string
            Section name
        section_area : float
            Cross-sectional area
        moment_inertia_zaxis : float
            Moment of inertia about centroidal local z-axis
        moment_inertia_yaxis : float
            Moment of inertial about centroidal local y-axis
        moment_inertia_polar : float
            Polar moment of inertia
        young_mod : float
            Young's modulus of the section material
        shear_mod : float
            Shear modulus of the section material

        """

        self._section_name = section_name
        self._area = section_area
        self._Izz = moment_inertia_zaxis
        self._Iyy = moment_inertia_yaxis
        self._modE = young_mod
        self._modG = shear_mod
        if moment_inertia_polar == None:
            self._J = self._Izz + self._Iyy
        else:
            self._J = moment_inertia_polar

    @property
    def section_name(self):
        return self._section_name

    @property
    def area(self):
        return self._area

    @property
    def Izz(self):
        "Based on local axis"
        return self._Izz

    @property
    def Iyy(self):
        "Based on local axis"
        return self._Iyy

    @property
    def J(self):
        """Polar moment of inertia"""
        return self._J

    @property
    def modE(self):
        """Young's modulus"""
        return self._modE

    @property
    def modG(self):
        """Shear modulus"""
        return self._modG

    def __repr__(self):
        return f"Section(name={self._section_name}, area={self._area}, Izz={self._Izz}, Iyy={self._Iyy}, E={self._modE}, G={self._modG})"


def gen_plt_data(beam, member):
    """Generate data for plotting at distances along the local x-axis.

    Parameters
    ----------

    beam : FEModel3d object
        Model object
    member : string
        Name of the member to generate the data

    Returns
    -------
    xs : 1D numpy array
        Locations along the location x-axis at which values are computed
    shear_vals : [] list of floats
        Values of shear at xs
    moment_vals : [] list of floats
        Values of bending moment at xs
    deflection_vals : [] list of floats
        Values of deflections at xs
    """

    xs = np.linspace(0, beam.Members[member].L(), 31)
    shear_vals = []
    moment_vals = []
    deflection_vals = []

    for x in xs:
        shear_vals.append(beam.Members[member].shear("Fy", x))
        moment_vals.append(beam.Members[member].moment("Mz", x))
        deflection_vals.append(beam.Members[member].deflection("dy", x))

    return xs, shear_vals, moment_vals, deflection_vals


def plt_combined(beam, fpath):
    """Plot shear force, bending moment and deflection in a single plot.

    Parameters
    ----------
    beam : FEModel3D object
        Beam FEModel3D object
    fpath : pathlib Path object or str
        Full path including filename of the figure to save the figure.
    """

    # figure
    figscale = 1.0
    figw = 3.487
    figh = 3.5 * figw / 1.618
    figsz = (figscale * figw, figscale * figh)
    fig, axs = plt.subplots(3, 1, figsize=figsz, sharex=True)

    # data for shear, moment and deflection diagrams
    xs = []
    shear_vals = []
    moment_vals = []
    deflection_vals = []

    for member in beam.Members:
        xss, shears, moments, deflections = gen_plt_data(beam, member)
        xs.append(xss)
        shear_vals += shears
        moment_vals += moments
        deflection_vals += deflections

    # fmt : off
    # distances along the beam, convert from local coordinates to global coordinates, must be a better way to do this
    xs = [list(np.diff(x, prepend=0)) for x in xs]
    xs = np.cumsum(np.concatenate(xs).ravel())
    # fmt : on

    # plots : shear
    axs[0].plot(xs, shear_vals)
    axs[0].plot(
        [xs[0], xs[0], xs[-1], xs[-1]],
        [shear_vals[0], 0, 0, shear_vals[-1]],
        color="k",
        lw=2.0,
    )
    axs[0].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[0].set_ylabel("Shear Force, $V$ [kip]")

    # moment
    axs[1].plot(xs, moment_vals)
    axs[1].plot(
        [xs[0], xs[0], xs[-1], xs[-1]],
        [moment_vals[0], 0, 0, moment_vals[-1]],
        color="k",
        lw=2.0,
    )
    axs[1].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[1].set_ylabel("Bending Moment, $M$ [kip-in]")

    # deflections
    axs[2].plot(xs, deflection_vals)
    axs[2].plot(
        [xs[0], xs[0], xs[-1], xs[-1]],
        [deflection_vals[0], 0, 0, deflection_vals[-1]],
        color="k",
        lw=2.0,
    )
    axs[2].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[2].set_xlabel("Distance, $x$ [in]")
    axs[2].set_ylabel("Deflection, $\delta$ [in]")

    fig.savefig(fpath, format="png", dpi=150, bbox_inches="tight")


def main():
    # paths and filenames, need to create these folder structure
    path_root = Path(r"../")
    path_plt = path_root.joinpath("plt/")
    path_report = path_root.joinpath("rpt/")
    fname_vtk = "example-01c-vtk.png"
    fname_fig = "example-01c.png"
    fname_rpt = "example-01c_rpt.pdf"
    ffigvtk = path_plt.joinpath(fname_vtk)
    frpt = path_report.joinpath(fname_rpt)
    ffig = path_plt.joinpath(fname_fig)

    # beam FEM model
    sbeam = FEModel3D()

    # define section: section_name, cross-section area, Iz, Iy, J, E, G
    # J=None, E=29000 ksi, G=11400 ksi are optional parameters
    section_1 = Section("section-1", 20, 150, 100)

    # Point load [kip]
    udl = -10
    # lengths of each span [ft]
    spans = [10, 10]
    # global x-coord for the nodes
    coord_nodes = [0] + [span * 12 for span in spans]
    coord_nodes = list(accumulate(coord_nodes))

    # add nodes (node_label, x, y, z)
    for idx, coord in enumerate(coord_nodes):
        print(f"Node: { 'n' + str(idx)}, Coords: ({coord}, 0, 0) ")
        sbeam.add_node("n" + str(idx), coord, 0, 0)

    # add element between the two nodes and define material
    # (element_label, nodei, nodej, E, G, Iy, Iz, J, A)
    for idx in range(len(spans)):
        print(
            f"Element: {'m' + str(idx)}, Node i: {'n' + str(idx)}, Node j: {'n' + str(idx+1)}"
        )
        sbeam.add_member(
            "m" + str(idx),
            "n" + str(idx),
            "n" + str(idx + 1),
            section_1.modE,
            section_1.modG,
            section_1.Iyy,
            section_1.Izz,
            section_1.J,
            section_1.area,
        )

    # simple support at the nodes(node_label, Dx, Dy, Dz, Rx, Ry, Rz)
    # if true the support is restrained in that direction,
    # D is displacement, Dx is True, then the support is restrained against displacement in the X direction,
    # R is rotation, if Rx is True, then the support is restrained
    # against rotation about the x axis
    sbeam.def_support("n0", True, True, True, True, False, False)
    for idx in range(len(spans)):
        print(f"Node: {'n' + str(idx + 1)}, Restraint: True, True, True")
        sbeam.def_support("n" + str(idx + 1), True, True, True)

    for idx, span in enumerate(spans):
        print(f"Element: {'m' + str(idx)}, 'Fy', {udl}, {udl}")
        sbeam.add_member_dist_load("m" + str(idx), "Fy", udl, udl)

    # analyze beam and check static checks
    sbeam.analyze(check_statics=True)

    # shear moment and deflection diagrams
    plt_combined(sbeam, fpath=ffig)

    if True:
        Visualization.render_model(
            sbeam,
            annotation_size=10,
            deformed_shape=True,
            deformed_scale=30,
            render_loads=True,
            combo_name="Combo 1",
            screenshot=str(ffigvtk),
        )

    if True:
        Reporting.create_report(
            sbeam,
            output_filepath=str(frpt),
            plates=False,
            plate_corner_forces=False,
            plate_center_forces=False,
            plate_corner_membrane=False,
            plate_center_membrane=False,
        )


if __name__ == "__main__":
    main()
