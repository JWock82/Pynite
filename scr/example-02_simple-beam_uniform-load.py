# Example 02: Simply supported beam with a uniform loading
# Units: inches, kips

import numpy as np
from pathlib import Path
from PyNite import FEModel3D
from PyNite import Reporting
import matplotlib.pyplot as plt

# Updated parameters for this plot
plt.style.use("../plt/figstylefile.mplstyle")


def plt_combined(feamodel, member_label=None, combo_name="Combo 1", fpath=None):
    """Plot shear force, bending moment and deflection in a single plot.

    Parameters
    ----------
    beam : FEModel3D object
        Beam FEModel3D object
    member_label : string
        Label/name of the element/member for which the diagrams are to be drawn
    combo_name : string
        Label/name of the load combinatio nfor which the diagrams are to be drawn
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
    xs = np.linspace(0, feamodel.Members[member_label].L(), 31)
    shear_vals = []
    moment_vals = []
    displacements = []
    for x in xs:
        shear_vals.append(
            feamodel.Members[member_label].shear("Fy", x, combo_name=combo_name)
        )
        moment_vals.append(
            feamodel.Members[member_label].moment("Mz", x, combo_name=combo_name)
        )
        displacements.append(
            feamodel.Members[member_label].deflection("dy", x, combo_name=combo_name)
        )

    axs[0].plot(xs, shear_vals)
    axs[0].plot(
        [xs[0], xs[0], xs[-1], xs[-1]],
        [shear_vals[0], 0, 0, shear_vals[-1]],
        color="k",
        lw=2.0,
    )
    axs[0].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[0].set_ylabel("Shear Force, $V$ [kip]")

    axs[1].plot(xs, moment_vals)
    axs[1].plot(
        [xs[0], xs[0], xs[-1], xs[-1]],
        [moment_vals[0], 0, 0, moment_vals[-1]],
        color="k",
        lw=2.0,
    )
    axs[1].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[1].set_ylabel("Bending Moment, $M$ [kip-in]")

    axs[2].plot(xs, displacements)
    axs[2].plot(
        [xs[0], xs[0], xs[-1], xs[-1]],
        [displacements[0], 0, 0, displacements[-1]],
        color="k",
        lw=2.0,
    )
    axs[2].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[2].set_xlabel("Distance, $x$ [in]")
    axs[2].set_ylabel("Displacement, $\delta$ [in]")

    fig.savefig(fpath, format="png", dpi=150, bbox_inches="tight")


def main():
    # paths and filenames
    path_root = Path(r"../")
    path_data = path_root.joinpath("data/")
    path_plt = path_root.joinpath("plt/")
    path_report = path_root.joinpath("rpt/")
    fname_vtk = "example-02-vtk.png"
    fname_fig = "example-02.png"
    fname_rpt = "example-02_rpt.pdf"
    ffigvtk = path_plt.joinpath(fname_vtk)
    frpt = path_report.joinpath(fname_rpt)
    ffig = path_plt.joinpath(fname_fig)
    # beam parameters
    length = 14
    # in inches
    length = length * 12
    # Young's modulus in ksi
    modE = 29000
    # shear modulus in ksi
    modG = 11400
    # moment of inertia [in^4]
    Iy = 100
    Iz = 150
    J = Iy + Iz
    # cross-section area [in^2]
    area = 20

    # uniform loading [lb/ft]
    uloading = 20
    # convert to kip/in
    uloading = uloading / 100 * 1 / 12

    # new finite element model
    sbeam = FEModel3D()

    # add nodes, L = 14 ft = 168 inches
    sbeam.add_node("n1", 0, 0, 0)
    sbeam.add_node("n2", 168, 0, 0)

    # define the beam
    sbeam.add_member("m1", "n1", "n2", modE, modG, Iy, Iz, J, area)

    # supports, restrain the beam in x, y, and z displacement
    sbeam.def_support("n1", True, True, True)
    sbeam.def_support("n2", True, True, True, True)

    # loading
    sbeam.add_member_dist_load("m1", "Fy", -uloading, -uloading, 0, 168)

    # analyze
    sbeam.analyze()

    # shear moment and deflection diagrams
    plt_combined(sbeam, member_label="m1", fpath=ffig)

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
