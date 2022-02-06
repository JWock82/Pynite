# PyNite Example
# Example: Simply supported beam with a point load
# Units: inches and kips, use consistent units

import numpy as np
from PyNite import FEModel3D
from PyNite import Visualization
from pathlib import Path
import matplotlib.pyplot as plt

# Updated parameters for this plot
plt.style.use("../plt/figstylefile.mplstyle")


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
    xs = np.linspace(0, beam.Members["M1"].L(), 31)
    shear_vals = []
    moment_vals = []
    displacements = []
    for x in xs:
        shear_vals.append(beam.Members["M1"].shear("Fy", x, combo_name="1.2D+1.6L"))
        moment_vals.append(beam.Members["M1"].moment("Mz", x, combo_name="1.2D+1.6L"))
        displacements.append(
            beam.Members["M1"].deflection("dy", x, combo_name="1.2D+1.6L")
        )

    axs[0].plot(xs, shear_vals)
    axs[0].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[0].set_ylabel("Shear Force, $V$ [kip]")

    axs[1].plot(xs, moment_vals)
    axs[1].grid(visible=True, which="major", axis="both", ls=":", alpha=0.5)
    axs[1].set_ylabel("Bending Moment, $M$ [kip-in]")

    axs[2].plot(xs, displacements)
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
    fname_vtk = "example-01-vtk.png"
    fname_fig = "example-01.png"
    fname_rpt = "example-01_rpt.pdf"
    ffigvtk = path_plt.joinpath(fname_vtk)
    frpt = path_report.joinpath(fname_rpt)
    ffig = path_plt.joinpath(fname_fig)

    # beam FEM model
    simple_beam = FEModel3D()

    # Beam span in [ft]
    span = 14
    # span in inch
    span = 14 * 12

    # Young's modulus, E ksi
    E = 29000
    # Shear modulus, G ksi
    G = 11400
    # Second moment of area, in^4
    Iy = 100
    Iz = 150
    J = 250
    # Cross section area of the beam, in^2
    beam_area = 20

    # add nodes (node_label, x, y, z)
    simple_beam.add_node("N1", 0, 0, 0)
    simple_beam.add_node("N2", span, 0, 0)

    # add element between the two nodes and define material
    # (element_label, nodei, nodej, E, G, Iy, Iz, J, A)
    simple_beam.add_member("M1", "N1", "N2", E, G, Iy, Iz, J, beam_area)

    # simple support at the nodes(node_label, Dx, Dy, Dz, Rx, Ry, Rz)
    # if true the support is restrained in that direction,
    # D is displacement, Dx is True, then the support is restrained against displacement in the X direction,
    # R is rotation, if Rx is True, then the support is restrained
    # against rotation about the x axis
    simple_beam.def_support("N1", True, True, True, True, False, False)
    simple_beam.def_support("N2", True, True, True, False, False, False)
    # add point load of 5 kips at mid span, (element label, loading_direction, load value, load location, load label)
    # the sign of the load value indicates the direction of loading
    # +ve in the same direction as Global Y direction
    simple_beam.add_member_pt_load("M1", "Fy", -5, span / 2, "D")
    # add live load, 8 kips of live load at mid span in the downward direction
    simple_beam.add_member_pt_load("M1", "Fy", -8, span / 2, "L")

    # add load combination
    simple_beam.add_load_combo("1.4D", {"D": 1.4})
    simple_beam.add_load_combo("1.2D+1.6L", {"D": 1.2, "L": 1.6})

    # analyze beam and check static checks
    simple_beam.analyze(check_statics=True)

    # shear moment and deflection diagrams
    plt_combined(simple_beam, ffig)

    # Visualization.render_model(
    #     simple_beam,
    #     annotation_size=10,
    #     deformed_shape=True,
    #     deformed_scale=30,
    #     render_loads=True,
    #     combo_name="1.2D+1.6L",
    #     screenshot=str(ffigvtk),
    # )

    # renderer = Visualization.Renderer(simple_beam)
    # renderer.screenshot(filepath=str(ffigvtk))

    # simple_beam.Members["M1"].plot_shear("Fy", "1.2D+1.6L")
    # simple_beam.Members["M1"].plot_moment("Mz", "1.2D+1.6L")
    # simple_beam.Members["M1"].plot_deflection("dy", "1.2D+1.6L")

    # from PyNite import Reporting

    # Reporting.create_report(
    #     simple_beam,
    #     output_filepath=str(frpt),
    #     plates=False,
    #     plate_corner_forces=False,
    #     plate_center_forces=False,
    #     plate_corner_membrane=False,
    #     plate_center_membrane=False,
    # )


if __name__ == "__main__":
    main()
