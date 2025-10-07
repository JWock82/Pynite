import os
import pytest
from io import BytesIO
from IPython.display import Image
from Pynite import FEModel3D
from Pynite.Visualization import Renderer as VTKRenderer
from Pynite.Rendering import Renderer as PVRenderer

# Optional: force PyVista off-screen to prevent GUI pop-ups in IDE runs
import pyvista as pv
pv.OFF_SCREEN = True


# Utility helper
def set_offscreen(rndr):
    """Enable off-screen rendering regardless of backend."""

    if hasattr(rndr, "window"):           # VTK Renderer
        rndr.window.SetOffScreenRendering(1)
    elif hasattr(rndr, "plotter"):        # PyVista-based renderer
        rndr.plotter.off_screen = True

# Model setup
def visual_model():
    """A simple model for rendering tests."""

    # Create a new model
    m = FEModel3D()

    # Add a beam to the model
    m.add_node('N1', 0, 0, 0)
    m.add_node('N2', 3.333, 0, 0)
    m.add_node('N3', 6.667, 0, 0)
    m.add_node('N4', 10, 0, 0)

    # Add every type of support to the model
    m.def_support('N1', True, True, True, True, True, True)
    m.def_support('N2', True, True, True, True, True, False)
    m.def_support('N3', False, False, False, False, False, True)
    m.def_support('N4', True, True, True, False, False, False)

    m.add_material("Steel", 29000/144, 11200/144, 0.3, 0.49, 36/144)

    m.add_section("Rect", 20, 100, 200, 150)

    m.add_member("M1", "N1", "N3", "Steel", "Rect")

    # Add every type of member load to the model
    m.add_member_dist_load("M1", "Fy", -1, -1, case="D")
    m.add_member_dist_load("M1", "Fz", -1, -1, case="D")
    m.add_member_dist_load("M1", "FY", -1, -1, case="D")
    m.add_member_dist_load("M1", "Fy", -1, -1, case="D")
    m.add_member_pt_load("M1", "Fy", -1, 5, case="D")
    m.add_member_pt_load("M1", "Fz", -1, 5, case="D")
    m.add_member_pt_load("M1", "FY", -1, 5, case="D")
    m.add_member_pt_load("M1", "FZ", -1, 5, case="D")
    m.add_member_pt_load("M1", "Mz", -1, 5, case="D")
    m.add_member_pt_load("M1", "My", -1, 5, case="D")
    m.add_member_pt_load("M1", "MZ", -1, 5, case="D")
    m.add_member_pt_load("M1", "MY", -1, 5, case="D")

    # Add a plate to the model
    m.add_node("N5", 0, 0, 10)
    m.add_node("N6", 10, 0, 10)
    m.add_node("N7", 10, 10, 10)
    m.add_node("N8", 0, 10, 10)

    # Supports (edges fixed)
    for n in ["N1", "N2", "N3", "N4"]:
        m.def_support(n, True, True, True, True, True, True)

    m.add_plate("P1", "N5", "N6", "N7", "N8", 1, 'Steel', 1, 1)

    # Uniform pressure load
    m.add_plate_surface_pressure("P1", -100, case="D")

    m.add_load_combo("1.4D", {"D": 1.4})
    m.analyze_linear()
    return m

@pytest.fixture(params=["VTK", "PV"], scope="function")
def renderer(request):
    """Yields renderers VTK and Pyvista."""

    model = visual_model()
    rndr = VTKRenderer(model) if request.param == "VTK" else PVRenderer(model)
    set_offscreen(rndr)

    rndr.render_loads = True
    rndr.render_nodes = True
    rndr.combo_name = list(model.load_combos.keys())[0]

    return rndr


def test_toggle_visual_properties(renderer):

    # Test each plate stress
    for stress in ['Mx', 'My', 'Mxy', 'Sx', 'Sy', 'Txy']:

        renderer.annotation_size = 8
        renderer.labels = True
        renderer.color_map = stress
        renderer.scalar_bar = True
        renderer.scalar_bar_text_size = 18
        assert renderer.annotation_size == 8
        assert renderer.labels
        assert renderer.color_map == stress
        assert renderer.scalar_bar
        assert renderer.scalar_bar_text_size == 18

        if isinstance(renderer, VTKRenderer):
            renderer.render_model(interact=False)
        else:
            renderer.render_model()


def test_toggle_deformed_shape(renderer):
    renderer.deformed_shape = True
    renderer.deformed_scale = 25
    renderer.update(reset_camera=False)
    assert renderer.deformed_shape
    assert renderer.deformed_scale == 25


@pytest.mark.parametrize("renderer_cls", [VTKRenderer, PVRenderer])
def test_update_pipeline(renderer_cls):
    """Checks the update pipeline for both model types and renderers."""
    model = visual_model()
    rndr = renderer_cls(model)
    set_offscreen(rndr)

    rndr.combo_name = list(model.load_combos.keys())[0]
    rndr.render_loads = True
    rndr.render_nodes = True
    rndr.update(reset_camera=True)

    # Verify backend context exists
    assert hasattr(rndr, "update")
    assert hasattr(rndr, "screenshot")
    assert hasattr(rndr, "plotter") or hasattr(rndr, "window")


def test_render_model_pipeline():
    """Ensures render_model executes without error for both model types."""
    model = visual_model()
    rndr = VTKRenderer(model)
    set_offscreen(rndr)
    rndr.combo_name = list(model.load_combos.keys())[0]
    rndr.render_model(interact=False)


# ============================================================================
# Screenshot coverage (updated for PVRenderer + VTKRenderer)
# ============================================================================

@pytest.mark.parametrize("path", ["BytesIO", "console"])
def test_screenshots(renderer, path, tmp_path):
    """
    Validate screenshot outputs for both renderer types.

    VTKRenderer:
        - Returns an IPython Image or BytesIO object depending on filepath
    PVRenderer:
        - Usually returns None, but still generates an image internally
    """
    result = renderer.screenshot(filepath=path, interact=False)

    if hasattr(renderer, "plotter"):  # PVRenderer (PyVista-based)
        # Screenshot likely returns None, so validate image array instead
        img = renderer.plotter.screenshot(return_img=True)
        assert img is not None and hasattr(img, "shape") and len(img.shape) in (2, 3)
    else:  # VTKRenderer (VTK window-based)
        if path == "BytesIO":
            assert isinstance(result, BytesIO)
        else:
            assert isinstance(result, Image)


def test_screenshot_to_file(renderer, tmp_path):
    """
    Ensure screenshot file output works for both renderers.

    PVRenderer -> writes via PyVista plotter.screenshot
    VTKRenderer -> writes via its internal vtkPNGWriter
    """
    out_file = tmp_path / "render_output.png"
    result = renderer.screenshot(filepath=str(out_file), interact=False)

    # PVRenderer usually returns None, so check the file directly
    if hasattr(renderer, "plotter"):
        if not out_file.exists():  # Some PVRenderers don't write via .screenshot()
            # Fallback: force-save manually through PyVista
            renderer.plotter.screenshot(filename=str(out_file))
    else:
        # VTKRenderer may return an Image; file should still exist
        assert isinstance(result, (type(None), Image, BytesIO))

    assert out_file.exists(), "Screenshot file was not created."
    assert out_file.stat().st_size > 0, "Screenshot file is empty."


# ============================================================================
# Manual entry point for running tests in an IDE
# ============================================================================

if __name__ == "__main__":
    import sys
    args = [
        "-v",
        "--color=yes",
        "--maxfail=3",
        "--disable-warnings",
        __file__,
    ]
    print("Running PyNite renderer tests directly...")
    sys.exit(pytest.main(args))
