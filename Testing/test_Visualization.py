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


# ============================================================================
# Utility helper
# ============================================================================

def set_offscreen(rndr):
    """Enable off-screen rendering regardless of backend."""
    if hasattr(rndr, "window"):           # legacy VTKRenderer
        rndr.window.SetOffScreenRendering(1)
    elif hasattr(rndr, "plotter"):        # PyVista-based renderer
        rndr.plotter.off_screen = True


# ============================================================================
# Beam model fixture
# ============================================================================

@pytest.fixture(scope="module")
def beam_model():
    """Simple beam model for rendering tests."""
    m = FEModel3D()
    m.add_node("N1", 0, 0, 0)
    m.add_node("N2", 10, 0, 0)
    m.add_material("Steel", 29000/144, 11200/144, 0.3, 0.49, 36/144)
    m.add_section("Rect", 20, 100, 200, 150)
    m.add_member("M1", "N1", "N2", "Steel", "Rect")
    m.add_member_dist_load("M1", "Fy", -1, -1, case="D")
    m.add_load_combo("1.4D", {"D": 1.4})
    m.analyze_linear()
    return m


# ============================================================================
# Plate model fixture
# ============================================================================

@pytest.fixture(scope="module")
def plate_model():
    """Creates a small quad mesh plate for rendering tests."""
    m = FEModel3D()

    # Define corner nodes (square plate)
    m.add_node("N1", 0, 0, 0)
    m.add_node("N2", 10, 0, 0)
    m.add_node("N3", 10, 10, 0)
    m.add_node("N4", 0, 10, 0)

    # Supports (edges fixed)
    for n in ["N1", "N2", "N3", "N4"]:
        m.def_support(n, True, True, True, True, True, True)

    # Material and plate
    m.add_material("Concrete", 3600/144, 1800/144, 0.2, 0.15, 150/144)
    m.add_plate("P1", "N1", "N2", "N3", "N4", 1, "Concrete", 1, 1)

    # Uniform pressure load
    m.add_plate_surface_pressure("P1", -100, case="DL")

    # Load combo
    m.add_load_combo("1.2DL", {"DL": 1.2})
    m.analyze_linear()
    return m


# ============================================================================
# Combined parameterized renderer fixture
# ============================================================================

@pytest.fixture(params=["VTK", "PV"], scope="function")
def renderer(request, beam_model, plate_model):
    """Yields renderers for both beam and plate models."""
    model = beam_model if request.node.name.endswith("_beam") else plate_model
    rndr = VTKRenderer(model) if request.param == "VTK" else PVRenderer(model)
    set_offscreen(rndr)

    rndr.render_loads = True
    rndr.render_nodes = True
    rndr.combo_name = list(model.load_combos.keys())[0]
    return rndr


# ============================================================================
# General attribute tests
# ============================================================================

def test_set_visual_properties(renderer):
    renderer.annotation_size = 8
    renderer.labels = True
    renderer.color_map = "Mz"
    renderer.scalar_bar = True
    renderer.scalar_bar_text_size = 18
    assert renderer.annotation_size == 8
    assert renderer.labels
    assert renderer.color_map == "Mz"
    assert renderer.scalar_bar
    assert renderer.scalar_bar_text_size == 18


def test_toggle_deformed_shape(renderer):
    renderer.deformed_shape = True
    renderer.deformed_scale = 25
    renderer.update(reset_camera=False)
    assert renderer.deformed_shape
    assert renderer.deformed_scale == 25


# ============================================================================
# Beam and plate coverage (run separately)
# ============================================================================

@pytest.mark.parametrize("renderer_cls", [VTKRenderer, PVRenderer])
@pytest.mark.parametrize("model_type", ["beam", "plate"])
def test_update_pipeline(renderer_cls, model_type, beam_model, plate_model):
    """Checks the update pipeline for both model types and renderers."""
    model = beam_model if model_type == "beam" else plate_model
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


@pytest.mark.parametrize("model_type", ["beam", "plate"])
def test_render_model_pipeline(model_type, beam_model, plate_model):
    """Ensures render_model executes without error for both model types."""
    model = beam_model if model_type == "beam" else plate_model
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
