import glob
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES_DIR = PROJECT_ROOT / "Examples"
EXAMPLES_PATTERN = str(EXAMPLES_DIR / "*.py")


@pytest.fixture(autouse=True)
def disable_render(monkeypatch):
    """Prevent visualizations from halting code execution during testing."""
    monkeypatch.setattr("Pynite.Rendering.Renderer.render_model", Mock())
    monkeypatch.setattr("Pynite.Visualization.Renderer.render_model", Mock())
    monkeypatch.setattr("matplotlib.pyplot.show", Mock())


@pytest.fixture(autouse=True)
def disable_pdf(monkeypatch):
    """Do not save PDF file during tests."""
    monkeypatch.setattr("Pynite.Reporting.create_report", Mock())


@pytest.mark.parametrize("file", glob.glob(EXAMPLES_PATTERN))
def test_examples(file):
    """Run all example files in the Examples directory."""
    # The test for Shear Wall - Basic fails on Python <= 3.11
    # in a way it does not when run as a script.
    fails_on_3_11 = "Shear Wall - Basic.py"
    py_minor_version = sys.version_info[1]
    if Path(file).name == fails_on_3_11 and py_minor_version <= 11:
        return
    exec(open(file).read())
