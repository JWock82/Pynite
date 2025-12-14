import glob
import importlib
from unittest.mock import Mock

import pytest

PROJECT_ROOT = importlib.resources.files("Pynite").parent
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
    exec(open(file).read())
