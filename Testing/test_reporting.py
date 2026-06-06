# -*- coding: utf-8 -*-
"""
Regression tests for the HTML reporting module (Pynite.Reporting).

Covers issue #268: generating a report for a model containing quads/plates
raised ``KeyError: 'Combo 1'`` because the report template passed the load
combination name into the elements' ``local`` argument instead of
``combo_name`` (leaving ``combo_name`` at its default ``'Combo 1'``, which was
never analyzed). A second, related template bug indexed the ``(3, 1)`` membrane
result as ``[n]`` instead of ``[n][0]``.
"""

import os
import tempfile
import unittest

from Pynite import FEModel3D

try:
    import jinja2  # noqa: F401

    from Pynite import Reporting

    HAS_JINJA = True
except ImportError:
    HAS_JINJA = False


def _build_model() -> FEModel3D:
    """A supported quad + plate under pressure, analyzed for a single custom
    load combination so the default ``'Combo 1'`` is never solved."""
    model = FEModel3D()
    model.add_material("Steel", 29000, 11200, 0.3, 0.284)

    # Quad nodes
    model.add_node("N1", 0, 0, 0)
    model.add_node("N2", 10, 0, 0)
    model.add_node("N3", 10, 10, 0)
    model.add_node("N4", 0, 10, 0)
    model.add_quad("Q1", "N1", "N2", "N3", "N4", 0.25, "Steel")

    # Plate nodes (offset so it is a separate element)
    model.add_node("N5", 10, 0, 0)
    model.add_node("N6", 20, 0, 0)
    model.add_node("N7", 20, 10, 0)
    model.add_node("N8", 10, 10, 0)
    model.add_plate("P1", "N5", "N6", "N7", "N8", 0.25, "Steel")

    for n in ("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"):
        model.def_support(n, True, True, True, True, True, True)

    model.add_quad_surface_pressure("Q1", 0.1, case="D")
    model.add_plate_surface_pressure("P1", 0.1, case="D")

    # Only a custom-named combo is analyzed; the default 'Combo 1' is not.
    model.add_load_combo("1.4D", {"D": 1.4})
    model.analyze_linear(check_statics=False)
    return model


@unittest.skipUnless(HAS_JINJA, "reporting extra (Jinja2) not installed")
class TestReporting(unittest.TestCase):
    def test_report_with_quads_and_plates(self) -> None:
        """create_report must not raise KeyError 'Combo 1' (issue #268)."""
        model = _build_model()
        out = os.path.join(tempfile.mkdtemp(), "report.html")

        # Must not raise (previously raised KeyError: 'Combo 1').
        Reporting.create_report(model, format="html", output_filepath=out)

        self.assertTrue(os.path.exists(out))
        html = open(out, encoding="utf-8").read()
        # The analyzed combo is reported, and the unanalyzed default is absent.
        self.assertIn("1.4D", html)
        self.assertNotIn("Combo 1", html)


if __name__ == "__main__":
    unittest.main()
