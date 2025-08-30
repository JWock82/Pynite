"""
reporting.py
------------
Generates HTML or PDF reports for Pynite FE models.

Key points:
- Uses Jinja2 to render an HTML template.
- Optionally converts that HTML into a PDF using pdfkit + wkhtmltopdf.
- Provides a single entry point: `create_report()`.
"""

# --- Standard library imports ---
import os
import platform
import shutil
import logging
from pathlib import Path
from typing import TYPE_CHECKING

# --- Third-party imports ---
from jinja2 import Environment, PackageLoader

if TYPE_CHECKING:
    from Pynite import FEModel3D


# =============================================================================
# Logging Setup
# =============================================================================

# Create a module-level logger
logger = logging.getLogger(__name__)

# Only configure logging if this file is run directly.
# When imported as a library, the caller controls logging configuration.
if __name__ == "__main__" and not logger.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s [%(name)s]: %(message)s"
    )


# =============================================================================
# Jinja2 Template Setup
# =============================================================================

# Determine the directory where this file lives
path = Path(__file__).parent

# Set up Jinja2 environment
# - PackageLoader tells Jinja2 where to find the HTML templates
env = Environment(
    loader=PackageLoader('Pynite', '.'),
)

# Load the main HTML report template
template = env.get_template('Report_Template.html')


# =============================================================================
# Helper: Find wkhtmltopdf
# =============================================================================

def get_wkhtmltopdf_path() -> str | None:
    """
    Try to locate wkhtmltopdf across platforms.
    Returns the full path if found, else None.
    """

    # 1. First check PATH (most common case)
    path = shutil.which("wkhtmltopdf")
    if path:
        return path

    # 2. Check platform-specific common install locations
    system = platform.system()

    if system == "Windows":
        candidates = [
            r"C:\Program Files\wkhtmltopdf\bin\wkhtmltopdf.exe",
            r"C:\Program Files (x86)\wkhtmltopdf\bin\wkhtmltopdf.exe",
        ]
    elif system == "Darwin":  # macOS
        candidates = [
            "/usr/local/bin/wkhtmltopdf",
            "/opt/homebrew/bin/wkhtmltopdf",
        ]
    elif system == "Linux":
        candidates = [
            "/usr/bin/wkhtmltopdf",
            "/usr/local/bin/wkhtmltopdf",
        ]
    else:
        candidates = []

    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate

    # None found
    return None


# =============================================================================
# PDFKit Setup
# =============================================================================

try:
    import pdfkit

    wkhtmltopdf_path = get_wkhtmltopdf_path()
    if wkhtmltopdf_path:
        config = pdfkit.configuration(executable=wkhtmltopdf_path)
        logger.info(f"wkhtmltopdf found at: {wkhtmltopdf_path}")
    else:
        raise FileNotFoundError("wkhtmltopdf not found on system. Install wkhtmltopdf.")
except ImportError:
    logger.error("pdfkit not installed. Run: pip install pdfkit")
    pdfkit = None
    config = None
except Exception as e:
    logger.error(f"PDFKit setup failed: {e}")
    pdfkit = None
    config = None


# =============================================================================
# Main Report Function
# =============================================================================

def create_report(model: FEModel3D,
                  output_filepath: Path | str = path / 'Pynite Report.pdf',
                  format: str = 'pdf',
                  **kwargs) -> None:
    """
    Creates a report for a given finite element model.

    :param model: The model to generate the report for.
    :param output_filepath: Filepath for the output file.
    :param format: Output format: 'pdf' or 'html'.
    :param kwargs: Report options (see below).

    Keyword Arguments:
        * node_table (bool): Include node data (default: True).
        * member_table (bool): Include member data (default: True).
        * member_releases (bool): Include member end releases (default: True).
        * plate_table (bool): Include plate/quad data (default: True).
        * node_reactions (bool): Include node reactions (default: True).
        * node_displacements (bool): Include node displacements (default: True).
        * member_end_forces (bool): Include member end forces (default: True).
        * member_internal_forces (bool): Include member internal forces (default: True).
        * plate_corner_forces (bool): Include plate corner out-of-plane forces (default: True).
        * plate_center_forces (bool): Include plate center out-of-plane forces (default: True).
        * plate_corner_membrane (bool): Include plate corner in-plane forces (default: True).
        * plate_center_membrane (bool): Include plate center in-plane forces (default: True).
    """

    # -------------------------------------------------------------------------
    # Default Settings
    # -------------------------------------------------------------------------
    defaults = {
        'node_table': True,
        'member_table': True,
        'member_releases': True,
        'plate_table': True,
        'node_reactions': True,
        'node_displacements': True,
        'member_end_forces': True,
        'member_internal_forces': True,
        'plate_corner_forces': True,
        'plate_center_forces': True,
        'plate_corner_membrane': True,
        'plate_center_membrane': True,
    }

    # Fill missing kwargs with defaults
    for key, val in defaults.items():
        kwargs.setdefault(key, val)

    # -------------------------------------------------------------------------
    # Add Model Data to Context
    # -------------------------------------------------------------------------
    kwargs['nodes'] = model.nodes.values()
    kwargs['members'] = model.members.values()
    kwargs['plates'] = model.plates.values()
    kwargs['quads'] = model.quads.values()
    kwargs['load_combos'] = model.load_combos.values()

    # -------------------------------------------------------------------------
    # Render HTML from Template
    # -------------------------------------------------------------------------
    try:
        HTML = template.render(**kwargs)
    except Exception as e:
        logger.error(f"Template rendering failed: {e}")
        raise RuntimeError("Failed to render report template.") from e

    # -------------------------------------------------------------------------
    # Export Report
    # -------------------------------------------------------------------------
    try:
        if format == 'pdf':
            if pdfkit is None or config is None:
                raise RuntimeError("pdfkit/wkhtmltopdf not available. Cannot create PDF report.")

            pdfkit.from_string(
                HTML,
                str(output_filepath),
                css=str(path / 'MainStyleSheet.css'),
                configuration=config
            )
            logger.info(f"PDF report generated at: {output_filepath}")

        elif format == 'html':
            with open(output_filepath, "w", encoding="utf-8") as file:
                file.write(HTML)
            logger.info(f"HTML report generated at: {output_filepath}")

        else:
            raise ValueError("Invalid format. Use 'pdf' or 'html'.")

    except Exception as e:
        logger.error(f"Report creation failed: {e}")
        raise


# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    # Example: Run a demo if this file is executed directly
    logger.info("This module is intended to be imported, not run directly.")
