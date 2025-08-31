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
from pathlib import Path
from typing import TYPE_CHECKING

# --- Third-party imports ---
from jinja2 import Environment, PackageLoader

# --- Pynite imports ---
from Pynite import FEModel3D

# --- Type checking ---
if TYPE_CHECKING:
    from Pynite import FEModel3D

# %%
# Jinja2 Template Setup

# Determine the directory where this file lives
path = Path(__file__).parent

# Set up Jinja2 environment
# - PackageLoader tells Jinja2 where to find the HTML templates
env = Environment(
    loader=PackageLoader('Pynite', '.'),
)

# Load the main HTML report template
template = env.get_template('Report_Template.html')

# %%
def create_report(model: FEModel3D,
                  output_filepath: Path | str = path / 'Pynite Report.pdf',
                  format: str = 'pdf',
                  log: bool = True,
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

    # Default Settings
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

    # Add model data to context
    kwargs['nodes'] = model.nodes.values()
    kwargs['members'] = model.members.values()
    kwargs['plates'] = model.plates.values()
    kwargs['quads'] = model.quads.values()
    kwargs['load_combos'] = model.load_combos.values()

    # Render HTML from Template
    HTML = template.render(**kwargs)

    # Check if a PDF file has been requested
    if format.upper() == 'PDF':

        # Import PDFKit
        try:
            import pdfkit
        except:
            raise ImportError('PDFKit is not installed. Install it using `pip install pdfkit`.')

        # Attempt to get the wkhtmltopdf path
        wkhtmltopdf_path = get_wkhtmltopdf_path()

        # Notify the user if wkthmltopdf is not found
        if wkhtmltopdf_path is None:
            raise Exception('Unable to locate wkhtmltopdf. If it is already installed add it to your system\'s PATH environmental variable so Pynite can find it.')

        # Set up a PDFKit configuration that uses the correct path to wkhtmltopdf
        config = pdfkit.configuration(wkhtmltopdf=wkhtmltopdf_path)

        # Generate the PDF report
        pdfkit.from_string(
            HTML,
            str(output_filepath),
            css=str(path / 'MainStyleSheet.css'),
            configuration=config
        )

        # Notify the user where the report was saved
        if log: print(f"- PDF report generated at: {output_filepath}")

    elif format.upper() == 'HTML':

        # Write directly to an HTML file
        with open(output_filepath, "w", encoding="utf-8") as file:
            file.write(HTML)

        # Notify the user of the report's location
        if log: print(f"- HTML report generated at: {output_filepath}")

    else:
        raise ValueError("Invalid format. Use 'pdf' or 'html'.")


def get_wkhtmltopdf_path(log: bool = True) -> str | None:
    """
    Tries to locate wkhtmltopdf (across any platform).
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
