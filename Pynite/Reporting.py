# Import libraries necessary for report printing
from __future__ import annotations # Allows more recent type hints features
from typing import TYPE_CHECKING

from jinja2 import Environment, PackageLoader
import pdfkit

if TYPE_CHECKING:
    from Pynite import FEModel3D

# Determine the filepath to the local Pynite installation
from pathlib import Path
path = Path(__file__).parent

# Set up the jinja2 template environment
env = Environment(
    loader=PackageLoader('Pynite', '.'),
)

# Get the report template
template = env.get_template('Report_Template.html')

def create_report(model: FEModel3D, output_filepath:Path=path/'./Pynite Report.pdf', **kwargs):
    """Creates a pdf report for a given finite element model.

    :param model: The model to generate the report for.
    :type model: ``FEModel3D``
    :param output_filepath: The filepath to send the report to. Defaults to 'Pynite Report.pdf' in your ``PYTHONPATH``
    :type output_filepath: ``str``, optional
    :param \\**kwargs: See below for a list of valid arguments.
    :Keyword Arguments:
        * *node_table* (``bool``) -- Set to ``True`` if you want node data included in the report. Defaults to ``True``.
        * *member_table* (``bool``) -- Set to ``True`` if you want member data included in the report. Defaults to ``True``.
        * *member_releases* (``bool``) -- Set to ``True`` if you want member end release data included in the report. Defaults to ``True``.
        * *plate_table* (``bool``) -- Set to ``True if you want plate/quad data included in the report. Defaults to ``True``.
        * *node_reactions* (``bool``) -- Set to ``True`` if you want node reactions included in the report. Defaults to ``True``
        * *node_displacements* (``bool``) -- Set to ``True`` if you want node displacement results included in the report. Defaults to ``True``.
        * *member_end_forces* (``bool``) -- Set to ``True`` if you want member end force results included in the report. Defaults to ``True``.
        * *member_internal_forces* (``bool``) -- Set to ``True`` if you want member internal force results included in the report. Defaults to ``True``
        * *plate_corner_forces* (``bool``) -- Set to ``True`` if you want plate/quad corner force results (out-of-plane/bending) included in the report. Defaults to ``True``.
        * *plate_center_forces* (``bool``) -- Set to ``True`` if you want plate/quad center force results (out-of-plane/bending) included in the report. Defaults to ``True``.
        * *plate_corner_membrane* (``bool``) -- Set to ``True`` if you want plate/quad corner membrane (in-plane) force results included in the report. Defaults to ``True``.
        * *plate_center_membrane* (``bool``) -- Set to ``True`` if you want plate/quad center membrane (in-plane) force results included in the report. Defaults to ``True``.
    """

    # Create default report settings
    if 'node_table' not in kwargs: kwargs['node_table'] = True
    if 'member_table' not in kwargs: kwargs['member_table'] = True
    if 'member_releases' not in kwargs: kwargs['member_releases'] = True
    if 'plate_table' not in kwargs: kwargs['plate_table'] = True
    if 'node_reactions' not in kwargs: kwargs['node_reactions'] = True
    if 'node_displacements' not in kwargs: kwargs['node_displacements'] = True
    if 'member_end_forces' not in kwargs: kwargs['member_end_forces'] = True
    if 'member_internal_forces' not in kwargs: kwargs['member_internal_forces'] = True
    if 'plate_corner_forces' not in kwargs: kwargs['plate_corner_forces'] = True
    if 'plate_center_forces' not in kwargs: kwargs['plate_center_forces'] = True
    if 'plate_corner_membrane' not in kwargs: kwargs['plate_corner_membrane'] = True
    if 'plate_center_membrane' not in kwargs: kwargs['plate_center_membrane'] = True

    # Pass the dictionaries to the report template
    kwargs['nodes'] = model.nodes.values()
    kwargs['members'] = model.members.values()
    kwargs['plates'] = model.plates.values()
    kwargs['quads'] = model.quads.values()
    kwargs['load_combos'] = model.load_combos.values()

    # Create the report HTML using jinja2
    HTML = template.render(**kwargs)

    # Convert the HTML to pdf format using PDFKit
    # Note that wkhtmltopdf must be installed on the system, and included on the system's PATH environment variable for PDFKit to work
    pdfkit.from_string(HTML, output_filepath, css=path / './MainStyleSheet.css')
    
    return
