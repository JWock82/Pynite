# Import libraries necessary for report printing
from jinja2 import Environment, PackageLoader
import pdfkit
import warnings

# Determine the filepath to the local PyNite installation
from pathlib import Path
path = Path(__file__).parent

# Set up the jinja2 template environment
env = Environment(
    loader=PackageLoader('PyNite', '.'),
)

# Get the report template
template = env.get_template('Report_Template.html')

def CreateReport(model, output_filepath=path / './PyNite Report.pdf', \
                 nodes=True, members=True, plates=True, member_releases=True, \
                 node_reactions=True, node_displacements=True, member_end_forces=True, member_internal_forces=True, \
                 plate_corner_forces=True, plate_center_forces=True, plate_corner_membrane=True, plate_center_membrane=True):
                 warnings.warn('`CreateReport` will be replaced with `create_report` in a future version of PyNite.', FutureWarning)
                 create_report(model, output_filepath, nodes, members, plates, member_releases, \
                               node_reactions, node_displacements, member_end_forces, member_internal_forces, \
                               plate_corner_forces, plate_center_forces, plate_corner_membrane, plate_center_membrane)

def create_report(model, output_filepath=path / './PyNite Report.pdf', \
                  nodes=True, members=True, plates=True, member_releases=True, \
                  node_reactions=True, node_displacements=True, member_end_forces=True, member_internal_forces=True, \
                  plate_corner_forces=True, plate_center_forces=True, plate_corner_membrane=True, plate_center_membrane=True):
    '''
    Creates a pdf report for a given finite element model.

    Parameters
    ----------
    model : FEModel3D
        The finite element model used to generate the report.
    output_filepath : string
        The filepath with filename to save the pdf report to.
    '''
    
    # Create the report HTML using jinja2
    HTML = template.render(nodes=model.Nodes.values(), members=model.Members.values(), plates=model.Plates.values(), quads=model.Quads.values(), load_combos=model.LoadCombos.keys(), \
                           node_table=nodes, member_table=members, plate_table=plates, member_releases=member_releases, \
                           node_reactions=node_reactions, node_displacements=node_displacements, member_end_forces=member_end_forces, member_internal_forces=member_internal_forces, \
                           plate_corner_forces=plate_corner_forces, plate_center_forces=plate_center_forces, plate_corner_membrane=plate_corner_membrane, plate_center_membrane=plate_center_membrane)

    # Convert the HTML to pdf format using PDFKit
    # Note that wkhtmltopdf must be installed on the system, and included on the system's PATH environment variable for PDFKit to work
    pdfkit.from_string(HTML, output_filepath, css=path / './MainStyleSheet.css')
    
    return
