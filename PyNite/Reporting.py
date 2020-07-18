# Import libraries necessary for report printing
from jinja2 import Environment, PackageLoader
import pdfkit

# Set up the jinja2 template environment
env = Environment(
    loader=PackageLoader('PyNite', '.'),
)

# Get the report template
template = env.get_template('Report_Template.html')

def CreateReport(model, output_filepath='.//PyNite Report.pdf'):
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
    HTML = template.render(nodes=model.Nodes, members=model.Members, plates=model.Plates)

    # Convert the HTML to pdf format using PDFKit
    # Note that wkhtmltopdf must be installed on the system, and included on the system's PATH environment variable for PDFKit to work
    pdfkit.from_string(HTML, output_filepath, css='.//PyNite//MainStyleSheet.css')
    
    return
