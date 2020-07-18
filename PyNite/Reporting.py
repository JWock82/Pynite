from jinja2 import Environment, PackageLoader
import pdfkit

env = Environment(
    loader=PackageLoader('PyNite', '.'),
)

# Get the report template
template = env.get_template('Report_Template.html')

def CreateReport(model, output_filepath='.//PyNite Report.pdf'):
    pdfkit.from_string(template.render(nodes=model.Nodes, members=model.Members), output_filepath, css='.//PyNite//MainStyleSheet.css')
    return
