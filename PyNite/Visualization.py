# Import the Visualization Toolkit
# You must be running a 64 bit version of Python for this to work
import vtk

## Create the visualization pipeline
# Create a data source
def RenderModel(model):

  lines = []
  mappers = []
  actors = []
  i = 0

  # Step through each member in the model
  for member in model.Members:

    # Create a new line data source
    lines.append(vtk.vtkLineSource())

    # Identify the member's i-node and j-node
    iNode = member.iNode
    jNode = member.jNode

    # Step through each node in the model and find the position of the
    # i-node and j-node
    for node in model.Nodes:

      # Check to see if the current node is the i-node
      if node.Name == iNode.Name:
        lines[i].SetPoint1(node.X, node.Y, node.Z)

      # Check to see if the current node is the j-node
      elif node.Name == jNode.Name:
        lines[i].SetPoint2(node.X, node.Y, node.Z)
    
    # Create a new mapper for the line
    # The mapper maps the data into graphics primitives
    mappers.append(vtk.vtkPolyDataMapper())
    mappers[i].SetInputConnection(lines[i].GetOutputPort())

    # Connect the mapper to an actor
    # The actor is concerned with how the model fits into the screen
    actors.append(vtk.vtkActor())
    actors[i].SetMapper(mappers[i])

    # Prepare for the next iteration/member
    i += 1

  window = vtk.vtkRenderWindow()
  # Sets the pixel width, length of the window.
  window.SetSize(500, 500)

  interactor = vtk.vtkRenderWindowInteractor()
  interactor.SetRenderWindow(window)

  renderer = vtk.vtkRenderer()
  window.AddRenderer(renderer)

  for actor in actors:
    renderer.AddActor(actor)

  # Setting the background to blue.
  renderer.SetBackground(0.1, 0.1, 0.4)

  window.Render()
  interactor.Start()