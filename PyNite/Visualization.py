# Import the Visualization Toolkit
# You must be running a 64 bit version of Python for this to work
import vtk

## Create the visualization pipeline
# Create a data source
def RenderModel(model):

  lines = []
  mappers = []
  actors = []

  memberLabels = []
  textActors = []
  textMappers = []

  NodeVisuals = []

  for node in model.Nodes:
    NodeVisuals.append(NodeVisual(node))

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
        Xi = node.X
        Yi = node.Y
        Zi = node.Z
        lines[i].SetPoint1(Xi, Yi, Zi)

      # Check to see if the current node is the j-node
      elif node.Name == jNode.Name:
        Xj = node.X
        Yj = node.Y
        Zj = node.Z
        lines[i].SetPoint2(Xj, Yj, Zj)
    
    # Create the text for the member label
    memberLabels.append(vtk.vtkVectorText())
    memberLabels[i].SetText(member.Name)

    # Create new mappers
    # The mapper maps the data into graphics primitives
    mappers.append(vtk.vtkPolyDataMapper())
    mappers[i].SetInputConnection(lines[i].GetOutputPort())

    textMappers.append(vtk.vtkPolyDataMapper())
    textMappers[i].SetInputConnection(memberLabels[i].GetOutputPort())

    # Connect the mappers to the actors
    # The actor is concerned with how the model fits into the screen
    actors.append(vtk.vtkActor())
    actors[i].SetMapper(mappers[i])

    textActors.append(vtk.vtkFollower())
    textActors[i].SetMapper(textMappers[i])
    textActors[i].SetScale(5, 5, 5)
    textActors[i].SetPosition((Xi+Xj)/2, (Yi+Yj)/2, (Zi+Zj)/2)

    # Prepare for the next iteration/member
    i += 1

  window = vtk.vtkRenderWindow()
  # Sets the pixel width, length of the window.
  window.SetSize(500, 500)

  # Set up the interactor
  # The interactor style determines how user interactions affect the view
  interactor = vtk.vtkRenderWindowInteractor()
  style = vtk.vtkInteractorStyleTrackballCamera() # The trackball camera style behaves a lot like most CAD programs
  interactor.SetInteractorStyle(style)
  interactor.SetRenderWindow(window)

  renderer = vtk.vtkRenderer()
  window.AddRenderer(renderer)

  # Add each member label's text actor
  for textActor in textActors:
    
    # Add the text actor
    renderer.AddActor(textActor)
    
    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    textActor.SetCamera(renderer.GetActiveCamera())
    
  for actor in actors:
    renderer.AddActor(actor)

  for visual in NodeVisuals:
    renderer.AddActor(visual.actor)
    renderer.AddActor(visual.lblActor)
    visual.lblActor.SetCamera(renderer.GetActiveCamera())

  # Setting the background to blue.
  renderer.SetBackground(0.1, 0.1, 0.4)

  # Reset the camera
  renderer.ResetCamera()

  window.Render()
  interactor.Start()

# Converts a node object into a node in the viewer
class NodeVisual():

  # Constructor
  def __init__(self, node):
    
    # Get the node's position
    self.X = node.X          # Global X coordinate
    self.Y = node.Y          # Global Y coordinate
    self.Z = node.Z          # Global Z coordinate

    # Generate a sphere for the node
    self.sphere = vtk.vtkSphereSource()
    self.sphere.SetCenter(self.X, self.Y, self.Z)
    self.sphere.SetRadius(3)

    # Set up a mapper for the node
    self.mapper = vtk.vtkPolyDataMapper()
    self.mapper.SetInputConnection(self.sphere.GetOutputPort())

    # Set up an actor for the node
    self.actor = vtk.vtkActor()
    self.actor.SetMapper(self.mapper)

    # Create the text for the node label
    self.label = vtk.vtkVectorText()
    self.label.SetText(node.Name)
    
    # Set up a mapper for the label
    self.lblMapper = vtk.vtkPolyDataMapper()
    self.lblMapper.SetInputConnection(self.label.GetOutputPort())

    # Set up an actor for the label
    self.lblActor = vtk.vtkFollower()
    self.lblActor.SetMapper(self.lblMapper)
    self.lblActor.SetScale(5, 5, 5)
    self.lblActor.SetPosition(self.X+4, self.Y+4, self.Z+4)


