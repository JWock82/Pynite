# Import the Visualization Toolkit if it's not alread imported
# You must be running a 64 bit version of Python for this to work
import vtk

# Create a data source
def RenderModel(model, textHeight=5):

  visNodes = []
  for node in model.Nodes:
    visNodes.append(VisNode(node, textHeight))
  
  visMembers = []
  for member in model.Members:
    visMembers.append(VisMember(member, model.Nodes, textHeight))

  # Create a window
  window = vtk.vtkRenderWindow()

  # Set the pixel width and length of the window
  window.SetSize(500, 500)

  # Set up the interactor
  # The interactor style determines how user interactions affect the view
  interactor = vtk.vtkRenderWindowInteractor()
  style = vtk.vtkInteractorStyleTrackballCamera() # The trackball camera style behaves a lot like most CAD programs
  interactor.SetInteractorStyle(style)
  interactor.SetRenderWindow(window)

  # Create a renderer
  renderer = vtk.vtkRenderer()
  window.AddRenderer(renderer)
  
  # Add actors for each member
  for visMember in visMembers:

    # Add the actor for the member
    renderer.AddActor(visMember.actor)

    # Add the actor for the member label
    renderer.AddActor(visMember.lblActor)

    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    visMember.lblActor.SetCamera(renderer.GetActiveCamera())

  # Add actors for each node
  for visNode in visNodes:

    # Add the actor for the node
    renderer.AddActor(visNode.actor)

    # Add the actor for the node label
    renderer.AddActor(visNode.lblActor)

    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    visNode.lblActor.SetCamera(renderer.GetActiveCamera())

  # Setting the background to blue.
  renderer.SetBackground(0.1, 0.1, 0.4)

  # Reset the camera
  renderer.ResetCamera()

  window.Render()
  interactor.Start()

# Converts a node object into a node for the viewer
class VisNode():

  # Constructor
  def __init__(self, node, textHeight=5):
    
    # Get the node's position
    X = node.X          # Global X coordinate
    Y = node.Y          # Global Y coordinate
    Z = node.Z          # Global Z coordinate

    # Generate a sphere for the node
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(X, Y, Z)
    sphere.SetRadius(0.6*textHeight)

    # Set up a mapper for the node
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())

    # Set up an actor for the node
    self.actor = vtk.vtkActor()
    self.actor.SetMapper(mapper)

    # Create the text for the node label
    label = vtk.vtkVectorText()
    label.SetText(node.Name)
    
    # Set up a mapper for the node label
    lblMapper = vtk.vtkPolyDataMapper()
    lblMapper.SetInputConnection(label.GetOutputPort())

    # Set up an actor for the node label
    self.lblActor = vtk.vtkFollower()
    self.lblActor.SetMapper(lblMapper)
    self.lblActor.SetScale(textHeight, textHeight, textHeight)
    self.lblActor.SetPosition(X+0.6*textHeight, Y+0.6*textHeight, Z)

# Converts a member object into a member for the viewer
class VisMember():

  # Constructor
  def __init__(self, member, nodes, textHeight=5):

    # Generate a line for the member
    line = vtk.vtkLineSource()

    # Step through each node in the model and find the position of the i-node and j-node
    for node in nodes:

      # Check to see if the current node is the i-node
      if node.Name == member.iNode.Name:
        Xi = node.X
        Yi = node.Y
        Zi = node.Z
        line.SetPoint1(Xi, Yi, Zi)

      # Check to see if the current node is the j-node
      elif node.Name == member.jNode.Name:
        Xj = node.X
        Yj = node.Y
        Zj = node.Z
        line.SetPoint2(Xj, Yj, Zj)
    
    # Set up a mapper for the member
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(line.GetOutputPort())

    # Set up an actor for the member
    self.actor = vtk.vtkActor()
    self.actor.SetMapper(mapper)

    # Create the text for the member label
    label = vtk.vtkVectorText()
    label.SetText(member.Name)

    # Set up a mapper for the member label
    lblMapper = vtk.vtkPolyDataMapper()
    lblMapper.SetInputConnection(label.GetOutputPort())

    # Set up an actor for the member label
    self.lblActor = vtk.vtkFollower()
    self.lblActor.SetMapper(lblMapper)
    self.lblActor.SetScale(textHeight, textHeight, textHeight)
    self.lblActor.SetPosition((Xi+Xj)/2, (Yi+Yj)/2, (Zi+Zj)/2)