# Import the Visualization Toolkit if it's not alread imported
# You must be running a 64 bit version of Python for this to work
import vtk
from numpy import array, empty, append
from math import isclose

# Create a data source
def RenderModel(model, textHeight=5):

  visNodes = []
  for node in model.Nodes:
    visNodes.append(VisNode(node, textHeight))

  visAuxNodes = []
  for auxnode in model.auxNodes:
    visAuxNodes.append(VisNode(auxnode, textHeight, color='red'))

  visMembers = []
  for member in model.Members:
    visMembers.append(VisMember(member, model.Nodes, textHeight))

  visPlates = []
  for plate in model.Plates:
    visPlates.append(VisPlate(plate, model.Nodes, textHeight))

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

  # Add actors for each plate
  for visPlate in visPlates:

    # Add the actors for the plate
    renderer.AddActor(visPlate.actor)

    # Add the actor for the plate label
    renderer.AddActor(visPlate.lblActor)

    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    visPlate.lblActor.SetCamera(renderer.GetActiveCamera())

  # Add actors for each auxnode
  for visAuxNode in visAuxNodes:

    # Add the actor for the node
    renderer.AddActor(visAuxNode.actor)

    # Add the actor for the node label
    renderer.AddActor(visAuxNode.lblActor)

    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    visAuxNode.lblActor.SetCamera(renderer.GetActiveCamera())

  # Setting the background to blue.
  renderer.SetBackground(0.1, 0.1, 0.4)

  # Reset the camera
  renderer.ResetCamera()

  window.Render()
  interactor.Start()

# Converts a node object into a node for the viewer
class VisNode():

  # Constructor
  def __init__(self, node, textHeight=5, color=None):
    
    # Create a 'vtkPolyData' object to represent the node
    polyData = vtk.vtkPolyData()

    # Create a 'vtkAppendPolyData' filter to append the node and it's supports together into a single dataset
    appendFilter = vtk.vtkAppendPolyData()

    # Get the node's position
    X = node.X # Global X coordinate
    Y = node.Y # Global Y coordinate
    Z = node.Z # Global Z coordinate

    # Generate a sphere for the node
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(X, Y, Z)
    sphere.SetRadius(0.6*textHeight)
    sphere.Update()

    polyData.ShallowCopy(sphere.GetOutput())
    appendFilter.AddInputData(polyData)
    
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
    self.lblActor.SetPosition(X + 0.6*textHeight, Y + 0.6*textHeight, Z)

    # Generate any supports that occur at the node
    # Check for a fixed suppport
    if node.SupportDX == True and node.SupportDY == True and node.SupportDZ == True \
      and node.SupportRX == True and node.SupportRY == True and node.SupportRZ == True:

      # Create the fixed support
      support = vtk.vtkCubeSource()
      support.SetCenter(node.X, node.Y, node.Z)
      support.SetXLength(textHeight*1.2)
      support.SetYLength(textHeight*1.2)
      support.SetZLength(textHeight*1.2)

      # Copy and append the support data to the append filter
      support.Update()
      polyData = vtk.vtkPolyData()
      polyData.ShallowCopy(support.GetOutput())
      appendFilter.AddInputData(polyData)
    
    # Check for a pinned support
    elif node.SupportDX == True and node.SupportDY == True and node.SupportDZ == True \
      and node.SupportRX == False and node.SupportRY == False and node.SupportRZ == False:
      
      # Create the pinned support
      support = vtk.vtkConeSource()
      support.SetCenter(node.X, node.Y-0.6*textHeight, node.Z)
      support.SetDirection((0, 1, 0))
      support.SetHeight(textHeight*1.2)
      support.SetRadius(textHeight*1.2)

      # Copy and append the support data to the append filter
      support.Update()
      polyData = vtk.vtkPolyData()
      polyData.ShallowCopy(support.GetOutput())
      appendFilter.AddInputData(polyData)
    
    # Other support conditions
    else:

      # Restrained against X translation
      if node.SupportDX == True:
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X-textHeight, node.Y, node.Z)
        support1.SetPoint2(node.X+textHeight, node.Y, node.Z)

        # Copy and append the support data to the append filter
        support1.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support1.GetOutput())
        appendFilter.AddInputData(polyData)

        support2 = vtk.vtkConeSource()
        support2.SetCenter(node.X-textHeight, node.Y, node.Z)
        support2.SetDirection((1, 0, 0))
        support2.SetHeight(textHeight*0.6)
        support2.SetRadius(textHeight*0.3)

        # Copy and append the support data to the append filter
        support2.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support2.GetOutput())
        appendFilter.AddInputData(polyData)

        support3 = vtk.vtkConeSource()
        support3.SetCenter(node.X+textHeight, node.Y, node.Z)
        support3.SetDirection((-1, 0, 0))
        support3.SetHeight(textHeight*0.6)
        support3.SetRadius(textHeight*0.3)

        # Copy and append the support data to the append filter
        support3.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support3.GetOutput())
        appendFilter.AddInputData(polyData)
      
      # Restrained against Y translation
      if node.SupportDY == True:
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y-textHeight, node.Z)
        support1.SetPoint2(node.X, node.Y+textHeight, node.Z)

        # Copy and append the support data to the append filter
        support1.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support1.GetOutput())
        appendFilter.AddInputData(polyData)

        support2 = vtk.vtkConeSource()
        support2.SetCenter(node.X, node.Y-textHeight, node.Z)
        support2.SetDirection((0, 1, 0))
        support2.SetHeight(textHeight*0.6)
        support2.SetRadius(textHeight*0.3)

        # Copy and append the support data to the append filter
        support2.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support2.GetOutput())
        appendFilter.AddInputData(polyData)

        support3 = vtk.vtkConeSource()
        support3.SetCenter(node.X, node.Y+textHeight, node.Z)
        support3.SetDirection((0, -1, 0))
        support3.SetHeight(textHeight*0.6)
        support3.SetRadius(textHeight*0.3)

        # Copy and append the support data to the append filter
        support3.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support3.GetOutput())
        appendFilter.AddInputData(polyData)
      
      # Restrained against Z translation
      if node.SupportDZ == True:
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y, node.Z-textHeight)
        support1.SetPoint2(node.X, node.Y, node.Z+textHeight)

        # Copy and append the support data to the append filter
        support1.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support1.GetOutput())
        appendFilter.AddInputData(polyData)

        support2 = vtk.vtkConeSource()
        support2.SetCenter(node.X, node.Y, node.Z-textHeight)
        support2.SetDirection((0, 0, 1))
        support2.SetHeight(textHeight*0.6)
        support2.SetRadius(textHeight*0.3)

        # Copy and append the support data to the append filter
        support2.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support2.GetOutput())
        appendFilter.AddInputData(polyData)

        support3 = vtk.vtkConeSource()
        support3.SetCenter(node.X, node.Y, node.Z+textHeight)
        support3.SetDirection((0, 0, -1))
        support3.SetHeight(textHeight*0.6)
        support3.SetRadius(textHeight*0.3)

        # Copy and append the support data to the append filter
        support3.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support3.GetOutput())
        appendFilter.AddInputData(polyData)

      # Restrained against rotation about the X-axis
      if node.SupportRX == True:
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X-1.6*textHeight, node.Y, node.Z)
        support1.SetPoint2(node.X+1.6*textHeight, node.Y, node.Z)

        # Copy and append the support data to the append filter
        support1.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support1.GetOutput())
        appendFilter.AddInputData(polyData)

        support2 = vtk.vtkCubeSource()
        support2.SetCenter(node.X-1.9*textHeight, node.Y, node.Z)
        support2.SetXLength(textHeight*0.6)
        support2.SetYLength(textHeight*0.6)
        support2.SetZLength(textHeight*0.6)

        # Copy and append the support data to the append filter
        support2.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support2.GetOutput())
        appendFilter.AddInputData(polyData)

        support3 = vtk.vtkCubeSource()
        support3.SetCenter(node.X+1.9*textHeight, node.Y, node.Z)
        support3.SetXLength(textHeight*0.6)
        support3.SetYLength(textHeight*0.6)
        support3.SetZLength(textHeight*0.6)

        # Copy and append the support data to the append filter
        support3.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support3.GetOutput())
        appendFilter.AddInputData(polyData)
      
      # Restrained against rotation about the Y-axis
      if node.SupportRY == True:
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y-1.6*textHeight, node.Z)
        support1.SetPoint2(node.X, node.Y+1.6*textHeight, node.Z)

        # Copy and append the support data to the append filter
        support1.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support1.GetOutput())
        appendFilter.AddInputData(polyData)

        support2 = vtk.vtkCubeSource()
        support2.SetCenter(node.X, node.Y-1.9*textHeight, node.Z)
        support2.SetXLength(textHeight*0.6)
        support2.SetYLength(textHeight*0.6)
        support2.SetZLength(textHeight*0.6)

        # Copy and append the support data to the append filter
        support2.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support2.GetOutput())
        appendFilter.AddInputData(polyData)

        support3 = vtk.vtkCubeSource()
        support3.SetCenter(node.X, node.Y+1.9*textHeight, node.Z)
        support3.SetXLength(textHeight*0.6)
        support3.SetYLength(textHeight*0.6)
        support3.SetZLength(textHeight*0.6)

        # Copy and append the support data to the append filter
        support3.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support3.GetOutput())
        appendFilter.AddInputData(polyData)
      
      # Restrained against rotation about the Z-axis
      if node.SupportRZ == True:
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y, node.Z-1.6*textHeight)
        support1.SetPoint2(node.X, node.Y, node.Z+1.6*textHeight)

        # Copy and append the support data to the append filter
        support1.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support1.GetOutput())
        appendFilter.AddInputData(polyData)

        support2 = vtk.vtkCubeSource()
        support2.SetCenter(node.X, node.Y, node.Z-1.9*textHeight)
        support2.SetXLength(textHeight*0.6)
        support2.SetYLength(textHeight*0.6)
        support2.SetZLength(textHeight*0.6)

        # Copy and append the support data to the append filter
        support2.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support2.GetOutput())
        appendFilter.AddInputData(polyData)

        support3 = vtk.vtkCubeSource()
        support3.SetCenter(node.X, node.Y, node.Z+1.9*textHeight)
        support3.SetXLength(textHeight*0.6)
        support3.SetYLength(textHeight*0.6)
        support3.SetZLength(textHeight*0.6)

        # Copy and append the support data to the append filter
        support3.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(support3.GetOutput())
        appendFilter.AddInputData(polyData)
    
    # Update the append filter
    appendFilter.Update()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(appendFilter.GetOutputPort())
    self.actor = vtk.vtkActor()

    # Add color to the actors
    if color == 'red':
      self.actor.GetProperty().SetColor(255, 0, 0) # Red
      self.lblActor.GetProperty().SetColor(255, 0, 0) # Red
    elif color == 'yellow':
      self.actor.GetProperty().SetColor(255, 255, 0) # Yellow
      self.lblActor.GetProperty().SetColor(255, 255, 0) # Yellow
    
    # Set the mapper for the node's actor
    self.actor.SetMapper(mapper)
               
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

# Converts a plate object into a plate for the viewer
class VisPlate():

  # Constructor
  def __init__(self, plate, nodes, textHeight=5):

    # Generate lines for the plate
    line1 = vtk.vtkLineSource()
    line2 = vtk.vtkLineSource()
    line3 = vtk.vtkLineSource()
    line4 = vtk.vtkLineSource()

    # Step through each node in the model and find the position of the i-node and j-node
    for node in nodes:

      # Check to see if the current node is the i-node
      if node.Name == plate.iNode.Name:
        Xi = node.X
        Yi = node.Y
        Zi = node.Z
        line1.SetPoint1(Xi, Yi, Zi)
        line4.SetPoint2(Xi, Yi, Zi)

      # Check to see if the current node is the j-node
      elif node.Name == plate.jNode.Name:
        Xj = node.X
        Yj = node.Y
        Zj = node.Z
        line1.SetPoint2(Xj, Yj, Zj)
        line2.SetPoint1(Xj, Yj, Zj)

      # Check to see if the current node is the m-node
      elif node.Name == plate.mNode.Name:
        Xm = node.X
        Ym = node.Y
        Zm = node.Z
        line2.SetPoint2(Xm, Ym, Zm)
        line3.SetPoint1(Xm, Ym, Zm)
      
      # Check to see if the current node is the n-node
      elif node.Name == plate.nNode.Name:
        Xn = node.X
        Yn = node.Y
        Zn = node.Z
        line3.SetPoint2(Xn, Yn, Zn)
        line4.SetPoint1(Xn, Yn, Zn)

    # Create a 'vtkAppendPolyData' filter to add all the lines to
    appendFilter = vtk.vtkAppendPolyData()

    line1.Update()
    polyData1 = vtk.vtkPolyData()
    polyData1.ShallowCopy(line1.GetOutput())
    appendFilter.AddInputData(polyData1)

    line2.Update()
    polyData2 = vtk.vtkPolyData()
    polyData2.ShallowCopy(line2.GetOutput())
    appendFilter.AddInputData(polyData2)

    line3.Update()
    polyData3 = vtk.vtkPolyData()
    polyData3.ShallowCopy(line3.GetOutput())
    appendFilter.AddInputData(polyData3)

    line4.Update()
    polyData4 = vtk.vtkPolyData()
    polyData4.ShallowCopy(line4.GetOutput())
    appendFilter.AddInputData(polyData4)
    
    # Set up a mapper
    mapper = vtk.vtkPolyDataMapper()
    appendFilter.Update()
    mapper.SetInputConnection(appendFilter.GetOutputPort())

    # Set up an actor
    self.actor = vtk.vtkActor()
    self.actor.SetMapper(mapper)

    # Create the text for the plate label
    label = vtk.vtkVectorText()
    label.SetText(plate.Name)

    # Set up a mapper for the plate label
    lblMapper = vtk.vtkPolyDataMapper()
    lblMapper.SetInputConnection(label.GetOutputPort())

    # Set up an actor for the plate label
    self.lblActor = vtk.vtkFollower()
    self.lblActor.SetMapper(lblMapper)
    self.lblActor.SetScale(textHeight, textHeight, textHeight)
    self.lblActor.SetPosition((Xi+Xj+Xm+Xn)/4, (Yi+Yj+Ym+Yn)/4, (Zi+Zj+Zm+Zn)/4)

# Converts a node object into a node in its deformed position for the viewer
class VisDeformedNode():
  
  def __init__(self, node, scale_factor, textHeight=5, combo_name='Default'):
  
    # Calculate the node's deformed position
    newX = node.X + scale_factor*(node.DX[combo_name])
    newY = node.Y + scale_factor*(node.DY[combo_name])
    newZ = node.Z + scale_factor*(node.DZ[combo_name])

    # Generate a sphere for the node
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(newX, newY, newZ)
    sphere.SetRadius(0.6*textHeight)

    # Set up a mapper for the node
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())

    # Set up an actor for the node
    self.actor = vtk.vtkActor()
    self.actor.GetProperty().SetColor(255, 255, 0) # Yellow
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
    self.lblActor.SetPosition(newX + 0.6*textHeight, newY + 0.6*textHeight, newZ)
    self.lblActor.GetProperty().SetColor(255, 255, 0) # Yellow

# Converts a member object into a member in its deformed position for the viewer
class VisDeformedMember():

  def __init__(self, member, nodes, scale_factor, textHeight=5, combo_name='Default'):

    L = member.L() # Member length
    T = member.T() # Member local transformation matrix

    cos_x = array([T[0,0:3]]) # Direction cosines of local x-axis
    cos_y = array([T[1,0:3]]) # Direction cosines of local y-axis
    cos_z = array([T[2,0:3]]) # Direction cosines of local z-axis

    # Find the initial position of the local i-node
    # Step through each node
    for node in nodes:
      
      # Check to see if the current node is the i-node
      if node.Name == member.iNode.Name:
        Xi = node.X
        Yi = node.Y
        Zi = node.Z

    # Calculate the local y-axis displacements at 20 points along the member's length
    DY_plot = empty((0,3))
    for i in range(20):

      # Displacements in local coordinates
      dy_tot = member.Deflection('dy', L/19*i, combo_name)

      # Magnified displacements in global coordinates
      DY_plot = append(DY_plot, dy_tot*cos_y*scale_factor, axis=0)

    # Calculate the local z-axis displacements at 20 points along the member's length
    DZ_plot = empty((0,3)) 
    for i in range(20):

      # Displacements in local coordinates
      dz_tot = member.Deflection('dz', L/19*i, combo_name)

      # Magnified displacements in global coordinates
      DZ_plot = append(DZ_plot, dz_tot*cos_z*scale_factor, axis=0)

    # Calculate the local x-axis displacements at 20 points along the member's length
    DX_plot = empty((0,3)) 
    for i in range(20):

      # Displacements in local coordinates
      dx_tot = [[Xi, Yi, Zi]] + (L/19*i + member.Deflection('dx', L/19*i, combo_name)*scale_factor)*cos_x
      
      # Magnified displacements in global coordinates
      DX_plot = append(DX_plot, dx_tot, axis=0)
    
    # Sum the component displacements to obtain overall displacement
    D_plot = (DY_plot + DZ_plot + DX_plot)

    # Generate vtk points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(D_plot))

    for i in range(len(D_plot)):
      points.SetPoint(i, D_plot[i, 0], D_plot[i, 1], D_plot[i, 2])

    # Generate vtk lines
    lines = vtk.vtkCellArray()
    lines.InsertNextCell(len(D_plot))

    for i in range(len(D_plot)):
      lines.InsertCellPoint(i)

    # Create a polyline from the defined points and lines
    polyline = vtk.vtkPolyData()
    polyline.SetPoints(points)
    polyline.SetLines(lines)

    # Set up a mapper
    polylineMapper = vtk.vtkPolyDataMapper()
    polylineMapper.SetInputData(polyline)
    polylineMapper.Update()

    # Set up an actor for the polyline
    self.polylineActor = vtk.vtkActor()
    self.polylineActor.SetMapper(polylineMapper)
    self.polylineActor.GetProperty().SetColor(255,255,0) # Yellow

def DeformedShape(model, scale_factor, textHeight=5, combo_name='Default'):

  visNodes = []
  for node in model.Nodes:
    visNodes.append(VisNode(node, textHeight))

  visAuxNodes = []
  for auxnode in model.auxNodes:
    visAuxNodes.append(VisNode(auxnode, textHeight, color='red'))

  visMembers = []
  for member in model.Members:
    visMembers.append(VisMember(member, model.Nodes, textHeight))

  visDeformedNodes = []
  for node in model.Nodes:
    visDeformedNodes.append(VisDeformedNode(node, scale_factor, textHeight))
  
  visDeformedMembers = []
  for member in model.Members:
    visDeformedMembers.append(VisDeformedMember(member, model.Nodes, scale_factor, textHeight))

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
  
  # Add actors for each aux node
  for visAuxNode in visAuxNodes:

    # Add the actor for the node
    renderer.AddActor(visAuxNode.actor)

    # Add the actor for the node label
    renderer.AddActor(visAuxNode.lblActor)

    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    visAuxNode.lblActor.SetCamera(renderer.GetActiveCamera())

  # Add actors for each deformed member
  for visDeformedMember in visDeformedMembers:

    # Add the actor for the deformed member
    renderer.AddActor(visDeformedMember.polylineActor)

  # Add actors for each deformed node
  for visDeformedNode in visDeformedNodes:

    # Add the actor for the node
    renderer.AddActor(visDeformedNode.actor)

    # Add the actor for the node label
    renderer.AddActor(visDeformedNode.lblActor)

    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    visDeformedNode.lblActor.SetCamera(renderer.GetActiveCamera())

  # Setting the background to blue.
  renderer.SetBackground(0.1, 0.1, 0.4)

  # Reset the camera
  renderer.ResetCamera()

  window.Render()
  interactor.Start()
    
#%%
# TODO: Finish the following code for load visualization at a future time

# # Creates an arrow for the viewer
# class VisArrow():

#   # Constructor
#   def __init__(self, Position, Direction, Length, ArrowHeadHeight, LabelText = None):
#     '''
#     Parameters
#     ----------
#     Position : tuple
#       A tuple of X, Y and Z coordinates for the point of the arrow: (X, Y, Z)
#     Direction : tuple
#       A tuple indicating the direction vector for the arrow: (i, j, k)
#     Length : number
#       The length of the arrow
#     ArrowHeadHeight : number
#       The height of the arrow head
#     LabelText : string
#       Text that will show up at the tail of the arrow. Defaults to None.
#     '''

#     # Create the arrow head
#     arrowHead = vtk.vtkConeSource()
#     arrowHead.SetCenter(Position)
#     arrowHead.SetDirection(Direction)
#     arrowHead.SetHeight(ArrowHeadHeight)
#     arrowHead.SetRadius(ArrowHeadHeight/4)

#     mapper.SetInputConnection(cone.GetOutputPort())
    

#     # Calculate a scalar to convert the direction to a unit vector
#     scalar = (Direction[0]**2 + Direction[1]**2 + Direction[2]**2)**0.5

#     # Calculate a unit vector in the direction of the arrow
#     unitVector = (Direction[0]/scalar, Direction[1]/scalar, Direction[2]/scalar)
    
#     # Create the tail
#     tail = vtk.vtkLineSource()
#     tail.SetPoint1(Position)
#     tail.SetPoint2((Position[0]-Length*unitVector[0], Position[1]-Length*unitVector[1], Position[2]-Length*unitVector[2]))
    
