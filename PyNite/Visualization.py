# Import the Visualization Toolkit if it's not alread imported
# You must be running a 64 bit version of Python for this to work
import vtk

# Create a data source
def RenderModel(model, textHeight=5):

  visNodes = []
  for node in model.Nodes:
    visNodes.append(VisNode(node, textHeight))

  visAuxNodes = []
  for auxnode in model.auxNodes:
    visAuxNodes.append(VisAuxNode(auxnode, textHeight))

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

    # Add the actors for the node supports
    for support in visNode.supportActors:
      renderer.AddActor(support)

  # Add actors for each plate
  for visPlate in visPlates:

    # Add the actors for the plate
    renderer.AddActor(visPlate.actor1)
    renderer.AddActor(visPlate.actor2)
    renderer.AddActor(visPlate.actor3)
    renderer.AddActor(visPlate.actor4)

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

    # Generate any supports that occur at the node
    supportMappers = []
    self.supportActors = []

    # Check for a fixed suppport
    if (node.SupportDX == True) and (node.SupportDY == True) and (node.SupportDZ == True) \
      and (node.SupportRX == True) and (node.SupportRY == True) and (node.SupportRZ == True):

      # Create the fixed support
      support = vtk.vtkCubeSource()
      support.SetCenter(node.X, node.Y, node.Z)
      support.SetXLength(textHeight*1.2)
      support.SetYLength(textHeight*1.2)
      support.SetZLength(textHeight*1.2)

      # Change the node's mapper to the fixed support
      mapper.SetInputConnection(support.GetOutputPort())
    
    # Check for a pinned support
    elif (node.SupportDX == True) and (node.SupportDY == True) and (node.SupportDZ == True) \
      and (node.SupportRX == False) and (node.SupportRY == False) and (node.SupportRZ == False):
      
      # Create the pinned support
      support = vtk.vtkConeSource()
      support.SetCenter(node.X, node.Y-0.6*textHeight, node.Z)
      support.SetDirection((0, 1, 0))
      support.SetHeight(textHeight*1.2)
      support.SetRadius(textHeight*1.2)

      # Change the node's mapper to the pinned support
      mapper.SetInputConnection(support.GetOutputPort())
    
    # Other support conditions
    else:

      # Restrained against X translation
      if (node.SupportDX == True):
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X-textHeight, node.Y, node.Z)
        support1.SetPoint2(node.X+textHeight, node.Y, node.Z)

        support2 = vtk.vtkConeSource()
        support2.SetCenter(node.X-textHeight, node.Y, node.Z)
        support2.SetDirection((1, 0, 0))
        support2.SetHeight(textHeight*0.6)
        support2.SetRadius(textHeight*0.3)

        support3 = vtk.vtkConeSource()
        support3.SetCenter(node.X+textHeight, node.Y, node.Z)
        support3.SetDirection((-1, 0, 0))
        support3.SetHeight(textHeight*0.6)
        support3.SetRadius(textHeight*0.3)

        # Set up mappers for the support
        supportMapper1 = vtk.vtkPolyDataMapper()
        supportMapper2 = vtk.vtkPolyDataMapper()
        supportMapper3 = vtk.vtkPolyDataMapper()
        supportMapper1.SetInputConnection(support1.GetOutputPort())
        supportMapper2.SetInputConnection(support2.GetOutputPort())
        supportMapper3.SetInputConnection(support3.GetOutputPort())
        supportMappers.append(supportMapper1)
        supportMappers.append(supportMapper2)
        supportMappers.append(supportMapper3)

        # Set up actors for the support
        supportActor1 = vtk.vtkActor()
        supportActor2 = vtk.vtkActor()
        supportActor3 = vtk.vtkActor()
        supportActor1.SetMapper(supportMapper1)
        supportActor2.SetMapper(supportMapper2)
        supportActor3.SetMapper(supportMapper3)
        self.supportActors.append(supportActor1)
        self.supportActors.append(supportActor2)
        self.supportActors.append(supportActor3)
      
      # Restrained against Y translation
      if (node.SupportDY == True):
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y-textHeight, node.Z)
        support1.SetPoint2(node.X, node.Y+textHeight, node.Z)

        support2 = vtk.vtkConeSource()
        support2.SetCenter(node.X, node.Y-textHeight, node.Z)
        support2.SetDirection((0, 1, 0))
        support2.SetHeight(textHeight*0.6)
        support2.SetRadius(textHeight*0.3)

        support3 = vtk.vtkConeSource()
        support3.SetCenter(node.X, node.Y+textHeight, node.Z)
        support3.SetDirection((0, -1, 0))
        support3.SetHeight(textHeight*0.6)
        support3.SetRadius(textHeight*0.3)

        # Set up mappers for the support
        supportMapper1 = vtk.vtkPolyDataMapper()
        supportMapper2 = vtk.vtkPolyDataMapper()
        supportMapper3 = vtk.vtkPolyDataMapper()
        supportMapper1.SetInputConnection(support1.GetOutputPort())
        supportMapper2.SetInputConnection(support2.GetOutputPort())
        supportMapper3.SetInputConnection(support3.GetOutputPort())
        supportMappers.append(supportMapper1)
        supportMappers.append(supportMapper2)
        supportMappers.append(supportMapper3)

        # Set up actors for the support
        supportActor1 = vtk.vtkActor()
        supportActor2 = vtk.vtkActor()
        supportActor3 = vtk.vtkActor()
        supportActor1.SetMapper(supportMapper1)
        supportActor2.SetMapper(supportMapper2)
        supportActor3.SetMapper(supportMapper3)
        self.supportActors.append(supportActor1)
        self.supportActors.append(supportActor2)
        self.supportActors.append(supportActor3)
      
      # Restrained against Z translation
      if (node.SupportDZ == True):
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y, node.Z-textHeight)
        support1.SetPoint2(node.X, node.Y, node.Z+textHeight)

        support2 = vtk.vtkConeSource()
        support2.SetCenter(node.X, node.Y, node.Z-textHeight)
        support2.SetDirection((0, 0, 1))
        support2.SetHeight(textHeight*0.6)
        support2.SetRadius(textHeight*0.3)

        support3 = vtk.vtkConeSource()
        support3.SetCenter(node.X, node.Y, node.Z+textHeight)
        support3.SetDirection((0, 0, -1))
        support3.SetHeight(textHeight*0.6)
        support3.SetRadius(textHeight*0.3)

        # Set up mappers for the support
        supportMapper1 = vtk.vtkPolyDataMapper()
        supportMapper2 = vtk.vtkPolyDataMapper()
        supportMapper3 = vtk.vtkPolyDataMapper()
        supportMapper1.SetInputConnection(support1.GetOutputPort())
        supportMapper2.SetInputConnection(support2.GetOutputPort())
        supportMapper3.SetInputConnection(support3.GetOutputPort())
        supportMappers.append(supportMapper1)
        supportMappers.append(supportMapper2)
        supportMappers.append(supportMapper3)

        # Set actors for the support
        supportActor1 = vtk.vtkActor()
        supportActor2 = vtk.vtkActor()
        supportActor3 = vtk.vtkActor()
        supportActor1.SetMapper(supportMapper1)
        supportActor2.SetMapper(supportMapper2)
        supportActor3.SetMapper(supportMapper3)
        self.supportActors.append(supportActor1)
        self.supportActors.append(supportActor2)
        self.supportActors.append(supportActor3)

      # Restrained against rotation about the X-axis
      if (node.SupportRX == True):
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X-1.6*textHeight, node.Y, node.Z)
        support1.SetPoint2(node.X+1.6*textHeight, node.Y, node.Z)

        support2 = vtk.vtkCubeSource()
        support2.SetCenter(node.X-1.9*textHeight, node.Y, node.Z)
        support2.SetXLength(textHeight*0.6)
        support2.SetYLength(textHeight*0.6)
        support2.SetZLength(textHeight*0.6)

        support3 = vtk.vtkCubeSource()
        support3.SetCenter(node.X+1.9*textHeight, node.Y, node.Z)
        support3.SetXLength(textHeight*0.6)
        support3.SetYLength(textHeight*0.6)
        support3.SetZLength(textHeight*0.6)

        # Set up mappers for the support
        supportMapper1 = vtk.vtkPolyDataMapper()
        supportMapper2 = vtk.vtkPolyDataMapper()
        supportMapper3 = vtk.vtkPolyDataMapper()
        supportMapper1.SetInputConnection(support1.GetOutputPort())
        supportMapper2.SetInputConnection(support2.GetOutputPort())
        supportMapper3.SetInputConnection(support3.GetOutputPort())
        supportMappers.append(supportMapper1)
        supportMappers.append(supportMapper2)
        supportMappers.append(supportMapper3)

        # Set up actors for the support
        supportActor1 = vtk.vtkActor()
        supportActor2 = vtk.vtkActor()
        supportActor3 = vtk.vtkActor()
        supportActor1.SetMapper(supportMapper1)
        supportActor2.SetMapper(supportMapper2)
        supportActor3.SetMapper(supportMapper3)
        self.supportActors.append(supportActor1)
        self.supportActors.append(supportActor2)
        self.supportActors.append(supportActor3)
      
      # Restrained against rotation about the Y-axis
      if (node.SupportRY == True):
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y-1.6*textHeight, node.Z)
        support1.SetPoint2(node.X, node.Y+1.6*textHeight, node.Z)

        support2 = vtk.vtkCubeSource()
        support2.SetCenter(node.X, node.Y-1.9*textHeight, node.Z)
        support2.SetXLength(textHeight*0.6)
        support2.SetYLength(textHeight*0.6)
        support2.SetZLength(textHeight*0.6)

        support3 = vtk.vtkCubeSource()
        support3.SetCenter(node.X, node.Y+1.9*textHeight, node.Z)
        support3.SetXLength(textHeight*0.6)
        support3.SetYLength(textHeight*0.6)
        support3.SetZLength(textHeight*0.6)

        # Set up mappers for the support
        supportMapper1 = vtk.vtkPolyDataMapper()
        supportMapper2 = vtk.vtkPolyDataMapper()
        supportMapper3 = vtk.vtkPolyDataMapper()
        supportMapper1.SetInputConnection(support1.GetOutputPort())
        supportMapper2.SetInputConnection(support2.GetOutputPort())
        supportMapper3.SetInputConnection(support3.GetOutputPort())
        supportMappers.append(supportMapper1)
        supportMappers.append(supportMapper2)
        supportMappers.append(supportMapper3)

        # Set up actors for the support
        supportActor1 = vtk.vtkActor()
        supportActor2 = vtk.vtkActor()
        supportActor3 = vtk.vtkActor()
        supportActor1.SetMapper(supportMapper1)
        supportActor2.SetMapper(supportMapper2)
        supportActor3.SetMapper(supportMapper3)
        self.supportActors.append(supportActor1)
        self.supportActors.append(supportActor2)
        self.supportActors.append(supportActor3)
      
      # Restrained against rotation about the Z-axis
      if (node.SupportRZ == True):
        
        # Create the support
        support1 = vtk.vtkLineSource()  # The line showing the support direction
        support1.SetPoint1(node.X, node.Y, node.Z-1.6*textHeight)
        support1.SetPoint2(node.X, node.Y, node.Z+1.6*textHeight)

        support2 = vtk.vtkCubeSource()
        support2.SetCenter(node.X, node.Y, node.Z-1.9*textHeight)
        support2.SetXLength(textHeight*0.6)
        support2.SetYLength(textHeight*0.6)
        support2.SetZLength(textHeight*0.6)

        support3 = vtk.vtkCubeSource()
        support3.SetCenter(node.X, node.Y, node.Z+1.9*textHeight)
        support3.SetXLength(textHeight*0.6)
        support3.SetYLength(textHeight*0.6)
        support3.SetZLength(textHeight*0.6)

        # Set up mappers for the support
        supportMapper1 = vtk.vtkPolyDataMapper()
        supportMapper2 = vtk.vtkPolyDataMapper()
        supportMapper3 = vtk.vtkPolyDataMapper()
        supportMapper1.SetInputConnection(support1.GetOutputPort())
        supportMapper2.SetInputConnection(support2.GetOutputPort())
        supportMapper3.SetInputConnection(support3.GetOutputPort())
        supportMappers.append(supportMapper1)
        supportMappers.append(supportMapper2)
        supportMappers.append(supportMapper3)

        # Set up actors for the support
        supportActor1 = vtk.vtkActor()
        supportActor2 = vtk.vtkActor()
        supportActor3 = vtk.vtkActor()
        supportActor1.SetMapper(supportMapper1)
        supportActor2.SetMapper(supportMapper2)
        supportActor3.SetMapper(supportMapper3)
        self.supportActors.append(supportActor1)
        self.supportActors.append(supportActor2)
        self.supportActors.append(supportActor3)

# Converts a auxNode object into a auxnode for the viewer
class VisAuxNode():
 
  # Constructor
  def __init__(self, auxnode, textHeight):
    
    # Get the node's position
    X = auxnode.X          # Global X coordinate
    Y = auxnode.Y          # Global Y coordinate
    Z = auxnode.Z          # Global Z coordinate

    # Generate a sphere for the node
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(X, Y, Z)
    sphere.SetRadius(0.6*textHeight)

    # Set up a mapper for the node
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())

    # Set up an actor for the node
    self.actor = vtk.vtkActor()
    self.actor.GetProperty().SetColor(255,0,0) #red
    self.actor.SetMapper(mapper)

    # Create the text for the node label
    label = vtk.vtkVectorText()
    label.SetText(auxnode.Name)
    
    # Set up a mapper for the node label
    lblMapper = vtk.vtkPolyDataMapper()
    lblMapper.SetInputConnection(label.GetOutputPort())

    # Set up an actor for the node label
    self.lblActor = vtk.vtkFollower()
    self.lblActor.SetMapper(lblMapper)
    self.lblActor.SetScale(textHeight, textHeight, textHeight)
    self.lblActor.SetPosition(X+0.6*textHeight, Y+0.6*textHeight, Z)
    self.lblActor.GetProperty().SetColor(255,0,0) #red
        
        
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

    # Set up mappers for the plate
    mapper1 = vtk.vtkPolyDataMapper()
    mapper2 = vtk.vtkPolyDataMapper()
    mapper3 = vtk.vtkPolyDataMapper()
    mapper4 = vtk.vtkPolyDataMapper()

    mapper1.SetInputConnection(line1.GetOutputPort())
    mapper2.SetInputConnection(line2.GetOutputPort())
    mapper3.SetInputConnection(line3.GetOutputPort())
    mapper4.SetInputConnection(line4.GetOutputPort())

    # Set up actors for the plate
    self.actor1 = vtk.vtkActor()
    self.actor2 = vtk.vtkActor()
    self.actor3 = vtk.vtkActor()
    self.actor4 = vtk.vtkActor()

    self.actor1.SetMapper(mapper1)
    self.actor2.SetMapper(mapper2)
    self.actor3.SetMapper(mapper3)
    self.actor4.SetMapper(mapper4)

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

class VisDeformedNode():

    def __init__(self, node, scalefactor, textHeight):

        #Node position in deformed configuration
        newX = node.X + scalefactor*(node.DX)
        newY = node.Y + scalefactor*(node.DY)
        newZ = node.Z + scalefactor*(node.DZ)

        # Generate a sphere for the node
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(newX, newY, newZ)
        sphere.SetRadius(0.6*textHeight)

        # Set up a mapper for the node
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())

        # Set up an actor for the node
        self.actor = vtk.vtkActor()
        self.actor.GetProperty().SetColor(255,255,0) #yellow
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
        self.lblActor.SetPosition(newX+0.6*textHeight, newY+0.6*textHeight, newZ)
        self.lblActor.GetProperty().SetColor(255,255,0)

class VisDeformedMember():

    def __init__(self, member, nodes, scalefactor, textHeight):


        L=member.L()
        T=member.T()
        cos_x=array([T[0,0:3]]) #direction cosines of local x-axis
        cos_y=array([T[1,0:3]]) #direction cosines of local y-axis
        cos_z=array([T[2,0:3]]) #direction cosines of local z-axis


        # Step through each node in the model and find the position of the i-node 
        for node in nodes:

            # Check to see if the current node is the i-node
            if node.Name == member.iNode.Name:
                newXi = (node.X) + (node.DX)
                newYi = (node.Y) + (node.DY)
                newZi = (node.Z) + (node.DZ)
                     
        dy_plot=empty((0,3))
        for i in range(20):

            dy_tot=member.Deflection("dy",L/ 19 * i)
            dy_plot=append(dy_plot,(dy_tot*cos_y)*scalefactor, axis=0)


        dz_plot=empty((0,3)) 
        for i in range(20):

            dz_tot=member.Deflection("dz",L/ 19 * i)
            dz_plot=append(dz_plot,(dz_tot*cos_z)*scalefactor, axis=0)


        x_plot=empty((0,3)) 
        for i in range(20):

            x= [[newXi, newYi, newZi]] + ((L/19 *i) + member.AxialTranslation(L/19 *i)*scalefactor)*cos_x
            x_plot=append(x_plot,x,axis=0)
        

        #Vector to plot        
        d_plot=(dy_plot+dz_plot+x_plot)

        # Generate points
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(len(d_plot))

        for i in range(len(d_plot)):
            points.SetPoint(i, d_plot[i,0], d_plot[i,1], d_plot[i,2])

        # Generate lines
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(len(d_plot))

        for i in range(len(d_plot)):
            lines.InsertCellPoint(i)

        polyline = vtk.vtkPolyData()
        polyline.SetPoints(points)
        polyline.SetLines(lines)

        polylineMapper = vtk.vtkPolyDataMapper()
        polylineMapper.SetInputData(polyline)
        polylineMapper.Update()

        self.polylineActor = vtk.vtkActor()
        self.polylineActor.SetMapper(polylineMapper)
        self.polylineActor.GetProperty().SetColor(255,255,0) #yellow


def DeformedShape(model, scalefactor, textHeight):

  visNodes = []
  for node in model.Nodes:

    visNodes.append(VisNode(node, textHeight))

  visAuxNodes = []
  for auxnode in model.auxNodes:
    visAuxNodes.append(VisAuxNode(auxnode, textHeight))

  visMembers = []
  for member in model.Members:
    visMembers.append(VisMember(member, model.Nodes, textHeight))

  visDeformedNodes = []
  for node in model.Nodes:
    visDeformedNodes.append(VisDeformedNode(node,scalefactor, textHeight))
  
  visDeformedMembers = []
  for member in model.Members:
    visDeformedMembers.append(VisDeformedMember(member, model.Nodes,scalefactor, textHeight))


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

    # Add the actor for the member
    renderer.AddActor(visDeformedMember.polylineActor)

  # Add actors for each deformed node
  for visDeformedNode in visDeformedNodes:

    # Add the actor for the node
    renderer.AddActor(visDeformedNode.actor)

    # Add the actor for the node label
    renderer.AddActor(visDeformedNode.lblActor)

    # Set the text to follow the camera as the user interacts
    # This next line will require us to reset the camera when we're done (below)
    visNode.lblActor.SetCamera(renderer.GetActiveCamera())

    # Add the actors for the node supports
    for support in visNode.supportActors:
      renderer.AddActor(support)

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
#     """
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
#     """

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
    
