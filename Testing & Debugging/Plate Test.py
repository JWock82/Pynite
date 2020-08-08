# A First Course in the Finite Element Method, 4th Edition
# Daryl L. Logan
# Example 12.1
# Units for this model are pounds and inches

# Import 'FEModel3D' from 'PyNite'
from PyNite import FEModel3D
from numpy import allclose

# Create a new model
plModel= FEModel3D()

plModel.AddNode('N1', 0, 0, 0)
plModel.AddNode('N2', 10, 0, 0)
plModel.AddNode('N3', 20, 0, 0)
plModel.AddNode('N4', 0, 10, 0)
plModel.AddNode('N5', 10, 10, 0)
plModel.AddNode('N6', 20, 10, 0)
plModel.AddNode('N7', 0, 20, 0)
plModel.AddNode('N8', 10, 20, 0)
plModel.AddNode('N9', 20, 20, 0)

plModel.AddPlate('P1', 'N1', 'N2', 'N5', 'N4', 0.1, 30000000, 0.3)
plModel.AddPlate('P2', 'N2', 'N3', 'N6', 'N5', 0.1, 30000000, 0.3)
plModel.AddPlate('P3', 'N4', 'N5', 'N8', 'N7', 0.1, 30000000, 0.3)
plModel.AddPlate('P4', 'N5', 'N6', 'N9', 'N8', 0.1, 30000000, 0.3)

plModel.AddNodeLoad('N5', 'FZ', -100)

plModel.DefineSupport('N1', True, True, True, True, True, True)
plModel.DefineSupport('N2', True, True, True, True, True, True)
plModel.DefineSupport('N3', True, True, True, True, True, True)
plModel.DefineSupport('N4', True, True, True, True, True, True)
plModel.DefineSupport('N6', True, True, True, True, True, True)
plModel.DefineSupport('N7', True, True, True, True, True, True)
plModel.DefineSupport('N8', True, True, True, True, True, True)
plModel.DefineSupport('N9', True, True, True, True, True, True)

plModel.DefineSupport('N5', True, True, False, False, False, True)

# Check to see if the stiffness matrix for each plate is symmetric
# print(allclose(plModel.Plates[0].K(), plModel.Plates[0].K().T))
# print(allclose(plModel.Plates[1].K(), plModel.Plates[1].K().T))
# print(allclose(plModel.Plates[2].K(), plModel.Plates[2].K().T))
# print(allclose(plModel.Plates[3].K(), plModel.Plates[3].K().T))

# Check to see if the global stiffness matrix is symmetric
# print(allclose(plModel.K(Renumber=True), plModel.K(Renumber=False).T))

plModel.Analyze(check_statics=True)

print('Displacement: ', plModel.GetNode('N5').DZ)