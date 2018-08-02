# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:04:56 2017

@author: D. Craig Brinck, SE
"""
# %%
from numpy import zeros
from PyNite.BeamSegment import BeamSegment
import PyNite.FixedEndReactions

# %%
# 3D frame element
class Member3D():

#%%
    # Constructor
    def __init__(self, Name, ID, iNode, jNode, E, G, Iy, Iz, J, A):
        
        self.Name = Name # A unique name for the member given by the user
        self.ID = ID # Unique index number for the member assigned by the program
        self.iNode = iNode # The element's i-node
        self.jNode = jNode # The element's j-node
        self.E = E # The modulus of elasticity of the element
        self.G = G # The shear modulus of the element
        self.Iy = Iy # The y-axis moment of inertia
        self.Iz = Iz # The z-axis moment of inertia
        self.J = J # The polar moment of inertia or torsional constant
        self.A = A # The cross-sectional area
        self.L = ((jNode.X-iNode.X)**2+(jNode.Y-iNode.Y)**2+(jNode.Z-iNode.Z)**2)**0.5 # Length of the element
        self.PtLoads = [] # A list of point loads applied to the element
        self.Moments = [] # A list of moments applied to the element
        self.DistLoads = [] # A list of linear distributed loads applied to the element
        self.SegmentsZ = [] # A list of mathematically continuous beam segments for z-bending
        self.SegmentsY = [] # A list of mathematically continuous beam segments for y-bending
        self.FER = zeros((12,1)) # The fixed end reaction vector
        self.D = zeros((12, 1))
        
#%%
    # Local stiffness matrix
    def k(self):
        
        E = self.E
        G = self.G
        Iy = self.Iy
        Iz = self.Iz
        J = self.J
        A = self.A
        L = self.L
        
        return matrix([[A*E/L, 0, 0 , 0, 0, 0, -A*E/L, 0, 0, 0, 0, 0],
                       [0, 12*E*Iz/L**3, 0, 0, 0, 6*E*Iz/L**2, 0, -12*E*Iz/L**3, 0, 0, 0, 6*E*Iz/L**2],
                       [0, 0, 12*E*Iy/L**3, 0, -6*E*Iy/L**2, 0, 0, 0, -12*E*Iy/L**3, 0, -6*E*Iy/L**2, 0],
                       [0, 0, 0, G*J/L, 0, 0, 0, 0, 0, -G*J/L, 0, 0],
                       [0, 0, -6*E*Iy/L**2, 0, 4*E*Iy/L, 0, 0, 0, 6*E*Iy/L**2, 0, 2*E*Iy/L, 0],
                       [0, 6*E*Iz/L**2, 0, 0, 0, 4*E*Iz/L, 0, -6*E*Iz/L**2, 0, 0, 0, 2*E*Iz/L],
                       [-A*E/L, 0, 0, 0, 0, 0, A*E/L, 0, 0, 0, 0, 0],
                       [0, -12*E*Iz/L**3, 0, 0, 0, -6*E*Iz/L**2, 0, 12*E*Iz/L**3, 0, 0, 0, -6*E*Iz/L**2],
                       [0, 0, -12*E*Iy/L**3, 0, 6*E*Iy/L**2, 0, 0, 0, 12*E*Iy/L**3, 0, 6*E*Iy/L**2, 0],
                       [0, 0, 0, -G*J/L, 0, 0, 0, 0, 0, G*J/L, 0, 0],
                       [0, 0, -6*E*Iy/L**2, 0, 2*E*Iy/L, 0, 0, 0, 6*E*Iy/L**2, 0, 4*E*Iy/L, 0],
                       [0, 6*E*Iz/L**2, 0, 0, 0, 2*E*Iz/L, 0, -6*E*Iz/L**2, 0, 0, 0, 4*E*Iz/L]])
    
#%%
    # Local fixed end reaction vector
    def fer(self):
        
        # Initialize the fixed end reaction vector
        ferVector = zeros(12,1)
        
        # Sum the fixed end reactions for the point loads
        for ptLoad in self.PtLoads:
            ferVector = ferVector + PyNite.FixedEndReactions.FER_PtLoad(ptLoad.P, ptLoad.x, ptLoad.L, ptLoad.Axis)
        
        # Sum the fixed end reactions for the concentrated moments
        for moment in self.Moments:
            ferVector = ferVector+PyNite.FixedEndReactions.FER_Moment(moment.m, moment.x, moment.L, moment.Axis)
        
        # Sum the fixed end reactions for the distributed loads
        for distLoad in self.DistLoads:
            ferVector = ferVector+PyNite.FixedEndReactions.FER_LinLoad(distLoad.w1, distLoad.w2, distLoad.x1, distLoad.x2, distLoad.L, distLoad.Axis)
        
        # Return the fixed end reaction vector        
        return ferVector
    
 #%%   
    # local end force vector
    def f(self):
        
        # Get the local stiffness matrix
        k = self.k
        
        # Get the local displacement vector

#%%
    # local displacement vector
    def d(self):
        
        # Calculate the local displacement vector
        return T*D
        
#%%  
    # Transformation matrix
    def T(self):
        
        x1 = self.iNode.X
        x2 = self.jNode.X
        y1 = self.iNode.Y
        y2 = self.jNode.Y
        z1 = self.iNode.Z
        z2 = self.jNode.Z
        L = self.L
        
        # Calculate direction cosines for the transformation matrix
        l = (x2-x1)/L
        m = (y2-y1)/L
        n = (z2-z1)/L
        D = (l**2+m**2)**0.5
        
        if l == 0 and m == 0 and n > 0:
            dirCos = matrix([[0, 0, 1],
                             [0, 1, 0],
                             [-1, 0, 0]])
        elif l == 0 and m == 0 and n < 0:
            dirCos = matrix([[0, 0, -1],
                             [0, 1, 0],
                             [1, 0, 0]])
        else:
            dirCos = matrix([[l, m, n],
                             [-m/D, l/D, 0],
                             [-l*n/D, -m*n/D, D]])
        
        # Build the transformation matrix
        transMatrix = zeros((12, 12))
        transMatrix[0:3, 0:3] = dirCos
        transMatrix[3:6, 3:6] = dirCos
        transMatrix[6:9, 6:9] = dirCos
        transMatrix[9:12, 9:12] = dirCos
        
        return transMatrix

#%%
    # Member global stiffness matrix
    def K(self):
        
        # Calculate and return the stiffness matrix in global coordinates
        return self.T.transpose()*self.k

#%% 
    # Global fixed end reaction vector
    def FER(self):
        
        # Get the local fixed end reaction vector
        ferVector = self.fer
        
        # Transform the vector into global coordinates
        transMatrix = self.T
        
        # Return the fixed end reaction vector
        return transMatrix*ferVector

#%%
    # Returns the global displacement vector
    def D(self):
        
        # Initialize the displacement vector
        Dvector = zeros(12,1)
        
        # Read in the global displacements
        Dvector.itemset((0, 0), self.iNode.DX)
        Dvector.itemset((1, 0), self.iNode.DY)

#%%
    # Adds a point load to the frame element
    def AddPointLoad(self, P, x, Axis):
        
        self.PtLoads.append([P, x, Axis])

#%%
    # Adds a concentrated moment to the frame element
    def AddMoment(self, M, x, Direction):
        
        self.Moments.append([M, x, Direction])

#%%
    def AddDistLoad(self, w1, w2, x1, x2, Axis, Type):
        
        """
        Adds a linear distributed load to the member
        
        Parameters
        ----------
        w1 : number
          Starting magnitude of the distributed load
        w2 : number
          Ending magnitude of the distributed load
        x1 : number
          Starting point of the distributed laod relative to the start of the member
        x2 : number
          Ending point of the distributed laod relative to the start of the member
        Axis : 1 or 2
          Use 1 for loads causing bending about the z-axis, or 2 for y-axis bending
        Type : 1 or 2
          Use 1 for transverse loading, and 2 for axial loading
        
        Return Value
        ------------
        Does not return a value. Simply adds the load to the model.
        
        Notes
        -----
        Any unit system may be used, as long as units are consistent with each other throughout the model.
        """
        
        self.DistLoads.append([w1, w2, x1, x2, Axis, Type])

#%%    
    # Divides the element up into mathematically continuous segments along each axis
    def Segment(self):
        
        # Get the members length and stiffness properties
        L = self.L
        E = self.E
        Iz = self.Iz
        Iy = self.Iy
        SegmentsZ = self.SegmentsZ
        SegmentsY = self.SegmentsY
        
        # Create a list of discontinuity locations
        disconts = [0, L]
        
        for load in self.PtLoads:
            disconts.append(load[1])
        
        for load in self.Moments:
            disconts.append(load[1])
        
        for load in self.DistLoads:
            disconts.append(load[2])
            disconts.append(load[3])
        
        # Sort the list and eliminate duplicate values
        disconts = sorted(set(disconts))
        
        # Clear out old data from any previous runs
        SegmentsZ.clear()
        SegmentsY.clear()
        
        # Create a list of mathematically continuous segments for each direction
        for index in range(len(disconts)):
            newSeg = BeamSegment()
            newSeg.x1 = disconts[index]
            newSeg.x2 = disconts[index+1]
            SegmentsZ.append(newSeg)
            SegmentsY.append(newSeg)
        
        # Get the member local end forces, local fixed end reactions, and local displacements
        f = self.f
        fer = self.fer
        d = self.d
        
        # Get the local deflections and calculate the local rotations at the start of the member
        m1z = f[5, 0]
        m2z = f[11, 0]
        m1y = f[4, 0]
        m2y = f[10, 0]
        fem1z = fer[5, 0]
        fem2z = fer[11, 0]
        fem1y = fer[4, 0]
        fem2y = fer[10, 0]
        delta1y = d[1, 0]
        delta2y = d[7, 0]
        delta1z = d[2, 0]
        delta2z = d[8, 0]
        SegmentsZ.Delta1 = delta1y
        SegmentsY.Delta1 = delta1z
        SegmentsZ.theta1 = 1/3*((m1z-fem1z)*L/(E*Iz)-(m2z-fem2z)*L/(2*E*Iz)+3*(delta2y-delta1y)/L)
        SegmentsY.theta1 = 1/3*((m1y-fem1y)*L/(E*Iy)-(m2y-fem2y)*L/(2*E*Iy)+3*(delta2z-delta1z)/L)
        
        # Add loads to each segment in the Z-direction
        for i in range(len(SegmentsZ)):
            
            # Get the starting point of the segment relative to the start of the beam
            x = SegmentsZ[i].x1
            
            # Initialize the distributed loads on the segment to zero
            SegmentsZ[i].w1 = 0
            SegmentsZ[i].w2 = 0
            
            # Initailize the rotation and displacement at the start of the segment
            if i > 0:
                SegmentsZ[i].theta1 = SegmentsZ[i-1].Rotation(SegmentsZ[i-1].Length, E*Iz)
                SegmentsZ[i].Delta1 = SegmentsZ[i-1].Deflection(SegmentsZ[i-1].Length, E*Iz)
            
            # Add the beam end forces to the segment
            SegmentsZ[i].p1 = f[0,0]
            SegmentsZ[i].V1 = f[1,0]
            SegmentsZ[i].M1 = f[3,0]+SegmentsZ[i].V1*x
            
            # Step through each distributed load on the member
            for j in range(len(self.DistLoads)):
                
                # Get the parameters for the distributed load
                w1 = self.DistLoads[j][0]
                w2 = self.DistLoads[j][1]
                x1 = self.DistLoads[j][2]
                x2 = self.DistLoads[j][3]
                Axis = self.DistLoads[j][4]
                Type = self.DistLoads[j][5]
                
                # Determine whether the load is on the segment
                if x >= x1 and x < x2 and Axis == 1:
                    # Break up the load and place it on the segment
                    if Type == 1:
                        SegmentsZ[i].w1 += (w2-w1)/(x2-x1)*(x-x1)+w1
                        SegmentsZ[i].w2 += (w2-w1)/(x2-x1)*(SegmentsZ[i].x2-x1)+w1
                    else:
                        SegmentsZ[i].p1 += (w2-w1)/(x2-x1)*(x-x1)+w1
                        SegmentsZ[i].p2 += (w2-w1)/(x2-x1)*(SegmentsZ[i].x2-x1)+w1
