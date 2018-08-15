# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:04:56 2017

@author: D. Craig Brinck, SE
"""
# %%
from numpy import zeros, matrix
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
        self.PtLoads = [] # A list of point loads applied to the element (Direction, P, x)
        self.DistLoads = [] # A list of linear distributed loads applied to the element (Direction, w1, w2, x1, x2)
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
        Dvector = zeros((12,1))
        
        # Read in the global displacements
        Dvector.itemset((0, 0), self.iNode.DX)
        Dvector.itemset((1, 0), self.iNode.DY)

#%%
    # Adds a concentrated moment to the frame element
    def AddMoment(self, M, x, Direction):
        
        self.Moments.append((M, x, Direction))

#%%    
    # Divides the element up into mathematically continuous segments along each axis
    def SegmentMember(self):
        
        # Get the member's length and stiffness properties
        L = self.L
        E = self.E
        Iz = self.Iz
        Iy = self.Iy
        SegmentsZ = self.SegmentsZ
        SegmentsY = self.SegmentsY
        
        # Create a list of discontinuity locations
        disconts = [0, L] # Member ends
        
        for load in self.PtLoads: 
            disconts.append(load[2]) # Point load locations
        
        for load in self.DistLoads: 
            disconts.append(load[3]) # Distributed load start locations
            disconts.append(load[4]) # Distributed load end locations
        
        # Sort the list and eliminate duplicate values
        disconts = sorted(set(disconts))
        
        # Clear out old data from any previous analyses
        SegmentsZ.clear()
        SegmentsY.clear()
        
        # Create a list of mathematically continuous segments for each direction
        for index in range(len(disconts)):
            
            # z-direction segments (bending about local z-axis)
            newSeg = BeamSegment()        # Create the new segment
            newSeg.x1 = disconts[index]   # Segment start location
            newSeg.x2 = disconts[index+1] # Segment end location
            newSeg.EI = E*Iz              # Segment flexural stiffness
            SegmentsZ.append(newSeg)      # Add the segment to the list
            
            # y-direction segments (bending about local y-axis)
            newSeg = BeamSegment()        # Create the new segment
            newSeg.x1 = disconts[index]   # Segment start location
            newSeg.x2 = disconts[index+1] # Segment end location
            newSeg.EI = E*Iy              # Segment flexural stiffness
            SegmentsY.append(newSeg)      # Add the segment to the list
        
        # Get the member local end forces, local fixed end reactions, and local displacements
        f = self.f     # Member local end force vector
        fer = self.fer # Member local fixed end reaction vector
        d = self.d     # Member local displacement vector
        
        # Get the local deflections and calculate the slope at the start of the member
        # Note that the slope may not be available directly from the local displacement vector if member end releases have been used,
        # so slope-deflection has been applied to solve for it
        m1z = f[5, 0]       # local z-axis moment at start of member
        m2z = f[11, 0]      # local z-axis moment at end of member
        m1y = f[4, 0]       # local y-axis moment at start of member
        m2y = f[10, 0]      # local y-axis moment at end of member
        fem1z = fer[5, 0]   # local z-axis fixed end moment at start of member
        fem2z = fer[11, 0]  # local z-axis fixed end moment at end of member
        fem1y = fer[4, 0]   # local y-axis fixed end moment at start of member
        fem2y = fer[10, 0]  # local y-axis fixed end moment at end of member
        delta1y = d[1, 0]   # local y displacement at start of member
        delta2y = d[7, 0]   # local y displacement at end of member
        delta1z = d[2, 0]   # local z displacement at start of member
        delta2z = d[8, 0]   # local z displacement at end of member
        SegmentsZ[0].Delta1 = delta1y
        SegmentsY[0].Delta1 = delta1z
        SegmentsZ[0].theta1 = 1/3*((m1z - fem1z)*L/(E*Iz) - (m2z - fem2z)*L/(2*E*Iz) + 3*(delta2y - delta1y)/L)
        SegmentsY[0].theta1 = 1/3*((m1y - fem1y)*L/(E*Iy) - (m2y - fem2y)*L/(2*E*Iy) + 3*(delta2z - delta1z)/L)
        
        # Add loads to each segment
        for i in range(len(SegmentsZ)):
            
            # Get the starting point of the segment
            x = SegmentsZ[i].x1
            
            # Initialize the distributed loads on the segment to zero
            SegmentsZ[i].w1 = 0
            SegmentsZ[i].w2 = 0
            SegmentsZ[i].p1 = 0
            SegmentsZ[i].p2 = 0
            SegmentsY[i].w1 = 0
            SegmentsY[i].w2 = 0
            SegmentsY[i].p1 = 0
            SegmentsY[i].p2 = 0
            
            # Initialize the slope and displacement at the start of the segment
            if i > 0: # The first segment has already been initialized
                SegmentsZ[i].theta1 = SegmentsZ[i-1].Slope(SegmentsZ[i-1].Length)
                SegmentsZ[i].Delta1 = SegmentsZ[i-1].Deflection(SegmentsZ[i-1].Length)
                SegmentsY[i].theta1 = SegmentsY[i-1].Slope(SegmentsY[i-1].Length)
                SegmentsY[i].Delta1 = SegmentsY[i-1].Deflection(SegmentsY[i-1].Length)
                
            # Add the effects of the beam end forces to the segment
            SegmentsZ[i].P1 = f[0,0]
            SegmentsZ[i].V1 = f[1,0]
            SegmentsZ[i].M1 = f[5,0] + f[1,0]*x
            SegmentsY[i].P1 = f[0,0]
            SegmentsY[i].V1 = f[2,0]
            SegmentsY[i].M1 = f[4,0] + f[2,0]*x
            
            # Add effects of point loads occuring prior to this segment
            for ptLoad in self.PtLoads:
                
                if ptLoad[2] <= x:
                    
                    if ptLoad[0] == "Fx":
                        SegmentsZ[i].P1 += ptLoad[1]
                    elif ptLoad[0] == "Fy":
                        SegmentsZ[i].V1 -= ptLoad[1]
                        SegmentsZ[i].M1 -= ptLoad[1]*(x-ptLoad[2])
                    elif ptLoad[0] == "Fz":
                        SegmentsY[i].V1 -= ptLoad[1]
                        SegmentsY[i].M1 -= ptLoad[1]*(x - ptLoad[2])
                    elif ptLoad[0] == "My":
                        SegmentsY[i].M1 -= ptLoad[1]
                    elif ptLoad[0] == "Mz":
                        SegmentsZ[i].M1 -= ptLoad[1]
            
            # Add distributed loads to the segment
            for distLoad in self.DistLoads:
                
                # Get the parameters for the distributed load
                Direction = distLoad[0]
                w1 = distLoad[1]
                w2 = distLoad[2]
                x1 = distLoad[3]
                x2 = distLoad[4]
            
                # Determine if the load affects the moment at the start of the segment
                if x1 < x:
                    
                    if Direction == "Fx":
                        
                        # Determine if the load ends after the start of the segment
                        if x2 > x:
                                                
                            # Break up the load and place it on the segment
                            # Note that 'w1' and 'w2' are really 'p1' and 'p2' here
                            SegmentsZ[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                            SegmentsZ[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1
                            SegmentsY[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                            SegmentsY[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1
                            
                            # Calculate the magnitude of the load at the start of the segment
                            w2 = w1+(w2-w1)/(x2-x1)*(x-x1)
                            x2 = x
                        
                        # Calculate the axial force at the start of the segment
                        SegmentsZ[i].P1 -= (w1 + w2)/2*(x2 - x1)
                        SegmentsY[i].P1 -= (w1 + w2)/2*(x2 - x1)
                    
                    elif Direction == "Fy":
                        
                        # Determine if the load ends after the start of the segment
                        if x2 > x:
                                                
                            # Break up the load and place it on the segment
                            SegmentsZ[i].w1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                            SegmentsZ[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1
                            
                            # Calculate the magnitude of the load at the start of the segment
                            # This will be used as the 'x2' value for the load prior to the start of the segment
                            w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                            x2 = x
                        
                        # Calculate the shear and moment at the start of the segment due to the load
                        SegmentsZ[i].V1 -= (w1 + w2)/2*(x2 - x1)
                        SegmentsZ[i].M1 -= (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6
                    
                    elif Direction == "Fz":
                        
                        # Determine if the load ends after the start of the segment
                        if x2 > x:
                                                
                            # Break up the load and place it on the segment
                            SegmentsY[i].w1 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x1 - x1) + w1
                            SegmentsY[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1
                            
                            # Calculate the magnitude of the load at the start of the segment
                            w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                            x2 = x
                        
                        # Calculate the shear and moment at the start of the segment due to the load
                        SegmentsY[i].V1 -= (w1 + w2)/2*(x2 - x1)
                        SegmentsY[i].M1 -= (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6