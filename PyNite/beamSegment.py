# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 20:52:31 2017

@author: D. Craig Brinck, SE
"""
# %%
# A mathematically continuous beam segment
class BeamSegment():
    
    # Constructor
    def __init__(self):
        self.x1 = None
        self.x2 = None
        self.w1 = None
        self.w2 = None
        self.p1 = None
        self.p2 = None
        self.V1 = None
        self.M1 = None
        self.P1 = None
        self.theta1 = None
        self.Delta1 = None

    # Returns the shear force at a location on the segment
    def Shear(self, x):
        V1 = self.V1
        w1 = self.w1
        w2 = self.w2
        L = self.L
        return V1-(w2-w1)/(2*L)*x**2-w1*x
    
    # Returns the moment at a location on the segment
    def Moment(self, x):
        V1 = self.V1
        M1 = self.M1
        w1 = self.w1
        w2 = self.w2
        L = self.L
        return M1+V1*x-(w2-w1)/(6*L)*x**3-w1*x**2/2    
    
    # Returns the axial force at a location on the segment
    def Axial(self, x):
        P1 = self.P1
        p1 = self.p1
        p2 = self.p2
        L = self.L
        return P1+(p2-p1)/(2*L)*x**2+p1*x

    #Returns the rotation at a location on the segment
    def Slope(self, x, EI):
        
        V1 = self.V1
        M1 = self.M1
        w1 = self.w1
        w2 = self.w2
        theta1 = self.theta1
        L = self.L
        
        return theta1+1/EI*(M1*x+V1/2*x**2-(w2-w1)/(24*L)*x**4-w1/6*x**3)

    # Returns the deflection at a location on the segment
    def Deflection(self, x, EI):
        
        V1 = self.V1
        M1 = self.M1
        w1 = self.w1
        w2 = self.w2
        theta1 = self.theta1
        Delta1 = self.Delta1
        L = self.L
        
        return Delta1+theta1*x+1/EI*(M1/2*x**2+V1/6*x**3-(w2-w1)/(120*L)*x**5-w1/24*x**4)

    # Returns the maximum shear in the segment
    def MaxShear(self):
        
        w1 = self.w1
        w2 = self.w2
        L = self.x2-self.x1
        
        # Determine possible locations of maximum shear
        if w1-w2 == 0:
            x1 = 0
        else:
            x1 = w1*L/(w1-w2)
    
        if x1 < 0 or x1 > L:
            x1 = 0
    
        x2 = 0
        x3 = L
    
        # Find the shear at each location of interest
        V1 = self.Shear(x1)
        V2 = self.Shear(x2)
        V3 = self.Shear(x3)
    
        # Return the maximum shear
        return max(V1, V2, V3)

    # Returns the minimum shear in the segment
    def MinShear(self):
        
        w1 = self.w1
        w2 = self.w2
        L = self.x2-self.x1
        
        # Determine possible locations of minimum shear
        if w1-w2 == 0:
            x1 = 0
        else:
            x1 = w1*L/(w1-w2)
    
        if x1 < 0 or x1 > L:
            x1 = 0
    
        x2 = 0
        x3 = L
    
        # Find the shear at each location of interest
        V1 = self.Shear(x1)
        V2 = self.Shear(x2)
        V3 = self.Shear(x3)
    
        # Return the minimum shear
        return min(V1, V2, V3)

    # Returns the maximum moment in the segment
    def MaxMoment(self):
        
        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.x2-self.x1
        
        # Find the quadratic equation parameters
        a = (w1-w2)/(2*L)
        b = -w1
        c = V1
    
        # Determine possible locations of maximum moment
        if a == 0:
            if b != 0:
                x1 = -c/b
            else:
                x1 = 0
            x2 = 0
        elif b**2-4*a*c < 0:
            x1 = 0
            x2 = 0
        else:
            x1 = (-b+(b**2-4*a*c)**0.5)/(2*a)
            x2 = (-b-(b**2-4*a*c)**0.5)/(2*a)
    
        x3 = 0
        x4 = L
    
        if x1 < 0 or x1 > L:
            x1 = 0
    
        if x2 < 0 or x2 > L:
            x2 = 0
    
        # Find the moment at each location of interest
        M1 = self.Moment(x1)
        M2 = self.Moment(x2)
        M3 = self.Moment(x3)
        M4 = self.Moment(x4)
    
        # Return the maximum moment
        return max(M1, M2, M3, M4)

    # Returns the minimum moment in the segment
    def MinMoment(self):
        
        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.x2-self.x1
        
        # Find the quadratic equation parameters
        a = (w1-w2)/(2*L)
        b = -w1
        c = V1
    
        # Determine possible locations of minimum moment
        if a == 0:
            if b != 0:
                x1 = -c/b
            else:
                x1 = 0
            x2 = 0
        elif b**2-4*a*c < 0:
            x1 = 0
            x2 = 0
        else:
            x1 = (-b+(b**2-4*a*c)**0.5)/(2*a)
            x2 = (-b-(b**2-4*a*c)**0.5)/(2*a)
    
        x3 = 0
        x4 = L
    
        if x1 < 0 or x1 > L:
            x1 = 0
    
        if x2 < 0 or x2 > L:
            x2 = 0
    
        # Find the moment at each location of interest
        M1 = self.Moment(x1)
        M2 = self.Moment(x2)
        M3 = self.Moment(x3)
        M4 = self.Moment(x4)
    
        # Return the minimum moment
        return min(M1, M2, M3, M4)

    # Returns the maximum axial force in the segment
    def MaxAxial(self):
        
        p1 = self.p1
        p2 = self.p2
        L = self.x2-self.x1
        
        # Determine possible locations of maximum axial force
        if p1-p2 != 0:
            x1 = L*p1/(p1-p2)
        else:
            x1 = 0
    
        if x1 < 0 or x1 > L:
            x1 = 0
    
        x2 = 0
        x3 = L
    
        # Find the axial force at each location of interest
        P1 = self.Axial(x1)
        P2 = self.Axial(x2)
        P3 = self.Axial(x3)
    
        # Return the maximum axial force
        return max(P1, P2, P3)

    # Returns the minimum axial force in the segment
    def MinAxial(self):
        
        p1 = self.p1
        p2 = self. p2
        L = self.x2-self.x1
        
        # Determine possible locations of minimum axial force
        if p1-p2 != 0:
            x1 = L*p1/(p1-p2)
        else:
            x1 = 0
    
        if x1<0 or x1>L:
            x1 = 0

        x2 = 0
        x3 = L
    
        # Find the axial force at each location of interest
        P1 = self.Axial(x1)
        P2 = self.Axial(x2)
        P3 = self.Axial(x3)
    
        # Return the minimum axial force
        return min(P1, P2, P3)

