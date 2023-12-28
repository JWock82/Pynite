# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 20:52:31 2017

@author: D. Craig Brinck, SE
"""
# %%
# A mathematically continuous beam segment
class BeamSegZ():
    """
    A mathematically continuous beam segment
    
    Properties
    ----------
    x1 : number
      The starting location of the segment relative to the start of the beam
    x2 : number
      The ending location of the segment relative to the start of the beam
    w1 : number
      The distributed load magnitude at the start of the segment
    w2 : number
      The distributed load magnitude at the end of the segment
    p1 : number
      The distributed axial load magnitude at the start of the segment
    p2 : number
      The distributed axial load magnitude at the end of the segment
    V1 : number
      The internal shear force at the start of the segment
    M1 : number
      The internal moment at the start of the segment
    P1 : number
      The internal axial force at the start of the segment
    T1 : number
      Torsional moment at start of segment
    theta1: number
      The slope (radians) at the start of the segment
    delta1: number
      The transverse displacement at the start of the segment
    delta_x1 : number
      The axial displacement at the start of the segment
    EI : number
      The flexural stiffness of the segment
    EA : number
      The axial stiffness of the segment
    
    Methods
    -------
    __init__()
      Constructor
    Length()
      Returns the length of the segment
    Shear(x)
      Returns the shear force at a location on the segment
    
    Notes
    -----
    Any unit system may be used as long as the units are consistent with each other
    
    """

    def __init__(self):
        """
        Constructor
        """
        
        self.x1 = None # Start location of beam segment (relative to start of beam)
        self.x2 = None # End location of beam segment (relative to start of beam)
        self.w1 = None # Linear distributed transverse load at start of segment
        self.w2 = None # Linear distributed transverse load at end of segment
        self.p1 = None # Linear distributed axial load at start of segment
        self.p2 = None # Linear distributed axial load at end of segment
        self.V1 = None # Internal shear force at start of segment
        self.M1 = None # Internal moment at start of segment
        self.P1 = None # Internal axial force at start of segment
        self.T1 = None # Torsional moment at start of segment
        self.theta1 = None # Slope at start of beam segment
        self.delta1 = None # Displacement at start of beam segment
        self.delta_x1 = None # Axial displacement at start of beam segment
        self.EI = None # Flexural stiffness of the beam segment
        self.EA = None # Axial stiffness of the beam segment

    # Returns the length of the segment
    def Length(self):
        """
        Returns the length of the segment
        """
        
        return self.x2 - self.x1

    # Returns the shear force at a location 'x' on the segment
    def Shear(self, x):
        
        V1 = self.V1
        w1 = self.w1
        w2 = self.w2
        L = self.Length()
        
        return V1 + w1*x + x**2*(-w1 + w2)/(2*L)

    # Returns the moment at a location on the segment
    def moment(self, x, P_delta=False):

        V1 = self.V1
        M1 = self.M1
        P1 = self.P1
        Px = self.axial(x)
        w1 = self.w1
        w2 = self.w2
        delta_1 = self.delta1
        delta_x = self.deflection(x)
        L = self.Length()
        
        M = M1 - V1*x - w1*x**2/2 - x**3*(-w1 + w2)/(6*L)

        # Include the P-little-delta moment if a P-Delta analysis was run
        if P_delta == True:
            M += P1*(delta_x - delta_1)
        
        # Return the computed moment
        return M

    # Returns the axial force at a location on the segment
    def axial(self, x):
        
        P1 = self.P1
        p1 = self.p1
        p2 = self.p2
        L = self.Length()
        
        return P1 + (p2 - p1)/(2*L)*x**2 + p1*x

    def Torsion(self):
        """
        Returns the torsional moment in the segment.
        """
        
        # The torsional moment is constant across the segment
        # This can be updated in the future for distributed torsional forces
        return self.T1
    
    def slope(self, x, P_delta=False):
        """Returns the slope of the elastic curve at any point `x` along the segment.

        :param x: Location (relative to start of segment) where slope is to be calculated.
        :type x: float
        :param P_delta: Indicates whether P-little-delta effects should be included in the calculation. Generally only used when a P-Delta analysis has been run. Defaults to False.
        :type P_delta: bool, optional
        :return: The slope of the elastic curve (radians) at location `x`.
        :rtype: float
        """
        
        V1 = self.V1
        M1 = self.M1
        P1 = self.P1
        w1 = self.w1
        w2 = self.w2
        delta_1 = self.delta1
        delta_x = self.deflection(x, P_delta)
        theta_1 = self.theta1
        L = self.Length()
        EI = self.EI

        theta_x = theta_1 - (-V1*x**2/2 - w1*x**3/6 + x*M1 + x**4*(w1 - w2)/(24*L))/EI
        theta_x = theta_1 - (-V1*x**2/2 - w1*x**3/6 + x*M1 + x**4*(w1 - w2)/(24*L))/EI

        if P_delta == True:
            theta_x -= x*(-P1*delta_1 + P1*delta_x)/EI
        
        # TODO: This is an old equation left for reference. Delete it after the new equation has been proved over time.
        # return theta_1 - (M1*x - V1*x**2/2 - w1*x**3/6 + x**4*(w1 - w2)/(24*L))/(EI)
        
        # Return the calculated slope
        return theta_x

    # Returns the deflection at a location on the segment
    def deflection(self, x, P_delta=False):
        
        V1 = self.V1
        M1 = self.M1
        P1 = self.P1
        w1 = self.w1
        w2 = self.w2
        theta_1 = self.theta1
        delta_1 = self.delta1
        d_delta = 1
        L = self.Length()
        EI = self.EI
        
        # Iteration is required to calculate P-little-delta effects
        if P_delta == True:

            # Initialize the deflection at `x` to match the deflection at the start of the segment
            delta_x = delta_1

            # Iterate until we reach a deflection convergence of 1%
            while d_delta > 0.01: 
                
                # Save the deflection value from the last iteration
                delta_last = delta_x

                # Compute the deflection
                # delta_x = delta1 + theta1*x + V1*x**3/(6*EI) + w1*x**4/(24*EI) + x**2*(-M1 + P1*delta1 - P1*delta_x)/(2*EI) + x**5*(-w1 + w2)/(120*EI*L)
                delta_x = delta_1 + theta_1*x + V1*x**3/(6*EI) + w1*x**4/(24*EI) + x**2*(-M1 + P1*delta_1 - P1*delta_x)/(2*EI) + x**5*(-w1 + w2)/(120*EI*L)
                
                # Check the change in deflection between iterations
                if delta_last != 0:
                    d_delta = abs(delta_x/delta_last - 1)
                else:
                    # Members with no relative deflection after at least one iteration need no further iterations
                    if delta_1 - delta_x == 0:
                        break
            
            # Return the calculated deflection
            return delta_x
        
        # Non-P-delta solutions are not iterative
        else:

            # Return the calcuated deflection
            # return delta1 + theta1*x - M1*x**2/(2*EI) + V1*x**3/(6*EI) + w1*x**4/(24*EI) + x**5*(-w1 + w2)/(120*EI*L)
            return delta_1 + theta_1*x + V1*x**3/(6*EI) + w1*x**4/(24*EI) + x**2*(-M1)/(2*EI) + x**5*(-w1 + w2)/(120*EI*L)

    def AxialDeflection(self, x):

        delta_x1 = self.delta_x1
        P1 = self.P1
        p1 = self.p1
        p2 = self.p2
        L = self.Length()
        EA = self.EA

        return delta_x1 - 1/EA*(P1*x + p1*x**2/2 + (p2 - p1)*x**3/(6*L))

    # Returns the maximum shear in the segment
    def max_shear(self):
        
        w1 = self.w1
        w2 = self.w2
        L = self.Length()
        
        # Determine possible locations of maximum shear
        if w1-w2 == 0:
            x1 = 0
        else:
            x1 = w1*L/(w1-w2)
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
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
    def min_shear(self):
        
        w1 = self.w1
        w2 = self.w2
        L = self.Length()
        
        # Determine possible locations of minimum shear
        if w1-w2 == 0:
            x1 = 0
        else:
            x1 = w1*L/(w1-w2)
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
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
    def max_moment(self, P_delta=False):
        
        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.Length()
        
        # Find the quadratic equation parameters
        a = -(w2-w1)/(2*L)
        b = -w1
        c = -V1
    
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
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
            x1 = 0
    
        if round(x2, 10) < 0 or round(x2, 10) > round(L, 10):
            x2 = 0
    
        # Find the moment at each location of interest
        M1 = self.moment(x1, P_delta)
        M2 = self.moment(x2, P_delta)
        M3 = self.moment(x3, P_delta)
        M4 = self.moment(x4, P_delta)
    
        # Return the maximum moment
        return max(M1, M2, M3, M4)

    # Returns the minimum moment in the segment
    def min_moment(self, P_delta=False):
        
        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.Length()
        
        # Find the quadratic equation parameters
        a = -(w2-w1)/(2*L)
        b = -w1
        c = -V1
    
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
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
            x1 = 0
    
        if round(x2, 10) < 0 or round(x2, 10) > round(L, 10):
            x2 = 0
    
        # Find the moment at each location of interest
        M1 = self.moment(x1, P_delta)
        M2 = self.moment(x2, P_delta)
        M3 = self.moment(x3, P_delta)
        M4 = self.moment(x4, P_delta)
    
        # Return the minimum moment
        return min(M1, M2, M3, M4)

    # Returns the maximum axial force in the segment
    def max_axial(self):
        
        p1 = self.p1
        p2 = self.p2
        L = self.Length()
        
        # Determine possible locations of maximum axial force
        if p1-p2 != 0:
            x1 = L*p1/(p1-p2)
        else:
            x1 = 0
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
            x1 = 0
    
        x2 = 0
        x3 = L
    
        # Find the axial force at each location of interest
        P1 = self.axial(x1)
        P2 = self.axial(x2)
        P3 = self.axial(x3)
    
        # Return the maximum axial force
        return max(P1, P2, P3)

    # Returns the minimum axial force in the segment
    def min_axial(self):
        
        p1 = self.p1
        p2 = self. p2
        L = self.Length()
        
        # Determine possible locations of minimum axial force
        if p1-p2 != 0:
            x1 = L*p1/(p1-p2)
        else:
            x1 = 0
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
            x1 = 0

        x2 = 0
        x3 = L
    
        # Find the axial force at each location of interest
        P1 = self.axial(x1)
        P2 = self.axial(x2)
        P3 = self.axial(x3)
    
        # Return the minimum axial force
        return min(P1, P2, P3)

    def MaxTorsion(self):
        """
        Returns the maximum torsional moment in the segment.
        """
    
        # Return the maximum torsional moment
        # Since the torsional moment is constant on the segment, the maximum torsional moment is T1
        # This can be updated in the future for distributed torsional forces
        return self.T1

    def MinTorsion(self):
        """
        Returns the minimum torsional moment in the segment.
        """

        # Return the minimum torsional moment
        # Since the torsional moment is constant on the segment, the minimum torsional moment is T1
        # This can be updated in the future for distributed torsional forces
        return self.T1
