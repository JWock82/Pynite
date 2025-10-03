from Pynite.BeamSegZ import BeamSegZ


# %%
class BeamSegY(BeamSegZ):

    # Returns the moment at a location on the segment
    def moment(self, x: float, P_delta: bool = False) -> float:
        '''
        Returns the moment at a location on the segment.

        Parameters
        ----------
        x : number
            Location (relative to start of segment) where moment is to be calculated
        '''

        V1 = self.V1
        M1 = self.M1
        P1 = self.P1
        w1 = self.w1
        w2 = self.w2
        L = self.Length()

        # M = M1 + V1*x + w1*x**2/2 + x**3*(-w1 + w2)/(6*L)
        M = -M1 - V1*x - w1*x**2/2 - x**3*(-w1 + w2)/(6*L)

        if P_delta == True:
            delta1 = self.delta1
            delta = self.deflection(x)
            M += P1*(delta - delta1)

        return M

    def slope(self, x: float, P_delta: bool = False) -> float:
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
        theta_1 = self.theta1

        L = self.Length()
        EI = self.EI

        if P_delta == True:
            delta_x = self.deflection(x, P_delta)
            delta_1 = self.delta1
            return theta_1 + (-V1*x**2/2 - w1*x**3/6 + x*(-M1 - P1*delta_1 + P1*delta_x) + x**4*(w1 - w2)/(24*L))/EI
        else:
            return theta_1 + (-V1*x**2/2 - w1*x**3/6 + x*(-M1) + x**4*(w1 - w2)/(24*L))/EI

    # Returns the deflection at a location on the segment
    def deflection(self, x: float, P_delta: bool = False) -> float:

        V1 = self.V1
        M1 = self.M1
        P1 = self.P1
        w1 = self.w1
        w2 = self.w2
        theta_1 = self.theta1
        delta_1 = self.delta1
        L = self.Length()
        EI = self.EI

        # Check if a P-delta solution is requested
        if P_delta == True:

            # Return the calculated deflection, amplified for P-delta effects
            return (delta_1 - theta_1*x + V1*x**3/(6*EI) + w1*x**4/(24*EI) - x**2*(-M1 - P1*delta_1)/(2*EI) - x**5*(w1 - w2)/(120*EI*L))/(1 + P1*x**2/(2*EI))

        else:

            # Return the calcuated deflection
            return delta_1 - theta_1*x + V1*x**3/(6*EI) + w1*x**4/(24*EI) - x**2*(-M1)/(2*EI) - x**5*(w1 - w2)/(120*EI*L)

    # Returns the maximum moment in the segment
    def max_moment(self, P_delta: bool = False) -> float:

        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.Length()

        # Find the quadratic equation parameters
        a = (w2-w1)/(2*L)
        b = w1
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

        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
            x1 = 0

        if round(x2, 10) < 0 or round(x2, 10) > round(L, 10):
            x2 = 0

        # Find the moment at each location of interest
        M1 = self.moment(x1)
        M2 = self.moment(x2)
        M3 = self.moment(x3)
        M4 = self.moment(x4)

        # Return the maximum moment
        return max(M1, M2, M3, M4)

    # Returns the minimum moment in the segment
    def min_moment(self, P_delta: bool = False) -> float:

        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.Length()

        # Find the quadratic equation parameters
        a = (w2-w1)/(2*L)
        b = w1
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
