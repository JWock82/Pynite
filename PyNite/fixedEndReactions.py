# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:58:03 2017

@author: D. Craig Brinck, SE
"""
# %%
import numpy as np

# %%
# Returns the fixed end reaction vector for a point load
def FER_PtLoad(P, x, L):
    
    # Define variables
    b = L - x
    
    # Create the fixed end reaction vector
    FER = np.zeros((6, 1))
    
    # Populate the fixed end reaction vector
    FER[1, 0] = P*b**2*(L+2*x)/L**3
    FER[2, 0] = P*x*b**2/L**2
    FER[4, 0] = P*x**2*(L+2*b)/L**3
    FER[5, 0] = -P*x**2*b/L**2
    
    return FER
    
# %%
# Returns the fixed end reaction vector for a moment
def FER_Moment(M, x, L):
    
    # Define variables
    b = L - x
    
    # Create the fixed end reaction vector
    FER = np.zeros((6, 1))
    
    # Populate the fixed end reaction vector
    FER[1, 0] = M*(x**2+b**2-4*x*b-L**2)/L**3
    FER[2, 0] = -M*b*(2*x-b)/L**2
    FER[4, 0] = -M*(x**2+b**2-4*x*b-L**2)/L**3
    FER[5, 0] = -M*x*(2*b-x)/L**2
    
    return FER
    
# %%
# Returns the fixed end reaction vector for a linear distributed load
def FER_LinLoad(w1, w2, x1, x2, L):
        
    # Create the fixed end reaction vector
    FER = np.zeros((6, 1))
    
    # Populate the fixed end reaction vector
    FER[1, 0] = -(x1 - x2)*(10*L**3*w1 + 10*L**3*w2 - 15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3)
    FER[2, 0] = -(x1 - x2)*(20*L**2*w1*x1 + 10*L**2*w1*x2 + 10*L**2*w2*x1 + 20*L**2*w2*x2 - 30*L*w1*x1**2 - 20*L*w1*x1*x2 - 10*L*w1*x2**2 - 10*L*w2*x1**2 - 20*L*w2*x1*x2 - 30*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2)
    FER[4, 0] = (x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3)
    FER[5, 0] = -(x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2)
    
    return FER

# %%
# Returns the fixed end reaction vector for an axial point load
def FER_AxialPtLoad(P, x, L):
    
    # Create the fixed end reaction vector
    FER = np.matrix((6, 1))
    
    # Populate the fixed end reaction vector
    FER[0, 0] = -P*(L-x)/L
    FER[3, 0] = -P*x/L
    
    return FER

# %%
# Returns the fixed end reaction vector for a distributed axial load
def FER_AxialLinLoad(p1, p2, x1, x2, L):
    
    # Create the fixed end reaction vector
    FER = np.zeros((6, 1))
    
    # Populate the fixed end reaction vector
    FER[0, 0] = 1/(6*L)*(x1-x2)*(3*L*p1+3*L*p2-2*p1*x1-p1*x2-p2*x1-2*p2*x2)
    FER[3, 0] = 1/(6*L)*(x1-x2)*(2*p1*x1+p1*x2+p2*x1+2*p2*x2)
    
    return FER
