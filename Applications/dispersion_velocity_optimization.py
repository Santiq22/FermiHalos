# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 09:30:31 2023

@author: RAR colaboration
"""

# ================================================= Packages ==================================================== #
from scipy.integrate import solve_ivp
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import os
import sys
import time
# =============================================================================================================== #

""" The approach here is to define a function that computes the dispersion velocity in the line of sight (LOS). 
This quantity is related with the mass distribution of the satellite and is a direct astrophysical observation. 
After solving the Jean equation, one can computes the mass density times the mean square radial velocity needed
to compute the LOS velocity dispersion. Given a dataset of LOS velocity dispersions, this code pretends to find
the best fit parameters that generates a model able to reproduce the observations. """

# =========================================== Loading the data set ============================================== #
""" You should load the data set as a N x 3 array, where N is the number of observations and 3 because we have the
LOS velocity dispersion observations, the radius of such observations and their corresponding uncertainty. The code 
will read the dataset from the Data folder, so you ougth to save the file there. """

# Here you should put as a string variable the name of the txt file containing the observations
observations_file = 'test.txt'
data = np.loadtxt(os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'Data', observations_file)))
# =============================================================================================================== #

# ============================================= Global constants ================================================ #
G = 4.3009e-6                                #  Newton's gravitational constant - kpc (km/s)^2 M_sun^-1
# =============================================================================================================== #

# =========================================== Functions definition ============================================== #
# 2D mass density, Plummer model
def I(R, r_half, M_0):
    # M_0: Total stellar mass
    return (M_0/(np.pi*r_half**2))*(1.0/(1.0 + (R/r_half)**2)**2)

# 3d mass density, Plummer model
def nu(x, r_half, M_0):
    # M_0: Total stellar mass
    return (3.0*M_0/(4.0*np.pi*r_half**3))*(1.0/(1.0 + (x/r_half)**2)**(5/2))

# Integrand of the 3D mass density times mean square velocity squared
def integrand_nu_v_r2(s, M, beta, r_half, M_0):
    return s**(2.0*beta - 2.0)*nu(s, r_half, M_0)*M(s)

# 3D mass density times mean square velocity squared
@np.vectorize
def nu_v_r2(r, M, beta, r_half, M_0):
    int = quad(integrand_nu_v_r2, r, np.inf, args=(M, beta, r_half, M_0))
    return G*r**(-2.0*beta)*int

# Integrand of squared line of sight velocity dispersion
def integrand_sigma_los2(r, R, M, beta, r_half, M_0):
    return (1.0 - beta*(R/r)**2)*((nu_v_r2(r, M, beta, r_half, M_0)*r)/(np.sqrt(r**2 - R**2)))

# Squared line of sight velocity dispersion
@np.vectorize
def sigma_los2(R, M, beta, r_half, M_0):
    int = quad(integrand_sigma_los2, R, np.inf, args=(R, M, beta, r_half, M_0))[0]
    return 2.0/I(R, r_half, M_0)*int
# =============================================================================================================== #

# ========================================== Optimization procedure ============================================= #

# =============================================================================================================== #