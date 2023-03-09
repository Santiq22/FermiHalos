# -*- coding: utf-8 -*-
"""
@author: Martín Mestre
RAR model.

Metric convention:
g_00 = e^(nu)
g_11 = -e^(lambda)
"""

from scipy.integrate import solve_ivp
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import os
import sys
import time


""" The approach here is to define a function that defines subfunctions needed to compute the right 
hand side of the RAR equations. Then, with a numerical solver, the equations are integrated. Before
the function gives back the output, the astrophysical quantities are reescaled so they have physical 
meaning. The arguments are the four free parameters of the model and four variables indicating if the
statistical mechanical variables and the density of the distribution have to be computed. """

def model(param, deg_var=False, cutoff_var=False, temp_var=False, den_var=False):
    # dark matter particle mass: param[0] - It has to be given in keV
    # degeneracy parameter: param[1]
    # cutoff parameter: param[2]
    # temperature parameter: param[3]
    if (param[0] <= 0.0 or param[3] <= 0.0):
        raise ValueError("The particle mass and the temperature parameter have to be non-zero positive values.")
        
    def fermi(eps, alpha_r, beta_r, eps_r):
        """
        Fermi distribution.

        It also depends on alpha(r), beta(r) and eps(r).
        It has no explicit dependence on r.
        """
        up = 1.0 - np.exp((eps - eps_r)/beta_r)
        down = 1.0 + np.exp((eps - alpha_r)/beta_r)
        return up/down

    def g_rho(eps, alpha_r, beta_r, eps_r):
        """ Integrand for the density. """
        return eps*eps*np.sqrt(eps**2 - 1.0)*fermi(eps, alpha_r, beta_r, eps_r)

    def g_P(eps, alpha_r, beta_r, eps_r):
        """ Integrand for the pressure. """
        return np.sqrt(eps**2 - 1.0)**3*fermi(eps, alpha_r, beta_r, eps_r)

    def density(n_step, alpha_r, beta_r, eps_r):
        """ Density. We handle the case where one component of the eps_r array is <= 1 or very near to 1. """
        if (eps_r > 1.0 + float(n_step)*machine_eps):
            eps = np.linspace(1.0, eps_r, n_step)
            return a*integrate.simps(g_rho(eps, alpha_r, beta_r, eps_r), eps, even='avg')
        else:
            return 0.0

    def pressure(n_step, alpha_r, beta_r, eps_r):
        """ Pressure. We handle the case where one component of the eps_r array is <= 1 or very near to 1. """
        if (eps_r > 1.0 + float(n_step)*machine_eps):
            eps = np.linspace(1.0, eps_r, n_step)
            return b*integrate.simps(g_P(eps, alpha_r, beta_r, eps_r), eps, even='avg')
        else:
            return 0.0

    def border_density(t, u):
        """ Border event function: density. """
        exponential = np.exp(-0.5*(u[1] - nu_0))
        alpha_r = alpha_0*exponential
        beta_r = beta_0*exponential
        eps_r = eps_0*exponential
        rho = density(n_eos, alpha_r, beta_r, eps_r)
        return rho*to_astro - 1.0e-10  # M_sun/pc^3

    border_density.terminal = True
    border_density.direction = -1

    def tov(t, u):
        """
        Tolman-Oppenheimer-Volkoff equation.

        t: independent variable, t = ln(r/R)
        u: 2D array whose components are:
        u[0] = z with z = ln(Psi) and Psi = (M(r)/M)*(R/r)
        u[1] = nu = metric potential
        where:
        M: mass scaling factor
        R: radius scaling factor
        """
        # Computing variables for density and pressure
        exponential = np.exp(-0.5*(u[1] - nu_0))
        alpha_r = alpha_0*exponential
        beta_r = beta_0*exponential
        eps_r = eps_0*exponential
        
        # This is important to avoid NaNs and to define the halo border.
        if (eps_r < 1.0):
            eps_r = 1.0

        # Computing density and pressure:
        rho = density(n_eos, alpha_r, beta_r, eps_r)
        P = pressure(n_eos, alpha_r, beta_r, eps_r)

        # Right hand side of the two differential equations
        d_z = np.exp(2.0*t - u[0])*rho/rho_rel - 1.0
        d_nu = (np.exp(u[0]) + np.exp(2.0*t)*P/(rho_rel*c*c))/(1.0 - np.exp(u[0]))

        return [d_z, d_nu]

    # Set physical constants
    kpc2pc = 1000.0
    cm2kpc = 1.0/3.08567758e21
    cm2pc = cm2kpc*kpc2pc
    g2Msun = 1.0/1.98847e33
    to_astro = g2Msun/(cm2pc**3)
    pi = np.pi
    c = 2.99792458e+10                          # Light speed - cgs
    G = 6.67430e-8                              #  Newton's gravitational constant - cgs
    h = 6.6260755e-27                           # Planck constant - cgs
    g = 2.0
    m_e = 9.1093837015e-28                      # Electron mass - cgs
    ener_e = 510.99895                          # Electron energy - keV
    ener_f = param[0]                           # Dark matter particle mass - keV
    m = (ener_f/ener_e)*m_e                     # g
    rho_rel = (g*m**4/h**3)*(pi*c*c)**1.5       # g/cm^3
    R = c/np.sqrt(8.0*pi*G*rho_rel)             # cm
    M = 4.0*pi*R**3*rho_rel                     # g
    a = 4.0*rho_rel/np.sqrt(pi)                 # g/cm^3
    b = a*c*c/3.0                               # g/(cm s^2)
    n_eos = 2**10 + 1
    tau = 1.0e-16
    min_r = 1.0e-16                             # Minimum value of galactic radius to integrate - kpc
    max_r = 1.0e3                               # Maximum value of galactic radius to integrate - kpc
    machine_eps = np.finfo(float).eps           # Machine epsilon

    # Set initial conditions
    psi_0 = 2.0*tau
    z_0 = np.log(psi_0)
    nu_0 = 2.0*tau
    theta_0 = param[1]                          # Degeneracy parameter
    W_0 = param[2]                              # Cutoff parameter
    beta_0 = param[3]                           # Temperature parameter
    alpha_0 = 1.0 + beta_0*theta_0
    eps_0 = 1.0 + beta_0*W_0
    rho_0 = density(n_eos, alpha_0, beta_0, eps_0) # Central density
    u_0 = [z_0, nu_0]
    t_0 = np.log(np.sqrt(6.0*tau*rho_rel/rho_0))
    t_f = np.log(max_r/cm2kpc/R)

    # Solving the TOV system
    sol = solve_ivp(tov, [t_0, t_f], u_0, method='LSODA', events=(border_density),
                    rtol=5.0e-12, atol=0.0)

    # Defining physical variables
    r = np.exp(sol.t)*R                         # Spherical radius - cm
    z = sol.y[0]
    nu = sol.y[1]
    psi = np.exp(z)
    mass = psi*M/R*r                            # Enclosed mass - g
    
    # Temperature variable
    exponential = np.exp(-0.5*(nu - nu_0))
    beta = beta_0*exponential
        
    # Shift the metric:
    nu_origin = 2.0*np.log(np.sqrt(1.0 - psi[-1])*beta[-1]/beta[0])
    nu = nu + nu_origin*np.ones(len(nu))

    # In astrophysical units
    r = r*cm2kpc                                # kpc
    mass = mass*g2Msun                          # M_sun
    bool_r = (r > min_r)
    r = r[bool_r]
    mass = mass[bool_r]
    nu = nu[bool_r]
    
    # Organizing output
    output = [r, mass, nu]
    
    if deg_var:
        W = (1.0 + beta_0*W_0 - np.exp(nu/2.0))/beta_0
        theta = theta_0 - W_0 + W
        output.append(theta)
        
    if cutoff_var:
        W = (1.0 + beta_0*W_0 - np.exp(nu/2.0))/beta_0
        output.append(W)
        
    if temp_var:
        beta = beta[bool_r]
        output.append(beta)
        
    if den_var:
        mass_spline = InterpolatedUnivariateSpline(r, mass, k=4)
        """Density function, afterwards."""
        deriv = mass_spline.derivative(1)
        output.append(deriv(r)/(4.0*np.pi*r*r))

    return output