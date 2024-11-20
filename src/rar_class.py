# -*- coding: utf-8 -*-
"""
@author: RAR collaboration
Extended RAR model.

Metric convention:
g_00 = e^(nu)
g_11 = -e^(lambda)
"""

# ================================================= Packages ==================================================== #
from model import model
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import fminbound
# =============================================================================================================== #

# --- Constants
G_u = 4.3009e-6          # Newton's constant - (km/s)^2*kpc/M_sun
c = 2.99792458e+5        # Light speed - km/s
k = 8.617333262e-8       # Boltzmann's constant - keV/K

class Rar:
    """ RAR mass distribution object.
    
        This kind of objects has to be instantiated as Rar(parameters, boolean_variables). These parameters are 
        the four RAR parameters, the boolean variables indicating if additional physical variables have to be 
        computed, and the parameters related with the integration of the equations. For more details see 
        (https://github.com/Santiq22/rar-model)."""
    
    def __init__(self, param, dens_func=False, nu_func=False, lambda_func=False, press_func=False, circ_vel_func=False,
                 accel_func=False, deg_var=False, cutoff_var=False, temp_var=False, chemical_func=False, cutoff_func=False,
                 temperature_func=False, log_dens_slope_func=False, core_func=False, plateau_func=False, maximum_r=1.0e3, 
                 relative_tolerance=5.0e-12, number_of_steps=2**10 + 1):
        
        # ======================================== Numerical instance attributes ======================================== #
        self.DM_mass = param[0]                              # Dark matter particle mass
        self.theta_0 = param[1]                              # Degeneracy parameter
        self.W_0 = param[2]                                  # Cut-off parameter
        self.beta_0 = param[3]                               # Temperature parameter
        self.maximum_r = maximum_r                           # Maximum radii of integration
        self.relative_tolerance = relative_tolerance         # Relative tolerance used to solve the equations
        self.number_of_steps = number_of_steps               # Number of steps of integration used in pressure and density
        # =============================================================================================================== #
        
        # ========================================= Boolean instance attributes ========================================= #
        self.dens_func = dens_func                           # Density
        self.nu_func = nu_func                               # Metric potential
        self.lambda_func = lambda_func                       # Lambda function
        self.press_func = press_func                         # Pressure
        self.circ_vel_func = circ_vel_func                   # General relativistic circular velocity
        self.accel_func = accel_func                         # Newtonian gravitational field
        self.deg_var = deg_var                               # Degeneracy variable
        self.cutoff_var = cutoff_var                         # Cutoff variable
        self.temp_var = temp_var                             # Temperature variable
        self.chemical_func = chemical_func                   # Chemical potential
        self.cutoff_func = cutoff_func                       # Cutoff function
        self.temperature_func = temperature_func             # Temperature function
        self.log_dens_slope_func = log_dens_slope_func       # Logarithmic slope of the density profile
        self.core_func = core_func                           # DM core
        self.plateau_func = plateau_func                     # DM plateau
        # =============================================================================================================== #
        
        # ========================================== Checks before integration ========================================== #
        # Checking if the DM particle mass or beta_0 are less than 0
        if (self.DM_mass <= 0.0 or self.beta_0 <= 0.0):
            raise ValueError("The particle mass and the temperature parameter have to be non-zero positive values.")
            
        # Checking if number_of_steps is greater than the suggested value
        if (self.number_of_steps < 2**10 + 1):
            raise ValueError("The number of steps of integration has to be greater than 2^10 + 1 to ensure precision. The value given is {}".format(self.number_of_steps))
        # =============================================================================================================== #
        
        # ========================================= Computation of the solutions ======================================== #    
        # Call model to solve the RAR equations. The function model returns ndarrys of shape (n,).
        if not (self.press_func or self.circ_vel_func or self.core_func or self.plateau_func):
            self.r, self.m, self.nu, self.temperature_variable = model(param, maximum_r=self.maximum_r, relative_tolerance=self.relative_tolerance, 
                                                                       number_of_steps=self.number_of_steps, press_func=False)
        else:
            self.r, self.m, self.nu, self.temperature_variable, self.P = model(param, maximum_r=self.maximum_r, relative_tolerance=self.relative_tolerance, 
                                                                               number_of_steps=self.number_of_steps, press_func=True)
        
        # Continous mass function. Allows easy computation of derivatives
        self.__mass_spline = InterpolatedUnivariateSpline(self.r, self.m, k=4)
        # =============================================================================================================== #
        
        # ===================================== Interpolation of optional variables ===================================== #
        # --- Density
        if (self.dens_func or self.log_dens_slope_func or self.plateau_func):
            # Continous density function. Allows easy computation of derivatives
            self.__density_spline = (self.__mass_spline).derivative(1)
        
        # --- Metric potential
        if self.nu_func:
            # Continous metric potential function. Allows easy computation of derivatives
            self.__nu_spline = InterpolatedUnivariateSpline(self.r, self.nu, k=4)
        
        # --- Pressure
        if (self.press_func or self.circ_vel_func or self.core_func or self.plateau_func):
            # Continous pressure function. Allows easy computation of derivatives
            self.__P_spline = InterpolatedUnivariateSpline(self.r, self.P, k=4)
            
        # --- Cutoff variable
        if self.cutoff_var:
            # Define the cutoff variable array
            self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp(self.nu/2.0))/self.beta_0
            
            # Continous cutoff variable function. Allows easy computation of derivatives
            self.__cutoff_variable_spline = InterpolatedUnivariateSpline(self.r, self.cutoff_variable, k=4)
        
        # --- Degeneracy variable
        if self.deg_var:
            if not self.cutoff_var:
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp(self.nu/2.0))/self.beta_0
            
                # Define the degeneracy variable array
                self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
            
                # Continous degeneracy variable function. Allows easy computation of derivatives
                self.__degeneracy_variable_spline = InterpolatedUnivariateSpline(self.r, self.degeneracy_variable, k=4)
            else:
                # Define the degeneracy variable array
                self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
            
                # Continous degeneracy variable function. Allows easy computation of derivatives
                self.__degeneracy_variable_spline = InterpolatedUnivariateSpline(self.r, self.degeneracy_variable, k=4)
            
        # --- Temperature variable
        if self.temp_var:
            # Continous temperature variable function. Allows easy computation of derivatives
            self.__temperature_variable_spline = InterpolatedUnivariateSpline(self.r, self.temperature_variable, k=4)
            
        # --- Cutoff function
        if self.cutoff_func:
            if not (self.cutoff_var or self.deg_var):
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp(self.nu/2.0))/self.beta_0
                
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the cutoff function array
                self.cutoff = k*self.cutoff_variable*self.temperature
                
                # Continous cutoff function. Allows easy computation of derivatives
                self.__e_c_spline = InterpolatedUnivariateSpline(self.r, self.cutoff, k=4)
            else:
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the cutoff function array
                self.cutoff = k*self.cutoff_variable*self.temperature
                
                # Continous cutoff function. Allows easy computation of derivatives
                self.__e_c_spline = InterpolatedUnivariateSpline(self.r, self.cutoff, k=4)
            
        # --- Chemical potential
        if self.chemical_func:
            if not (self.cutoff_var or self.deg_var or self.cutoff_func):
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp(self.nu/2.0))/self.beta_0
                
                # Define the degeneracy variable array
                self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
                
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the chemical potential array
                self.chemical_potential = k*self.degeneracy_variable*self.temperature
                
                # Continous chemical potential function. Allows easy computation of derivatives
                self.__mu_spline = InterpolatedUnivariateSpline(self.r, self.chemical_potential, k=4)
            elif not (self.deg_var):
                # Define the degeneracy variable array
                self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
                
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the chemical potential array
                self.chemical_potential = k*self.degeneracy_variable*self.temperature
                
                # Continous chemical potential function. Allows easy computation of derivatives
                self.__mu_spline = InterpolatedUnivariateSpline(self.r, self.chemical_potential, k=4)
            else:
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the chemical potential array
                self.chemical_potential = k*self.degeneracy_variable*self.temperature
                
                # Continous chemical potential function. Allows easy computation of derivatives
                self.__mu_spline = InterpolatedUnivariateSpline(self.r, self.chemical_potential, k=4)
        
        # --- Temperature
        if self.temperature_func:
            # Define the temperature array
            self.temperature = self.DM_mass*self.temperature_variable/k
            
            # Continous temperature function. Allows easy computation of derivatives
            self.__T_spline = InterpolatedUnivariateSpline(self.r, self.temperature, k=4)
            
        # --- Logarithmic slope of the density profile
        if self.log_dens_slope_func:
            # Continous second derivative function. Allows easy computation of derivatives
            self.__second_derivative_of_mass = (self.__mass_spline).derivative(2)
        # =============================================================================================================== #

    # ============================================ Instance private methods ============================================= #
    def __mass(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.__mass_spline(r), self.__mass_spline(r_max))
        
    def __density(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.__density_spline(r)/(4.0*np.pi*r*r), 0.0)
        
    def __lambda_potential(self, r):
        return -np.log(1.0 - 2.0*G_u*self.__mass(r)/(c*c*r))
    
    def __pressure(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.__P_spline(r), 0.0)
    
    def __dnu_dr(self, r):
        return 1.0/r*((8.0*np.pi*G_u/c**4*self.__pressure(r)*r**2 + 1.0)/(1.0 - 2.0*G_u*self.__mass(r)/(c**2*r)) - 1.0)
    
    def __circular_velocity(self, r):
        return np.sqrt(0.5*c*c*r*self.__dnu_dr(r))
    # ==================================================================================================================== #
    
    # ============================================== Instance public methods ============================================= #
    def mass(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.__mass_spline(r), self.__mass_spline(r_max))
    
    def density(self, r):
        if not self.dens_func:
            raise NameError("The 'density' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__density_spline(r)/(4.0*np.pi*r*r), 0.0)
    
    def metric_potential(self, r):
        if not self.nu_func:
            raise NameError("The 'metric_potential' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__nu_spline(r), -self.__lambda_potential(r))
            
    def lambda_potential(self, r):
        if not self.lambda_func:
            raise NameError("The 'lambda_potential' method is not defined.")
        else:
            return -np.log(1.0 - 2.0*G_u*self.__mass(r)/(c*c*r))
        
    def pressure(self, r):
        if not self.press_func:
            raise NameError("The 'pressure' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__P_spline(r), 0.0)
    
    def circular_velocity(self, r):
        if not self.circ_vel_func:
            raise NameError("The 'circular_velocity' method is not defined.")
        else:
            return np.sqrt(0.5*c*c*r*self.__dnu_dr(r))

    def acceleration(self, x, y, z):
        if not self.accel_func:
            raise NameError("The 'acceleration' method is not defined.")
        else:
            r = np.sqrt(x*x + y*y + z*z)
            return -G_u*self.__mass(r)/(r**3)*np.array([x, y, z])
        
    def theta(self, r):
        if not self.deg_var:
            raise NameError("The 'theta' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__degeneracy_variable_spline(r), 0.0)
        
    def W(self, r):
        if not self.cutoff_var:
            raise NameError("The 'W' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__cutoff_variable_spline(r), 0.0)
        
    def beta(self, r):
        if not self.temp_var:
            raise NameError("The 'beta' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__temperature_variable_spline(r), 0.0)
        
    def mu(self, r):
        if not self.chemical_func:
            raise NameError("The 'mu' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__mu_spline(r), 0.0)
        
    def e_c(self, r):
        if not self.cutoff_func:
            raise NameError("The 'e_c' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__e_c_spline(r), 0.0)
        
    def T(self, r):
        if not self.temperature_func:
            raise NameError("The 'T' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__T_spline(r), 0.0)
        
    def logarithmic_density_slope(self, r):
        if not self.log_dens_slope_func:
            raise NameError("The 'logarithmic_density_slope' method is not defined.")
        else:
            return 2.0 - 1.0/(4.0*np.pi*r*self.__density(r))*self.__second_derivative_of_mass(r)
    
    def core(self):
        if not self.core_func:
            raise NameError("The 'core' method is not defined.")
        else:
            v = self.__circular_velocity(self.r)
            """ Old way of doing it:
            from scipy.signal import argrelextrema
            arg_max = argrelextrema(v, np.greater)
            r_cand = self.r[arg_max[0][0]]"""
            arg_max = np.argmax(v)
            r_cand = self.r[arg_max]
            bounds = np.array([r_cand*0.5, r_cand*1.5])
            r_core = fminbound(lambda r : -self.__circular_velocity(r), bounds[0], bounds[1], xtol=0.5e-12, maxfun=1000)
            m_core = self.__mass(r_core)
            return r_core, m_core
        
    def plateau(self):
        if not self.plateau_func:
            raise NameError("The 'plateau' method is not defined.")
        else:
            v = self.__circular_velocity(self.r)
            arg_max = np.argmax(v)
            r_new = self.r[arg_max:]
            v_new = self.__circular_velocity(r_new)
            arg_min = np.argmin(v_new)
            r_cand = r_new[arg_min]
            bounds = np.array([r_cand*0.5, r_cand*1.5])
            r_plateau = fminbound(lambda r : self.__circular_velocity(r), bounds[0], bounds[1], xtol=0.5e-12, maxfun=1000)
            rho_plateau = self.__density(r_plateau)
            return r_plateau, rho_plateau
    # ==================================================================================================================== #