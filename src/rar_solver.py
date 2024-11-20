# -*- coding: utf-8 -*-
"""
@author: RAR collaboration
Extended RAR model.

Metric convention:
g_00 = e^(nu)
g_11 = -e^(lambda)
"""

# ================================================= Packages ==================================================== #
import numpy as np
from scipy.integrate import solve_ivp
import scipy.integrate as integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import optimize
# =============================================================================================================== #

""" The approach here is to define a function that defines subfunctions needed to compute the right 
hand side of the RAR equations. Then, with a numerical solver, the equations are integrated. Before
the function gives back the output, the astrophysical quantities are reescaled so they have physical 
meaning. The arguments are the four free parameters of the model, the maximum radii of integration, 
the relative tolerance used by the solver and the number of steps used in the computation of the 
density and pressure as integrals. """

def model(param, maximum_r, relative_tolerance, number_of_steps):
    # dark matter particle mass: param[0] - It has to be given in keV
    # degeneracy parameter: param[1]
    # cutoff parameter: param[2]
    # temperature parameter: param[3]
        
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
    
    def _pressure(n_step, alpha_r, beta_r, eps_r):
        """ Pressure. We handle the case where one component of the eps_r array is <= 1 or very near to 1. """
        # Array of pressures, with zero values.
        p = np.zeros(len(alpha_r))
        # Array of booleans. True are the components satisfying the condition.
        eps_r_bool = eps_r > 1.0 + float(n_step)*machine_eps
        # Indexes of True values.
        true_indexes = np.where(eps_r_bool == True)[0]
        # Filter the values where the condition is satisfyed.
        alpha_r_true, beta_r_true, eps_r_true = alpha_r[true_indexes], beta_r[true_indexes], eps_r[true_indexes]
        # This will generate an array of shape n_step x len(eps_r_true).
        eps = np.linspace(1.0, eps_r_true, n_step)
        #p[true_indexes] = b*quad(g_P, eps[0,:], eps[-1,:], args=(alpha_r_true, beta_r_true, eps_r_true))[0]
        p[true_indexes] = b*integrate.simps(g_P(eps, alpha_r_true, beta_r_true, eps_r_true), eps, axis=0, even='avg')
        return p

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
    kpc2cm = 3.08567758e21
    cm2pc = cm2kpc*kpc2pc
    g2Msun = 1.0/1.98847e33
    Msun2g = 1.98847e33
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
    n_eos = number_of_steps
    tau = 1.0e-16
    min_r = 1.0e-16                             # Minimum value of galactic radius to integrate - kpc
    max_r = maximum_r                           # Maximum value of galactic radius to integrate - kpc
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
    t_0 = 0.5*np.log(6.0*tau*rho_rel/rho_0)
    t_f = np.log(max_r/cm2kpc/R)

    # Solving the TOV system
    sol = solve_ivp(tov, [t_0, t_f], u_0, method='LSODA', events=(border_density),
                    rtol=relative_tolerance, atol=0.0)

    # Defining physical variables
    r = np.exp(sol.t)*R                         # Spherical radius - cm
    z = sol.y[0]
    nu = sol.y[1]
    psi = np.exp(z)
    mass = psi*M/R*r                            # Enclosed mass - g
    
    # Temperature variable and pressure
    exponential = np.exp(-0.5*(nu - nu_0))
    beta = beta_0*exponential
    alpha = alpha_0*exponential
    eps = eps_0*exponential
    P = _pressure(n_eos, alpha, beta, eps)
    
    # Shift the metric:
    nu_origin = 2.0*np.log(np.sqrt(1.0 - psi[-1])*beta[-1]/beta[0])
    nu = nu + nu_origin

    # In astrophysical units
    r = r*cm2kpc                                # kpc
    mass = mass*g2Msun                          # M_sun
    bool_r = (r > min_r)
    r = r[bool_r]
    mass = mass[bool_r]
    nu = nu[bool_r]
    beta = beta[bool_r]
    P = P[bool_r]/Msun2g*kpc2cm                 # Msun/(kpc s^2)

    return r, mass, nu, beta, P

# Constants
G_u = 4.3009e-6          # (km/s)^2*kpc/M_sun
c = 2.99792458e+5        # Light speed - km/s
k = 8.617333262e-8       # Boltzmann constant - keV/K

class Rar():
    """ RAR mass distribution object.
    
        This kind of objects has to be instantiated as Rar(parameters, boolean_variables). These parameters are 
        the four RAR parameters, the boolean variables indicating if additional physical variables have to be 
        computed, and the parameters related with the integration of the equations. For more details see 
        (https://github.com/Santiq22/rar-model)."""
    
    def __init__(self, param, dens_func=False, nu_func=False, lambda_func=False, press_func=False, circ_vel_func=False,
                 accel_func=False, deg_var=False, cutoff_var=False, temp_var=False, chemical_func=False, cutoff_func=False,
                 temperature_func=False, log_dens_slope_func=False, core_func=False, plateau_func=False, maximum_r=1.0e3, 
                 relative_tolerance=5.0e-12, number_of_steps=2**10 + 1):
        
        # Individual parameters
        self.DM_mass, self.theta_0, self.W_0, self.beta_0 = param[0], param[1], param[2], param[3]
        
        # Checking if the DM particle mass or beta_0 are less than 0
        if (self.DM_mass <= 0.0 or self.beta_0 <= 0.0):
            raise ValueError("The particle mass and the temperature parameter have to be non-zero positive values.")
            
        # Checking if number_of_steps is greater than the suggested value
        if (number_of_steps < 2**10 + 1):
            raise ValueError("The number of steps of integration has to be greater than 2^10 + 1 to ensure precision. The value given is {}".format(number_of_steps))
            
        # Call model to solve the RAR equations. The function model returns ndarrys of shape (n,).
        self.r, self.m, self.nu, self.temperature_variable, self.P = model(param, maximum_r, relative_tolerance, number_of_steps)
        
        # Continous mass function. Allows easy computation of derivatives
        self.mass_spline = InterpolatedUnivariateSpline(self.r, self.m, k=4)
        
        # ============================================= Boolean attributes ============================================== #
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
        
        # ===================================== Interpolation of optional variables ===================================== #
        # --- Density
        if (self.dens_func or self.log_dens_slope_func or self.plateau_func):
            # Continous density function. Allows easy computation of derivatives
            self.density_spline = (self.mass_spline).derivative(1)
        
        # --- Metric potential
        if self.nu_func:
            # Continous metric potential function. Allows easy computation of derivatives
            self.nu_spline = InterpolatedUnivariateSpline(self.r, self.nu, k=4)
        
        # --- Pressure
        if (self.press_func or self.circ_vel_func or self.core_func or self.plateau_func):
            # Continous pressure function. Allows easy computation of derivatives
            self.P_spline = InterpolatedUnivariateSpline(self.r, self.P, k=4)
        
        # --- Degeneracy variable
        if self.deg_var:
            # Define the cutoff variable array
            self.cutoff_variable = (1.0 + self.beta_0*self.W_0)/self.beta_0*(1.0 - (1.0 - 2.0*G_u*self.m[-1]/(c**2*self.r[-1]))**(-0.5)*np.exp(self.nu/2.0))
            
            # Define the degeneracy variable array
            self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
            
            # Continous degeneracy variable function. Allows easy computation of derivatives
            self.degeneracy_variable_spline = InterpolatedUnivariateSpline(self.r, self.degeneracy_variable, k=4)
            
        # --- Cutoff variable
        if self.cutoff_var:
            # Define the cutoff variable array
            self.cutoff_variable = (1.0 + self.beta_0*self.W_0)/self.beta_0*(1.0 - (1.0 - 2.0*G_u*self.m[-1]/(c**2*self.r[-1]))**(-0.5)*np.exp(self.nu/2.0))
            
            # Continous cutoff variable function. Allows easy computation of derivatives
            self.cutoff_variable_spline = InterpolatedUnivariateSpline(self.r, self.cutoff_variable, k=4)
            
        # --- Temperature variable
        if self.temp_var:
            # Continous temperature variable function. Allows easy computation of derivatives
            self.temperature_variable_spline = InterpolatedUnivariateSpline(self.r, self.temperature_variable, k=4)
            
        # --- Chemical potential
        if self.chemical_func:
            if not self.deg_var:
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0)/self.beta_0*(1.0 - (1.0 - 2.0*G_u*self.m[-1]/(c**2*self.r[-1]))**(-0.5)*np.exp(self.nu/2.0))
                
                # Define the degeneracy variable array
                self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
                
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the chemical potential array
                self.chemical_potential = k*self.degeneracy_variable*self.temperature
                
                # Continous chemical potential function. Allows easy computation of derivatives
                self.mu_spline = InterpolatedUnivariateSpline(self.r, self.chemical_potential, k=4)
            else:
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the chemical potential array
                self.chemical_potential = k*self.degeneracy_variable*self.temperature
                
                # Continous chemical potential function. Allows easy computation of derivatives
                self.mu_spline = InterpolatedUnivariateSpline(self.r, self.chemical_potential, k=4)
                
        # --- Cutoff function
        if self.cutoff_func:
            if not self.deg_var:
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0)/self.beta_0*(1.0 - (1.0 - 2.0*G_u*self.m[-1]/(c**2*self.r[-1]))**(-0.5)*np.exp(self.nu/2.0))
                
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the cutoff function array
                self.cutoff = k*self.cutoff_variable*self.temperature
                
                # Continous cutoff function. Allows easy computation of derivatives
                self.e_c_spline = InterpolatedUnivariateSpline(self.r, self.cutoff, k=4)
            else:
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
                # Define the cutoff function array
                self.cutoff = k*self.cutoff_variable*self.temperature
                
                # Continous cutoff function. Allows easy computation of derivatives
                self.e_c_spline = InterpolatedUnivariateSpline(self.r, self.cutoff, k=4)
        
        # --- Temperature
        if self.temperature_func:
            # Define the temperature array
            self.temperature = self.DM_mass*self.temperature_variable/k
            
            # Continous temperature function. Allows easy computation of derivatives
            self.T_spline = InterpolatedUnivariateSpline(self.r, self.temperature, k=4)
            
        # --- Logarithmic slope of the density profile
        if self.log_dens_slope_func:
            # Continous second derivative function. Allows easy computation of derivatives
            self.second_derivative_of_mass = (self.mass_spline).derivative(2)
        # =============================================================================================================== #
        
    # =============================================== Static methods ================================================ #
    @staticmethod
    def _mass(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.mass_spline(r), self.mass_spline(r_max))
        
    @staticmethod
    def _density(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.density_spline(r)/(4.0*np.pi*r*r), 0.0)
        
    @staticmethod
    def _lambda_potential(self, r):
        return -np.log(1.0 - 2.0*G_u*self._mass(self, r)/(c*c*r))
    
    @staticmethod
    def _pressure(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.P_spline(r), 0.0)
    
    @staticmethod
    def _dnu_dr(self, r):
        return 1.0/r*((8.0*np.pi*G_u/c**4*self._pressure(self, r)*r**2 + 1.0)/(1.0 - 2.0*G_u*self._mass(self, r)/(c**2*r)) - 1.0)
    
    @staticmethod
    def _circular_velocity(self, r):
        return np.sqrt(0.5*c*c*r*self._dnu_dr(self, r))
    # =============================================================================================================== #
    
    # ============================================== Instance methods =============================================== #
    def mass(self, r):
        r_max = self.r[-1]
        return np.where(r < r_max, self.mass_spline(r), self.mass_spline(r_max))
    
    def density(self, r):
        if not self.dens_func:
            raise NameError("The 'density' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.density_spline(r)/(4.0*np.pi*r*r), 0.0)
    
    def metric_potential(self, r):
        if not self.nu_func:
            raise NameError("The 'metric_potential' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.nu_spline(r), -self._lambda_potential(self, r))
            
    def lambda_potential(self, r):
        if not self.lambda_func:
            raise NameError("The 'lambda_potential' method is not defined.")
        else:
            return -np.log(1.0 - 2.0*G_u*self._mass(self, r)/(c*c*r))
        
    def pressure(self, r):
        if not self.press_func:
            raise NameError("The 'pressure' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.P_spline(r), 0.0)
    
    def circular_velocity(self, r):
        if not self.circ_vel_func:
            raise NameError("The 'circular_velocity' method is not defined.")
        else:
            return np.sqrt(0.5*c*c*r*self._dnu_dr(self, r))

    def acceleration(self, x, y, z):
        if not self.accel_func:
            raise NameError("The 'acceleration' method is not defined.")
        else:
            r = np.sqrt(x*x + y*y + z*z)
            return -G_u*self._mass(self, r)/(r**3)*np.array([x, y, z])
        
    def theta(self, r):
        if not self.deg_var:
            raise NameError("The 'theta' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.degeneracy_variable_spline(r), 0.0)
        
    def W(self, r):
        if not self.cutoff_var:
            raise NameError("The 'W' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.cutoff_variable_spline(r), 0.0)
        
    def beta(self, r):
        if not self.temp_var:
            raise NameError("The 'beta' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.temperature_variable_spline(r), 0.0)
        
    def mu(self, r):
        if not self.chemical_func:
            raise NameError("The 'mu' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.mu_spline(r), 0.0)
        
    def e_c(self, r):
        if not self.cutoff_func:
            raise NameError("The 'e_c' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.e_c_spline(r), 0.0)
        
    def T(self, r):
        if not self.temperature_func:
            raise NameError("The 'T' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.T_spline(r), 0.0)
        
    def logarithmic_density_slope(self, r):
        if not self.log_dens_slope_func:
            raise NameError("The 'logarithmic_density_slope' method is not defined.")
        else:
            return 2.0 - 1.0/(4.0*np.pi*r*self._density(self, r))*self.second_derivative_of_mass(r)
    
    def core(self):
        if not self.core_func:
            raise NameError("The 'core' method is not defined.")
        else:
            v = self._circular_velocity(self, self.r)
            """ Old way of doing it:
            from scipy.signal import argrelextrema
            arg_max = argrelextrema(v, np.greater)
            r_cand = self.r[arg_max[0][0]]"""
            arg_max = np.argmax(v)
            r_cand = self.r[arg_max]
            bounds = np.array([r_cand*0.5, r_cand*1.5])
            r_core = optimize.fminbound(lambda r : -self._circular_velocity(self, r), bounds[0],
                                        bounds[1], xtol=0.5e-12, maxfun=1000)
            m_core = self._mass(self, r_core)
            return r_core, m_core
        
    def plateau(self):
        if not self.plateau_func:
            raise NameError("The 'plateau' method is not defined.")
        else:
            v = self._circular_velocity(self, self.r)
            arg_max = np.argmax(v)
            r_new = self.r[arg_max:]
            v_new = self._circular_velocity(self, r_new)
            arg_min = np.argmin(v_new)
            r_cand = r_new[arg_min]
            bounds = np.array([r_cand*0.5, r_cand*1.5])
            r_plateau = optimize.fminbound(lambda r : self._circular_velocity(self, r), bounds[0],
                                        bounds[1], xtol=0.5e-12, maxfun=1000)
            rho_plateau = self._density(self, r_plateau)
            return r_plateau, rho_plateau
    # =============================================================================================================== #