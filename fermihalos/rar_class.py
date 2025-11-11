# -*- coding: utf-8 -*-
"""
@author: RAR collaboration
Extended RAR model.

Metric convention:
g_00 = e^(nu)
g_11 = -e^(lambda)
"""

# ================================================= Packages ==================================================== #
from fermihalos.model import model
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import fminbound
from typing import Tuple
# =============================================================================================================== #

# --- Constants
G_u = 4.3009e-6          # Newton's constant - (km/s)^2*kpc/M_sun
c = 2.99792458e+5        # Light speed - km/s
k = 8.617333262e-8       # Boltzmann's constant - keV/K

class Rar:
    """ 
    Extended RAR mass distribution object.
    
    This class instantiates an object of the ``Rar`` class, which represents a fermionic dark matter halo 
    following the extended RAR model. The constructor of this class will integrate the differential equations 
    defining the halo, based on the input parameters stated below.
    
    Parameters
    ----------
    param: ndarray
        Array containing the dark matter particle mass, the degeneracy parameter, the cutoff variable and the
        temperature variable, respectively. 
    dens_func: bool, optional, default=False
        Flag allowing the computation of the density profile.
    nu_func: bool, optional, default=False
        Flag allowing the computation of the metric potential.
    particles_func: bool, optional, default=False
        Flag allowing the computation of the enclosed particle number.
    lambda_func: bool, optional, default=False
        Flag allowing the computation of the lambda potential.
    press_func: bool, optional, default=False
        Flag allowing the computation of the pressure profile.
    n_func: bool, optional, default=False
        Flag allowing the computation of the particle number density.
    circ_vel_func: bool, optional, default=False
        Flag allowing the computation of the circular velocity profile.
    accel_func: bool, optional, default=False
        Flag allowing the computation of the acceleration field.
    deg_var: bool, optional, default=False
        Flag allowing the computation of the degeneracy variable.
    cutoff_var: bool, optional, default=False
        Flag allowing the computation of the cutoff variable.
    temp_var: bool, optional, default=False
        Flag allowing the computation of the temperature variable.
    chemical_func: bool, optional, default=False
        Flag allowing the computation of the chemical potential function.
    cutoff_func: bool, optional, default=False
        Flag allowing the computation of the cutoff function.
    temperature_func: bool, optional, default=False
        Flag allowing the computation of the temperature function.
    log_dens_slope_func: bool, optional, default=False
        Flag allowing the computation of the logarithmic density slope function.
    core_func: bool, optional, default=False
        Flag allowing the computation of the core function.
    plateau_func: bool, optional, default=False
        Flag allowing the computation of the plateau function.
    maximum_r: float, optional, default=1.0e3
        Maximum radius of integration in kpc.
    relative_tolerance: float, optional, default=5.0e-12
        Relative tolerance used by the integrator to solve the equations.
    number_of_steps: int, optional, default=2\*\*10 + 1
        Number of steps used to integrate the density, pressure, and particle number density used to compute the 
        right-hand side of the differential equations. We strongly suggest that the value of `number_of_steps` 
        is greater than the minimum value 2**10 + 1 to ensure precision at the time of computing the solutions.
        
    Attributes
    ----------
    DM_mass: float
        Float variable representing the dark matter particle mass in :math:`keV/c^{2}`.
    theta_0: float
        Float variable representing the degeneracy parameter :math:`\\theta_{0}` of the system.
    W_0: float
        Float variable representing the cutoff parameter :math:`W_{0}` of the system.
    beta_0: float
        Float variable representing the temperature parameter :math:`\\beta_{0}` of the system.
    dens_func: bool
        Boolean variable that enables the computation of the density profile of the distribution.
    nu_func: bool
        Boolean variable that enables the computation of the metric potential.
    particles_func: bool
        Boolean variable that enables the computation of the enclosed particle number.
    lambda_func: bool
        Boolean variable that enables the computation of the lambda potential.
    press_func: bool
        Boolean variable that enables the computation of the pressure profile.
    n_func: bool
        Boolean variable that enables the computation of the particle number density.
    circ_vel_func: bool
        Boolean variable that enables the computation of the circular velocity profile.
    accel_func: bool
        Boolean variable that enables the computation of the Newtonian gravitational field exerted by the dark matter halo.
    deg_var: bool
        Boolean variable that enables the computation of the degeneracy variable.
    cutoff_var: bool
        Boolean variable that enables the computation of the cutoff variable.
    temp_var: bool
        Boolean variable that enables the computation of the temperature variable.
    chemical_func: bool
        Boolean variable that enables the computation of the chemical potential.
    cutoff_func: bool
        Boolean variable that enables the computation of the cutoff energy function.
    temperature_func: bool
        Boolean variable that enables the computation of the temperature function.
    log_dens_slope_func: bool
        Boolean variable that enables the computation of the logarithmic density slope function.
    core_func: bool
        Boolean variable that enables the computation of the radii of the dark matter core and its mass.
    plateau_func: bool
        Boolean variable that enables the computation of the radii of the dark matter plateau and its density.
    maximum_r: float
        Float variable representing the maximum radius of integration in :math:`kpc`.
    relative_tolerance: float
        Float variable representing the relative tolerance used by the integrator to solve the equations.
    number_of_steps: int
        Integer variable representing the number of steps used to integrate the density, pressure, and particle 
        number density used to compute the right-hand side of the differential equations. We strongly suggest 
        that the value of ``number_of_steps`` is greater than the minimum value :math:`2^{10} + 1` to ensure 
        precision at the time of computing the solutions.
    r: float ndarray of shape [n,]
        Array of the radius where the solution was computed in :math:`kpc`.
    m: float ndarray of shape [n,]
        Array of enclosed masses at the radius given in ``r``. The units are :math:`M_{\odot}`.
    nu: float ndarray of shape [n,]
        Array of metric potentials (dimensionless) at the radius given in ``r``.
    N: float ndarray of shape [n,]
        Array of enclosed particles number at the radius given in ``r``. 
    nu_0: float
        Value of the metric potential at the center of the distribution, :math:`\\nu_{0}`.
    P: float ndarray of shape [n,]
        Array of pressures at the radius given in ``r``. The units are :math:`M_{\odot}/(kpc\ s^{2})`. It is 
        only available if ``press_func`` is ``True``.
    n: float ndarray of shape [n,]
        Array of particles number densities at the radius given in ``r``. The units are [:math:`kpc^{-3}`]. 
        It is only available if ``n_func`` is ``True``.
    degeneracy_variable: float ndarray of shape [n,]
        Array of values of the degeneracy variable (dimensionless) at the radius given in ``r``. It is only 
        available if ``deg_var`` or ``chemical_func`` is ``True``.
    cutoff_variable: float ndarray of shape [n,]
        Array of values of the cutoff variable (dimensionless) at the radius given in ``r``. It is only 
        available if ``deg_var``, ``cutoff_var``, ``chemical_func`` or ``cutoff_func`` is ``True``.
    temperature_variable: float ndarray of shape [n,]
        Array of values of the temperature variable (dimensionless) at the radius given in ``r``.
    chemical_potential: float ndarray of shape [n,]
        Array of values of the chemical potential at the radius given in ``r``. The units are :math:`keV`. 
        It is only available if ``chemical_func`` is ``True``.
    cutoff: float ndarray of shape [n,]
        Array of values of the cutoff energy function at the radius given in ``r``. The units are :math:`keV`.
        It is only available if ``cutoff_func`` is ``True``.
    temperature: float ndarray of shape [n,]
        Array of values of the temperature function at the radius given in ``r``. The units are :math:`K`.
        it is only available if ``chemical_func``, ``cutoff_func`` or ``temperature_func`` is ``True``.
        
    Notes
    -----
    For more details see `README <https://github.com/Santiq22/FermiHalos/README.md>`_.
    """
    
    def __init__(self, param: np.ndarray, 
                 dens_func: bool = False, 
                 nu_func: bool = False,
                 particles_func: bool = False, 
                 lambda_func: bool = False, 
                 press_func: bool = False,
                 n_func: bool = False,
                 circ_vel_func: bool = False,
                 accel_func: bool = False, 
                 deg_var: bool = False, 
                 cutoff_var: bool = False, 
                 temp_var: bool = False, 
                 chemical_func: bool = False, 
                 cutoff_func: bool = False,
                 temperature_func: bool = False, 
                 log_dens_slope_func: bool = False, 
                 core_func: bool = False, 
                 plateau_func: bool = False, 
                 maximum_r: float = 1.0e3, 
                 relative_tolerance: float = 5.0e-12, 
                 number_of_steps: int = 2**10 + 1):
        
        # ======================================== Numerical instance attributes ======================================== #
        self.DM_mass = param[0]                              # Dark matter particle mass - keV/c^2
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
        self.particles_func = particles_func                 # Particle number
        self.lambda_func = lambda_func                       # Lambda function
        self.press_func = press_func                         # Pressure
        self.n_func = n_func                                 # Particle number density
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
            raise ValueError(f"The number of steps of integration has to be greater than 2^10 + 1 to ensure precision. The value given is {self.number_of_steps}")
        # =============================================================================================================== #
        
        # ========================================= Computation of the solutions ======================================== #    
        # Call model to solve the RAR equations. The function model returns ndarrys of shape (n,).
        self.r, self.m, self.nu, self.N, self.temperature_variable, self.nu_0, self.P, self.n = model(
            param, 
            maximum_r = self.maximum_r, 
            relative_tolerance = self.relative_tolerance, 
            number_of_steps = self.number_of_steps, 
            press_func = (self.press_func or self.circ_vel_func or self.core_func or self.plateau_func), 
            n_func = self.n_func)
        
        # Continous mass function. Allows easy computation of derivatives
        self.__mass_spline = InterpolatedUnivariateSpline(self.r, self.m, k = 4)
        # =============================================================================================================== #
        
        # ===================================== Interpolation of optional variables ===================================== #
        # --- Density
        if (self.dens_func or self.log_dens_slope_func or self.plateau_func):
            # Continous density function. Allows easy computation of derivatives
            self.__density_spline = (self.__mass_spline).derivative(1)
        
        # --- Metric potential
        if self.nu_func:
            # Continous metric potential function. Allows easy computation of derivatives
            self.__nu_spline = InterpolatedUnivariateSpline(self.r, self.nu, k = 4)
            
        # --- Particle number
        if self.particles_func:
            # Continous particle number function. Allows easy computation of derivatives
            self.__N_spline = InterpolatedUnivariateSpline(self.r, self.N, k = 4)
        
        # --- Pressure
        if (self.press_func or self.circ_vel_func or self.core_func or self.plateau_func):
            # Continous pressure function. Allows easy computation of derivatives
            self.__P_spline = InterpolatedUnivariateSpline(self.r, self.P, k = 4)
            
        # --- Particle number density
        if self.n_func:
            # Continous particle number density function. Allows easy computation of derivatives
            self.__n_spline = InterpolatedUnivariateSpline(self.r, self.n, k = 4)
            
        # --- Cutoff variable
        if self.cutoff_var:
            # Define the cutoff variable array
            self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp((self.nu - self.nu_0)/2.0))/self.beta_0
            
            # Continous cutoff variable function. Allows easy computation of derivatives
            self.__cutoff_variable_spline = InterpolatedUnivariateSpline(self.r, self.cutoff_variable, k = 4)
        
        # --- Degeneracy variable
        if self.deg_var:
            # Check if the cut-off variable was computed before
            if not self.cutoff_var:
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp((self.nu - self.nu_0)/2.0))/self.beta_0
                
            # Define the degeneracy variable array
            self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
            
            # Continous degeneracy variable function. Allows easy computation of derivatives
            self.__degeneracy_variable_spline = InterpolatedUnivariateSpline(self.r, self.degeneracy_variable, k = 4)
            
        # --- Temperature variable
        if self.temp_var:
            # Continous temperature variable function. Allows easy computation of derivatives
            self.__temperature_variable_spline = InterpolatedUnivariateSpline(self.r, self.temperature_variable, k = 4)
            
        # --- Cutoff function
        if self.cutoff_func:
            # Define the temperature array
            self.temperature = self.DM_mass*self.temperature_variable/k
            
            # Check if the cut-off variable was computed before
            if not (self.cutoff_var or self.deg_var):
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp((self.nu - self.nu_0)/2.0))/self.beta_0
                
            # Define the cutoff function array
            self.cutoff = k*self.cutoff_variable*self.temperature
                
            # Continous cutoff function. Allows easy computation of derivatives
            self.__e_c_spline = InterpolatedUnivariateSpline(self.r, self.cutoff, k = 4)
            
        # --- Chemical potential
        if self.chemical_func:
            # Check if the cut-off variable was computed before
            if not (self.cutoff_var or self.deg_var or self.cutoff_func):
                # Define the cutoff variable array
                self.cutoff_variable = (1.0 + self.beta_0*self.W_0 - np.exp((self.nu - self.nu_0)/2.0))/self.beta_0
                
            # Check if the degeneracy variable was computed before
            if not self.deg_var:
                # Define the degeneracy variable array
                self.degeneracy_variable = self.theta_0 - self.W_0 + self.cutoff_variable
                
            # Check if the temperature function was computed before
            if not self.cutoff_func:
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
                
            # Define the chemical potential array
            self.chemical_potential = k*self.degeneracy_variable*self.temperature
                
            # Continous chemical potential function. Allows easy computation of derivatives
            self.__mu_spline = InterpolatedUnivariateSpline(self.r, self.chemical_potential, k = 4)
        
        # --- Temperature
        if self.temperature_func:
            # Check if the temperature function was computed before
            if not self.cutoff_func:
                # Define the temperature array
                self.temperature = self.DM_mass*self.temperature_variable/k
            
            # Continous temperature function. Allows easy computation of derivatives
            self.__T_spline = InterpolatedUnivariateSpline(self.r, self.temperature, k = 4)
            
        # --- Logarithmic slope of the density profile
        if self.log_dens_slope_func:
            # Continous second derivative function. Allows easy computation of derivatives
            self.__second_derivative_of_mass = (self.__mass_spline).derivative(2)
        # =============================================================================================================== #
        
    def __repr__(self):
        return f'<{type(self).__name__} object at 0x{id(self):x}>'

    # ============================================ Instance private methods ============================================= #
    def __mass(self, r: float | np.ndarray) -> np.ndarray:
        r_max = self.r[-1]
        return np.where(r < r_max, self.__mass_spline(r), self.__mass_spline(r_max))
        
    def __density(self, r: float | np.ndarray) -> np.ndarray:
        r_max = self.r[-1]
        return np.where(r < r_max, self.__density_spline(r)/(4.0*np.pi*r*r), 0.0)
        
    def __lambda_potential(self, r: float | np.ndarray) -> np.ndarray:
        return -np.log(1.0 - 2.0*G_u*self.__mass(r)/(c*c*r))
    
    def __pressure(self, r: float | np.ndarray) -> np.ndarray:
        r_max = self.r[-1]
        return np.where(r < r_max, self.__P_spline(r), 0.0)
    
    def __dnu_dr(self, r: float | np.ndarray) -> np.ndarray:
        return 1.0/r*((8.0*np.pi*G_u/c**4*self.__pressure(r)*r**2 + 1.0)/(1.0 - 2.0*G_u*self.__mass(r)/(c**2*r)) - 1.0)
    
    def __circular_velocity(self, r: float | np.ndarray) -> np.ndarray:
        return np.sqrt(0.5*c*c*r*self.__dnu_dr(r))
    # ==================================================================================================================== #
    
    # ============================================== Instance public methods ============================================= #
    def mass(self, 
             r: float | np.ndarray) -> np.ndarray:
        """
        .. _mass-function:
        
        Enclosed mass profile in :math:`M_{\odot}`.
        
        This function computes the enclosed mass profile of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`. For 
        :math:`r \geq R`, being :math:`R` the radial size of the halo, this function takes a constant value that is the mass 
        of the self-gravtiating system.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the enclosed mass
            
        Returns
        -------
        out: ndarray
            Enclosed mass at r
            
        Notes
        -----
        Setting:
        
        .. math::
        
            M(r) = \int_{0}^{r}4\pi r'^{\ 2}\\rho(r')dr',
            
        where :math:`\\rho(r)` is the :ref:`density <density-function>` of the distribution, the method is defined as a piecewise function of 2 parts, whose expression is:
        
        .. math::
        
            mass(r)=
                \\begin{cases}
                    M(r) & \\text{if } r < R,\\\\
                    M(r_{\mathrm{max}}) & \\text{if } r \geq R.
                \\end{cases}
        """
        
        r_max = self.r[-1]
        
        return np.where(r < r_max, self.__mass_spline(r), self.__mass_spline(r_max))
    
    def density(self, 
                r: float | np.ndarray) -> np.ndarray:
        """
        .. _density-function:
        
        Density profile in :math:`M_{\odot}/kpc^{3}`.
        
        This function computes the density profile of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`.
        For :math:`r \geq R`, being :math:`R` the radial size of the halo, this function takes a constant value that is 0.
        This is due to the cut-off behaviour that naturaly has the model at the edge of the system.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the density
            
        Returns
        -------
        out: ndarray
            Density at r
            
        Raises
        ------
        NameError
            If ``self.dens_func`` is ``False``
            
        Notes
        -----
        Setting:
        
        .. math::
            
            \\rho(r) = \\frac{2m}{h^{3}}\int f_{\mathrm{RAR}}\\left(\epsilon(\\vec{p}),r\\right) \\left(1+\\frac{\epsilon(\\vec{p})}{mc^{2}} \\right) d^{3}p,
        
        where :math:`c` is the speed of light, :math:`h` is the Planck's constant, :math:`m` is the dark matter particle mass, 
        :math:`f_{\mathrm{RAR}}` is the coarse-grained distribution function characterising the fermionic model, 
        and :math:`\epsilon` is the relativistic energy of a dark matter particle without its rest mass, the method is defined as a 
        piecewise function of 2 parts, whose expression is:
        
        .. math::
        
            density(r)=
                \\begin{cases}
                    \\rho(r) & \\text{if } r < R,\\\\
                    0 & \\text{if } r \geq R.
                \\end{cases}
        """
        
        if not self.dens_func:
            raise NameError("The 'density' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__density_spline(r)/(4.0*np.pi*r*r), 0.0)
    
    def metric_potential(self, 
                         r: float | np.ndarray) -> np.ndarray:
        """
        .. _metric-potential-function:
        
        Metric potential (dimensionless). 
        
        This function computes the metric potential of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`. 
        For :math:`r \geq R`, being :math:`R` the radial size of the halo, this function takes the form of 
        :math:`-\lambda(r)` the :ref:`lambda potential <lambda-potential-function>`, as imposed by the continuity 
        condition of the metric potential.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the metric potential
            
        Returns
        -------
        out: ndarray
            Metric potential at r
            
        Raises
        ------
        NameError
            If ``self.nu_func`` is ``False``
            
        Notes
        -----
        This function has no analytical expresion. Instead, it is the solution of the differential equation:
        
        .. math::
        
            \\frac{d\\nu(r)}{dr} = \\frac{1}{r}\\left[\\left(1-\\frac{2GM(r)}{c^{2}r}\\right)^{-1}\\left(\\frac{8\pi G}{c^{4}}P(r)r^{2}+1\\right)-1\\right],
        
        where :math:`M(r)` is the :ref:`enclosed mass <mass-function>` and :math:`P(r)` is the :ref:`pressure <pressure-function>` 
        of the system.
        """
        
        if not self.nu_func:
            raise NameError("The 'metric_potential' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__nu_spline(r), -self.__lambda_potential(r))
        
    def particle_number(self, 
                        r: float | np.ndarray) -> np.ndarray:
        """
        Enclosed particle number (dimensionless).
        
        This function computes the enclosed particle number of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`. 
        For :math:`r \geq R`, being :math:`R` the radial size of the halo, this function takes a constant value that is the 
        number of particles composing the system.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the enclosed particle number
            
        Returns
        -------
        out: ndarray
            Enclosed particle number at r
            
        Raises
        ------
        NameError
            If ``self.particles_func`` is ``False``
            
        Notes
        -----
        Defining:
        
        .. math::
        
            N(r) = \int_{0}^{r}4\pi r'^{\ 2}e^{\\lambda(r')/2}n(r')dr',
         
        being :math:`\\lambda(r)` the :ref:`lambda potential <lambda-potential-function>` and :math:`n(r)` the 
        :ref:`particle number density <particle-number-density-function>`, then the expresion for the particle 
        number method is:
        
        .. math::
        
            particle\_number(r)=
                \\begin{cases}
                    N(r) & \\text{if } r < R,\\\\
                    N(r_{\mathrm{max}}) & \\text{if } r \geq R.
                \end{cases}
        """
        
        if not self.particles_func:
            raise NameError("The 'particle_number' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__N_spline(r), self.__N_spline(r_max))
            
    def lambda_potential(self, 
                         r: float | np.ndarray) -> np.ndarray:
        """
        .. _lambda-potential-function:
        
        Lambda potential (dimensionless). 
        
        This function computes the lambda potential of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`. 
        For :math:`r \geq R`, being :math:`R` the radial size of the halo, this function depending on the enclosed mass 
        declines as :math:`\propto 1/r`.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the lambda potential
            
        Returns
        -------
        out: ndarray
            Lambda potential at r
            
        Raises
        ------
        NameError
            If ``self.lambda_func`` is ``False``
            
        Notes
        -----
        The expression defining this function is:
        
        .. math::
        
            lambda\_potential(r)=
                \\begin{cases}
                    -\mathrm{ln}\\left[1 - \\frac{2GM(r)}{c^{2}r}\\right] & \\text{if } r < R,\\\\
                    -\mathrm{ln}\\left[1 - \\frac{2GM(r_{\\textrm{max}})}{c^{2}r}\\right] & \\text{if } r \geq R,
                \\end{cases}
                
        where :math:`M(r)` is the :ref:`enclosed mass <mass-function>` and :math:`P(r)` is the :ref:`pressure <pressure-function>`
        of the fermionic distribution.
        """
        
        if not self.lambda_func:
            raise NameError("The 'lambda_potential' method is not defined.")
        else:
            return -np.log(1.0 - 2.0*G_u*self.__mass(r)/(c*c*r))
        
    def pressure(self, 
                 r: float | np.ndarray) -> np.ndarray:
        """
        .. _pressure-function:
        
        Pressure profile in :math:`M_{\odot}/(kpc\ s^{2})`.
        
        This function computes the pressure profile of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`. 
        For :math:`r \geq R`, being R the radial size of the halo, this function takes a constant value that is 0, as 
        expected outside a distribution, since there is no presence of matter exerting pressure.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the pressure profile
            
        Returns
        -------
        out: ndarray
            Pressure profile at r
            
        Raises
        ------
        NameError
            If ``self.press_func`` is ``False``
            
        Notes
        -----
        Setting:
        
        .. math::
        
            P(r) = \\frac{2}{3h^{3}}\int f_{\mathrm{RAR}}\\left(\epsilon(\\vec{p}),r\\right)\epsilon(\\vec{p}) \\frac{1+\epsilon(\\vec{p})/2mc^{2}}{1+\epsilon(\\vec{p})/mc^{2}}\ d^{3}p,
        
        where :math:`c` is the speed of light, :math:`h` is the Planck's constant, :math:`m` is the dark matter particle mass, 
        :math:`f_{\mathrm{RAR}}` is the coarse-grained distribution function characterising the fermionic model, 
        and :math:`\epsilon` is the relativistic energy of a dark matter particle without its rest mass, the 
        method is defined as a piecewise function of 2 parts whose expresion is:
        
        .. math::
        
            pressure(r)=
                \\begin{cases}
                    P(r) & \\text{if } r < R,\\\\
                    0 & \\text{if } r \geq R.
                \end{cases}
        """
        
        if not self.press_func:
            raise NameError("The 'pressure' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__P_spline(r), 0.0)
        
    def n(self, 
          r: float | np.ndarray) -> np.ndarray:
        """
        .. _particle-number-density-function:
        
        Particle number density in :math:`kpc^{-3}`.
        
        This function computes the particle number density of the fermionic dark matter distribution which 
        is the solution of the differential equations for the given free parameters of the model. It is 
        defined for every :math:`r > 0`. For :math:`r \geq R`, being :math:`R` the radial size of the halo, 
        this function takes a constant value that is 0, as expected outside the system where there are no
        presence of particles.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the particle number density
            
        Returns
        -------
        out: ndarray
            Particle number density at r
            
        Raises
        ------
        NameError
            If ``self.n_func`` is ``False``
            
        Notes
        -----
        The expression defining this quantity is:
        
        .. math::
        
            n(r)=
                \\begin{cases}
                    \\frac{1}{h^{3}}\int f_{\mathrm{RAR}}(\epsilon(\\vec{p}), r)d^{3}p & \\text{if } r < R,\\\\
                    0 & \\text{if } r \geq R,
                \\end{cases}

        where :math:`h` is the Planck's constant, :math:`f_{\mathrm{RAR}}` is the coarse-grained distribution function that characterizes 
        the fermionic model, and :math:`\epsilon` is the relativistic energy of a dark matter particle without its rest mass.
        """
        
        if not self.n_func:
            raise NameError("The 'n' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__n_spline(r), 0.0)
    
    def circular_velocity(self, 
                          r: float | np.ndarray) -> np.ndarray:
        """
        General relativistic expression of the circular velocity profile in :math:`km/s`.
        
        This function computes the circular velocity profile of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the circular velocity profile
            
        Returns
        -------
        out: ndarray
            Circular velocity profile density at r
            
        Raises
        ------
        NameError
            If ``self.circ_vel_func`` is ``False``
            
        Notes
        -----
        The expression defining this function is:
        
        .. math::
        
            circular\_velocity(r)=
                \\begin{cases}
                    \sqrt{\\frac{1}{2}c^{2}\\left[\\left(\\frac{8\pi G}{c^{4}}P(r)r^{2} + 1\\right)\\left(1 - \\frac{2GM(r)}{c^{2}r}\\right)^{-1} - 1\\right]} & \\text{if } r < R,\\\\
                    \sqrt{\\frac{1}{2}c^{2}\\left[\\left(\\frac{8\pi G}{c^{4}}P(r)r^{2} + 1\\right)\\left(1 - \\frac{2GM(r_{\mathrm{max}})}{c^{2}r}\\right)^{-1} - 1\\right]} & \\text{if } r \geq R,
                \end{cases}
                
        where :math:`c` is the speed of light, :math:`G` is the Newtonian gravitational constant, :math:`P(r)` is the 
        :ref:`pressure <pressure-function>` of the system, and :math:`M(r)` its :ref:`enclosed mass <mass-function>`.
        """
        
        if not self.circ_vel_func:
            raise NameError("The 'circular_velocity' method is not defined.")
        else:
            return np.sqrt(0.5*c*c*r*self.__dnu_dr(r))

    def acceleration(self, 
                     x: float, 
                     y: float, 
                     z: float) -> np.ndarray:
        """
        Newtonian gravitational field in :math:`(km/s)^{2}\ kpc^{-1}`. 
        
        This function computes the acceleration field of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for a given position 
        vector :math:`(x, y, z)`, supposing the origin of the reference system is in the center of the distribution.
        
        Parameters
        ----------
        x: float
            x coordinate where to compute the acceleration field
        y: float
            y coordinate where to compute the acceleration field
        z: float
            z coordinate where to compute the acceleration field
            
        Returns
        -------
        out: ndarray of shape (3,)
            Acceleration field for a given position vector :math:`(x, y, z)`
            
        Raises
        ------
        NameError
            If ``self.accel_func`` is ``False``
            
        Notes
        -----
        The expression defining the gravitational field is:
        
        .. math::

            acceleration(x, y, z)=
                \\begin{cases}
                    -\\frac{GM(r)}{r^{3}}\\vec{r} & \\text{if } r < R,\\\\
                    -\\frac{GM(r_{\mathrm{max}})}{r^{3}}\\vec{r} & \\text{if } r \geq R,
                \end{cases}
                
        where :math:`\\vec{r} = (x, y, z)`, :math:`r = ||\\vec{r}||`, :math:`G` is the Newtonian gravitational constant, and
        :math:`M(r)` is the :ref:`enclosed mass <mass-function>` of the system.
        """
        
        if not self.accel_func:
            raise NameError("The 'acceleration' method is not defined.")
        else:
            r = np.sqrt(x*x + y*y + z*z)
            return -G_u*self.__mass(r)/(r**3)*np.array([x, y, z])
        
    def theta(self, 
              r: float | np.ndarray) -> np.ndarray:
        """
        Degeneracy variable (dimensionless).
        
        This function computes the degeneracy variable of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`. 
        For :math:`r \geq R`, being :math:`R` the radial size of the halo, this function takes a constant value that is 0.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the degeneracy variable
            
        Returns
        -------
        out: ndarray
            Degeneracy variable at r
            
        Raises
        ------
        NameError
            If ``self.deg_var`` is ``False``
            
        Notes
        -----
        This function has the form:
        
        .. math::
        
            theta(r)=
                \\begin{cases}
                    \\theta_{0} - W_{0} + W(r) & \\text{if } r < R,\\\\
                    0 & \\text{if } r \geq R,
                \\end{cases}
        
        where, :math:`\\theta_{0}` is the degeneracy parameter, :math:`W_{0}` is the cut-off parameter, and 
        :math:`W(r)` is the :ref:`cut-off variable <W-function>.
        """
        
        if not self.deg_var:
            raise NameError("The 'theta' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__degeneracy_variable_spline(r), 0.0)
        
    def W(self, 
          r: float | np.ndarray) -> np.ndarray:
        """
        .. _W-function:
        
        Cutoff variable (dimensionless).
        
        This function computes the cutoff variable of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every :math:`r > 0`. 
        For :math:`r \geq R`, being :math:`R` the radial size of the halo, this function takes a constant value that is 0.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the cutoff variable
            
        Returns
        -------
        out: ndarray
            Cutoff variable at r
            
        Raises
        ------
        NameError
            If ``self.cutoff_var`` is ``False``
            
        Notes
        -----
        This function has the form:
        
        .. math::
        
            W(r) = 
                \\begin{cases}
                    \\frac{1 + \\beta_{0}W_{0} - e^{(\\nu(r) - \\nu_{0})/2}}{\\beta_{0}} & \\text{if } r < R,\\\\
                    0 & \\text{if } r \geq R,
                \\end{cases}
                
        where :math:`\\nu(r)` is the :ref:`metric potential <metric-potential-function>`, :math:`W_{0}` is the cut-off
        parameter, :math:`\\beta_{0}` is the temperature parameter, and :math:`\\nu_{0}` is the value of the metric 
        potential at the origin of the distribution fulfilling the continuity condition with the Schwarzschild metric.
        """
        
        if not self.cutoff_var:
            raise NameError("The 'W' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__cutoff_variable_spline(r), 0.0)
        
    def beta(self, 
             r: float | np.ndarray) -> np.ndarray:
        """
        Temperature variable. 
        
        This function computes the temperature variable of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every r > 0. For r >= R, 
        being R the radial size of the halo, this function takes a constant value that is 0.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the temperature variable
            
        Returns
        -------
        ndarray
            Temperature variable at r
            
        Raises
        ------
        NameError
            If self.temp_var is False
        """
        
        if not self.temp_var:
            raise NameError("The 'beta' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__temperature_variable_spline(r), 0.0)
        
    def mu(self, 
           r: float | np.ndarray) -> np.ndarray:
        """
        Chemical potential function. 
        
        This function computes the chemical potential function of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every r > 0. For r >= R, 
        being R the radial size of the halo, this function takes a constant value that is 0.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the chemical potential function
            
        Returns
        -------
        ndarray
            Chemical potential function at r
            
        Raises
        ------
        NameError
            If self.chemical_func is False
        """
        
        if not self.chemical_func:
            raise NameError("The 'mu' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__mu_spline(r), 0.0)
        
    def e_c(self, 
            r: float | np.ndarray) -> np.ndarray:
        """
        Cutoff energy function. 
        
        This function computes the cutoff energy function of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every r > 0. For r >= R, 
        being R the radial size of the halo, this function takes a constant value that is 0.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the cutoff energy function
            
        Returns
        -------
        ndarray
            Cutoff energy function at r
            
        Raises
        ------
        NameError
            If self.cutoff_func is False
        """
        
        if not self.cutoff_func:
            raise NameError("The 'e_c' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__e_c_spline(r), 0.0)
        
    def T(self, 
          r: float | np.ndarray) -> np.ndarray:
        """
        Temperature function. 
        
        This function computes the temperature function of the fermionic dark matter distribution which is the solution 
        of the differential equations for the given free parameters of the model. It is defined for every r > 0. For r >= R, 
        being R the radial size of the halo, this function takes a constant value that is 0.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the temperature function
            
        Returns
        -------
        ndarray
            Temperature function at r
            
        Raises
        ------
        NameError
            If self.temperature_func is False
        """
        
        if not self.temperature_func:
            raise NameError("The 'T' method is not defined.")
        else:
            r_max = self.r[-1]
            return np.where(r < r_max, self.__T_spline(r), 0.0)
        
    def logarithmic_density_slope(self, 
                                  r: float | np.ndarray) -> np.ndarray:
        """
        Logarithmic density slope function. 
        
        This function computes the logarithmic density slope function of the fermionic dark matter distribution which is the 
        solution of the differential equations for the given free parameters of the model. It is defined for every r > 0. 
        Warning: for r >= R, being R the radial size of the halo, this function takes value arbitrarily large.
        
        Parameters
        ----------
        r: float or ndarray
            Radius where to compute the logarithmic density slope function
            
        Returns
        -------
        ndarray
            Logarithmic density slope function at r
            
        Raises
        ------
        NameError
            If self.log_dens_slope_func is False
        """
        
        if not self.log_dens_slope_func:
            raise NameError("The 'logarithmic_density_slope' method is not defined.")
        else:
            return 2.0 - 1.0/(4.0*np.pi*r*self.__density(r))*self.__second_derivative_of_mass(r)
    
    def core(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Core function. 
        
        This function computes the radius and mass of the core of the fermionic dark matter distribution which is the 
        solution of the differential equations for the given free parameters of the model. It does so by using an 
        optimization algorithm.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        tuple(ndarray, ndarray)
            Core radius and core mass, respectively
            
        Raises
        ------
        NameError
            If self.core_func is False
        """
        
        if not self.core_func:
            raise NameError("The 'core' method is not defined.")
        else:
            v = self.__circular_velocity(self.r)
            # Old way of doing it:
            # from scipy.signal import argrelextrema
            # arg_max = argrelextrema(v, np.greater)
            # r_cand = self.r[arg_max[0][0]]"""
            arg_max = np.argmax(v)
            r_cand = self.r[arg_max]
            bounds = np.array([r_cand*0.5, r_cand*1.5])
            r_core = fminbound(lambda r : -self.__circular_velocity(r), 
                               bounds[0], 
                               bounds[1], 
                               xtol = 0.5e-12, 
                               maxfun = 1000)
            m_core = self.__mass(r_core)
            
            return (r_core, m_core)
        
    def plateau(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Plateau function. 
        
        This function computes the radius and density of the plateau of the fermionic dark matter distribution which is the 
        solution of the differential equations for the given free parameters of the model. It does so by using an 
        optimization algorithm.
        
        Parameters
        ----------
        None
            
        Returns
        -------
        tuple(ndarray, ndarray)
            Radius and density of the plateau, respectively
            
        Raises
        ------
        NameError
            If self.plateau_func is False
        """
        
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
            r_plateau = fminbound(lambda r : self.__circular_velocity(r), 
                                  bounds[0], 
                                  bounds[1], 
                                  xtol = 0.5e-12, 
                                  maxfun = 1000)
            rho_plateau = self.__density(r_plateau)
            
            return (r_plateau, rho_plateau)
    # ==================================================================================================================== #
    
if __name__ == '__main__':
    p = np.array([200.0, 35.0, 69.0, 1.0e-4])
    halo = Rar(p)