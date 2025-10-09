# Extended RAR model
Repository containing the programs needed to solve the RAR model's equations.

## Files description

### `rar_class.py`:

The code defines a class whose name is `Rar`. It has to be instantiated as:

```python
halo_object = Rar(parameters, dens_func=False, nu_func=False, lambda_func=False, press_func=False, 
                 circ_vel_func=False, accel_func=False, deg_var=False, cutoff_var=False, 
                 temp_var=False, chemical_func=False, cutoff_func=False, temperature_func=False, 
                 log_dens_slope_func=False, core_func=False, plateau_func=False, maximum_r=1.0e3, 
                 relative_tolerance=5.0e-12, number_of_steps=2**10 + 1)
```

where the `parameters` variable is a numpy array object of shape (4,), whose components are (in this order): the dark matter particle mass in keV, the degeneracy parameter, the cutoff parameter, and the temperature parameter (the last three adimensional). The boolean variables are used as flags to compute astrophysical and statistical mechanical variables. To do so, change `False` to `True`. See Rar's attributes section to further details on the variables passed to the Rar class. 

Once the class `Rar` is instantiated, it automatically calls `model`, a function that solves the RAR model equations. This function is called as:

```python
model(parameters, maximum_r, relative_tolerance, number_of_steps, press_func)
```
where `parameters` is the array used in the instance of the class.

The `model` function defines several subfunctions needed to compute the right-hand side of the TOV equations. These subfunctions include a Fermi-Dirac-like distribution function, two integrands for computing the density and pressure, and two functions for computing the density and pressure themselves.

The TOV equations are solved using the `solve_ivp` function from the *scipy.integrate module*. The right-hand side of the TOV equations is computed using the function called `tov`. The solution of the TOV equations is then re-scaled to obtain the physical quantities, including the radius, enclosed mass, metric potential, and pressure of the dark matter halo. The optional attributes of the `Rar` class then compute the other physical variables using the outputs of the `model` function.

#### Rar's attributes

All the following attributes are *boolean instance attributes* and *set up instance attributes*.

- `dens_func`: Boolean variable that enables the computation of the density profile of the distribution. The default value is `False`.
- `nu_func`: Boolean variable that enables the computation of the metric potential. The default value is `False`.
- `lambda_func`: Boolean variable that enables the computation of the lambda potential. The default value is `False`.
- `press_func`: Boolean variable that enables the computation of the pressure profile. The default value is `False`.
- `circ_vel_func`: Boolean variable that enables the computation of the circular velocity profile. The default value is `False`.
- `accel_func`: Boolean variable that enables the computation of the Newtonian gravitational field exerted by the dark matter halo. The default value is `False`.
- `deg_var`: Boolean variable that enables the computation of the degeneracy variable. The default value is `False`.
- `cutoff_var`: Boolean variable that enables the computation of the cutoff variable. The default value is `False`.
- `temp_var`: Boolean variable that enables the computation of the temperature variable. The default value is `False`.
- `chemical_func`: Boolean variable that enables the computation of the chemical potential. The default value is `False`.
- `cutoff_func`: Boolean variable that enables the computation of the cutoff energy function. The default value is `False`.
- `temperature_func`: Boolean variable that enables the computation of the temperature function. The default value is `False`.
- `log_dens_slope_func`: Boolean variable that enables the computation of the logarithmic density slope function. The default value is `False`.
- `core_func`: Boolean variable that enables the computation of the radii of the dark matter core and its mass. The default value is `False`.
- `plateau_func`: Boolean variable that enables the computation of the radii of the dark matter plateau and its density. The default value is `False`.
- `maximum_r`: _(float)_ Maximum radius of integration in $kpc$.
- `relative_tolerance`: _(float)_ Relative tolerance used by the integrator to solve the equations.
- `number_of_steps`: _(int)_ Number of steps used to integrate the density and pressure used to compute the right-hand side of the differential equations. We strongly suggest that the value of `number_of_steps` is greater than the minimum value $`2^{10} + 1`$ to ensure precision at the time of computing the solutions.

In addition, there are some instance attributes representing physical quantities, which are:

- `DM_mass` [$`keV`$]: Dark matter particle mass.
- `theta_0`: Degeneracy parameter $\theta_{0}$ of the system.
- `W_0`: Cutoff parameter $W_{0}$ of the system.
- `beta_0`: Temperature parameter $\beta_{0}$ of the system.
- `r` [$`kpc`$]: Array of the radius where the solution was computed. It is a numpy ndarray of shape (n,). Available by default.
- `m` [$`M_{\odot}`$]: Array of enclosed masses at the radius given in `r`. It is a numpy ndarray of shape (n,). Available by default.
- `nu`: Array of the metric potentials (dimensionless) at the radius given in `r`. It is a numpy ndarray of shape (n,). Available by default.
- `nu_0`: Value of the metric potential at the center of the distribution. Available by default.
- `P` [$`M_{\odot}/(kpc\ s^{2})`$]: Array of pressures at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `press_func` is `True`.
- `degeneracy_variable`: Array of values of the degeneracy variable at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `deg_var` or `chemical_func` is `True`.
- `cutoff_variable`: Array of values of the cutoff variable at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `deg_var`, `cutoff_var`, `chemical_func` or `cutoff_func` is `True`.
- `temperature_variable`: Array of values of the temperature variable at the radius given in `r`. It is a numpy ndarray of shape (n,). Available by default.
- `chemical_potential` [$`keV`$]: Array of values of the chemical potential at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `chemical_func` is `True`.
- `cutoff` [$`keV`$]: Array of values of the cutoff energy function at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `cutoff_func` is `True`.
- `temperature` [$`K`$]: Array of values of the temperature function at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `chemical_func`, `cutoff_func` or `temperature_func` is `True`.

#### Rar's methods

The only one method that is computed by default is the enclosed mass of the dark matter distribution. To enable the computation of other astrophysical and statistical mechanical variables just change the boolean attributes to `True` while instantiating the object. The methods defined in the class `Rar` are:

- `mass` [$`M_{\odot}`$]: Enclosed mass of the distribution defined for all spherical radii. It is defined as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  mass(r)=
      \begin{cases}
          M(r) & \text{if } r < r_{\mathrm{max}},\\
          M(r_{\textrm{max}}) & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
$`M(r)`$ is the integrated enclosed mass. It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `density` [$`M_{\odot}/kpc^{3}`$]: Mass density of the distribution defined for all spherical radii. It is computed when `dens_func=True`. It is defined as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  density(r)=
      \begin{cases}
          \frac{\mathrm{d}M(r)}{\mathrm{d}r}\frac{1}{4\pi r^{2}} & \text{if } r < r_{\mathrm{max}},\\
          0 & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `metric_potential`: Metric potential function (dimensionless) of the fermionic distribution. It is computed when `nu_func=True`. It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `lambda_potential`: Lambda function (dimensionless) used in the definition of the spacetime metric. It is computed when `lambda_func=True`. It is defined for all the spherical radius as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  lambda\_potential(r)=
      \begin{cases}
          -\mathrm{ln}\left[1 - \frac{2GM(r)}{c^{2}r}\right] & \text{if } r < r_{\mathrm{max}},\\
          -\mathrm{ln}\left[1 - \frac{2GM(r_{\textrm{max}})}{c^{2}r}\right] & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `pressure` [$`M_{\odot}/(kpc\ s^{2})`$]: Pressure of the fermionic distribution. It is computed when `press_func=True`. It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `circular_velocity` [$`km/s`$]: General relativistic expression of the circular velocity of the mass distribution defined for all spherical radii. It is computed when `circ_vel_func=True`. It is defined as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  circular\_velocity(r)=
      \begin{cases}
          \sqrt{\frac{1}{2}c^{2}\left[\left(\frac{8\pi G}{c^{4}}P(r)r^{2} + 1\right)\left(1 - \frac{2GM(r)}{c^{2}r}\right)^{-1} - 1\right]} & \text{if } r < r_{\mathrm{max}},\\
          \sqrt{\frac{1}{2}c^{2}\left[\left(\frac{8\pi G}{c^{4}}P(r)r^{2} + 1\right)\left(1 - \frac{2GM(r_{\textrm{max}})}{c^{2}r}\right)^{-1} - 1\right]} & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `acceleration` [$`(km/s)^{2}kpc^{-1}`$]: Newtonian gravitational field of the distribution, defined for all the spherical radius. It is computed when `accel_func=True`. It takes as input three floating numbers, the cartesian coordinates of the reference system. They have to be given as three separate variables. It is defined as a piecewise function of 2 parts, whose expression is (supposing the center of mass of the distribution is in the origin):
```math
\begin{equation}
  acceleration(x, y, z)=
      \begin{cases}
          -\frac{GM(r)}{r^{3}}\vec{r} & \text{if } r < r_{\mathrm{max}},\\
          -\frac{GM(r_{\textrm{max}})}{r^{3}}\vec{r} & \text{if } r \geq r_{\mathrm{max}},
      \end{cases}
\end{equation}
```
where $\vec{r} = (x, y, z)$ and $r = ||\vec{r}||$. It returns a numpy ndarray of shape (3,). It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `theta`: Degeneracy variable (dimensionless). It is computed when `deg_var=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `W`: Cutoff variable (dimensionless). It is computed when `cutoff_var=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `beta`: Temperature variable (dimensionless). It is computed when `temp_var=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `mu` [$`keV`$]: Chemical potential. It is computed when `chemical_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `e_c` [$`keV`$]: Cutoff energy function. It is computed when `cutoff_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `T` [$`K`$]: Temperature function. It is computed when `temperature_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`.
- `logarithmic_density_slope`: Logarithmic slope of the density profile (dimensionless). This dimensionless quantity is computed when `log_dens_slope_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `scipy.interpolate.InterpolatedUnivariateSpline`. It is computed up to the maximum radii of integration, which can be computed by doing `halo_object.r[-1]`. Its mathematical expression is:
```math
\begin{equation}
  \gamma(r)= -\frac{d\mathrm{ln}\rho(r)}{d\mathrm{ln}r} = 2 - \frac{1}{4\pi r\rho(r)}\frac{d^{2}M(r)}{dr^{2}}
\end{equation}
```
- `core`: It is computed when `core_func=True`. When this function with no argument is called, it returns the *core radius* $r_{\mathrm{core}}$ (in $kpc$) and the *mass of the core* (in $M_{\odot}$) of the distribution, computed as $M(r_{\mathrm{core}})$. They are returned as two separate floating numbers in a tuple ($r_{\mathrm{core}}$, $`M(r_{\mathrm{core}})`$).
- `plateau`: It is computed when `plateau_func=True`. When this function with no argument is called, it returns the *plateau radius* $r_{\mathrm{plateau}}$ (in $kpc$) and the *density of the plateau* (in $M_{\odot}/kpc^{3}$) of the distribution, computed as $\rho(r_{\mathrm{plateau}})$. They are returned as two separate floating numbers in a tuple ($r_{\mathrm{plateau}}$, $`\rho(r_{\mathrm{plateau}})`$).

## Dependencies

The dependencies of the program are:
- Numpy
- Scipy

## Downloading the program

#### Using GitHub

To download the program, just clone the repository to the directory where you would like to have it saved. For more details on how to clone a GitHub repository see [Cloning a repository](https://docs.github.com/es/repositories/creating-and-managing-repositories/cloning-a-repository).

#### Using Conda

To install this programm using conda as a package manager you can do:

```python
conda ...
```

#### Using PyPI

To install this package from the PyPI repositories just write in a terminal:

```python
pip ...
```

## How to use it?

To use the program, it has to be instantiated a `Rar` object as indicated in the box above. To do so, the `Rar` class has to be imported previously as:

```python
from rar_class import Rar
```

Take into account that the file where it will be the importing of the `Rar` class has to be in the same directory as the `rar_class.py` file. Otherwise, it has to be changed the path to where `rar_class.py` is saved, going from the directory where the object is going to be instantiated. This is easily done by writing the following lines in the code where it is going to be the instance:

```python
import sys
import os

# Real path from the root directory to where rar_class.py is
path_rar = os.path.realpath(os.path.join(os.path.dirname(__file__), 'path/', 'to/', 'rar_class.py folder'))

# Insert the path location of rar_class.py into the list of folders giving by sys.path
sys.path.insert(0, path_rar)

# Import the class
from rar_class import Rar
```

## Example of use

There are plenty of aplications where a RAR dark matter halo can play a role (see the main paper for further details). Here we show a classical application of galactic dynamics, the fitting of a rotation curve.

```python
# =============================== Packages =================================== #
import numpy as np
from scipy.optimize import differential_evolution
import os
from rar_class import Rar

# Module containing two classes representing an arbitrary disk and an arbitrary bulge
from mass_dist import ExponentialDisk, deVaucouleurs
# ============================================================================ #

# ============================= Observations ================================= #
""" The observations were taken from Sofue, PASJ, 67, 4, 2015, 75 """
# ---  radius (kpc), velocity (km/s), standard deviation (km/s)
observations = np.loadtxt('M31_grand_rotation_curve.txt')
radius = observations[:,0]
velocities = observations[:,1]
deviations = observations[:,2]

# The first value of deviations is 0 on the list
deviations[0] = deviations[0] + 4.0          # km/s

# Maximum radii to work with
r_gal_max = 250.0                            # kpc
radius = radius[radius < r_gal_max]
max_idx = len(radius)
velocities = velocities[:max_idx]
deviations = deviations[:max_idx]
# ============================================================================ #

# =========================== Initial parameters ============================= #
# - RAR seed to do the fitting
m_DM = 70.0                                  # keV

# Degeneracy parameter
theta_0 = 36.0

# Cut-off parameter
W_0 = 61.0

# Temperature parameter
beta_0 = 4.3e-3

# Mass of the core and its error
m_core = 1.0e8                               # M_Sun
m_core_err = 0.1*m_core                      # M_Sun

# - de Vaucouleurs bulge parameters - Best fits of Sofue 2015
M_b = 1.7e10                                 # M_sun
a_b = 1.35                                   # kpc
kappa = 7.6695
eta = 22.665

# - Exponential disk parameters - Best fits of Sofue 2015
M_d = 1.8e11                                 # M_sun
a_d = 5.28                                   # kpc

# Initial parameter vector
p_0 = np.array([theta_0, W_0, beta_0, M_b, a_b, M_d, a_d])
# ============================================================================ #

# ========================= Bounds for parameters ============================ #
""" Limits to constrain the parameter space and hence ease the best-fitting procedure """
bounds = ((35.5, 37.5),
          (60.0, 62.0),
          (4.0e-3, 4.7e-3),
          (0.5*M_b, 1.5*M_b),
          (a_b - 0.5, a_b + 0.5),
          (1.6e11, 1.9e11),
          (4.6, 6.0))
# ============================================================================ #

# ====================== Reduced chi squared function ======================== #
def chi_2_red(p):
    # Bulge object
    bulge = deVaucouleurs(np.array([p[3], p[4], kappa, eta]))

    # Disk object
    disk = ExponentialDisk(p[-2:])

    # Dark matter halo object
    halo = Rar(np.array([m_DM, p[0], p[1], p[2]]), circ_vel_func=True, core_func=True)

    # Total circular velocity
    v = np.sqrt(bulge.circular_velocity(radius)**2 + disk.circular_velocity(radius)**2 + 
                halo.circular_velocity(radius)**2)

    # Reduced chi squared function
    chi = (np.sum((v - velocities)**2/deviations**2) + (halo.core()[1] - m_core)**2/m_core_err**2)/(len(velocities) - len(p))

    return chi
# ============================================================================ #

# ============================ Fitting function ============================== #
def fit(cores):
    opt = differential_evolution(chi_2_red, bounds, strategy='best1bin', maxiter=150, 
                                 popsize=75, recombination=0.4, mutation=(0.2, 0.5),
                                 tol=1.0e-10, atol=0.0, disp=True, polish=True, 
                                 x0=p_0, seed=10, workers=cores)
    
    solution = opt.x
    np.savetxt('best_fit_parameters_for_M31_70kev.txt', solution)
    return solution
# ============================================================================ #

# =========================== Fitting procedure ============================== #
# Number of CPU cores to run differential_evolution in parallel
cores = 10
best_fit_parameters = fit(cores)
print("Fitting procedure completed\n")
print("The best fit parameters are: ", best_fit_parameters)
# ============================================================================ #
```

This code will generate a dark matter mass distribution with a radial profile as follows:

![density_profile](figures/density_profile.png)

and a total circular velocity profile as:

![circular_velocity_profile](figures/circular_velocity_profile.png)

## License


This code is subject to MIT License. In addition, it is asked the following conditions in order to use the code:

- The potential works published in scientific journals that have used this code have to cite this collaboration. The official publication of the code can be found [here]().
- In case any issue or bug is found, please, report it as an issue on the GitHub page of the repository. This way, we can work to solve it as soon as possible.
