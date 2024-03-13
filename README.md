# Extended RAR model
Repository containing the programs needed to solve the RAR model's equations.

## Files description

### `rar_solver.py`:

The code defines a class whose name is `Rar`. It has to be instantiated as:

```python
halo_object = Rar(parameters, dens_func=False, nu_func=False, lambda_func=False, press_func=False, 
                 circ_vel_func=False, accel_func=False, deg_var=False, cutoff_var=False, 
                 temp_var=False, chemical_func=False, cutoff_func=False, temperature_func=False, 
                 log_dens_slope_func=False, core_func=False, plateau_func=False, maximum_r=1.0e3, 
                 relative_tolerance=5.0e-12, number_of_steps=2**10 + 1)
```

where the `parameters` variable is a numpy array object of shape (4,), whose components are (in this order): the dark matter particle mass in keV, the degeneracy parameter, the cutoff parameter, and the temperature parameter (the last three adimensional). `dens_func`, `nu_func`, `lambda_func`, `press_func`, `circ_vel_func`, `accel_func`, `deg_var`, `cutoff_var`, `temp_var`, `chemical_func`, `cutoff_func`, `temperature_func` and `core_func` are default boolean variables whose values are those indicated above. They are used as flags to compute astrophysical and statistical mechanical variables. To do so, change `False` to `True`. Regarding the other three variables, `maximum_r` is the maximum radii of integration in [$`kpc`$], `relative_tolerance` is the relative tolerance used by the integrator to solve the equations, and `number_of_steps` is the number of points used to integrate the density and pressure used to compute the right-hand side of the differential equations. They are float variables whose default values are indicated in the box above. We strongly suggest that the value of `number_of_steps` is greater than the minimum value $`2^{10} + 1`$ to ensure precision at the time of computing the solutions. Once the class `Rar` is instantiated, it automatically calls `model`, a function that solves the RAR model equations. This function is called as:

```python
model(parameters, maximum_r, relative_tolerance, number_of_steps)
```
where `parameters` is the array used in the instance of the class.

The `model` function defines several subfunctions needed to compute the right-hand side of the TOV equations. These subfunctions include a Fermi-Dirac-like distribution function, two integrands for computing the density and pressure, and two functions for computing the density and pressure themselves.

The TOV equations are solved using the `solve_ivp` function from the *scipy.integrate module*. The right-hand side of the TOV equations is computed using the function called `tov`. The solution of the TOV equations is then re-scaled to obtain the physical quantities, including the radius, enclosed mass, metric potential, and pressure of the dark matter halo. The optional attributes of the `Rar` class then compute the other physical variables using the outputs of the `model` function.

#### Rar's attributes

All the following attributes are *boolean instance attributes*.

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

In addition, there are some instance attributes representing physical quantities, which are:

- `DM_mass`: Dark matter particle mass in [$`keV`$].
- `theta_0`: Degeneracy parameter $\theta_{0}$ of the system.
- `W_0`: Cutoff parameter $W_{0}$ of the system.
- `beta_0`: Temperature parameter $\beta_{0}$ of the system.
- `r` [$`kpc`$]: Array of the radius where the solution was computed. It is a numpy ndarray of shape (n,). Available by default.
- `m` [$`M_{\odot}`$]: Array of enclosed mass at the radius given in `r`. It is a numpy ndarray of shape (n,). Available by default.
- `nu`: Array of the metric potential at the radius given in `r`. It is a numpy ndarray of shape (n,). Available by default.
- `P` [$`M_{\odot}/(kpc\ s^{2})`$]: Array of pressure at the radius given in `r`. It is a numpy ndarray of shape (n,). Available by default.
- `degeneracy_variable`: Array of values of the degeneracy variable at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `deg_var` is `True`.
- `cutoff_variable`: Array of values of the cutoff variable at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `deg_var` or `cutoff_var` is `True`.
- `temperature_variable`: Array of values of the temperature variable at the radius given in `r`. It is a numpy ndarray of shape (n,). Available by default.
- `chemical_potential` [$`keV`$]: Array of values of the chemical potential at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `chemical_func` is `True`.
- `cutoff` [$`keV`$]: Array of values of the cutoff energy function at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `cutoff_func` is `True`.
- `temperature` [$`K`$]: Array of values of the temperature function at the radius given in `r`. It is a numpy ndarray of shape (n,). Only available if `chemical_func`, `cutoff_func` or `temperature_func` is `True`.

#### Rar's methods

The only method that is computed by default is the enclosed mass of the dark matter distribution. To enable the computation of other astrophysical and statistical mechanical variables just change the boolean attributes to `True` while instantiating the object. The methods defined in the class `Rar` are:

- `mass` [$`M_{\odot}`$]: Enclosed mass of the distribution defined for all spherical radii. It is defined as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  mass(r)=
      \begin{cases}
          M(r) & \text{if } r < r_{\mathrm{max}}\\
          M(r_{\textrm{max}}) & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
$`M(r)`$ is the integrated enclosed mass. It takes a number or a numpy ndarray of shape (n,) as input and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `density` [$`M_{\odot}/kpc^{3}`$]: Mass density of the distribution defined for all spherical radii. It is computed when `dens_func=True`. It is defined as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  density(r)=
      \begin{cases}
          \frac{\mathrm{d}M(r)}{\mathrm{d}r}\frac{1}{4\pi r^{2}} & \text{if } r < r_{\mathrm{max}}\\
          0 & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `metric_potential`: Metric potential function of the fermionic distribution. It is computed when `nu_func=True`. It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `lambda_potential`: Lambda function used in the definition of the spacetime metric. It is computed when `lambda_func=True`. It is defined for all the spherical radius as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  lambda\_potential(r)=
      \begin{cases}
          -\mathrm{ln}\left[1 - \frac{2GM(r)}{c^{2}r}\right] & \text{if } r < r_{\mathrm{max}}\\
          -\mathrm{ln}\left[1 - \frac{2GM(r_{\textrm{max}})}{c^{2}r}\right] & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `pressure` [$`M_{\odot}/(kpc\ s^{2})`$]: Pressure of the fermionic distribution. It is computed when `press_func=True`. It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `circular_velocity` [$`km/s`$]: General relativistic expression of the circular velocity of the mass distribution defined for all spherical radii. It is computed when `circ_vel_func=True`. It is defined as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  circular\_velocity(r)=
      \begin{cases}
          \sqrt{\frac{1}{2}c^{2}\left[\left(\frac{8\pi G}{c^{4}}P(r)r^{2} + 1\right)\left(1 - \frac{2GM(r)}{c^{2}r}\right)^{-1} - 1\right]} & \text{if } r < r_{\mathrm{max}}\\
          \sqrt{\frac{1}{2}c^{2}\left[\left(\frac{8\pi G}{c^{4}}P(r)r^{2} + 1\right)\left(1 - \frac{2GM(r_{\textrm{max}})}{c^{2}r}\right)^{-1} - 1\right]} & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `acceleration` [$`(km/s)^{2}kpc^{-1}`$]: Newtonian gravitational field of the distribution, defined for all the spherical radius. It is computed when `accel_func=True`. It takes as input three floating numbers, the cartesian coordinates of the reference system. They have to be given as three separate variables. It is defined as a piecewise function of 2 parts, whose expression is:
```math
\begin{equation}
  acceleration(x, y, z)=
      \begin{cases}
          -\frac{GM(r)}{r^{3}}\vec{r} & \text{if } r < r_{\mathrm{max}}\\
          -\frac{GM(r_{\textrm{max}})}{r^{3}}\vec{r} & \text{if } r \geq r_{\mathrm{max}},
      \end{cases}
\end{equation}
```
where $\vec{r} = (x, y, z)$ and $r = ||\vec{r}||$. It returns a numpy ndarray of shape (3,). It is an object of the class `InterpolatedUnivariateSpline`.
- `theta`: Degeneracy variable. It is computed when `deg_var=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `W`: Cutoff variable. It is computed when `cutoff_var=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `beta`: Temperature variable. It is computed when `temp_var=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `mu` [$`keV`$]: Chemical potential. It is computed when `chemical_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `e_c` [$`keV`$]: Cutoff energy function. It is computed when `cutoff_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `T` [$`K`$]: Temperature function. It is computed when `temperature_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `logarithmic_density_slope`: Logarithmic slope of the density profile. This dimensionless quantity is computed when `log_dens_slope_func=True`. This function takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`. It is computed up to the maximum radii of integration, which can be computed by doing `halo_object.r[-1]`.
```math
\begin{equation}
  \gamma(r)= -\frac{d\mathrm{ln}\rho(r)}{d\mathrm{ln}r} = 2 - \frac{1}{4\pi r\rho(r)}\frac{d^{2}M(r)}{dr^{2}}
\end{equation}
```
- `core`: It is computed when `core_func=True`. When this function with no argument is called, it returns the *core radius* (in [$`kpc`$]) and the *mass of the core* (in [$`M_{\odot}`$]) of the distribution, as two separate floating numbers, in this order.
- `plateau`: It is computed when `plateau_func=True`. When this function with no argument is called, it returns the *plateau radius* (in [$`kpc`$]) and the *density of the plateau* (in [$`M_{\odot}kpc^{-3}`$]) of the distribution, as two separate floating numbers, in this order.

## Dependencies

The dependencies of the program are:
- Numpy
- Scipy

## Downloading the program

To download the program, just clone the repository to the directory where you would like to have it saved. For more details on how to clone a GitHub repository see [Cloning a repository](https://docs.github.com/es/repositories/creating-and-managing-repositories/cloning-a-repository).

## How to use it?

To use the program, it has to be instantiated a `Rar` object as indicated in the box above. To do so, the `Rar` class has to be imported previously as:

```python
from rar_solver import Rar
```

Take into account that the file where it will be the importing of the `Rar` class has to be in the same directory as the `rar_solver.py` file. Otherwise, it has to be changed the path to where `rar_solver.py` is saved, going from the directory where the object is going to be instantiated. This is easily done by writing the following lines in the code where it is going to be the instance:

```python
import sys

path_rar = os.path.realpath(os.path.join(os.path.dirname(__file__), 'path/', 'to/', 'rar_solver.py folder'))
sys.path.insert(0, path_rar)
from rar_solver import Rar
```

## License

This code is subject to MIT License. In addition, it is asked the following conditions in order to use the code:

- The potential works published in scientific journals that have used this code have to cite this collaboration. The official publication of the code can be found [here]().
- In case any issue or bug is found, please, report it as an issue on the GitHub page of the repository. This way, we can work to solve it as soon as possible.
