# Extended RAR model
Repository containing the programs needed to solve the RAR model's equations.

## Files description

### `rar_solver.py`:

The code defines a class whose name is `Rar`. It has to be instantiated as:

```python
halo_object = Rar(parameters, deg_var=False, cutoff_var=False, temp_var=False)
```

where the `parameters` variable is a numpy array object of shape (4,), whose components are (in this order): the dark matter particle mass in keV, the degeneracy parameter, the cutoff parameter, and the temperature parameter (the last three adimensional). `deg_var`, `cutoff_var` and `temp_var` are default boolean variables whose values are those indicated above. They are used in case the statistical mechanical variables need to be computed. To do so, just change `False` by `True`. Once the class `Rar` is instantiated, it calls `model`, a function that solves the RAR model equations. This function takes the `parameters` variable as input.

The `model` function defines several subfunctions needed to compute the right-hand side of the TOV equations. These subfunctions include a Fermi-Dirac like distribution function, two integrands for computing the density and pressure, and two functions for computing the density and pressure themselves.

The TOV equations are solved using the `solve_ivp` function from the *scipy.integrate module*. The right-hand side of the TOV equations are computed using the function called `tov`. The solution of the TOV equations is then re-scaled to obtain the physical quantities, including the mass, radius an metric potential of the dark matter halo.

#### Rar's attributes

All the following attributes are *instance attributes*.

- `parameter`: Array containing the 4 RAR's parameters, whose order is: $m$, $\theta_{0}$, $W_{0}$ and $\beta_{0}$. It is a numpy ndarray of shape (4,).
- `r` [$`kpc`$]: Array of radius where the solution was computed. It is a numpy ndarray of shape (n,).
- `m` [$`M_{\odot}`$]: Array of enclosed mass at the radius given in `r`. It is a numpy ndarray of shape (n,).
- `nu`: Array of the metric potential at the radius given in `r`. It is a numpy ndarray of shape (n,).
- `degeneracy`: Array of values of the degeneracy variable at the radius given in `r`. It is a numpy ndarray of shape (n,).
- `deg_var`: Boolean variable that enables the computation of the degeneracy variable. Default value is `False`.
- `cutoff`: Array of values of the cutoff variable at the radius given in `r`. It is a numpy ndarray of shape (n,).
- `cutoff_var`: Boolean variable that enables the computation of the cutoff variable. Default value is `False`.
- `temperature`: Array of values of the temperature variable at the radius given in `r`. It is a numpy ndarray of shape (n,).
- `temp_var`: Boolean variable that enables the computation of the temperature variable. Default value is `False`.

#### Rar's methods

- `theta`: Degeneracy variable. This functions takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `W`: Cutoff variable. This functions takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `beta`: Temperature variable. This functions takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `mass` [$`M_{\odot}`$]: Enclosed mass of the distribution defined for all spherical radii. It is defined as a piecewise function of 2 parts, whose expresion is:
```math
\begin{equation}
  mass(r)=
      \begin{cases}
          M(r) & \text{if } r < r_{\mathrm{max}}\\
          M(r_{\textrm{max}}) & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
$`M(r)`$ is the integrated enclosed mass. It takes a number or a numpy ndarray as input of shape (n,) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `density` [$`M_{\odot}/kpc^{3}`$]: Mass density of the distribution defined for all spherical radii. It is defined as a piecewise function of 2 parts, whose expresion is:
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
- `metric_potential`: Metric potential function of the fermionic distribution. It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `lambda_potential`: Lambda function used in the definition of the spacetime metric. It is defined for all the spherical radius as a piecewise function of 2 parts, whose expresion is:
```math
\begin{equation}
  \lambda(r)=
      \begin{cases}
          -\mathrm{ln}\left[1 - \frac{2GM(r)}{c^{2}r}\right] & \text{if } r < r_{\mathrm{max}}\\
          -\mathrm{ln}\left[1 - \frac{2GM(r_{\textrm{max}})}{c^{2}r}\right] & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `circular_velocity` [$`km/s`$]: General relativistic expresion of the circular velocity of the mass distribution defined for all spherical radius. It is defined as a piecewise function of 2 parts, whose expresion is:
```math
\begin{equation}
  \lambda(r)=
      \begin{cases}
           & \text{if } r < r_{\mathrm{max}}\\
           & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation} 
```
It takes a number or a numpy ndarray of shape (n,) as input (spherical radius) and returns a value or a numpy ndarray of shape (n,), respectively. It is an object of the class `InterpolatedUnivariateSpline`.
- `acceleration` [$`(km/s)^{2}kpc^{-1}`$]: Newtonian gravitational field of the distribution, defined for all the spherical radius. It takes as input three variables, the cartesian coordinates of the reference system. They have to be pased as a numpy ndarray of shape (3,). It is defined as a piecewise function of 2 parts, whose expresion is:
```math
\begin{equation}
  acceleration(\vec{r})=
      \begin{cases}
          -\frac{GM(r)}{r^{3}}\vec{r} & \text{if } r < r_{\mathrm{max}}\\
          -\frac{GM(r_{\textrm{max}})}{r^{3}}\vec{r} & \text{if } r \geq r_{\mathrm{max}}.
      \end{cases}
\end{equation}
```
It returns a numpy ndarray of shape (3,). It is an object of the class `InterpolatedUnivariateSpline`.
- `core`: This function, with no argument, when it is called returns the *core radius* (in $`kpc`$) and the *mass of the core* (in $`M_{\odot}`$). It returns this two values as float numbers.

## Dependencies

The dependencies of the programm are:
- Numpy
- Scipy

## Downloading the programm

To download the programm, just clone the repository to the directory where you would like to have it saved. For more details in how to clone a GitHub repository see [Cloning a repository](https://docs.github.com/es/repositories/creating-and-managing-repositories/cloning-a-repository).

## How to use it?

To use the programm, it has to be instantiated a `Rar` object as indicated in the box above. To do so, it has to be imported the `Rar` class before as:

```python
from rar_solver import Rar
```

Take into account that the file where it will be the importing of the `Rar` class has to be in the same directory than the `rar_solver.py` file. Otherwise, it has to be changed the path to where `rar_solver.py` is saved, going from the directory where the object is going to be instantiated. This is easyly done by writing the following lines in the code where it is going to be the instance:

```python
import sys

path_rar = os.path.realpath(os.path.join(os.path.dirname(__file__), 'path/', 'to/', 'rar_solver.py folder'))
sys.path.insert(0, path_rar)
from rar_solver import Rar
```

## License

This code is subject to the GNU General Public License v3.0. In addition, it is asked the following conditions in order to use the code:

- The potential works published in scientific journals that have used this code have to cite this collaboration. The official publication of the code can be found [here]().
- In case any issue or bug is found, please, report it as a issue on the GitHub page of the repository. This way, we can work to solve it as soon as possible.
