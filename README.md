# Extended RAR model
Repository containing the programs needed to solve the RAR model's equations and their applications.

## Files description

### `rar_solver.py`:

The code defines a class whose name is `Rar`. It has to be instantiated as:

```python
halo_object = Rar(parameters, deg_var=False, cutoff_var=False, temp_var=False)
```

where the `parameters` variable is a numpy array object of shape (4,), whose components are (in this order): the dark matter particle mass in keV, the degeneracy parameter, the cutoff parameter, and the temperature parameter (the last three adimensional). `deg_var`, `cutoff_var` and `temp_var` are default boolean variables whose values are those indicated above. They are used in case the statistical mechanical variables and the density of the distribution need to be computed. To do so, just change `False` by `True`. The class `Rar` calls `model`, a function that solves the RAR model equations. This function takes the `parameters` variable as input.

The `model` function defines several subfunctions needed to compute the right-hand side of the TOV equations. These subfunctions include a Fermi-Dirac like distribution function, two integrands for computing the density and pressure, and two functions for computing the density and pressure themselves.

The TOV equations are solved using the `solve_ivp` function from the scipy.integrate module. The right-hand side of the TOV equations are computed using the function called `tov`. The solution of the TOV equations is then re-scaled to obtain the physical quantities, including the mass, radius an metric potential of the dark matter halo.

## Dependencies

The dependencies of the programm are:
- Numpy
- Scipy

## Downloading the programm

To download the programm, just clone the repository to the directory where you would like to have it saved.

## How to use it?

To use the programm, it has to be instantiated a `Rar` object as indicated in the box above. To do so, it has to be imported before as:

```python
from rar_solver import Rar
```

Take into account that the file where it will be the importing of the `Rar` class has to be in the same directory than the `rar_solver.py` file. Otherwise, it has to be changed the path to where `rar_solver.py` is saved, going from the directory where the object has to be imported. This is easyly done as:

```python
import sys

path_rar = os.path.realpath(os.path.join(os.path.dirname(__file__), 'path/', 'to/', 'rar_solver.py folder/'))
sys.path.insert(0, path_rar)
from rar_solver import Rar
```
