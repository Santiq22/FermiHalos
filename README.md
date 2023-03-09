# Extended RAR model
Repository containing the programs needed to solve the RAR model's equations and their applications.

## Files description

### `rar_solver.py':

The code defines a function called `model' that solves the RAR model equations. The model function takes four parameters as inputs (in this order): the particle mass of the dark matter particle in keV, the degeneracy parameter, the cutoff parameter, and the temperature parameter (the last three adimensional). These parameter has to be passed as numpy array objects. It also takes four optional boolean variables (`deg_var', `cutoff_var', `temp_var', and `den_var') indicating whether the statistical mechanical variables and the density of the distribution need to be computed.

The code defines several subfunctions needed to compute the right-hand side of the TOV equations. These subfunctions include a Fermi-Dirac distribution function, two integrands for computing the density and pressure, and two functions for computing the density and pressure themselves.

The TOV equations are solved using the solve_ivp function from the scipy.integrate module. The right-hand side of the TOV equations are computed using the function called `tov'. The solution of the TOV equations is then re-scaled to obtain the physical quantities, including the mass, radius an metric potential of the dark matter halo. Finally, the model function returns a list containing these physical quantities, as well as the statistical mechanical variables and the density of the distribution if requested.

### `dwarfs.py'

rar-model is free software, distributed under the GNU General Public License v3.0. This implies that you may freely distribute and copy the software. You may also modify it as you wish, and distribute these modified versions as long as you indicate prominently any changes you made in the original code, and as long as you leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details. Note that the authors retain their copyright on the code.
