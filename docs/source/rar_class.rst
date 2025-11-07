Rar class
=========

This class instantiates an object of the ``Rar`` class, which represents a fermionic dark matter halo following the extended RAR model. The constructor of this class will integrate the differential equations defining the halo, based on the input free parameters. An object of this class will have different attributes and methods, as stated below.

Rar's attributes
---------------

All the following attributes are *boolean instance attributes* and *set up instance attributes*.

- ``dens_func``: Boolean variable that enables the computation of the density profile of the distribution. The default value is ``False``.
- ``nu_func``: Boolean variable that enables the computation of the metric potential. The default value is ``False``.
- ``particles_func``: Boolean variable that enables the computation of the enclosed particle number. The default value is ``False``.
- ``lambda_func``: Boolean variable that enables the computation of the lambda potential. The default value is ``False``.
- ``press_func``: Boolean variable that enables the computation of the pressure profile. The default value is ``False``.
- ``n_func``: Boolean variable that enables the computation of the particle number density. The default value is ``False``.
- ``circ_vel_func``: Boolean variable that enables the computation of the circular velocity profile. The default value is ``False``.
- ``accel_func``: Boolean variable that enables the computation of the Newtonian gravitational field exerted by the dark matter halo. The default value is ``False``.
- ``deg_var``: Boolean variable that enables the computation of the degeneracy variable. The default value is ``False``.
- ``cutoff_var``: Boolean variable that enables the computation of the cutoff variable. The default value is ``False``.
- ``temp_var``: Boolean variable that enables the computation of the temperature variable. The default value is ``False``.
- ``chemical_func``: Boolean variable that enables the computation of the chemical potential. The default value is ``False``.
- ``cutoff_func``: Boolean variable that enables the computation of the cutoff energy function. The default value is ``False``.
- ``temperature_func``: Boolean variable that enables the computation of the temperature function. The default value is ``False``.
- ``log_dens_slope_func``: Boolean variable that enables the computation of the logarithmic density slope function. The default value is ``False``.
- ``core_func``: Boolean variable that enables the computation of the radii of the dark matter core and its mass. The default value is ``False``.
- ``plateau_func``: Boolean variable that enables the computation of the radii of the dark matter plateau and its density. The default value is ``False``.
- ``maximum_r``: (*float*) Maximum radius of integration in :math:`kpc`.
- ``relative_tolerance``: (*float*) Relative tolerance used by the integrator to solve the equations.
- ``number_of_steps``: *(int)* Number of steps used to integrate the density and pressure used to compute the right-hand side of the differential equations. We strongly suggest that the value of ``number_of_steps`` is greater than the minimum value :math:`2^{10} + 1` to ensure precision at the time of computing the solutions.

In addition, there are some instance attributes representing physical quantities, which are:

- ``DM_mass`` :math:`keV/c^{2}`: Dark matter particle mass.
- ``theta_0``: Degeneracy parameter :math:`\theta_{0}` of the system.
- ``W_0``: Cutoff parameter :math:`W_{0}` of the system.
- ``beta_0``: Temperature parameter :math:`\beta_{0}` of the system.
- ``r`` :math:`kpc`: Array of the radius where the solution was computed. It is a numpy ndarray of shape (n,). Available by default.
- ``m`` :math:`M_{\odot}`: Array of enclosed masses at the radius given in ``r``. It is a numpy ndarray of shape (n,). Available by default.
- ``nu``: Array of metric potentials (dimensionless) at the radius given in ``r``. It is a numpy ndarray of shape (n,). Available by default.
- ``N``: Array of enclosed particle numbers at the radius given in ``r``. It is a numpy ndarray of shape (n,). Available by default. 
- ``nu_0``: Value of the metric potential at the center of the distribution, :math:`\nu_{0}`. Available by default.
- ``P`` :math:`M_{\odot}/(kpc\ s^{2})`: Array of pressures at the radius given in ``r``. It is a numpy ndarray of shape (n,). Only available if ``press_func`` is ``True``.
- ``n`` :math:`kpc^{-3}`: Array of particle number densities at the radius given in ``r``. It is a numpy ndarray of shape (n,). Only available if ``n_func`` is ``True``.
- ``degeneracy_variable``: Array of values of the degeneracy variable at the radius given in ``r``. It is a numpy ndarray of shape (n,). Only available if ``deg_var`` or ``chemical_func`` is ``True``.
- ``cutoff_variable``: Array of values of the cutoff variable at the radius given in ``r``. It is a numpy ndarray of shape (n,). Only available if ``deg_var``, ``cutoff_var``, ``chemical_func`` or ``cutoff_func`` is ``True``.
- ``temperature_variable``: Array of values of the temperature variable at the radius given in ``r``. It is a numpy ndarray of shape (n,). Available by default.
- ``chemical_potential`` :math:`keV`: Array of values of the chemical potential at the radius given in ``r``. It is a numpy ndarray of shape (n,). Only available if ``chemical_func`` is ``True``.
- ``cutoff`` :math:`keV`: Array of values of the cutoff energy function at the radius given in ``r``. It is a numpy ndarray of shape (n,). Only available if ``cutoff_func`` is ``True``.
- ``temperature`` :math:`K`: Array of values of the temperature function at the radius given in ``r``. It is a numpy ndarray of shape (n,). Only available if ``chemical_func``, ``cutoff_func`` or ``temperature_func`` is ``True``.