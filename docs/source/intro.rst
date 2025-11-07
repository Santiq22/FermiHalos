Introduction
============

The code defines a class whose name is ``Rar``. It has to be instantiated as:

.. code-block:: python

    >>> halo_object = Rar(parameters, dens_func=False, nu_func=False, particles_func=False, lambda_func=False, press_func=False, n_func=False, circ_vel_func=False, accel_func=False, deg_var=False, cutoff_var=False, temp_var=False, chemical_func=False, cutoff_func=False, temperature_func=False, log_dens_slope_func=False, core_func=False, plateau_func=False, maximum_r=1.0e3, relative_tolerance=5.0e-12, number_of_steps=2**10 + 1)

where the ``parameters`` variable is a numpy array object of shape (4,), whose components are (in this order): the dark matter particle mass in :math:`keV/c^{2}`, the degeneracy parameter, the cutoff parameter, and the temperature parameter (the last three are dimensionless). The boolean variables are used as flags to compute astrophysical and statistical mechanical variables. To do so, change ``False`` to ``True``. See Rar's attributes section to further details on the variables passed to the Rar class.

Inside ``Rar`` class
--------------------

Once the class ``Rar`` is instantiated, it automatically calls ``model``, a function that solves the RAR model equations. This function is called as:

.. code-block:: python

    >>> model(parameters, maximum_r, relative_tolerance, number_of_steps, press_func, n_func)

where ``parameters`` is the array used in the instance of the class.

The ``model`` function defines several subfunctions needed to compute the right-hand side of the Tolman-Oppenheimer-Volkoff (TOV) equations. These subfunctions include a Fermi-Dirac-like distribution function, three integrands for computing the density, pressure, and particle number density, and three functions for computing the density, pressure, and particle number density themselves.

The TOV equations are solved using the ``solve_ivp`` function from the ``scipy.integrate`` module. The right-hand side of the TOV equations is computed using the function called ``tov``. The solution of the TOV equations is then re-scaled to obtain the physical quantities, including the radius, enclosed mass, metric potential, pressure, particle number density, the metric potential at the origin of the distribution, the temperature variable, and the enclosed particle number of the dark matter halo. The optional attributes of the ``Rar`` class then enable the computation of the other physical variables using the outputs of the ``model`` function.