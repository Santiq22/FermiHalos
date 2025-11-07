Quick start
===========

To use the program, it has to be instantiated a ``Rar`` object as indicated at the beggining of the documentation. To do so, the ``Rar`` class has to be imported previously as:

.. code-block:: python

    >>> from fermihalos import Rar

Then, by setting up the 4 free parameters of the model and the boolean flags to allow for the computation of different astrophysical quantities, the differential equations are integrated during the instance of the ``Rar`` object and the methods are ready to use.

.. code-block:: python
    
    >>> import numpy as np
    >>> p = np.array([150.0, 39.0, 70.0, 1.0e-5])   # [m, theta_0, W_0, beta_0]
    >>> halo = Rar(p, dens_func=True)
    >>> halo.density(0.15)                          # [r] = kpc
    >>> array(3.10583169e+11)                       # [density] = M_sun/kpc^3