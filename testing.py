# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:32:25 2023

@author: Santiago Collazo
"""

import matplotlib.pyplot as plt
import numpy as np
from rar_solver import Rar

# ----------- Plot features -----------
# Properties to decorate the plots.
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'serif'   
plt.rcParams['font.sans-serif'] = 'New Century Schoolbook' # 'Times', 'Liberation Serif', 'Times New Roman'
#plt.rcParams['font.serif'] = ['Helvetica']
plt.rcParams['font.size'] = 13
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.edgecolor'] = 'k'
plt.rcParams['legend.markerscale'] = 7
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.top'] = False
plt.rcParams['ytick.right'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width']= 0.5
plt.rcParams['xtick.major.size']= 5.0
plt.rcParams['xtick.minor.width']= 0.5
plt.rcParams['xtick.minor.size']= 3.0
plt.rcParams['ytick.major.width']= 0.5
plt.rcParams['ytick.major.size']= 5.0
plt.rcParams['ytick.minor.width']= 0.5
plt.rcParams['ytick.minor.size']= 3.0

h = Rar(np.array([5.480880070125579806e+01,
                  3.780867927802387385e+01,
                  6.644915273597560201e+01,
                  1.113903337971913738e-05]), dens_var=True)
r = np.logspace(-10, 2, 10**5)


fig, ax = plt.subplots(1, 1, figsize=(6,6), dpi=380)
plt.xscale('log')
plt.yscale('log')
ax.plot(r, h.density(r), lw=1, color='#91430e', label='Best fit model')
ax.set_xlim(1e-10, 1e2)
ax.set_ylim(1e2, 1e28)
ax.set_xlabel("r [kpc]")
ax.set_ylabel(r'$\rho(r)\ [M_{\odot}/kpc^{3}]$')
ax.legend()
plt.show()
