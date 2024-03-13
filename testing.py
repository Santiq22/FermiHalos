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

param_pwr = np.array([20.0,
                      3.441729863769217701e+01,
                      2.0*3.441729863769217701e+01,
                      1.274212705635805845e-07])

param_bec = np.array([56.0,
                      37.765595,
                      66.34067,
                      1.1977342e-05])

h_pwr = Rar(param_pwr, dens_func=True)
h_bec = Rar(param_bec, dens_func=True)

r = np.logspace(-10, np.log10(h_pwr.r[-1]), 10**5, endpoint=False)

fig, ax = plt.subplots(1, 1, figsize=(6,6), dpi=380)

ax.plot(r, r*r*h_pwr.density(r), lw=2, color='#91430e', label='RAR 2')
ax.plot(r, r*r*h_bec.density(r), lw=2, color='black', label='Becerra et al.')
ax.set_xlim(5.0e-3, 60)
plt.xscale('log')
plt.yscale('log')
ax.set_xlabel("r [kpc]")
ax.set_ylabel(r'$r^{2}\rho(r)\ [M_{\odot}/kpc]$')
#ax.set_ylabel(r'$T(r)$ [K]')
#ax.set_ylabel(r'$\beta(r)$')
ax.legend()
plt.show()