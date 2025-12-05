"""
File: density_plots.py
Created on 2025-12-05 12:22:15
Author: Santiago Collazo
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Real path from the root directory to where rar_class.py is
#path_rar = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'fermihalos'))

# Insert the path location of rar_class.py into the list of folders giving by sys.path
#sys.path.insert(0, path_rar)
# from rar_class import Rar
from fermihalos import Rar

# ======================================== Plot features ======================================== #
# Properties to decorate the plots.
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'serif'   
plt.rcParams['font.sans-serif'] = 'New Century Schoolbook' # 'Times', 'Liberation Serif', 'Times New Roman'
#plt.rcParams['font.serif'] = ['Helvetica']
plt.rcParams['font.size'] = 17
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
# =============================================================================================== #

# ========================================= Parameters ========================================== #
theta_0 = 37.765595
W_0 = 66.34067
beta_0 = 1.1977342e-05
m_DM = 56.0
c = 9.7156140203724e-12                                  # kpc/s
machine_eps = np.finfo(float).eps                        # Machine epsilon
integration_methods = ['LSODA', 'BDF', 'DOP853', 'Radau']
# =============================================================================================== #

# Creating figure and axis
fig, ax = plt.subplots(1, 1, figsize = (6, 6), dpi = 380)

for method in integration_methods:
    # ==================================== Solving the RAR model ==================================== #
    halo = Rar(np.array([m_DM, theta_0, W_0, beta_0]), 
            dens_func = True,
            int_method=method)

    r = np.logspace(np.log10(halo.r[0]), np.log10(halo.r[-1]), 10**6, endpoint = False)

    # Density and pressure factor
    rho = halo.density(r)
    # =============================================================================================== #
    
    # ============================================ Plot ============================================= #
    ax.plot(r, rho, lw = 2, label=method)
    # =============================================================================================== #
    
plt.xscale('log')
plt.yscale('log')
#ax.set_xlim(1.0e-10, 80.0)
ax.set_xlim(1.0e-6, 1.0e-5)
ax.set_ylim(1.0e12, 1.0e18)
ax.set_xlabel("r [kpc]")
ax.set_ylabel(r"$\rho(r)$ $[M_{\odot}/kpc^{3}]$")
ax.legend(loc = 'upper right')
fig.savefig('../figures/density_profile_diffetent_methods_zoom_in.png', bbox_inches = 'tight')
# =============================================================================================== #