"""
File: plots.py
Created on 2025-09-12 11:38:02
Author: Santiago Collazo
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Real path from the root directory to where rar_class.py is
path_rar = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'fermihalos'))

# Insert the path location of rar_class.py into the list of folders giving by sys.path
sys.path.insert(0, path_rar)
from rar_class import Rar

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

# ====================================== Extra parameters ======================================= #
machine_eps = np.finfo(float).eps                      # Machine epsilon
# =============================================================================================== #

# ==================================== Solving the RAR model ==================================== #
"""beta_0 = 1.113903337971913738e-05
theta_0 = 3.780867927802387385e+01
W_0 = 6.644915273597560201e+01
m_DM = 5.480880070125579806e+01                          # keV"""
theta_0 = 37.765595
W_0 = 66.34067
beta_0 = 1.1977342e-05
m_DM = 56.0
c = 9.7156140203724e-12                                  # kpc/s

halo = Rar(np.array([m_DM, theta_0, W_0, beta_0]), 
           nu_func = True, 
           dens_func = True,
           press_func = True,
           core_func = True,
           plateau_func = True,
           number_of_steps = 20000)

r = np.logspace(np.log10(halo.r[0]), np.log10(halo.r[-1]), 10**7, endpoint = False)

# Gradient of metric potential - Accessing using the mangled name
dnu_dr = halo._Rar__dnu_dr(r)

# Density and pressure factor
rho_P = halo.density(r)*c*c + halo.pressure(r)

# Gradient of pressure - Accessing using the mangled name
dP_dr = (halo._Rar__P_spline).derivative(1)
dP_dr = dP_dr(r)
# =============================================================================================== #

# ============================================ Plot ============================================= #
fig, ax = plt.subplots(1, 1, figsize = (6, 6), dpi = 380)
ax.plot(r, dP_dr/rho_P*halo.core()[0], lw = 2.0, label = 'P', color = "#064980")
ax.plot(r, 0.5*dnu_dr*halo.core()[0], lw = 2.0, label = r'$\nu$', color = "#367A12")
ax.plot(r, 0.5*dnu_dr*halo.core()[0] + dP_dr/rho_P*halo.core()[0], 
        lw = 1.0, color = "#3D3F3C")
ax.hlines(y = min(dP_dr/rho_P*halo.core()[0]), xmin = 1.0e-10, xmax = 80.0, lw = 1.0, color = 'grey')
ax.hlines(y = -min(dP_dr/rho_P*halo.core()[0]), xmin = 1.0e-10, xmax = 80.0, lw = 1.0, color = 'grey')
ax.axvline(halo.core()[0], lw = 3, color = 'black', label = r'$r_{\mathrm{core}}$')
ax.axvline(halo.plateau()[0], lw = 3, color = 'violet', ls = '--', label = r'$r_{\mathrm{plateau}}$')
ax.axvline(halo.r[-1], lw = 3, color = 'darkblue', ls = '-.', label = r'$r_{\mathrm{max}}$')
plt.xscale('log')
ax.set_xlim(1.0e-10, 80.0)
#ax.set_xlim(halo.core()[0], 80.0)
ax.set_ylim(-1.0e-6, 1.0e-6)
#ax.set_ylim(-1 - 9.0*machine_eps, -1 + 9.0*machine_eps)
#ax.set_yticks([-8.0*machine_eps, -6.0*machine_eps, -4.0*machine_eps, -2.0*machine_eps, 0.0, 
#               2.0*machine_eps, 4.0*machine_eps, 6.0*machine_eps, 8.0*machine_eps])
#ax.set_yticklabels([r'-8$\epsilon$', r'-6$\epsilon$', r'-4$\epsilon$', r'-2$\epsilon$', '0', 
#                    r'2$\epsilon$', r'4$\epsilon$', r'6$\epsilon$', r'8$\epsilon$'])
ax.set_xlabel("r [kpc]")
ax.legend(loc = 'upper left')
fig.savefig('../figures/grad_P_vs_grad_nu.png', bbox_inches = 'tight')
# =============================================================================================== #