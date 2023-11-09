# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:32:25 2023

@author: Santiago Collazo
"""

import matplotlib.pyplot as plt
import numpy as np
from rar_solver import Rar

# ============================= Plot features ================================= #
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
# ============================================================================= #

# ========================= Solving the RAR model ============================= #
beta_0 = 1.113903337971913738e-05
theta_0 = 3.780867927802387385e+01
W_0 = 6.644915273597560201e+01
m_DM = 5.480880070125579806e+01   # keV
G_u = 4.3009e-6          # (km/s)^2*kpc/M_sun
c = 2.99792458e+5        # Light speed - km/s

h = Rar(np.array([m_DM,
                  theta_0,
                  W_0,
                  beta_0]), nu_func=True, temp_var=True)
r = np.logspace(np.log10(h.r[0]), np.log10(h.r[-1]), 10**6, endpoint=False)

nu_0 = 2.0*np.log(np.sqrt(1.0 - 2.0*G_u*h.m[-1]/(c**2*h.r[-1]))*1.0/(1.0 + beta_0*W_0))
# ============================================================================= #

# ================================== Plot ===================================== #
fig, ax = plt.subplots(1, 1, figsize=(6,6), dpi=380)
plt.xscale('log')
ax.plot(r, np.exp((h.metric_potential(r) - nu_0)/2)*h.beta(r)/beta_0, lw=3, ls=':', color='#91430e', label=r'$\frac{\beta(r)}{\beta_{0}}\mathrm{e}^{(\nu(r) - \nu_{0})/2}$')
ax.axhline(1, lw=1, color='black')
ax.axvline(h.r[-1], lw=1, color='darkblue', ls='-.', label=r'$r_{\mathrm{max}}$')
ax.set_xlim(1.0e-10, 80)
ax.set_ylim(0.99975, 1.00150)
ax.set_xlabel("r [kpc]")
ax.legend()
plt.show()
# ============================================================================= #