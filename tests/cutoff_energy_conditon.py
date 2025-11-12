# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:32:25 2023
@author: Santiago Collazo
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Real path from the root directory to where rar_class.py is
#path_rar = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'fermihalos'))

# Insert the path location of rar_class.py into the list of folders giving by sys.path
#sys.path.insert(0, path_rar)
#from rar_class import Rar
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

# ====================================== Extra parameters ======================================= #
machine_eps = np.finfo(float).eps                      # Machine epsilon
# =============================================================================================== #

# ==================================== Solving the RAR model ==================================== #
beta_0 = 1.113903337971913738e-05
theta_0 = 3.780867927802387385e+01
W_0 = 6.644915273597560201e+01
m_DM = 5.480880070125579806e+01                        # keV

halo = Rar(np.array([m_DM, theta_0, W_0, beta_0]), 
           nu_func = True, 
           cutoff_func = True,
           core_func = True,
           plateau_func = True)

r = np.logspace(np.log10(halo.r[0]), np.log10(halo.r[-1]), 10**6, endpoint = False)

nu_0 = halo.nu_0
# =============================================================================================== #

# ============================================ Plot ============================================= #
fig, ax = plt.subplots(1, 1, figsize = (6, 6), dpi = 380)
ax.plot(r, np.exp((halo.metric_potential(r) - nu_0)/2.0)*(halo.e_c(r) + m_DM)/(m_DM*(1.0 + beta_0*W_0)) - 1.0, 
        lw = 0.5, ls = ':', color = '#91430e') 
        #label = r'$\frac{\mathrm{e}^{(\nu(r) - \nu_{0})/2}(\epsilon_{c}(r) + mc^{2})}{mc^{2}(1 + \beta_{0}W_{0})} - 1$')
ax.axvline(halo.core()[0], lw = 3, color = 'black', label = r'$r_{\mathrm{core}}$')
ax.axvline(halo.plateau()[0], lw = 3, color = 'violet', ls = '--', label = r'$r_{\mathrm{plateau}}$')
ax.axvline(halo.r[-1], lw = 3, color = 'darkblue', ls = '-.', label = r'$r_{\mathrm{max}}$')
plt.xscale('log')
ax.set_xlim(1.0e-10, 80)
ax.set_ylim(-7.0*machine_eps, 7.0*machine_eps)
ax.set_yticks([-6.0*machine_eps, -4.0*machine_eps, -2.0*machine_eps, 0.0, 
               2.0*machine_eps, 4.0*machine_eps, 6.0*machine_eps])
ax.set_yticklabels([r'-6$\epsilon$', r'-4$\epsilon$', r'-2$\epsilon$', '0', 
                    r'2$\epsilon$', r'4$\epsilon$', r'6$\epsilon$'])
ax.set_xlabel("r [kpc]")
ax.legend(loc = 'upper left')
#fig.savefig('../figures/cutoff_energy_condition.png', bbox_inches = 'tight')
# =============================================================================================== #