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
plt.rcParams['font.size'] = 10
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

h = Rar(np.array([5.000000000000000000e+01,
                  3.770201401012075593e+01,
                  6.628801944438806970e+01,
                  8.712527791070693930e-06]), circ_vel_var=True)
r = np.logspace(-7, 2, 10**5)

fig, ax = plt.subplots(1, 1, figsize=(6,6), dpi=320)
plt.xscale('log')
plt.yscale('log')
ax.plot(r, h.circular_velocity_RG(r), lw=1, color='black', label='GR')
ax.plot(r, h.circular_velocity(r), lw=1, color='darkblue', label='Newtonian')
ax.legend()
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(6,6), dpi=320)
plt.xscale('log')
plt.yscale('log')
ax.plot(r, h.pressure(r), lw=1, color='black')
ax.legend()
plt.show()
