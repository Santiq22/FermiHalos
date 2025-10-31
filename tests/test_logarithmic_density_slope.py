"""
File: test_logarithmic_density_slope.py
Created on 2025-10-30 11:23:17
Author: Santiago Collazo
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
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

# ======================================== Burkert model ======================================== #
def rho_b(r: np.ndarray | float, 
          rho_0: float, 
          r_0: float) -> np.ndarray | float:
    """ Burkert density profile """
    
    x = r/r_0
    return rho_0/(x + 1.0)/(x**2 + 1.0)

class Burkert:
    def __init__(self, 
                 rho_0: float, 
                 r_0: float):
        self.rho_0 = rho_0
        self.r_0 = r_0
        
    def density(self, 
                r: np.ndarray | float) -> np.ndarray | float:
        """ Density profile """
        
        return rho_b(r, self.rho_0, self.r_0)
    
    def logarithmic_density_slope(self, 
                                  r: np.ndarray) -> np.ndarray:
        """ Logarithmic density slope """
        
        x = r/self.r_0
        prefactor = 8.0*np.pi*r*self.rho_0
        simple = (x + 1.0)**(-1)
        double = (x**2 + 1.0)**(-1)
        d2M_dr2 = prefactor*simple*double*(-r/2.0/self.r_0*simple \
            - x**2*double + 1)
        
        return 2.0 - 1.0/(4.0*np.pi*r*self.density(r))*d2M_dr2
# =============================================================================================== #

# ================================== Solving the Burkert model ================================== #
rho_0 = 2.197416380230734870e+07       # M_sun/kpc^3
r_0 = 1.008698047532039155e+01         # kpc
halo_b = Burkert(rho_0=rho_0, r_0=r_0)
# =============================================================================================== #

# ==================================== Solving the RAR model ==================================== #
beta_0 = 1.208047426367269421e-05
theta_0 = 3.741551266351441996e+01
W_0 = 6.574489281155653941e+01
m_DM = 56.0                            # keV

halo = Rar(np.array([m_DM, theta_0, W_0, beta_0]),
           log_dens_slope_func = True, 
           core_func = True,
           plateau_func = True)

r = np.logspace(np.log10(halo.r[0]), np.log10(halo.r[-1]), 10**6, endpoint = False)
# =============================================================================================== #

# ============================================ Plot ============================================= #
fig, ax = plt.subplots(1, 1, figsize = (6, 6), dpi = 380)
ax.plot(r, halo.logarithmic_density_slope(r), lw = 2.0, 
        color = '#91430e', label = 'RAR 2')
ax.plot(r, halo_b.logarithmic_density_slope(r), lw = 2.0, ls = '-.', 
        color = "#3d700f", label = 'Burkert')
ax.axvline(halo.core()[0], lw = 2, color = 'black', label = r'$r_{\mathrm{core}}$')
ax.axvline(halo.plateau()[0], lw = 2, color = 'violet', ls = '--', label = r'$r_{\mathrm{plateau}}$')
plt.xscale('log')
ax.set_xlim(1.0e-10, halo.r[-1])
ax.set_xlabel('r [kpc]')
plt.yscale('log')
#ax.set_ylim(-15.0, 100000.0)
ax.set_ylabel('Logarithmic density slope')
ax.legend()
#fig.savefig('../figures/logarithmic_density_slope.png', bbox_inches = 'tight')
# =============================================================================================== #