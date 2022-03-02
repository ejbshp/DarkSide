#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 2021

@author: EB

Plotting neutrino spec and the ratio of the old spec to the new spec - All in Recoil energy
"""


import numpy as np
from matplotlib import pyplot as plt

# loading in SM neutrino data from theorists - New SM
er , spec = np.loadtxt('argon_spec.txt', delimiter=' ')
# convert er from GeV to KeV
er = er * 1e6

# Data from RH before conversion to PE
old_er , old_spec = np.loadtxt('rateTOTAr_old_spec_for_comparison.txt', delimiter='\t', unpack=True)

ex_er , ex_spec = np.loadtxt('approx_of_recoil_E_from_NRfromCEVNSdoc.csv', delimiter=',', unpack=True)

# =============================================================================
# plotting spectra together
# =============================================================================


f=plt.figure(figsize=(12,8))
plt.loglog(er,spec, label='New')
plt.loglog(old_er, old_spec, label='Old')

plt.loglog(ex_er, ex_spec, label='Approx of expected spec')
plt.ylabel(r'$\frac{dR}{dE_R}$ $[(keV\ ton\ yr)^{-1}]$', size=26)
plt.xlabel(r'$E_R$ $[keV]$', size=26)
plt.legend(fontsize=18,frameon=False,loc='upper right')


# =============================================================================
# plotting ratio
# =============================================================================

xvals = np.geomspace(0.01,100,5000)

# interpolate new specs so we can ratio
nyvals = np.interp(xvals,er,spec)
oyvals = np.interp(xvals,old_er,old_spec)

ratio = nyvals/oyvals

f2=plt.figure(figsize=(12,8))
plt.plot(xvals,ratio,'-o')
plt.ylabel(r'Ratio of rate ', size=26)
plt.xlabel(r'$E_R$ $[keV]$', size=26)

