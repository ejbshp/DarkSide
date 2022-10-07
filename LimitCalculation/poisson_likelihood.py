#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
DarkSide20k_neutrino_limit.py

Author: EB
Date: Nov 2021

Calculating the likelihood for a single bin

Edit June 22: Using new sm and scaling background

'''


### imports ###
import numpy as np
from scipy import optimize
from scipy.stats import poisson


firstbin = 4 #!!! usually 4
lastbin =  15 #!!! usually 50

# =============================================================================
# Poisson likelihood calculation - works only for one bin
# =============================================================================

def poisson_likelihood_limit (counts,bkgrd,observed) :
   '''
    Gives the 90% CL for a single bin using a poisson likelihood method
    bkgrd: number of background events
    observed: number of events observed (bkgrd + sig)
    counts: number of signal counts expected for SM (c=1)
   '''
   if counts==0:
       return np.inf
   else:
        # Sum m=0 to m=Nobvs P(m|Nb,theta) - 0.1 = 0 find root to get Nobvs. Nobvs - Nb = Nsig
        # mu is the shape factor, second parameter for poisson.cdf is our inital guess for the root
        # CDF of X (random variable) evaluated at x is the P(X<=x) so finding roots of CDF(Nobs) - 0.1 = 0 gives our 90% CL
       Nobs = float(optimize.root(lambda mu: poisson.cdf(observed, mu) - 0.1,observed).x)
       N_90 = Nobs - bkgrd
       return N_90/counts

scale = 2.5

ar38bkgrd = np.sum(np.loadtxt('Data_files/ds20k-39Ar_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1])*scale

gammabkgrd = np.sum(np.loadtxt('Data_files/ds20k-gamma_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1])*scale
# neutrino signal
sm_file = '/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt'
cennssig = np.sum(np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,1]*100) # multiply by exposure

prediction = ar38bkgrd + gammabkgrd
measured_values = prediction

# we are predicting the background
print("Poisson likelihood limit:")
print(poisson_likelihood_limit(cennssig, prediction, measured_values))


