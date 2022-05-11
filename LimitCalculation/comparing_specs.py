#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparing specs

"""

import numpy as np
from matplotlib import pyplot as plt
from Func_DS20K_vlim import get_limit, single_get_limit

# =============================================================================
# Loading in data
# =============================================================================

firstbin = 0
lastbin = 60


oldcenns = np.loadtxt('Data_files/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]/100 # ds20k exposure
bins = np.loadtxt('Data_files/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,0]
mycenns = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A.txt',delimiter=' ')[firstbin:lastbin,1]
mybins = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A.txt',delimiter=' ')[firstbin:lastbin,0]

# comparing total number of events
print(sum(oldcenns))
print(sum(mycenns))



# plotting
f=plt.figure(figsize=(10,8))

plt.plot(mybins, mycenns, '-+',markersize=15, label='A my_cevns (new spec)', color='firebrick')
plt.plot(bins, oldcenns, '-+',markersize=15, label='RH (old spec)', color='royalblue')

plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Events per tyr',fontsize=26) 
plt.yscale('log')
plt.legend(fontsize=18,frameon=False,loc='upper right')


# =============================================================================
# Limit and single limit - uses old spec
# =============================================================================

print('Limit ', get_limit(4, 50, 0.15))

print('Single limit ', single_get_limit(0, 10, 0.15))