#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 2021

@author: EB

Plotting new vs old neutrino spectrum - both PE
"""

from matplotlib import pyplot as plt
import numpy as np

firstbin = 0
lastbin = 50

# loading data
ds20k_cevns = np.loadtxt('ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]
bins = np.loadtxt('ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,0]
new_cevns = np.loadtxt('argon_spec_PE.txt',delimiter=' ')[firstbin:lastbin,1]
new_bins = np.loadtxt('argon_spec_PE.txt',delimiter=' ')[firstbin:lastbin,0]



# lew SM data energy vs rate
er , spec = np.loadtxt('argon_spec.txt', delimiter=' ')
# convert er from GeV to KeV
er = er * 1e6

# convert old cevns into events/tyr - dividing by ds20k exposure in tyr
old_cevns = ds20k_cevns / 100

# plotting
f=plt.figure(figsize=(10,8))
plt.plot(new_bins,new_cevns, label='New CEvNS', color='firebrick')
plt.plot(bins,old_cevns, label='Old CEvNS')
plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Events per tyr',fontsize=26) 
plt.yscale('log')
plt.legend(fontsize=18,frameon=False,loc=9)
