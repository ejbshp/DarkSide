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
ds20k_cevns = np.loadtxt('data/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]
bins = np.loadtxt('data/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,0]
new_cevns = np.loadtxt('data/argon_spec_PE_multw.txt',delimiter=' ')[firstbin:lastbin,1]
new_bins = np.loadtxt('data/argon_spec_PE_multw.txt',delimiter=' ')[firstbin:lastbin,0]
old_cevns = np.loadtxt('data/old_argon_spec_PE_multw.txt',delimiter=' ')[firstbin:lastbin,1]
old_bins = np.loadtxt('data/old_argon_spec_PE_multw.txt',delimiter=' ')[firstbin:lastbin,0]


# convert old cevns into events/tyr - dividing by ds20k exposure in tyr
ds20k_cevns = ds20k_cevns / 100

# plotting
f=plt.figure(figsize=(10,8))
plt.plot(new_bins,new_cevns, label='New CEvNS', color='firebrick')
plt.plot(bins,ds20k_cevns, label='Old CEvNS already in PE')
plt.plot(old_bins,old_cevns, label='Old CEvNS using my conversion')
plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Events per tyr',fontsize=26) 
plt.yscale('log')
plt.legend(fontsize=18,frameon=False,loc=9)
