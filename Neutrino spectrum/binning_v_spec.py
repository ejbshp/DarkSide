#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 2021

author: EB

Plotting the methods of binning the new neutrino spectra.

"""
import numpy as np
import pandas as pd
import bisect
from matplotlib import pyplot as plt

# loading in SM neutrino data from theorists - New SM
er , spec = np.loadtxt('argon_spec.txt', delimiter=' ')
# convert er from GeV to KeV
er = er * 1e6

# max energy value is 100
# want to bin the data in 0.1keV wide bins 
# integrate the values eg bettween 0 and 0.1keV
# 5000 data points so am not interpolating - as x values not linear there could
# be some loss of shape in the higher keVs but it should be ok
# response function is binned in 0.1keV wide bins until 19.9keV where the bins increase to 0.5keV

# =============================================================================
# Binning the data - using trapz - binning the data into the same size bins as the response function
# =============================================================================
# loadinng in bin ends 
col_list = ['energy_bin_end_kev']
df = pd.read_csv('remake_ds20k_ne_response_nr.csv', usecols=col_list)
max_e = df.energy_bin_end_kev.to_list()


binned_rate = []
max_e_check = []
n_check = 0
last_ind = 0

for e in max_e:
    # index of last value below 0.1
    ind = bisect.bisect_left(er, e) -1
    x = er[last_ind:ind]
    y = spec[last_ind:ind]
    rate = np.trapz(y,x)
    # add to list
    max_e_check.append(x)
    binned_rate.append(rate)
    n_check += len(x)
    last_ind = ind


# =============================================================================
# Binning the data - using trapz - using each x data point as a bin
# =============================================================================
# loadinng in bin ends 
col_list = ['energy_bin_end_kev']
df = pd.read_csv('remake_ds20k_ne_response_nr.csv', usecols=col_list)
max_e = df.energy_bin_end_kev.to_list()


binned_rate2 = []


for j in range(len(er)-1): 
    rate = np.trapz([spec[j],spec[j+1]],[er[j],er[j+1]])
    # add to list
    binned_rate2.append(rate)

# =============================================================================
# Binning the data - width
# =============================================================================
# approximate the bin width as the distance between bins

rate_arr = []
width_arr = []
# getting PE count for the data given
for j in range(len(er)-1): # can't calculate diff for last point
    rate = spec[j]
    # multiply by the width of the energy bin - approx as distance between each point
    width = er[j+1]-er[j]
    rate = width * rate
    rate_arr.append(rate)
    width_arr.append(width)
    
    
print(bisect.bisect_left(er, 25e-3))
print(np.sum(rate_arr[498:]))
    
# =============================================================================
# Plotting
# =============================================================================
f=plt.figure(figsize=(12,8))
plt.yscale('log')
plt.xscale('log')
plt.plot(er[0:4999],rate_arr ,'o' , color='firebrick')
plt.plot(er[0:4999],binned_rate2)

plt.ylabel(r'Rate', size=26)
plt.xlabel(r'$E_R$ $[keV]$', size=26)
plt.legend(fontsize=18,frameon=False,loc='upper right')


# same bins as response function
f2=plt.figure(figsize=(12,8))
plt.bar(max_e, binned_rate, width=0.1,alpha=0.7, log=True)
plt.plot(max_e,binned_rate,'o')
plt.yscale('log')
plt.xscale('log')



