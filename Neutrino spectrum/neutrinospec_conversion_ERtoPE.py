#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 2021

author: EB

Converts neutrino spectrum data from recoil energy into Photoelectrons and then writes to file

"""
import numpy as np
import pandas as pd
import random
import bisect

# loading in SM neutrino data from theorists - New SM
#er , spec = np.loadtxt('argon_spec.txt', delimiter=' ')
old_er , old_spec = np.loadtxt('rateTOTAr_old_spec_for_comparison.txt', delimiter='\t', unpack=True)

# convert er from GeV to KeV
er = er * 1e6


# importing energy bins
col_list = ['energy_bin_start_kev','energy_bin_end_kev']
df = pd.read_csv('remake_ds20k_ne_response_nr.csv', usecols=col_list)
ebin_start = df.energy_bin_start_kev.to_list()
ebin_end = df.energy_bin_end_kev.to_list()

# importing electron bin probabilities
df2 = pd.read_csv('remake_ds20k_ne_response_nr.csv', usecols=range(3,203))
ne_probs = df2.values.tolist()

'''
# =============================================================================
# Converting er to No of PE
# =============================================================================

# empty array for np count
pe_count = np.zeros(len(ne_probs[0]))


# getting PE count for the data given
for i in range(len(ebins)):
    # get rate at that energy by interpolating data
    n = np.interp(ebins[i],er,spec)
    # multiplying by bin width
    index = next(x for x, val in enumerate(er) if val > ebins[i])
    width = er[index] - er[index-1]
    n = n * width
    # sampling n times
    while n > 0:
        # sampling from our distribution for that energy bin
        ne = random.choices(list(range(0,200)), ne_probs[i])
        # adding to count
        pe_count[ne.pop()] += 1
        n -= 1
   '''     
   
#%%
# =============================================================================
# Converting er to No of PE
# =============================================================================

# empty array for np count
pe_count = np.zeros(len(ne_probs[0]))

# getting PE count for the data given
for j in range(len(er)-1): # can't calculate diff for last point
    energy = er[j]
    rate = spec[j]
    # multiply by the width of the energy bin - approx as distance between each point
    width = er[j+1]-er[j]
    rate = width * rate
    # if energy is less than 0.1 no PE produced - so can skip sampling
    if energy < 0.1:
        pe_count[0] += rate
    else:
        # loop though to find energy bin
        for i in range(len(ebin_end)):
            emin = ebin_start[i]
            emax = ebin_end[i]
            if energy >= emin and energy < emax:
                # number of times we sample = rate for that energy bin
                while rate > 0:  
                    # sampling from our distribution for that energy bin
                    ne = random.choices(list(range(0,200)), ne_probs[i])
                    # adding to count
                    pe_count[ne.pop()] += 1
                    rate -= 1
                break
            




#%%

col_list = ['energy_bin_end_kev']
df = pd.read_csv('remake_ds20k_ne_response_nr.csv', usecols=col_list)
max_e = df.energy_bin_end_kev.to_list()


binned_rate = []
max_e_check = []
n_check = 0
last_ind = 0

pe_count = np.zeros(len(ne_probs[0]))

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
    
for i in range(len(binned_rate)):
    rate = binned_rate[i]
    while rate > 0:  
        # sampling from our distribution for that energy bin
        ne = random.choices(list(range(0,200)), ne_probs[i])
        # adding to count
        pe_count[ne.pop()] += 1
        rate -= 1
    





#%%







# Ne bins
bins = np.arange(len(pe_count))

# =============================================================================
# Writing data to file
# =============================================================================

file = open('argon_spec_PE.txt', 'w')
for n in bins:
    file.write(str(n) + ' ' + str(pe_count[n]) + '\n')
file.close()

diff = []
for i in range(len(er)-1):
    diff.append(er[i+1]-er[i])
    
    
    

    

