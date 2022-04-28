#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
my_cevns.py

March 2 2022 - Copy of neutrinospec_conversion_ERtoPE - cutting out a binning method.
March 4 2022 - Editing to change the bins, bins in the response function are half electron.
April 7 2022 - Changing upper limit to use last value in the dist - removing upper threshold

author: EB

Converts neutrino spectrum data from recoil energy into Photoelectrons and then writes to file

"""

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import bisect



# =============================================================================
# Importing data
# =============================================================================

# loading in SM neutrino data from Andrew
er , spec = np.loadtxt('data/argon_spec.txt', delimiter=' ')
# convert er from GeV to KeV
er = er * 1e6


# importing energy bins - from response func
col_list = ['energy_bin_start_kev','energy_bin_end_kev','energy_kev']
df = pd.read_csv('data/remake_ds20k_ne_response_nr.csv', usecols=col_list)
ebin_start = df.energy_bin_start_kev.to_list()
ebin_end = df.energy_bin_end_kev.to_list()
ebin_centre = df.energy_kev.to_list()

# importing electron bin probabilities - from response func
df2 = pd.read_csv('data/remake_ds20k_ne_response_nr.csv', usecols=range(3,203))
s2_probs_list = df2.values.tolist()

# importing binning information - from binning info file - describes how bins relate to NPE
col_list = ['linear_center_ne','start_ne']
df = pd.read_csv('data/ds20k_s2_binning_info.csv', usecols=col_list)
s2_pe_bin_centres = df.linear_center_ne.to_list()
s2_start_pe = df.start_ne.to_list()


# =============================================================================
# Converting er to No of PE - multiplying by width
# =============================================================================

# empty list to add pe response to - list of number of PE produced in an event
pe_list = []

# define variables
# multiply to sample more - smooth things out
mult = 5000
lower_threshold = 0.1 #kev
# who is the input spec from
data_source = 'A' # andrew

# to get last probs in response map
ran = False

# getting PE count for the data given - loop though every energy
for j in range(len(er)-1): # can't calculate diff for last point
    
    # get energy and rate at that point
    energy = er[j]
    rate = spec[j]
    
    # multiply by the width of the energy bin - approx as distance between each point
    width = er[j+1]-er[j]
    rate = width * rate * mult # mult to sample more - divide by this later
    
    # if energy is less than 0.1 no PE produced - so can skip sampling - round rate to int
    if energy > lower_threshold: # not recording values under threshold
        
        # find index of energy bin - note doesn't count from zero
        index = bisect.bisect_left(ebin_start, energy) - 1
        
        # get response func probabilities for that energy bin
        # if not empty normalise so all the values add up to 1
        
        # getting last value in reponse map to use for the higher energy values
        if sum(s2_probs_list[index]) == 0 and ran == False:
            ran = True
            # get last bin probs and normalise
            norm = 1.0 / sum(s2_probs_list[index-1])
            last_bin_prob = [i * norm for i in s2_probs_list[index-1]]
            s2_bin_probs = last_bin_prob
        elif sum(s2_probs_list[index]) == 0: # if the last probs has already been found
            s2_bin_probs = last_bin_prob
        else: # if the probs are not zero
            normaliser = 1.0 / sum(s2_probs_list[index])
            s2_bin_probs = [i * normaliser for i in s2_probs_list[index]] # just multipying the whole list by norm
            
        # sample from the probs_list
        # number of times we sample = rate for that energy bin
        while rate > 1:
            
            rate -= 1
            # sampling from our distribution for that energy bin
            s2_bin_no = np.random.choice( len(s2_bin_probs), p=s2_bin_probs, size=None)
            # get number of PE and scale within bin width
            pe_no = s2_pe_bin_centres[s2_bin_no] * energy / ebin_centre[index]
            # add to list
            pe_list.append(pe_no)

            
# =============================================================================
# Binning the pe_list data - and dividing by mult
# =============================================================================
# getting histogram
pe_no, bin_edges = np.histogram(pe_list, s2_start_pe)

# dividing by mult
pe_no = pe_no / mult



# rebinning into 1 PE wide bins

my_cevns = []
for i in range(int(len(pe_no)/2)) :
    my_cevns.append((pe_no[i*2]+pe_no[i*2+1])) # joint adjacent bins as they are 0.5 e wide ( to make them 1e- wide)
my_cevns = np.array(my_cevns)

# getting new bin centres
bins = bin_edges[1::2]

# =============================================================================
# Writing data to file
# =============================================================================
filename = 'output_my_cevns/PE_argon_SM_' + data_source + '.txt'
file = open(filename, 'w')

for n in range(len(bins)-1):
    file.write(str(bins[n]) + ' ' + str(my_cevns[n]) + '\n')
    
file.close()



#%%

# =============================================================================
# Plotting
# =============================================================================


firstbin = 0
lastbin = 60

cevns = np.loadtxt('output_my_cevns/PE_argon_SM_A.txt',delimiter=' ')[firstbin:lastbin,1]
bins = np.loadtxt('output_my_cevns/PE_argon_SM_A.txt',delimiter=' ')[firstbin:lastbin,0]


# loading data - for comparison
ds20k_cevns = np.loadtxt('data/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]
ds20k_bins = np.loadtxt('data/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,0]


# convert old cevns into events/tyr - dividing by ds20k exposure in tyr
ds20k_cevns = ds20k_cevns / 100


# plotting
f=plt.figure(figsize=(10,8))

plt.plot(bins, cevns, '-+',markersize=15, label='Using my_cevns', color='firebrick')
plt.plot(ds20k_bins,ds20k_cevns, '-+',markersize=15, label='RH spec in PE', color='royalblue')

plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Events per tyr',fontsize=26) 
plt.yscale('log')
plt.legend(fontsize=18,frameon=False,loc='upper right')


