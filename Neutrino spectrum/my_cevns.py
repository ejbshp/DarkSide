#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
my_cevns.py

March 2 2022 - Copy of neutrinospec_conversion_ERtoPE - cutting out a binning method.
March 4 2022 - Editing to change the bins, bins in the response function are half electron.

author: EB

Converts neutrino spectrum data from recoil energy into Photoelectrons and then writes to file

"""

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import bisect
import uproot


# =============================================================================
# Importing data
# =============================================================================

# loading in SM neutrino data from Andrew
er , spec = np.loadtxt('data/argon_spec.txt', delimiter=' ')
# convert er from GeV to KeV
er = er * 1e6

####### remove this - using paolo's spec overwriting andrew's spec
# loading in old SM - already in 
er , spec = np.loadtxt('data/rateTOTAr_old_spec_for_comparison.txt', delimiter='\t', unpack=True)


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

# multiply to sample more - smooth things out
mult = 1000

# getting PE count for the data given - loop though every energy
for j in range(len(er)-1): # can't calculate diff for last point
    
    # get energy and rate at that point
    energy = er[j]
    rate = spec[j]
    
    # multiply by the width of the energy bin - approx as distance between each point
    width = er[j+1]-er[j]
    rate = width * rate * mult # mult to sample more - divide by this later
    
    # if energy is less than 0.1 no PE produced - so can skip sampling - round rate to int
    if energy < 0.1 or energy >5.0: # adding in greater than 5 as in cenns.py
    
        # add 0 PE to the pe_list as many times as you would sample
        l = [0.0] * int(rate)
        pe_list.extend(l)
        
    else:
        
        # find index of energy bin - note doesn't count from zero
        index = bisect.bisect_left(ebin_start, energy) - 1
        
        # get response func probabilities for that energy bin
        # if not empty normalise so all the values add up to 1
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


#%%

# rebinning into 1 PE wide bins

cenns1 = []
for i in range(int(len(pe_no)/2)) :
    cenns1.append((pe_no[i*2]+pe_no[i*2+1])) # joint adjacent bins as they are 0.5 e wide ( to make them 1e- wide)
cenns1 = np.array(cenns1)[0:50]




# =============================================================================
# Plotting
# =============================================================================


firstbin = 0
lastbin = 50

# loading data - for comparison
ds20k_cevns = np.loadtxt('data/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]
bins = np.loadtxt('data/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,0]


# convert old cevns into events/tyr - dividing by ds20k exposure in tyr
ds20k_cevns = ds20k_cevns / 100


# get data from cenns.py root files

# not in darkside 50 - from RH a simulation they did
cenns = 'data/spectra_cenns.root'
cenns_data = uproot.open(cenns)
hist = cenns_data['hNR_per_tonne_yr']
cenns_events_pertonneyr = np.array(hist.values(),dtype='float64')
cennsbkgrd = []
for i in range(int(len(cenns_events_pertonneyr)/2)) :
    cennsbkgrd.append((cenns_events_pertonneyr[i*2]+cenns_events_pertonneyr[i*2+1])) # joint adjacent bins as they are 0.5 e wide ( to make them 1e- wide)
cennsbkgrd = np.array(cennsbkgrd)[0:50]





# plotting
f=plt.figure(figsize=(10,8))
plt.plot(bins,ds20k_cevns, label='Old CEvNS already in PE')
plt.plot(bins, cenns1, label='Using s2 binning info + 1PE wide bins')
plt.plot(bins, cennsbkgrd, label='from cenns.py')

#plt.grid()
#plt.xlim(0,50)
plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Events per tyr',fontsize=26) 
plt.yscale('log')
plt.legend(fontsize=18,frameon=False,loc=9)

