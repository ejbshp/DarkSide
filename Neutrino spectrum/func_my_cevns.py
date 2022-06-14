#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
func_my_cevns.py

author: EB

Function that onverts neutrino spectrum data from recoil energy into Photoelectron

May 12 2022: Function now returns stat error of the cevns sig

"""

import numpy as np
import pandas as pd
import bisect

# =============================================================================
# Importing data
# =============================================================================

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



def spectope(er, spec):
    '''
    Converts neutrino spectrum to be in terms of photoelectrons
    using the detector response
    

    Parameters
    ----------
    er : Array of recoil energy in Gev
    spec : Rate

    Returns
    -------
    my_cevns: The events per tyr
    bins: the binning information (PE)
    my_cevns_err: relative stat error on my_cevns - when plotting multiply by the spec

    '''
    # empty list to add pe response to - list of number of PE produced in an event
    pe_list = []
    
    # define variables
    # multiply to sample more - smooth things out
    mult = 1000
    lower_threshold = 0.1 #kev
    # to get last probs in response map
    ran = False
    
    # convert er from GeV to KeV
    er = er * 1e6
    
    # getting PE count for the data given - loop though every energy
    for j in range(len(er)-1): # can't calculate diff for last point
        print(j)
        # get energy and rate at that point
        energy = er[j]
        rate = spec[j]
        
        # multiply by the width of the energy bin - approx as distance between each point
        width = er[j+1]-er[j]
        rate = width * rate * mult # mult to sample more - divide by this later
        
        # if energy is less than 0.1 no PE produced - so can skip sampling - round rate to int
        if energy > lower_threshold:
            
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
    
    # making sure that bins and my_cevns are the same length
    bins = bins[0:len(my_cevns)]
        
    # =============================================================================
    # Stat errors - 1/sqrt(N)
    # =============================================================================
    my_cevns_err = 1 / np.sqrt(my_cevns*mult)
    
    return my_cevns, bins, my_cevns_err


