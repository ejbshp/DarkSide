#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
18th May 2022

EB

Looping through all of the grid values to get upper limits for the spec

"""



#TODO copy func_DS20k_vlim and make a callable function - however the cevns signal will need to be passed to each func
# do not declare it as a global variable
# use the grid_values.txt to loop through

import numpy as np
from func_limit_calc import get_limit

firstbin = 4
lastbin = 50

bkgrd_err = 0.15

ds20k_exposure_tonneyear = 100

# standard model for testing
cennssig_sm = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[firstbin:lastbin,1]*ds20k_exposure_tonneyear
rel_cenns_err_sm = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[firstbin:lastbin,2]

# unpacking pairs of values I want to use for the grid
g_values = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/grid_values.txt',delimiter=' ')[:,0]
m_values = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/grid_values.txt',delimiter=' ')[:,1]

#%%


# empty list for upper limits
limit_list = []
# for sanity
gs = []
ms = []

for i in range(len(g_values)):
    
    g = g_values[i]
    m = m_values[i]
    
    spec_file = '/Users/user/DarkSide/Neutrino spectrum/output_grid/mutau_pe_spec_' + str(g) + '_' + str(m) + '.txt'
    
    #loading data
    cennssig = np.loadtxt(spec_file ,delimiter=' ')[firstbin:lastbin,1]*ds20k_exposure_tonneyear
    rel_cenns_err = np.loadtxt(spec_file ,delimiter=' ')[firstbin:lastbin,2]
    # removing infinite or nans from errirs
    rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value
    
    # get limit
    upper_limit = get_limit(cennssig, rel_cenns_err, firstbin, lastbin, bkgrd_err)
    
    limit_list.append(upper_limit)
    gs.append(g)
    ms.append(m)
    
    

#%% Save to file

filename = 'grid_limits_' + str(firstbin) + '_' + str(lastbin) + '_' + str(bkgrd_err) + '_' + str(ds20k_exposure_tonneyear)  + '.txt'
file = open(filename, 'w')

for n in range(len(g_values)):
    file.write(str(g_values[n]) + ' ' + str(m_values[n]) + ' ' + str(limit_list[n]) + '\n')
    
file.close()
    
