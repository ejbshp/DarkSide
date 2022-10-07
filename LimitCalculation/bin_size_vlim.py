#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
bin_size_vlim.py

Author: EB
Date: Dec 2021

Calculating the limit for varying bin sizes, using function. Usinsg asypmtotic gaussian multibin calc and singlebin.

Edit July 2022 - new SM signal
'''

from func_limit_calc import get_limit_var, single_get_limit_var
from matplotlib import pyplot as plt
import numpy as np

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import beepy


# =============================================================================
# Calculate scaling factor for different values of firstbin and lastbin
# =============================================================================

per_err = 0.15
lastbin=50

# SM CEvNS - loading data
sm_file = '/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt'
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[0:50,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[0:50,2]
# removing infinite or nans from errors
sm_rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value


def get_values(threshold, lastbin,per_err):
    '''
    Returns array of scaling factor values for an increasing number of bins from the threshold

    '''
    limit_arr = []
    slimit_arr= []
    lastbin_count = threshold + 1
    while lastbin_count < lastbin + 1:
        
        limit = get_limit_var(sm_cenns[threshold:lastbin_count], sm_rel_cenns_err[threshold:lastbin_count], threshold, lastbin_count, bkgrd_err=per_err)
        slimit = single_get_limit_var(sm_cenns[threshold:lastbin_count], sm_rel_cenns_err[threshold:lastbin_count], threshold, lastbin_count, bkgrd_err=per_err)

        limit_arr.append(limit)
        slimit_arr.append(slimit)
        lastbin_count += 1
    # bins
    bins = np.arange(threshold,lastbin)
    return limit_arr, slimit_arr, bins


#%% 0e threshold

limits0, slimits0, bins0 = get_values(0,lastbin,per_err) # multi


beepy.beep(sound=5)

#%% 1e threshold

limits1, slimits1, bins1 = get_values(1,lastbin,per_err)

beepy.beep(sound=5)

#%% 2e threshold

limits2, slimits2, bins2 = get_values(2,lastbin,per_err)

beepy.beep(sound=5)
#%% 3e threshold

limits3, slimits3, bins3 = get_values(3,lastbin,per_err)

beepy.beep(sound=5)
#%% 4e threshold

limits4, slimits4, bins4 = get_values(4,lastbin,per_err)

beepy.beep(sound=5)

#%% Plotting limit for different values of xmax

fig,ax=plt.subplots(2,figsize=(16,16))



# multibin
ax[0].plot(bins0, limits0, '-o', linewidth=2, color='blue')
ax[0].plot(bins1, limits1, '-o', linewidth=2, color='red')
ax[0].plot(bins2, limits2, '-o', linewidth=2, color='green')
ax[0].plot(bins3, limits3, '-o', linewidth=2, color='orange')
ax[0].plot(bins4, limits4, '-o', linewidth=2, color = 'purple')


# singlebin
ax[1].plot(bins0, slimits0, '-D', linewidth=2, color='blue')
ax[1].plot(bins1, slimits1, '-D', linewidth=2, color='red')
ax[1].plot(bins2, slimits2, '-D', linewidth=2, color='green')
ax[1].plot(bins3, slimits3, '-D', linewidth=2, color='orange')
ax[1].plot(bins4, slimits4, '-D', linewidth=2, color = 'purple')

# y label
fig.text(0.06, 0.5, r'Scaling factor on SM CE$\nu$NS signal', ha='center', va='center', rotation='vertical', fontsize=26)
plt.xlabel(r'Upper electron threshold',fontsize=26)

# y log
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_title('Multiple bins', fontsize=20)
ax[1].set_title('Single bin', fontsize=20)

# legend
zeroe = mpatches.Patch(color='blue', label=r'$0e$')
onee = mpatches.Patch(color='red', label=r'$1e$')
twoe = mpatches.Patch(color='green', label=r'$2e$')
threee = mpatches.Patch(color='orange', label=r'$3e$')
foure = mpatches.Patch(color='purple', label=r'$4e$')

legend_elements = [zeroe, onee, twoe, threee, foure]


legend = ax[0].legend(handles=legend_elements, title='Lower electron threshold',fontsize=18,frameon=False,loc='upper right')
legend.get_title().set_fontsize('20') 

plt.savefig('/Users/user/Documents/Plots/updated_sm_studies/bin_size_new_sm.eps', format='eps')

plt.savefig('/Users/user/Documents/Plots/updated_sm_studies/bin_size_new_sm.jpg', format='jpg')
            


