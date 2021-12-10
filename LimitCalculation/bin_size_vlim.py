#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
bin_size_vlim.py

Author: EB
Date: Dec 2021

Calculating the limit for varying bin sizes, using function. Usinsg asypmtotic gaussian multibin calc and singlebin.

'''

from Func_DS20K_vlim import get_limit, single_get_limit
from matplotlib import pyplot as plt
import numpy as np

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import beepy


# =============================================================================
# Calculate scaling factor for different values of firstbin and lastbin
# =============================================================================
# 15% uncertainty on backgrounds
per_err = 0.15
lastbin = 51

def get_values(threshold, lastbin,per_err):
    '''
    Returns array of scaling factor values for an increasing number of bins from the threshold

    '''
    limit_arr = []
    slimit_arr= []
    lastbin_count = threshold + 1
    while lastbin_count < lastbin + 1:
        limit = get_limit(threshold, lastbin_count, per_err)
        slimit = single_get_limit(threshold, lastbin_count, per_err)
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
fig.text(0.06, 0.5, 'c', ha='center', va='center', rotation='vertical', fontsize=26)
plt.xlabel(r'$x_{max}$',fontsize=26)
# y log
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_title('Multiple bins', fontsize=20)
ax[1].set_title('Single bin', fontsize=20)

# legend
zeroe = mpatches.Patch(color='blue', label=r'$x_{min} = 0e$')
onee = mpatches.Patch(color='red', label=r'$x_{min} = 1e$')
twoe = mpatches.Patch(color='green', label=r'$x_{min} = 2e$')
threee = mpatches.Patch(color='orange', label=r'$x_{min} = 3e$')
foure = mpatches.Patch(color='purple', label=r'$x_{min} = 4e$')

legend_elements = [zeroe, onee, twoe, threee, foure]

ax[0].legend(handles=legend_elements,fontsize=18,frameon=False,loc='upper right')

