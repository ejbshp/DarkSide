#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
bin_size_vlim.py

Author: EB
Date: Dec 2021

Calculating the limit for different background uncertainty

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
lastbin = 51
err_arr = np.geomspace(0.01,0.5,50)

def get_values(threshold, lastbin,per_err_arr):
    '''
    Returns array of scaling factor values for an array of different background errors.

    '''
    limit_arr = []
    slimit_arr= []
    for err in per_err_arr:
        print(err)
        limit = get_limit(threshold, lastbin, err)
        slimit = single_get_limit(threshold, lastbin, err)
        limit_arr.append(limit)
        slimit_arr.append(slimit)
    # bins
    bins = np.arange(threshold,lastbin)
    return limit_arr, slimit_arr, bins

   
#%% 0e threshold

limits0, slimits0, bins0 = get_values(0,lastbin,err_arr)


beepy.beep(sound=5)

#%% 1e threshold

limits1, slimits1, bins1 = get_values(1,lastbin,err_arr)

beepy.beep(sound=5)

#%% 2e threshold

limits2, slimits2, bins2 = get_values(2,lastbin,err_arr)

beepy.beep(sound=5)
#%% 3e threshold

limits3, slimits3, bins3 = get_values(3,lastbin,err_arr)

beepy.beep(sound=5)
#%% 4e threshold

limits4, slimits4, bins4 = get_values(4,lastbin,err_arr)

beepy.beep(sound=5)

#%% Plotting limit for different values of xmax

fig,ax=plt.subplots(1,figsize=(20,14))

# multibin
ax.plot(err_arr, limits0, '-o', linewidth=2, color='blue')
ax.plot(err_arr, limits1, '-o', linewidth=2, color='red')
ax.plot(err_arr, limits2, '-o', linewidth=2, color='green')
ax.plot(err_arr, limits3, '-o', linewidth=2, color='orange')
ax.plot(err_arr, limits4, '-o', linewidth=2, color = 'purple')


plt.ylabel('c', fontsize=26)
plt.xlabel(r'Background uncertainty',fontsize=20)
# y log
ax.set_xscale('log')


# legend
zeroe = mpatches.Patch(color='blue', label=r'$x_{min} = 0e$')
onee = mpatches.Patch(color='red', label=r'$x_{min} = 1e$')
twoe = mpatches.Patch(color='green', label=r'$x_{min} = 2e$')
threee = mpatches.Patch(color='orange', label=r'$x_{min} = 3e$')
foure = mpatches.Patch(color='purple', label=r'$x_{min} = 4e$')

legend_elements = [zeroe, onee, twoe, threee, foure]

ax.legend(handles=legend_elements,fontsize=18,frameon=False,loc='upper left')


