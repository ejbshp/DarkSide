#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
bin_size_vlim.py

Author: EB
Date: Dec 2021

Calculating the limit for different background uncertainty

Edit: June 2022 - Using updated SM Spec
                - Adding in scaling background size

'''

from func_limit_calc import get_limit_var
from Func_DS20K_vlim import get_limit, single_get_limit
from matplotlib import pyplot as plt
import numpy as np

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker

import beepy

#%% Import data and def func

lastbin =50

err_arr = np.geomspace(0.01,0.5,50)

# SM CEvNS
sm_file = '/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt'



#%% 0e threshold

firstbin = 0
#loading data
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,2]
# removing infinite or nans from errirs
sm_rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value



limits0 = []
for err in err_arr:
    print(err)
    limit = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err)
    print(limit)
    limits0.append(limit)
   
#%% 1e threshold

firstbin = 1
#loading data
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,2]
# removing infinite or nans from errirs
sm_rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value


limits1 = []
for err in err_arr:
    print(err)
    limit = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err)
    print(limit)
    limits1.append(limit)
    


#%% 2e threshold

firstbin = 2
#loading data
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,2]
# removing infinite or nans from errirs
sm_rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value


limits2 = []
for err in err_arr:
    print(err)
    limit = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err)
    print(limit)
    limits2.append(limit)
    
#%%
    
# #!!! quick fix
# lim = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=0.20)
# print(lim)
# limits2[15] = lim

#%% 3e threshold

firstbin = 3
#loading data
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,2]
# removing infinite or nans from errirs
sm_rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value


limits3 = []
for err in err_arr:
    print(err)
    limit = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err)
    print(limit)
    limits3.append(limit)
#%% 4e threshold

firstbin = 4
#loading data
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,2]
# removing infinite or nans from errirs
sm_rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value


limits4 = []
for err in err_arr:
    print(err)
    limit = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err)
    print(limit)
    limits4.append(limit)
#%% Plotting limit for different values of xmax

fig,ax=plt.subplots(1,figsize=(20,14))

# multibin
ax.plot(err_arr, limits0, '-o', linewidth=2, color='blue')
ax.plot(err_arr, limits1, '-o', linewidth=2, color='red')
ax.plot(err_arr, limits2, '-o', linewidth=2, color='green')
ax.plot(err_arr, limits3, '-o', linewidth=2, color='orange')
ax.plot(err_arr, limits4, '-o', linewidth=2, color = 'purple')


plt.ylabel(r'Scaling factor on SM CE$\nu$NS signal', fontsize=26)
plt.xlabel(r'Fractional uncertainty on the radioactive backgrounds',fontsize=26)
# y log
ax.set_xscale('log')
ax.set_ylim(0.0,0.3)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks( [0.01,0.02,0.05,0.1,0.25,0.5])
#ax.set_yticklabels([0.0,0.05,0.1,0.15,0.2,0.25,0.3],fontsize=20)
# legend
zeroe = mpatches.Patch(color='blue', label=r'$0$')
onee = mpatches.Patch(color='red', label=r'$1$')
twoe = mpatches.Patch(color='green', label=r'$2$')
threee = mpatches.Patch(color='orange', label=r'$3$')
foure = mpatches.Patch(color='purple', label=r'$4$')

legend_elements = [zeroe, onee, twoe, threee, foure]

legend = ax.legend(handles=legend_elements,fontsize=18,frameon=False,loc='upper left', title='Electron Threshold')
legend.get_title().set_fontsize('22') 

plt.savefig('/Users/user/Documents/Plots/updated_sm_studies/bkgrderr_scaling.eps', format='eps')

