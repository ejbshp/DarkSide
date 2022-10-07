#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
bkgrd_scaling_vlim.py

Author: EB
Date: June 2022

Scaling the total radioactive background

Edit: July 2022 Scaling backgrounds individually

'''

from func_limit_calc import get_limit_var
from Func_DS20K_vlim import get_limit, single_get_limit
from matplotlib import pyplot as plt
import numpy as np

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import beepy

#%% Import data and def func

lastbin = 50
err = 0.15
# scale_list = np.geomspace(0.85,1.75,10)
scale_list = np.arange(0.5,2.75,0.1)
# SM CEvNS
sm_file = '/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt'



#%% 4e threshold

firstbin = 4
#loading data
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[firstbin:lastbin,2]
# removing infinite or nans from errirs
sm_rel_cenns_err = np.nan_to_num(rel_cenns_err, nan=1.0, posinf=1.0, neginf=1.0) # replae with 1 so errors will be as large as the value


total_limits = []
ar_limits = []
gam_limits = []

for scale in scale_list:
    print('Scale ',scale)
    
    # total background
    limit = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err, ar_scale= scale, gam_scale= scale)
    total_limits.append(limit)
    
    # ar background
    limit2 = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err, ar_scale= scale, gam_scale=1)
    ar_limits.append(limit2)
    
    # gam background
    limit3 = get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err, ar_scale= 1, gam_scale=scale)
    gam_limits.append(limit3)

#%% test

# print(get_limit_var(sm_cenns, sm_rel_cenns_err, firstbin, lastbin, bkgrd_err=err, ar_scale= 1, gam_scale=2.4))
    

#%% Plotting limit for different values of xmax

fig,ax=plt.subplots(1,figsize=(10,8))

ax.plot(scale_list, total_limits, '-o', linewidth=2, color='purple')



plt.axvline(x=1,linestyle='dashed',color='k',linewidth=2,label='Expected background')
plt.ylabel(r'Scaling factor on SM CE$\nu$NS signal', fontsize=20)
plt.xlabel(r'Multiplicative factor on the total radioactive background',fontsize=20)
plt.legend(fontsize=12,frameon=False,loc='upper center')

plt.savefig('/Users/user/Documents/Plots/updated_sm_studies/bkgrd_scaling_4e_0.15.eps', format='eps')

