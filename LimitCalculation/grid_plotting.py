#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
May 2022

EB

Plotting grid of g, m and limits
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import glob
from scipy import interpolate


# unpacking values - lower e thresh _ upper e thresh _ bkgrd err _ exposure tyr
g_values = np.loadtxt('more_grid_limits_4_50_0.15_100.txt',delimiter=' ')[:,0]
m_values = np.loadtxt('more_grid_limits_4_50_0.15_100.txt',delimiter=' ')[:,1]
limits = np.loadtxt('more_grid_limits_4_50_0.15_100.txt',delimiter=' ')[:,2]


#%% trying to get into 3d matrix

# remove repeats from m and g lists
y_val = list(dict.fromkeys(g_values))
x_val = list(dict.fromkeys(m_values))
y_val.sort()
x_val.sort()

# make empty matrix for values 
limit_matrix = np.zeros((len(y_val),len(x_val)))

# make the values nan - as we haven't calculated every value in this matrix
limit_matrix[:] = np.nan

# filling matrix - values I don't have will be left as 0
for i in range(len(g_values)):
    
    g = g_values[i]
    m = m_values[i]
    c = limits[i]
    
    # find indicies for matrix values
    y_ind = y_val.index(g)
    x_ind = x_val.index(m)
    
    # add to matrix
    limit_matrix[y_ind, x_ind] = c
    

#%% Plotting matrix

f, ax = plt.subplots(figsize=(20,20))

for i, txt in enumerate(limits):
    if limits[i] >= 0.01: ax.annotate('{:.2f}'.format(txt), (m_values[i], g_values[i]))
    elif limits[i] >= 0.0025: ax.annotate('{:.3f}'.format(txt), (m_values[i], g_values[i]))
    elif limits[i] >= 0.0001: ax.annotate('{:.4f}'.format(txt), (m_values[i], g_values[i]), va = 'bottom', ha = 'right')
# add empty values to end of lists as pcolor drops the last row
x_val.append(1)
y_val.append(1)
X,Y=np.meshgrid(x_val,y_val) 

cmap = plt.cm.magma.copy()

cmap.set_bad(color='lightgrey')
im = plt.pcolormesh(X,Y,limit_matrix, cmap=cmap)

plt.xscale('log')
plt.yscale('log')

plt.ylabel('$g_x$ [GeV]', size=16)
plt.xlabel('$M_A$ [GeV]', size=16)

plt.xlim(1.5e-3,1)
plt.ylim(1.25e-4,1e-2)

cbar = plt.colorbar()

cbar.set_label(r'$c$', size=20)
plt.show()


#%% Plotting as scatter plot


fig, ax = plt.subplots()

# for i, txt in enumerate(limits):
#     if limits[i] > 0.01: ax.annotate('{:.2f}'.format(txt), (m_values[i], g_values[i]))
#     elif limits[i] > 0.001: ax.annotate('{:.3f}'.format(txt), (m_values[i], g_values[i]))
#     elif limits[i] > 0.0001: ax.annotate('{:.4f}'.format(txt), (m_values[i], g_values[i]))

# plotting the diff with gx and M
plt.scatter(m_values, g_values, c=limits, vmin = 0.1, vmax = 1, cmap = plt.get_cmap('plasma'))
plt.ylabel('$g_x$ [GeV]', size=16)
plt.xlabel('$M_A$ [GeV]', size=16)
plt.xscale('log')
plt.yscale('log')


cbar = plt.colorbar()
cbar.set_label(r'$c$', size=20)
plt.show()


#%% Interpolating data for all of our datapoints

# Loading all the data into a pandas dataframe

LMuTaupaths = glob.glob('/Users/user/DarkSide/Neutrino spectrum/input/*')

df_MuTau = pd.DataFrame()

for p in LMuTaupaths:
    filepath = p + '/spectra_Ar.npy'
    ERs, g_x, m_A, spec=np.load(filepath, allow_pickle=True)
    for i in range(len(m_A)):
        temp_entry = {'ERs': ERs, 'g_x':g_x, 'm_A':m_A[i], 'spec':spec[i]}
        #print(temp_entry)
        df_MuTau = df_MuTau.append(temp_entry, ignore_index=True)
        
# add column with the interpolated values   
# interpolating function
f = interpolate.interp2d(m_values, g_values, limits)
df_MuTau['interp_lim'] = df_MuTau.apply(lambda x: f(x.m_A,g_x), axis=1)


#%% Plotting scatter plot of interpolated values

fig, ax = plt.subplots()

for i, txt in enumerate(limits):
    if limits[i] > 0.01: ax.annotate('{:.2f}'.format(txt), (m_values[i], g_values[i]))
    elif limits[i] > 0.001: ax.annotate('{:.3f}'.format(txt), (m_values[i], g_values[i]))
    elif limits[i] > 0.0001: ax.annotate('{:.4f}'.format(txt), (m_values[i], g_values[i]))

# plotting interpolated values for full colour
plt.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['interp_lim'].values, vmin=0, cmap = plt.get_cmap('plasma'))


# plotting the calculated values as points
plt.plot(m_values, g_values, 'o')
plt.ylabel('$g_x$ [GeV]', size=16)
plt.xlabel('$M_A$ [GeV]', size=16)
plt.xscale('log')
plt.yscale('log')


cbar = plt.colorbar()
cbar.set_label(r'$c$', size=20)
plt.show()