#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
May 2022

EB

Plotting grid of g, m and limits
"""

import numpy as np
from matplotlib import pyplot as plt


# unpacking values - lower e thresh _ upper e thresh _ bkgrd err _ exposure tyr
g_values = np.loadtxt('grid_limits_4_50_0.15_100.txt',delimiter=' ')[:,0]
m_values = np.loadtxt('grid_limits_4_50_0.15_100.txt',delimiter=' ')[:,1]
limits = np.loadtxt('grid_limits_4_50_0.15_100.txt',delimiter=' ')[:,2]


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

for i, txt in enumerate(limits):
    if limits[i] > 0.01: ax.annotate('{:.2f}'.format(txt), (m_values[i], g_values[i]))
    elif limits[i] > 0.001: ax.annotate('{:.3f}'.format(txt), (m_values[i], g_values[i]))
    elif limits[i] > 0.0001: ax.annotate('{:.4f}'.format(txt), (m_values[i], g_values[i]))

# plotting the diff with gx and M
plt.scatter(m_values, g_values, c=limits, cmap = plt.get_cmap('plasma'))
plt.ylabel('$g_x$ [GeV]', size=16)
plt.xlabel('$M_A$ [GeV]', size=16)
plt.xscale('log')
plt.yscale('log')


cbar = plt.colorbar()
cbar.set_label(r'$c$', size=20)
plt.show()
