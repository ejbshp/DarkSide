#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aug 2022

Importing files from grid jobs and plotting specs


"""

import numpy as np
import pandas as pd
import glob
from matplotlib import pyplot as plt


# import sm
sm_cevns = np.loadtxt('output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[0:50,1]
sm_bins = np.loadtxt('output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[0:50,0]
sm_err = np.loadtxt('output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[0:50,2]


#%% Loop through all files in grid_jobs_output to load all the grid proccessed
# pe specs into a dataframe

# Loading all the data into a pandas dataframe
df_pe = pd.DataFrame()



g_paths = glob.glob('grid_jobs_output/*')

for g_file in g_paths:
    # loop through all directories in grid_jobs_output
   
    #g_file = 'output_g_0.0001'
    
    # get value of g
    indx = g_file.rfind('_g_')
    g_val = g_file[indx+3:]
    g_x = float(g_val)
    
    pepaths = glob.glob(g_file + '/*')
    
    for p in pepaths:
        # loop through contents of g_file
        
        # get value of g
        indx = p.rfind('_m_')
        m_val = p[indx+3:-4]
        m_a = float(m_val)
        
        # load data
        bins = np.loadtxt(p,delimiter=' ')[:,0]
        cev = np.loadtxt(p,delimiter=' ')[:,1]
        err = np.loadtxt(p,delimiter=' ')[:,2]
    
        temp_entry = { 'g_x':g_x, 'm_a':m_a, 'bins':bins, 'cev':cev, 'err':err}
        
        # add to dataframe
        df_pe = df_pe.append(temp_entry, ignore_index=True)
    

print(len(df_pe.index))


#%% save dataframe as a csv

df_pe.to_csv('bsm_grid_pe_specs_dataframe.csv')

#%% plotting to check df_pe
fig=plt.figure(3,figsize=(10,8))

# plotting the diff with gx and M
plt.scatter(df_pe['m_a'].values, df_pe['g_x'].values, c=df_pe.index,  vmax = 2, cmap = plt.get_cmap('Spectral'))
plt.ylabel('$g_x$ [GeV]', size=16)
plt.xlabel('$M_A$ [GeV]', size=16)
plt.xscale('log')
plt.yscale('log')

plt.show()
