#!/usr/bin/env python
# coding: utf-8

"""
limit_grid.py

Author: EB

Aug 2022

Extract chosen ER spec files in order to process them on the grid

"""


#%% Imports and defining funcs


import numpy as np
import pandas as pd
import glob
import os
import csv

def getindex(g, m):
    '''
    func to help find rows in the dataframe which have the values of g and m 
    we are looking for
    '''
    # selecting the firt row that meets the condtions
    # allowing 5% either side
    mult = 1.05
    mult2 = 0.95
    if g < 1e-3: mult = 1.1; mult2 = 0.99

    index = df_MuTau[(df_MuTau['g_x'] < mult*g) & (df_MuTau['g_x'] > mult2*g) & (df_MuTau['m_A']< 1.05*m) & (df_MuTau['m_A'] > 0.95*m)].index[0]
    return index



#%% LOADING IN DATA

# Loading all the data into a pandas dataframe

LMuTaupaths = glob.glob('input/*')

df_MuTau = pd.DataFrame()

for p in LMuTaupaths:
    filepath = p + '/spectra_Ar.npy'
    ERs, g_x, m_A, spec=np.load(filepath, allow_pickle=True)
    for i in range(len(m_A)):
        temp_entry = {'ERs': ERs, 'g_x':g_x, 'm_A':m_A[i], 'spec':spec[i]}
        #print(temp_entry)
        df_MuTau = df_MuTau.append(temp_entry, ignore_index=True)


#%% Grid values
#Picking values of g_x and M_A for the grid

# # g values are between 1e-4 and 1e-2 - getting log distributed array of values
g_values = np.logspace(np.log10(1e-4),np.log10(1e-2),num=100)
# # m values are between 1e-3 and 1 - getting log distributed array of values
m_values = np.logspace(np.log10(1e-3),np.log10(1),num=100)

# save lists as strings to 3 dp
for g in g_values:
    # write to list of .3 g_values
    file = open('grid_jobs_gvalues.txt', 'a')
    file.write(f'{g:.3}' + '\n')
    file.close()
    
for m in m_values:
    # write to list of .3 g_values
    file = open('grid_jobs_mvalues.txt', 'a')
    file.write(f'{m:.3}' + '\n')
    file.close()
    

#%% Loop through values


# loop through g values
for g in g_values:
    
    # create directroy
    g_dir_name = 'g_' + f'{g:.3}'
    parent_dir = 'grid_jobs_input/'
    path = os.path.join(parent_dir, g_dir_name)
    os.mkdir(path)
    
    
    # # write to list of .3 g_values
    # file = open('atest_grid_jobs_gvalues.txt', 'a')
    # file.write(f'{g:.3}' + '\n')
    # file.close()
    
    for m in m_values:
    
        # get filename
        er_filename = 'mutau_er_spec_g_' + f'{g:.3}' + '_m_' + f'{m:.3}' + '.txt'

        # find the index corresponding to these values
        index = getindex(g, m)
        
        # get the er_spec and energies
        spec = df_MuTau.spec[index]
        er = df_MuTau.ERs[index]
        
        filename = 'grid_jobs_input/' + g_dir_name + '/' + er_filename
        
        # saving to file
        file = open(filename, 'w')
        writer = csv.writer(file, delimiter=' ')
        writer.writerows(zip(er,spec))
        file.close()
        
