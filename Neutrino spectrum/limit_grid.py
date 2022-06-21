#!/usr/bin/env python
# coding: utf-8

"""
limit_grid.py

Author: EB

17 May 2022

processing data for the limit grid

Converting specs of selected values to pe specs and writing to file

files are saved to Neutrino spectrum/output_grid

Edit: June 22 using diff to choose points to convert

"""


# In[40]: Imports and defining funcs


import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import matplotlib
import bisect
from func_my_cevns import spectope


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



# In[4]: LOADING IN DATA


# load standard model data - recoil energy spec
sm_er , sm_spec = np.loadtxt('data/argon_spec.txt', delimiter=' ')

# sm data - pe spec
mycenns = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[:,1]
mybins = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[:,0]
myerr = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[:,2]

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


# adding in diff

sm_index = bisect.bisect_left(sm_er, 0.1e-6)

index = bisect.bisect_left(df_MuTau.iloc[1].ERs, 0.1e-6)
df_MuTau['diff'] = df_MuTau.apply(lambda x: np.trapz(x.spec[index:], x.ERs[index:])/ np.trapz(sm_spec[sm_index:], sm_er[sm_index:]), axis=1)

diffs = df_MuTau['diff']


# #%% First value - takes a long time

# index = getindex(0.003, 0.002)
# # get the er_spec and energies
# spec = df_MuTau.spec[index]
# er = df_MuTau.ERs[index]
# # convert spec to pe
# cev, bins, err = spectope(er, spec)

# # saving to file just in case
# filename = 'first_value_pe_spec.txt'
# file = open(filename, 'w')

# for n in range(len(bins)-1):
#     file.write(str(bins[n]) + ' ' + str(cev[n]) + ' ' + str(err[n]) + '\n')
    
# file.close()

# # check that the values have been added
# print( df_MuTau[ (df_MuTau.pe_spec.notnull()) ] )

#%% Grid values
#Picking values of g_x and M_A for the grid

# # g values are between 1e-4 and 1e-2 - getting log distributed array of values
g_values = np.logspace(np.log10(1e-4),np.log10(1e-2),num=100)
# # m values are between 1e-3 and 1 - getting log distributed array of values
m_values = np.logspace(np.log10(1e-3),np.log10(1),num=100)


# # unpacking pairs of values I want to use for the grid
# g_values = np.loadtxt('grid_values.txt',delimiter=' ')[:,0]
# m_values = np.loadtxt('grid_values.txt',delimiter=' ')[:,1]



#%% Get values which haven't been processed
# index of values that I sampled less for for the sake of time
less_list = []
#TODO add if statement so I don't redo sapling for data I already have
import os.path
count = 0
# loop through all the values and add to the data frame
for i in range(len(g_values)):
    count += 1
    # get pair from list
    g = g_values[i]
    m = m_values[i]
    
    # get filename
    filename = 'output_grid/mutau_pe_spec_' + str(g) + '_' + str(m) + '.txt'

    # check if this has been processed before
    if os.path.exists(filename) == False:
        # file exists
        print('no file with this name ')
        # add index to list fmi
        less_list.append(i)
    else:
        print(' ')
        

#%% Loop through values
# index of values that I sampled less for for the sake of time
less_list = []
#TODO add if statement so I don't redo sapling for data I already have
import os.path
count = 0
# loop through all the values and add to the data frame
for i in range(len(g_values)):
    count += 1
    # get pair from list
    g = g_values[i]
    m = m_values[i]
    
    # get diff for these values of g and m
    df_index = getindex(g, m)
    # skip if diff isn't interesting
    if df_MuTau.loc[df_index, 'diff'] > 2.0 or df_MuTau.loc[df_index, 'diff'] < 1.1:
        continue
    
    # get filename
    filename = 'output_grid/mutau_pe_spec_' + f'{g:.3}' + '_' + f'{m:.3}' + '.txt'

    # check if this has been processed before
    if os.path.exists(filename) == False:

        # find the index corresponding to these values
        index = getindex(g, m)
        
        # get the er_spec and energies
        spec = df_MuTau.spec[index]
        er = df_MuTau.ERs[index]
        
        # convert spec to pe
        cev, bins, err = spectope(er, spec)
        
        # saving to file
        file = open(filename, 'w')
        for n in range(len(bins)-1):
            file.write(str(bins[n]) + ' ' + str(cev[n]) + ' ' + str(err[n]) + '\n')
        file.close()
        
        # add g and m values to grid_values file
        file = open('grid_values.txt', 'a')
        file.write(f'{g:.3}' + ' ' + f'{m:.3}' + '\n')
        file.close()
        
        # add index to list fmi
        less_list.append(index)
        




