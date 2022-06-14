#!/usr/bin/env python
# coding: utf-8

"""
limit_grid.py

Author: EB

17 May 2022

processing data for the limit grid

Converting specs of selected values to pe specs and writing to file

files are saved to Neutrino spectrum/output_grid



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

#%% Picking values of g_x and M_A for the grid

# # g values are between 1e-4 and 1e-2 - getting log distributed array of values
# g_values = np.logspace(np.log10(1e-4),np.log10(1e-2),num=3)
# # m values are between 1e-3 and 1 - getting log distributed array of values
# m_values = np.logspace(np.log10(1e-3),np.log10(1),num=3)


#%% Getting PE specs for the values that have been selected

print(df_MuTau)

# new empty column for values
df_MuTau['pe_spec'] = pd.NaT
df_MuTau['pe_bins'] = pd.NaT
df_MuTau['pe_err'] = pd.NaT

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


#%% First value

index = getindex(0.003, 0.002)
# get the er_spec and energies
spec = df_MuTau.spec[index]
er = df_MuTau.ERs[index]
# convert spec to pe
cev, bins, err = spectope(er, spec)

#%% saving to file just in case
filename = 'first_value_pe_spec.txt'
file = open(filename, 'w')

for n in range(len(bins)-1):
    file.write(str(bins[n]) + ' ' + str(cev[n]) + ' ' + str(err[n]) + '\n')
    
file.close()


#%% add to data frame

df_MuTau.iat[index,4] = cev
df_MuTau.iat[index,5] = bins
df_MuTau.iat[index,6] = err

# check that the values have been added
print( df_MuTau[ (df_MuTau.pe_spec.notnull()) ] )

#%% Grid values

# unpacking pairs of values I want to use for the grid
g_values = np.loadtxt('grid_values.txt',delimiter=' ')[:,0]
m_values = np.loadtxt('grid_values.txt',delimiter=' ')[:,1]



# #%%

# # plotting to see if its correct
# fig=plt.figure(3,figsize=(10,8))

# # plotting the diff with gx and M
# plt.plot(m_values, g_values,'o')
# plt.ylabel('$g_x$ [GeV]', size=16)
# plt.xlabel('$M_A$ [GeV]', size=16)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim(1e-3,1)
# plt.ylim(1e-4,1e-2)
# plt.show()


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
    
    print(g, m)
    
    # get filename
    filename = 'output_grid/mutau_pe_spec_' + str(g) + '_' + str(m) + '.txt'

    # check if this has been processed before
    if os.path.exists(filename) == False:

        # find the index corresponding to these values
        index = getindex(g, m)
        
        # get the er_spec and energies
        spec = df_MuTau.spec[index]
        er = df_MuTau.ERs[index]
        
        # convert spec to pe
        cev, bins, err = spectope(er, spec)
        
        # saving to file just in case
        file = open(filename, 'w')
        for n in range(len(bins)-1):
            file.write(str(bins[n]) + ' ' + str(cev[n]) + ' ' + str(err[n]) + '\n')
        file.close()
        
        # add to data frame
        df_MuTau.iat[index,4] = cev
        df_MuTau.iat[index,5] = bins
        df_MuTau.iat[index,6] = err
        
        # add index to list fmi
        less_list.append(index)
        

        








