#!/usr/bin/env python
# coding: utf-8

"""
lmultau_plotting.py

Plotting neutrino specs in kev

NOTE: moved to using a juypter notebook

"""


# In[40]:


import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import matplotlib
import bisect

from func_my_cevns import spectope

# In[4]: LOADING IN DATA


# load standard model data
sm_er , sm_spec = np.loadtxt('data/argon_spec.txt', delimiter=' ')

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


#%% plotting one spec
fig=plt.figure(1,figsize=(10,8))

plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)
plt.loglog(df_MuTau[(df_MuTau['g_x']==0.01) & (df_MuTau['m_A']==1.0)]['ERs'].values[0]*1e6, df_MuTau[(df_MuTau['g_x']==0.01) & (df_MuTau['m_A']==1.0)]['spec'].values[0])
plt.loglog(sm_er*1e6,sm_spec)


# In[11]: GETTING THE RATIO OF BSM TO SM


# Adding diff column to data frame
# using traps for the bsm specs then divding by the value for the sm - gives an indication of the difference between sma nd bsm
# find index above threshold for both and compare above 0.1 kev as this is what will make the difference in the PE plot
sm_index = bisect.bisect_left(sm_spec, 0.1e-6)
index = bisect.bisect_left(df_MuTau.iloc[1].ERs, 0.1e-6)
df_MuTau['diff'] = df_MuTau.apply(lambda x: np.trapz(x.spec[index:], x.ERs[index:])/ np.trapz(sm_spec[sm_index:], sm_er[sm_index:]), axis=1)



# In[54]: PLOTTING SM VS SOMETHING ELSE
fig=plt.figure(2,figsize=(10,8))

plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)
plt.loglog(df_MuTau.iloc[1998].ERs*1e6, df_MuTau.iloc[19998].spec)
plt.loglog(sm_er*1e6, sm_spec) # sm only from 1e-2 whereas other spec are from 1e-3 - won't make a difference for PE spec!


# In[50]: PLOTTING THE DIFF COLOUR SCATTER PLOT
fig=plt.figure(3,figsize=(10,8))

# plotting the diff with gx and M
f=plt.figure(figsize=(15,10))
plt.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['diff'].values, norm=matplotlib.colors.LogNorm(), vmax = 2)
plt.ylabel('$g_x$ [GeV]')
plt.xlabel('$m_A$ [GeV]')
plt.xscale('log')
plt.yscale('log')
cbar = plt.colorbar()
cbar.set_label('BSM/SM', rotation=270)
plt.show()

#%% PLOTTING THE INDEX

fig=plt.figure(0,figsize=(10,8))

# plotting the diff with gx and M
plt.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau.index)
plt.ylabel('$g_x$ [GeV]')
plt.xlabel('$m_A$ [GeV]')
plt.xscale('log')
plt.yscale('log')
cbar = plt.colorbar()
cbar.set_label('index', rotation=270)
plt.show()




# In[29]: PLOTTING SM VS A VARIETY - USING INDEX


# Picking some points to plot - plot here and do for PE
f=plt.figure(4,figsize=(10,8))
plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)

plt.loglog(df_MuTau.iloc[1998].ERs*1e6, df_MuTau.iloc[19998].spec, label = 'g_x = ' + str(df_MuTau.iloc[1998].g_x) + ', m_A = ' + str(df_MuTau.iloc[1998].m_A))

plt.loglog(df_MuTau.iloc[100].ERs*1e6, df_MuTau.iloc[100].spec, label = 'g_x = ' + str(df_MuTau.iloc[100].g_x) + ', m_A = ' + str(df_MuTau.iloc[100].m_A))

plt.loglog(df_MuTau.iloc[19995].ERs*1e6, df_MuTau.iloc[19995].spec, label = 'g_x = ' + str(df_MuTau.iloc[19995].g_x) + ', m_A = ' + str(df_MuTau.iloc[19995].m_A))

plt.loglog(sm_er*1e6, sm_spec, label = 'SM') # sm only from 1e-2 whereas other spec are from 1e-3 - won't make a difference for PE spec!
plt.legend(fontsize=12,frameon=False,loc='upper right')




# In[ ]: SELECTING VALUES
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
     
# selecting values - starred points
m = 1e-2
a = getindex(2e-3,m)
b = getindex(4e-4,m)
c = getindex(1.5e-4,m)
m = 1e-1
d = getindex(2e-4, m)
m = 1
e = getindex(2e-3, m)
f = getindex(4e-4, m)

print(a,b,c,d,e,f)

#%%
# plotting

f=plt.figure(5,figsize=(10,8))
plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)

f = 1800

plt.loglog(df_MuTau.iloc[a].ERs*1e6, df_MuTau.iloc[a].spec, label = 'g_x = ' + str(df_MuTau.iloc[a].g_x) + ', m_A = ' + str(df_MuTau.iloc[a].m_A))
plt.loglog(df_MuTau.iloc[b].ERs*1e6, df_MuTau.iloc[b].spec, label = 'g_x = ' + str(df_MuTau.iloc[b].g_x) + ', m_A = ' + str(df_MuTau.iloc[b].m_A))
plt.loglog(df_MuTau.iloc[c].ERs*1e6, df_MuTau.iloc[c].spec, label = 'g_x = ' + str(df_MuTau.iloc[c].g_x) + ', m_A = ' + str(df_MuTau.iloc[c].m_A))
plt.loglog(df_MuTau.iloc[d].ERs*1e6, df_MuTau.iloc[d].spec, label = 'g_x = ' + str(df_MuTau.iloc[d].g_x) + ', m_A = ' + str(df_MuTau.iloc[d].m_A))
plt.loglog(df_MuTau.iloc[e].ERs*1e6, df_MuTau.iloc[e].spec, label = 'g_x = ' + str(df_MuTau.iloc[e].g_x) + ', m_A = ' + str(df_MuTau.iloc[e].m_A))
plt.loglog(df_MuTau.iloc[f].ERs*1e6, df_MuTau.iloc[f].spec, label = 'g_x = ' + str(df_MuTau.iloc[f].g_x) + ', m_A = ' + str(df_MuTau.iloc[f].m_A))


plt.legend(fontsize=12,frameon=False,loc='upper right')


# In[ ]: CONVERTING TO PE


# USING MY CONVERSION TP PE
acev, abins = spectope(df_MuTau.iloc[a].ERs, df_MuTau.iloc[a].spec)
bcev, bbins = spectope(df_MuTau.iloc[b].ERs, df_MuTau.iloc[b].spec)
ccev, cbins = spectope(df_MuTau.iloc[c].ERs, df_MuTau.iloc[c].spec)
dcev, dbins = spectope(df_MuTau.iloc[d].ERs, df_MuTau.iloc[d].spec)
ecev, ebins = spectope(df_MuTau.iloc[e].ERs, df_MuTau.iloc[e].spec)
fcev, fbins = spectope(df_MuTau.iloc[1800].ERs, df_MuTau.iloc[1800].spec) # weird problem with f


#%%

# plotting
f=plt.figure(6,figsize=(10,8))


# =============================================================================
# NOT WORKING!!
# =============================================================================
plt.loglog(abins, acev, '-+',markersize=15, label = 'g_x = ' + str(df_MuTau.iloc[a].g_x) + ', m_A = ' + str(df_MuTau.iloc[a].m_A))
plt.loglog(bbins, bcev, '-+',markersize=15, label = 'g_x = ' + str(df_MuTau.iloc[b].g_x) + ', m_A = ' + str(df_MuTau.iloc[b].m_A))
plt.loglog(cbins, ccev, '-+',markersize=15, label = 'g_x = ' + str(df_MuTau.iloc[c].g_x) + ', m_A = ' + str(df_MuTau.iloc[c].m_A))
plt.loglog(dbins, dcev, '-+',markersize=15, label = 'g_x = ' + str(df_MuTau.iloc[d].g_x) + ', m_A = ' + str(df_MuTau.iloc[d].m_A))
plt.loglog(ebins, ecev, '-+',markersize=15, label = 'g_x = ' + str(df_MuTau.iloc[e].g_x) + ', m_A = ' + str(df_MuTau.iloc[e].m_A))
plt.loglog(fbins, fcev, '-+',markersize=15, label = 'g_x = ' + str(df_MuTau.iloc[1800].g_x) + ', m_A = ' + str(df_MuTau.iloc[1800].m_A)) # weird problem with f


plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Events per tyr',fontsize=26) 
plt.yscale('log')
plt.legend(fontsize=18,frameon=False,loc='upper right')

