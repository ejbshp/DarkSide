#!/usr/bin/env python
# coding: utf-8

"""
lmultau_plotting.py

Plotting neutrino specs in kev

NOTE: moved to using a juypter notebook

"""


# In[40]: Imports and defininf funcs


import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import matplotlib
import bisect
from matplotlib.cm import ScalarMappable

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

def getratio(index):
    '''
    returns ratio of bsm/sm for a given row ion df
    
    index: index of row in df

    '''
    # points at which we want to sample
    x_sample = sm_er # sm spec starts at a higher ER
    ratio = np.interp(x_sample, df_MuTau.iloc[index].ERs, df_MuTau.iloc[index].spec) / sm_spec
    
    return ratio

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
sm_index = bisect.bisect_left(sm_er, 0.1e-6)

index = bisect.bisect_left(df_MuTau.iloc[1].ERs, 0.1e-6)
df_MuTau['diff'] = df_MuTau.apply(lambda x: np.trapz(x.spec[index:], x.ERs[index:])/ np.trapz(sm_spec[sm_index:], sm_er[sm_index:]), axis=1)

diffs = df_MuTau['diff']

# In[54]: PLOTTING SM VS SOMETHING ELSE
fig=plt.figure(2,figsize=(10,8))

plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)
plt.loglog(df_MuTau.iloc[1998].ERs*1e6, df_MuTau.iloc[19998].spec)
plt.loglog(sm_er*1e6, sm_spec) # sm only from 1e-2 whereas other spec are from 1e-3 - won't make a difference for PE spec!


# In[50]: PLOTTING THE DIFF COLOUR SCATTER PLOT
fig=plt.figure(3,figsize=(10,8))

# plotting the diff with gx and M
plt.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['diff'].values,  vmax = 2, cmap = plt.get_cmap('Spectral'))
plt.ylabel('$g_x$ [GeV]', size=16)
plt.xlabel('$M_A$ [GeV]', size=16)
plt.xscale('log')
plt.yscale('log')

fmt = matplotlib.ticker.FormatStrFormatter('%1.1f')
cbar = plt.colorbar(format = fmt)
cbar.set_label(r'$\frac{(BSM + SM )}{SM}$', size=20)
plt.show()

#%% PLOTTING THE INDEX

fig=plt.figure(4,figsize=(10,8))

# plotting the diff with gx and M
plt.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau.index)
plt.ylabel('$g_x$ [GeV]')
plt.xlabel('$M_A$ [GeV]')
plt.xscale('log')
plt.yscale('log')
cbar = plt.colorbar()
cbar.set_label('Index')
plt.show()




#%% M1 PLOTTING FOR FIXED MASS

fixed_mass = 4e-3 # GeV

a = getindex(1e-4,fixed_mass)
b = getindex(2e-4,fixed_mass)
c = getindex(8e-4,fixed_mass)
d = getindex(1e-3, fixed_mass)
e = getindex(5e-3, fixed_mass)
g = getindex(1e-2, fixed_mass)


# plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25,10))
fig.suptitle('Fixed mass: m_A =' + str(fixed_mass) + ' [GeV]', fontsize = 20)

# Er spec for the 6 points
ax1.set_xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
ax1.set_ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)


ax1.loglog(df_MuTau.iloc[a].ERs*1e6, df_MuTau.iloc[a].spec, label = '$g_x=$' + str(df_MuTau.iloc[a].g_x), color ='C0' )
ax1.loglog(df_MuTau.iloc[b].ERs*1e6, df_MuTau.iloc[b].spec, label = '$g_x=$' + str(df_MuTau.iloc[b].g_x), color ='C1' )
ax1.loglog(df_MuTau.iloc[c].ERs*1e6, df_MuTau.iloc[c].spec, label = '$g_x=$' + str(df_MuTau.iloc[c].g_x), color ='C2' )
ax1.loglog(df_MuTau.iloc[d].ERs*1e6, df_MuTau.iloc[d].spec, label = '$g_x=$' + str(df_MuTau.iloc[d].g_x), color ='C3' )
ax1.loglog(df_MuTau.iloc[e].ERs*1e6, df_MuTau.iloc[e].spec, label = '$g_x=$' + str(df_MuTau.iloc[e].g_x), color ='C4' )
ax1.loglog(df_MuTau.iloc[g].ERs*1e6, df_MuTau.iloc[g].spec, label = '$g_x=$' + str(df_MuTau.iloc[g].g_x), color ='C5' )

ax1.legend(fontsize=12,frameon=False,loc='lower left')

# Color graph of all points
cmap = plt.get_cmap('Spectral')
ax2.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['diff'].values, vmax = 2, cmap=cmap)
ax2.set_ylabel('$g_x$ [GeV]', size=16)
ax2.set_xlabel('$M_A$ [GeV]', size=16)
ax2.set_xscale('log')
ax2.set_yscale('log')

fmt = matplotlib.ticker.FormatStrFormatter('%1.1f')
sm =  ScalarMappable(norm=matplotlib.colors.Normalize(vmin = 1, vmax = 2), cmap=cmap)
cbar = fig.colorbar(sm, ax=ax2, format=fmt)
cbar.set_label(r'$\frac{(BSM + SM )}{SM}$', size=20)




# plotting points on colour graph
ax2.plot(fixed_mass, df_MuTau.iloc[a].g_x, marker='o', markersize=10, color ='C0' )
ax2.plot(fixed_mass, df_MuTau.iloc[b].g_x, marker='o', markersize=10, color ='C1' )
ax2.plot(fixed_mass, df_MuTau.iloc[c].g_x, marker='o', markersize=10, color ='C2' )
ax2.plot(fixed_mass, df_MuTau.iloc[d].g_x, marker='o', markersize=10, color ='C3' )
ax2.plot(fixed_mass, df_MuTau.iloc[e].g_x, marker='o', markersize=10, color ='C4' )
ax2.plot(fixed_mass, df_MuTau.iloc[g].g_x, marker='o', markersize=10, color ='C5' )

plt.show()


#%% M2 PLOTTING BSM/SM vs ER


    
f=plt.figure(50,figsize=(10,8))

# different lengths sampled at different points

# line at 1
plt.axhline(y=1, color='r', linestyle='--')

#plt.ylim(0, 10)

plt.plot(sm_er*1e6, getratio(a) , label = '$g_x=$' + str(df_MuTau.iloc[a].g_x), color ='C0')
plt.plot(sm_er*1e6, getratio(b) , label = '$g_x=$' + str(df_MuTau.iloc[b].g_x), color ='C1')
plt.plot(sm_er*1e6, getratio(c) , label = '$g_x=$' + str(df_MuTau.iloc[c].g_x), color ='C2')
plt.plot(sm_er*1e6, getratio(d) , label = '$g_x=$' + str(df_MuTau.iloc[d].g_x), color ='C3')
plt.plot(sm_er*1e6, getratio(e) , label = '$g_x=$' + str(df_MuTau.iloc[e].g_x), color ='C4')
plt.plot(sm_er*1e6, getratio(g) , label = '$g_x=$' + str(df_MuTau.iloc[g].g_x), color ='C5')

plt.title('Fixed mass: m_A = ' + str(fixed_mass) + ' [GeV]')
plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\frac{(BSM + SM )}{SM}$', size=20)

plt.xscale('log')
plt.yscale('log')

plt.legend(fontsize=12,frameon=False,loc='upper right')

#%% G1 PLOTTING FOR FIXED COUPLING

fixed_g = 3e-4 # GeV

a = getindex(fixed_g, 1e-3)
b = getindex(fixed_g, 5e-3)
c = getindex(fixed_g, 1e-2)
d = getindex(fixed_g, 5e-2)
e = getindex(fixed_g, 1e-1)
g = getindex(fixed_g, 1)


# plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25,10))
fig.suptitle('Fixed coupling: g_x =' + str(fixed_g) + ' [GeV]', fontsize = 20)

# Er spec for the 6 points
ax1.set_xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
ax1.set_ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)


ax1.loglog(df_MuTau.iloc[a].ERs*1e6, df_MuTau.iloc[a].spec, label = ' m_A = ' + str(df_MuTau.iloc[a].m_A), color ='C0' )
ax1.loglog(df_MuTau.iloc[b].ERs*1e6, df_MuTau.iloc[b].spec, label = ' m_A = ' + str(df_MuTau.iloc[b].m_A), color ='C1' )
ax1.loglog(df_MuTau.iloc[c].ERs*1e6, df_MuTau.iloc[c].spec, label = ' m_A = ' + str(df_MuTau.iloc[c].m_A), color ='C2' )
ax1.loglog(df_MuTau.iloc[d].ERs*1e6, df_MuTau.iloc[d].spec, label = ' m_A = ' + str(df_MuTau.iloc[d].m_A), color ='C3' )
ax1.loglog(df_MuTau.iloc[e].ERs*1e6, df_MuTau.iloc[e].spec, label = ' m_A = ' + str(df_MuTau.iloc[e].m_A), color ='C4' )
ax1.loglog(df_MuTau.iloc[g].ERs*1e6, df_MuTau.iloc[g].spec, label = ' m_A = ' + str(df_MuTau.iloc[g].m_A), color ='C5' )

ax1.legend(fontsize=12,frameon=False,loc='lower left')

# Color graph of all points
cmap = plt.get_cmap('Spectral')
ax2.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['diff'].values, vmax = 2, cmap=cmap)
ax2.set_ylabel('$g_x$ [GeV]', size=16)
ax2.set_xlabel('$M_A$ [GeV]', size=16)
ax2.set_xscale('log')
ax2.set_yscale('log')

fmt = matplotlib.ticker.FormatStrFormatter('%1.1f')
sm =  ScalarMappable(norm=matplotlib.colors.Normalize(vmin = 1, vmax = 2), cmap=cmap)
cbar = fig.colorbar(sm, ax=ax2, format=fmt)
cbar.set_label(r'$\frac{(BSM + SM )}{SM}$', size=20)

# plotting points on colour graph
ax2.plot(df_MuTau.iloc[a].m_A, df_MuTau.iloc[a].g_x, marker='o', markersize=10, color ='C0' )
ax2.plot(df_MuTau.iloc[b].m_A, df_MuTau.iloc[b].g_x, marker='o', markersize=10, color ='C1' )
ax2.plot(df_MuTau.iloc[c].m_A, df_MuTau.iloc[c].g_x, marker='o', markersize=10, color ='C2' )
ax2.plot(df_MuTau.iloc[d].m_A, df_MuTau.iloc[d].g_x, marker='o', markersize=10, color ='C3' )
ax2.plot(df_MuTau.iloc[e].m_A, df_MuTau.iloc[e].g_x, marker='o', markersize=10, color ='C4' )
ax2.plot(df_MuTau.iloc[g].m_A, df_MuTau.iloc[g].g_x, marker='o', markersize=10, color ='C5' )

plt.show()


#%% G2 Plotting ratio for fixed coupling
f=plt.figure(51,figsize=(10,8))
# line at 1
plt.axhline(y=1, color='r', linestyle='--')



plt.plot(sm_er*1e6, getratio(a) , label = ' m_A = ' + str(df_MuTau.iloc[a].m_A), color ='C0')
plt.plot(sm_er*1e6, getratio(b) , label = ' m_A = ' + str(df_MuTau.iloc[b].m_A), color ='C1')
plt.plot(sm_er*1e6, getratio(c) , label = ' m_A = ' + str(df_MuTau.iloc[c].m_A), color ='C2')
plt.plot(sm_er*1e6, getratio(d) , label = ' m_A = ' + str(df_MuTau.iloc[d].m_A), color ='C3')
plt.plot(sm_er*1e6, getratio(e) , label = ' m_A = ' + str(df_MuTau.iloc[e].m_A), color ='C4')
plt.plot(sm_er*1e6, getratio(g) , label = ' m_A = ' + str(df_MuTau.iloc[g].m_A), color ='C5')

plt.title('Fixed coupling: g_x = ' + str(fixed_g) + ' [GeV]')
plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\frac{(BSM + SM )}{SM}$', size=20)

plt.xscale('log')
plt.yscale('log')

plt.legend(fontsize=12,frameon=False,loc='upper right')





# In[ ]: SELECTING VALUES  - BENCHMARK POINTS
     
# selecting values - starred points
m = 1e-2
a = getindex(2e-3,m)
b = getindex(4e-4,m)
c = getindex(1.5e-4,m)
m = 1e-1
d = getindex(2e-4, m)
m = 1
e = getindex(2e-3, m)
g = getindex(4e-4, m)

print(a,b,c,d,e,g)

#%%

# plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25,10))
fig.suptitle('Benchmark Points', size=20)

# Er spec for the 6 points
ax1.set_xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
ax1.set_ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)


ax1.loglog(df_MuTau.iloc[a].ERs*1e6, df_MuTau.iloc[a].spec, label = '$g_x=$' + str(df_MuTau.iloc[a].g_x) + '$M_A=$' + str(df_MuTau.iloc[a].m_A), color ='C0' )
ax1.loglog(df_MuTau.iloc[b].ERs*1e6, df_MuTau.iloc[b].spec, label = '$g_x=$' + str(df_MuTau.iloc[b].g_x) + '$M_A=$' + str(df_MuTau.iloc[b].m_A), color ='C1' )
ax1.loglog(df_MuTau.iloc[c].ERs*1e6, df_MuTau.iloc[c].spec, label = '$g_x=$' + str(df_MuTau.iloc[c].g_x) + '$M_A=$' + str(df_MuTau.iloc[c].m_A), color ='C2' )
ax1.loglog(df_MuTau.iloc[d].ERs*1e6, df_MuTau.iloc[d].spec, label = '$g_x=$' + str(df_MuTau.iloc[d].g_x) + '$M_A=$' + str(df_MuTau.iloc[d].m_A), color ='C3' )
ax1.loglog(df_MuTau.iloc[e].ERs*1e6, df_MuTau.iloc[e].spec, label = '$g_x=$' + str(df_MuTau.iloc[e].g_x) + '$M_A=$' + str(df_MuTau.iloc[e].m_A), color ='C4' )
ax1.loglog(df_MuTau.iloc[g].ERs*1e6, df_MuTau.iloc[g].spec, label = '$g_x=$' + str(df_MuTau.iloc[g].g_x) + '$M_A=$' + str(df_MuTau.iloc[g].m_A), color ='C5' )

ax1.legend(fontsize=12,frameon=False,loc='lower left')

# Color graph of all points
cmap = plt.get_cmap('Spectral')
ax2.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['diff'].values, vmax = 2, cmap=cmap)
ax2.set_ylabel('$g_x$ [GeV]', size=16)
ax2.set_xlabel('$M_A$ [GeV]', size=16)
ax2.set_xscale('log')
ax2.set_yscale('log')

fmt = matplotlib.ticker.FormatStrFormatter('%1.1f')
sm =  ScalarMappable(norm=matplotlib.colors.Normalize(vmin = 1, vmax = 2), cmap=cmap)
cbar = fig.colorbar(sm, ax=ax2, format=fmt)
cbar.set_label(r'$\frac{(BSM + SM )}{SM}$', size=20)

# plotting points on colour graph
ax2.plot(df_MuTau.iloc[a].m_A, df_MuTau.iloc[a].g_x, marker='o', markersize=10, color ='C0' )
ax2.plot(df_MuTau.iloc[b].m_A, df_MuTau.iloc[b].g_x, marker='o', markersize=10, color ='C1' )
ax2.plot(df_MuTau.iloc[c].m_A, df_MuTau.iloc[c].g_x, marker='o', markersize=10, color ='C2' )
ax2.plot(df_MuTau.iloc[d].m_A, df_MuTau.iloc[d].g_x, marker='o', markersize=10, color ='C3' )
ax2.plot(df_MuTau.iloc[e].m_A, df_MuTau.iloc[e].g_x, marker='o', markersize=10, color ='C4' )
ax2.plot(df_MuTau.iloc[g].m_A, df_MuTau.iloc[g].g_x, marker='o', markersize=10, color ='C5' )

plt.show()


#%% RATIO OF BENCHMARK POINTS

f=plt.figure(51,figsize=(10,8))
# line at 1
plt.axhline(y=1, color='r', linestyle='--')



plt.plot(sm_er*1e6, getratio(a) , label = '$g_x=$' + str(df_MuTau.iloc[a].g_x) + '$M_A=$' + str(df_MuTau.iloc[a].m_A), color ='C0' )
plt.plot(sm_er*1e6, getratio(b) , label = '$g_x=$' + str(df_MuTau.iloc[b].g_x) + '$M_A=$' + str(df_MuTau.iloc[b].m_A), color ='C1' )
plt.plot(sm_er*1e6, getratio(c) , label = '$g_x=$' + str(df_MuTau.iloc[c].g_x) + '$M_A=$' + str(df_MuTau.iloc[c].m_A), color ='C2' )
plt.plot(sm_er*1e6, getratio(d) , label = '$g_x=$' + str(df_MuTau.iloc[d].g_x) + '$M_A=$' + str(df_MuTau.iloc[d].m_A), color ='C3' )
plt.plot(sm_er*1e6, getratio(e) , label = '$g_x=$' + str(df_MuTau.iloc[e].g_x) + '$M_A=$' + str(df_MuTau.iloc[e].m_A), color ='C4' )
plt.plot(sm_er*1e6, getratio(g) , label = '$g_x=$' + str(df_MuTau.iloc[g].g_x) + '$M_A=$' + str(df_MuTau.iloc[g].m_A), color ='C5' )

plt.title('Ratio: Benchmark points', size=20)
plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\frac{(BSM + SM )}{SM}$', size=20)

plt.xscale('log')
plt.yscale('log')

plt.legend(fontsize=12,frameon=False,loc='center left')





# In[ ]: CONVERTING TO PE a

acev, abins, aerr = spectope(df_MuTau.iloc[a].ERs, df_MuTau.iloc[a].spec)




# In[ ]: CONVERTING TO PE b

bcev, bbins, berr = spectope(df_MuTau.iloc[b].ERs, df_MuTau.iloc[b].spec)



# In[ ]: CONVERTING TO PE c

ccev, cbins, cerr = spectope(df_MuTau.iloc[c].ERs, df_MuTau.iloc[c].spec)


# In[ ]: CONVERTING TO PE d

dcev, dbins, derr = spectope(df_MuTau.iloc[d].ERs, df_MuTau.iloc[d].spec)


# In[ ]: CONVERTING TO PE e

ecev, ebins, eerr = spectope(df_MuTau.iloc[e].ERs, df_MuTau.iloc[e].spec)


# In[ ]: CONVERTING TO PE g

gcev, gbins, gerr = spectope(df_MuTau.iloc[g].ERs, df_MuTau.iloc[g].spec)


#%% PE PLOTS BENCHMARK

# plotting

f, ax = plt.subplots(figsize=(15,10))

ax.set_yscale('log')
ax.set_xlabel('Number of electrons',fontsize=26)
ax.set_ylabel('Events per tyr',fontsize=26) 

ax.errorbar(abins, acev,yerr=aerr*acev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3, label = r'$\frac{N_{Events}}{\sqrt{N_{Samples}}}$')
ax.plot(abins, acev, label = f'$g_x$ =  {df_MuTau.iloc[a].g_x:.4f}, $M_A$ = {df_MuTau.iloc[a].m_A:.4f}')
ax.errorbar(bbins, bcev,yerr=berr*bcev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(bbins, bcev, label = f'$g_x$ =  {df_MuTau.iloc[b].g_x:.4f}, $M_A$ = {df_MuTau.iloc[b].m_A:.4f}')
ax.errorbar(cbins, ccev,yerr=cerr*ccev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(cbins, ccev, label = f'$g_x$ =  {df_MuTau.iloc[c].g_x:.4f}, $M_A$ = {df_MuTau.iloc[c].m_A:.4f}')
ax.errorbar(dbins, dcev,yerr=derr*dcev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(dbins, dcev, label = f'$g_x$ =  {df_MuTau.iloc[d].g_x:.4f}, $M_A$ = {df_MuTau.iloc[d].m_A:.4f}')
ax.errorbar(ebins, ecev,yerr=eerr*ecev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(ebins, ecev, label = f'$g_x$ =  {df_MuTau.iloc[e].g_x:.4f}, $M_A$ = {df_MuTau.iloc[e].m_A:.4f}')
ax.errorbar(gbins, gcev,yerr=gerr*gcev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(gbins, gcev, label = f'$g_x$ =  {df_MuTau.iloc[g].g_x:.4f}, $M_A$ = {df_MuTau.iloc[g].m_A:.4f}')

ax.set_xlim(0,65)

ax.set_title('Benchmark Points', size=20)

ax.legend(fontsize=18,frameon=False,loc='upper right')



#%% Comparing to SM spec in PE

mycenns = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[:,1]
mybins = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[:,0]
myerr = np.loadtxt('/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt',delimiter=' ')[:,2]

#%% Plotting ratio of pe with sm


f, ax = plt.subplots(figsize=(15,10))

ax.set_yscale('log')
ax.set_xlabel('Number of electrons',fontsize=26)
ax.set_ylabel(r'$\frac{(BSM + SM )}{SM}$', size=26)

ax.errorbar(abins, acev/mycenns,yerr=np.sqrt(aerr**2+myerr**2)*acev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3, label = r'$\frac{N_{Events}}{\sqrt{N_{Samples}}}$')
ax.plot(abins, acev/mycenns, '-',markersize=15, label = f'$g_x$ =  {df_MuTau.iloc[a].g_x:.4f}, $M_A$ = {df_MuTau.iloc[a].m_A:.4f}')
ax.errorbar(bbins, bcev/mycenns,yerr=np.sqrt(berr**2+myerr**2)*bcev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(bbins, bcev/mycenns, '-',markersize=15, label = f'$g_x$ =  {df_MuTau.iloc[b].g_x:.4f}, $M_A$ = {df_MuTau.iloc[b].m_A:.4f}')
ax.errorbar(cbins, ccev/mycenns,yerr=np.sqrt(cerr**2+myerr**2)*ccev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(cbins, ccev/mycenns, '-',markersize=15, label = f'$g_x$ =  {df_MuTau.iloc[c].g_x:.4f}, $M_A$ = {df_MuTau.iloc[c].m_A:.4f}')
ax.errorbar(dbins, dcev/mycenns,yerr=np.sqrt(derr**2+myerr**2)*dcev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(dbins, dcev/mycenns, '-',markersize=15, label = f'$g_x$ =  {df_MuTau.iloc[d].g_x:.4f}, $M_A$ = {df_MuTau.iloc[d].m_A:.4f}')
ax.errorbar(ebins, ecev/mycenns,yerr=np.sqrt(eerr**2+myerr**2)*ecev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(ebins, ecev/mycenns, '-',markersize=15, label = f'$g_x$ =  {df_MuTau.iloc[e].g_x:.4f}, $M_A$ = {df_MuTau.iloc[e].m_A:.4f}')
ax.errorbar(gbins, gcev/mycenns,yerr=np.sqrt(gerr**2+myerr**2)*gcev, fmt='o', capsize=3, color='k',linewidth=2, markersize=3)
ax.plot(gbins, gcev/mycenns, '-',markersize=15, label = f'$g_x$ =  {df_MuTau.iloc[g].g_x:.4f}, $M_A$ = {df_MuTau.iloc[g].m_A:.4f}')

ax.set_xlim(0,50)

ax.set_title('Ratio: Benchmark Points', size=20)

ax.legend(fontsize=18,frameon=False,loc='upper right')


#%% P1 PLOTTING ER SPEC AND COLOUR PLOT SIDE BY SIDE - CHOOSE 6 POINTS

# select points

a = getindex(1e-4,1e-3)
b = getindex(2e-4,4e-3)
c = getindex(8e-4,2e-2)
d = getindex(4e-3, 1e-1)
e = getindex(1e-2, 0.4)
g = getindex(4e-4, 1e-2)


# plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25,10))

# Er spec for the 6 points
ax1.set_xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
ax1.set_ylabel(r'$\rm{d}R/\rm{d}E_R\,\,\,\left[\rm{keV}^{-1}\,\,\rm{ton}^{-1}\,\,\rm{yr}^{-1}\right]$', size=16)

ax1.loglog(sm_er*1e6, sm_spec, '--' ,label=r'SM CE$\nu$NS', color='C6')
ax1.loglog(df_MuTau.iloc[a].ERs*1e6, df_MuTau.iloc[a].spec, label = f'$g_x$ =  {df_MuTau.iloc[a].g_x:.4f}, $M_A$ = {df_MuTau.iloc[a].m_A:.4f}', color ='C0' )
ax1.loglog(df_MuTau.iloc[b].ERs*1e6, df_MuTau.iloc[b].spec, label = f'$g_x$ =  {df_MuTau.iloc[b].g_x:.4f}, $M_A$ = {df_MuTau.iloc[b].m_A:.4f}', color ='C1' )
ax1.loglog(df_MuTau.iloc[c].ERs*1e6, df_MuTau.iloc[c].spec, label = f'$g_x$ =  {df_MuTau.iloc[c].g_x:.4f}, $M_A$ = {df_MuTau.iloc[c].m_A:.4f}', color ='C2' )
ax1.loglog(df_MuTau.iloc[d].ERs*1e6, df_MuTau.iloc[d].spec, label = f'$g_x$ =  {df_MuTau.iloc[d].g_x:.4f}, $M_A$ = {df_MuTau.iloc[d].m_A:.4f}', color ='C3' )
ax1.loglog(df_MuTau.iloc[e].ERs*1e6, df_MuTau.iloc[e].spec, label = f'$g_x$ =  {df_MuTau.iloc[e].g_x:.4f}, $M_A$ = {df_MuTau.iloc[e].m_A:.4f}', color ='C4' )
ax1.loglog(df_MuTau.iloc[g].ERs*1e6, df_MuTau.iloc[g].spec, label = f'$g_x$ =  {df_MuTau.iloc[g].g_x:.4f}, $M_A$ = {df_MuTau.iloc[g].m_A:.4f}', color ='C5' )

ax1.set_xlim(1e-2,1e2)
ax1.legend(fontsize=12,frameon=False,loc='lower left')

# Color graph of all points
cmap = plt.get_cmap('Spectral')
ax2.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['diff'].values, vmax = 2, cmap=cmap)
ax2.set_ylabel('$g_x$ [GeV]', size=16)
ax2.set_xlabel('$M_A$ [GeV]', size=16)
ax2.set_xscale('log')
ax2.set_yscale('log')

fmt = matplotlib.ticker.FormatStrFormatter('%1.1f')
sm =  ScalarMappable(norm=matplotlib.colors.Normalize(vmin = 1, vmax = 2), cmap=cmap)
cbar = fig.colorbar(sm, ax=ax2, format=fmt)
cbar.set_label(r'$\frac{(BSM + SM )}{SM}$', size=20)




# plotting points on colour graph
ax2.plot(df_MuTau.iloc[a].m_A, df_MuTau.iloc[a].g_x, marker='o', markersize=10, color ='C0' )
ax2.plot(df_MuTau.iloc[b].m_A, df_MuTau.iloc[b].g_x, marker='o', markersize=10, color ='C1' )
ax2.plot(df_MuTau.iloc[c].m_A, df_MuTau.iloc[c].g_x, marker='o', markersize=10, color ='C2' )
ax2.plot(df_MuTau.iloc[d].m_A, df_MuTau.iloc[d].g_x, marker='o', markersize=10, color ='C3' )
ax2.plot(df_MuTau.iloc[e].m_A, df_MuTau.iloc[e].g_x, marker='o', markersize=10, color ='C4' )
ax2.plot(df_MuTau.iloc[g].m_A, df_MuTau.iloc[g].g_x, marker='o', markersize=10, color ='C5' )

plt.show()

#%% P2 RATIO PLOT OF THE 6 POINTS
    
f=plt.figure(50,figsize=(10,8))

# different lengths sampled at different points

# line at 1
plt.axhline(y=1, color='r', linestyle='--')

#plt.ylim(0, 10)

plt.plot(sm_er*1e6, getratio(a) , label = '$g_x=$' + str(df_MuTau.iloc[a].g_x), color ='C0')
plt.plot(sm_er*1e6, getratio(b) , label = '$g_x=$' + str(df_MuTau.iloc[b].g_x), color ='C1')
plt.plot(sm_er*1e6, getratio(c) , label = '$g_x=$' + str(df_MuTau.iloc[c].g_x), color ='C2')
plt.plot(sm_er*1e6, getratio(d) , label = '$g_x=$' + str(df_MuTau.iloc[d].g_x), color ='C3')
plt.plot(sm_er*1e6, getratio(e) , label = '$g_x=$' + str(df_MuTau.iloc[e].g_x), color ='C4')
plt.plot(sm_er*1e6, getratio(g) , label = '$g_x=$' + str(df_MuTau.iloc[g].g_x), color ='C5')

plt.xlabel(r'$E_R\,\,\,\left[\rm{keV}\right]$', size=16)
plt.ylabel(r'$\frac{(BSM + SM )}{SM}$', size=20)

plt.xscale('log')
plt.yscale('log')

plt.legend(fontsize=12,frameon=False,loc='upper right')


#%% INVESTIGATING THE TAIL

# when using np.trapz and comparing to the SM spec integration no ratio is below 1, meaning that there doesn't seem to be any destructive effects
# however this doesn't take into account the possibility of major shape difference - meaning that in one area the spec dips below the SM
# but in another the BSM additions are so great that it cancels this dip out in the integral over the entire spec.
# Many of the plots below show major increases in rate at the lower energy regions which could be masking changes in the tail.
# Below I wll investigate changes in only the tail region.



# Adding taildiff column to data frame
# using traps for the bsm specs then divding by the value for the sm - gives an indication of the difference between sma nd bsm
# find index above threshold for both and compare above a given energy


sm_index = bisect.bisect_left(sm_er, 6.25e-6)
index = bisect.bisect_left(df_MuTau.iloc[1].ERs, 6.25e-6)
df_MuTau['taildiff'] = df_MuTau.apply(lambda x: np.trapz(x.spec[index:], x.ERs[index:])/ np.trapz(sm_spec[sm_index:], sm_er[sm_index:]), axis=1)

taildiffs = df_MuTau['taildiff']

fig=plt.figure(13,figsize=(10,8))

# plotting the diff with gx and M
plt.scatter(df_MuTau['m_A'].values, df_MuTau['g_x'].values, c=df_MuTau['taildiff'].values, vmax=2, cmap = plt.get_cmap('Spectral'))
plt.ylabel('$g_x$ [GeV]', size=16)
plt.xlabel('$M_A$ [GeV]', size=16)
plt.title('Tail only: $E_R$ > 6.25keV')
plt.xscale('log')
plt.yscale('log')

fmt = matplotlib.ticker.FormatStrFormatter('%1.1f')
cbar = plt.colorbar(format = fmt)
cbar.set_label(r'$\frac{(BSM + SM )}{SM}$', size=20)
plt.show()




