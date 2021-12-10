#!/usr/bin/env python
# coding: utf-8



import numpy as np
from   matplotlib import pyplot as plt


ds20k_exposure = 365*1000*100
ds20k_exposure_tonneyear = 100
xsection_data = 1e-40
ds20k_firstbin = 4
ds20k_lastbin = 50

# =============================================================================
# Retrieve data from files (already adjusted for ds20k)
# =============================================================================

ar39bkgrd = np.loadtxt('ds20k-39Ar_bkgrd.dat',delimiter=' ')[:,1]

gammabkgrd = np.loadtxt('ds20k-gamma_bkgrd.dat',delimiter=' ')[:,1]

cennssig = np.loadtxt('ds20k-cenns_bkgrd.dat',delimiter=' ')[:,1]

# getting bins
bins = (np.loadtxt('ds20k-39Ar_bkgrd.dat',delimiter=' '))[:,0]

# combining backgrounds
totalbkgrd = ar39bkgrd + gammabkgrd

# assuming that we observe only background
ds20k_observed = totalbkgrd

# =============================================================================
# Plotting data
# =============================================================================

f=plt.figure(figsize=(10,8))

plt.errorbar(bins,ds20k_observed,yerr=(ds20k_observed)**0.5, fmt='o', capsize=3, color='k',label='Predicted observed',linewidth=2)
plt.bar(bins,gammabkgrd,width=1,color='darkorchid',log=True,alpha=0.7,label='SiPM gamma background') #/widths
plt.bar(bins,ar39bkgrd,width=1,color='deepskyblue',log=True,bottom=gammabkgrd,alpha=0.7,label='39Ar background')
plt.plot(bins,cennssig, linewidth=3, color='firebrick',label=r'CE$\nu$NS spectrum')
plt.plot(bins,0.2*cennssig, linewidth=3, color='firebrick',label=r'CE$\nu$NS spectrum c=0.2')

plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label='4e$^-$ analysis threshold')


plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Number of events',fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)

plt.yscale('log')
plt.ylim(3e2,5e7)

plt.xlim(0,50)

plt.legend(fontsize=18,frameon=False,loc=9)
plt.tight_layout()

