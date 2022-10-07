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

ar39bkgrd = np.loadtxt('Data_files/ds20k-39Ar_bkgrd.dat',delimiter=' ')[0:50,1]

gammabkgrd = np.loadtxt('Data_files/ds20k-gamma_bkgrd.dat',delimiter=' ')[0:50,1]

#cennssig = np.loadtxt('Data_files/ds20k-cenns_bkgrd.dat',delimiter=' ')[:,1]
sm_file = '/Users/user/DarkSide/Neutrino spectrum/output_my_cevns/PE_argon_SM_A_with_err.txt'
sm_cenns = np.loadtxt(sm_file ,delimiter=' ')[0:50,1]*100 # multiply by exposure
rel_cenns_err = np.loadtxt(sm_file ,delimiter=' ')[0:50,2]

# getting bins
bins = (np.loadtxt('Data_files/ds20k-39Ar_bkgrd.dat',delimiter=' '))[:,0]

# combining backgrounds
totalbkgrd = ar39bkgrd + gammabkgrd

# assuming that we observe only background
ds20k_observed = totalbkgrd

# =============================================================================
# Plotting data
# =============================================================================

f=plt.figure(figsize=(10,8))

#plt.errorbar(bins,ds20k_observed,yerr=(ds20k_observed)**0.5, fmt='o', capsize=3, color='k',label='Predicted observed',linewidth=2)
plt.bar(bins,gammabkgrd,width=1,color='firebrick',log=True,alpha=0.7,label='SiPM gamma background') #/widths
plt.bar(bins,sm_cenns, width=1, color='darkorchid',log=True,bottom=gammabkgrd, alpha=0.7, label=r'CE$\nu$NS signal')
plt.bar(bins,ar39bkgrd,width=1,color='deepskyblue',log=True,bottom=gammabkgrd+sm_cenns,alpha=0.7,label=r'$^{39}Ar$ background')

plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label=r'4 $e^{-}$ analysis threshold')

# plotting smaller cenns signal on top
#plt.plot(bins,0.2*cennssig, linewidth=3, color='firebrick',label=r'CE$\nu$NS spectrum c=0.2')

# electron analysis threshold
#plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label='4e$^-$ analysis threshold')


plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Number of events',fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)

plt.yscale('log')
plt.ylim(3e2,1e7)

plt.xlim(0,50)

plt.legend(fontsize=18,frameon=False,loc=1)
plt.tight_layout()

plt.savefig('/Users/user/Documents/Plots/backgrounds_and_signal.eps', format='eps')
plt.savefig('/Users/user/Documents/Plots/backgrounds_and_signal.png', format='png')


#%% plotting cevns alone with error - divide by 100 for per tyr


# plotting
f=plt.figure(figsize=(10,8))

plt.errorbar(bins, sm_cenns/100,yerr=(sm_cenns/100)*rel_cenns_err, fmt='o', capsize=3, color='k',linewidth=2, label = r'$\frac{N_{Events}}{\sqrt{N_{Samples}}}$')
plt.bar(bins,sm_cenns/100,width=1, color='slateblue', log=True, alpha=0.5, label = r'SM CE$\nu$NS')

plt.xlim(0, 50)

plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel(r'Events $[tyr]^{-1}$',fontsize=26) 
plt.yscale('log')
plt.legend(fontsize=18,frameon=False,loc='upper right')

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


plt.savefig('/Users/user/Documents/Plots/cevns_signal_with_err.eps', format='eps')
plt.savefig('/Users/user/Documents/Plots/cevns_signal_with_err.png', format='png')

#%% Plotting 0.25 cevns


f=plt.figure(figsize=(10,8))

#plt.errorbar(bins,ds20k_observed,yerr=(ds20k_observed)**0.5, fmt='o', capsize=3, color='k',label='Predicted observed',linewidth=2)
plt.bar(bins,gammabkgrd,width=1,color='firebrick',log=True,alpha=0.7,label='SiPM gamma background') #/widths
plt.bar(bins,0.25*sm_cenns, width=1, color='darkorchid',log=True,bottom=gammabkgrd, alpha=0.7, label=r'CE$\nu$NS signal multiplied by 0.25')
plt.bar(bins,ar39bkgrd,width=1,color='deepskyblue',log=True,bottom=gammabkgrd+0.25*sm_cenns,alpha=0.7,label=r'$^{39}Ar$ background')

plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label=r'4 $e^{-}$ analysis threshold')

# plotting smaller cenns signal on top
#plt.plot(bins,0.2*cennssig, linewidth=3, color='firebrick',label=r'CE$\nu$NS spectrum c=0.2')

# electron analysis threshold
#plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label='4e$^-$ analysis threshold')


plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Number of events',fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)

plt.yscale('log')
plt.ylim(3e2,1e7)

plt.xlim(0,50)

plt.legend(fontsize=18,frameon=False,loc=1)
plt.tight_layout()


plt.savefig('/Users/user/Documents/Plots/0.25cevns_signal_with_err.eps', format='eps')
plt.savefig('/Users/user/Documents/Plots/0.25cevns_signal_with_err.png', format='png')

#%% Plotting 0.25 cevns ratio


f=plt.figure(figsize=(10,8))

total_sig = gammabkgrd + ar39bkgrd + sm_cenns
total_sig_25 = gammabkgrd + ar39bkgrd + sm_cenns*0.25

plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label=r'4 $e^{-}$ analysis threshold')

plt.plot(total_sig_25/total_sig, linewidth=3, label=r'$\frac{0.25CE \nu NS+bkgrd}{CE \nu NS+bkgrd}$')


# plotting smaller cenns signal on top
#plt.plot(bins,0.2*cennssig, linewidth=3, color='firebrick',label=r'CE$\nu$NS spectrum c=0.2')

# electron analysis threshold
#plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label='4e$^-$ analysis threshold')


plt.xlabel(r'Number of electrons',fontsize=26)
plt.ylabel(r'Ratio of total number of events',fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)

plt.ylim(0.96,1.01)
# plt.yscale('log')
plt.xlim(0,50)

plt.legend(fontsize=18,frameon=False,loc=1)
plt.tight_layout()


plt.savefig('/Users/user/Documents/Plots/0.25cevns_signal_ratio.eps', format='eps')
plt.savefig('/Users/user/Documents/Plots/0.25cevns_signal_ratio.png', format='png')


