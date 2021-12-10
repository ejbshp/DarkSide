#!/usr/bin/env python
# coding: utf-8



import numpy as np
import uproot
from   matplotlib import pyplot as plt


ds20k_exposure = 365*1000*100
ds20k_exposure_tonneyear = 100
xsection_data = 1e-40
ds20k_firstbin = 4
ds20k_lastbin = 50
expected=True
correlated=True


# get data from text files for backgrounds

ar39bkgrdfile = './DS50/input_spectra/ds20k-39Ar_85Kr-spec-binned.dat'
ar39bkgrd = (0.7/2.7)*ds20k_exposure*(np.loadtxt(ar39bkgrdfile, delimiter=" "))[:,1] # 0.7/2.7 is neglecting Kr  darkside 50 sims -> 20k
bins = (np.loadtxt(ar39bkgrdfile, delimiter=" "))[:,0]

gammabkgrdfile = './DS50/input_spectra/ds20k-gammaPDMs-binned.dat'
gammabkgrd = (2.51e-4)*(22.1)*ds20k_exposure*(np.loadtxt(gammabkgrdfile, delimiter=" "))[:,1] # 2.5 geometry factor, 22.1 more sensors but each sensor less radioactive darkside 50 sims -> 20k

# not in darkside 50 - from RH a simulation they did
cenns = './DS50/input_spectra/spectra_cenns.root'
cenns_data = uproot.open(cenns)
hist = cenns_data['hNR_per_tonne_yr']
cenns_events_pertonneyr = np.array(hist.values(),dtype='float64')*ds20k_exposure_tonneyear # chaning units
cennsbkgrd = []
for i in range(int(len(cenns_events_pertonneyr)/2)) :
    cennsbkgrd.append((cenns_events_pertonneyr[i*2]+cenns_events_pertonneyr[i*2+1])) # joint adjacent bins as they are 0.5 e wide ( to make them 1e- wide)
cennsbkgrd = np.array(cennsbkgrd)[0:50]

bins = bins+0.5

totalbkgrd = ar39bkgrd + gammabkgrd + cennsbkgrd

ds20k_observed = totalbkgrd




i=19 # 1.08 GeV
#filename = './DS50/spectra_database_with_DSX_response/spectra_ds20k_15Feb_Ar_O1_%i.root'%(i)
filename = './DS50/Ar_O1/spectra_ds20k_Ar_O1_%i.root'%(i)

dmdata = uproot.open(filename)
hist = dmdata["hSum%i"%i]
dmevents_perc = np.array(hist.values(),dtype='float64')*ds20k_exposure

BSM_prediction_per_c = []
#print('dm before',dmevents_perc)
for i in range(int(len(dmevents_perc)/2)) :
    BSM_prediction_per_c.append((dmevents_perc[i*2]+dmevents_perc[i*2+1])) # joint adjacent bins as they are 0.5 e wide ( to make them 1e- wide)

print('dm',BSM_prediction_per_c)

i=19
hist = dmdata["hNR%i"%i]
dmevents_perc2 = np.array(hist.values(),dtype='float64')*ds20k_exposure

BSM_prediction_per_c2 = []
#print('dm before',dmevents_perc2)
for i in range(int(len(dmevents_perc2)/2)) :
    BSM_prediction_per_c2.append((dmevents_perc2[i*2]+dmevents_perc2[i*2+1])) # joint adjacent bins as they are 0.5 e wide ( to make them 1e- wide)

print('dm no mig',BSM_prediction_per_c2)
BSM_prediction_per_c=np.array(BSM_prediction_per_c[0:50]) # threshold to 50
BSM_prediction_per_c2=np.array(BSM_prediction_per_c2[0:50]) # threshold to 50



#pick a specific dark matter mass
# i=4 is 0.1GeV, i=14 is 0.48GeV, i=19 is 1.08GeV, i=23 is 2GeV, i=29 is 5.3 GeV
f=plt.figure(figsize=(10,8))


plt.errorbar(bins,ds20k_observed,yerr=(ds20k_observed)**0.5, fmt='o', capsize=3, color='k',label='Predicted observed [DS20k]',linewidth=2)
plt.bar(bins,gammabkgrd,width=1,color='orangered',log=True,alpha=0.7,label='SiPM gamma background') #/widths
plt.bar(bins,cennsbkgrd,width=1,color='tab:purple',log=True,alpha=0.7,label='CEvNS background',bottom=gammabkgrd)

plt.bar(bins,ar39bkgrd,width=1,color='limegreen',log=True,bottom=cennsbkgrd+gammabkgrd,alpha=0.7,label='39Ar background')
#plt.plot(binmid,er+cevns+cathode,'x',color='tab:purple',markersize=8,label='SM only')

plt.plot(bins,BSM_prediction_per_c2*100,color='tab:red',label='DM Spectrum, $\mathcal{O}_1$ 1.1 GeV',linewidth=3)

plt.plot(bins,BSM_prediction_per_c*100,'--',color='tab:red',label='DM Spectrum with Migdal, $\mathcal{O}_1$ 1.1 GeV',linewidth=3)

plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label='4 e- analysis threshold')

#plt.plot(binmid,er+cevns+cathode+BSM_prediction*100,'x',color='tab:red',markersize=8,label='SM + DM (0.5GeV,1e-38cm2)')
plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Number of events',fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.xscale('log')
plt.yscale('log')
plt.ylim(3e2,5e8)
#plt.ylim(1e-6,10)
plt.xlim(0,50)
#plt.axvline(x=3,linestyle='dashed',color='k',linewidth=2)
plt.legend(fontsize=18,frameon=False,loc=9)
plt.tight_layout()
#f.savefig('ds20kspectra.png',facecolor='w',transparent=False)
#f.savefig('ds20kspectra.pdf',facecolor='w',transparent=False)


