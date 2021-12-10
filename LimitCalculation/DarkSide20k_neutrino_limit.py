#!/usr/bin/env python
# coding: utf-8

'''
DarkSide20k_neutrino_limit.py

Author: EB
Date: Oct 2021

**Binned profile likelihood ratio analysis for DS20k using a asymptotic gaussian likelihood**

Adapted from DarkSide20k_FullProjectedLimit. Now using neutrinos as signal not background, also we are
not using any DM signal.

'''


### imports ###
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from iminuit import Minuit

### defining variables ###
# use a 100tonneyear exposure
ds20k_exposure = 365*1000*100 # exposure in kgday
ds20k_exposure_tonneyear = 100 # exposure in tyr - 5yr run with fiducial mass of 20t

# the data is in number of electrons. use a 4e- threshold
firstbin = 0 #!!! usually 4
lastbin =  50 #!!! usually 50

# 15% uncertainty on backgrounds
per_err = 0.15

#!!! multiplying cevns signal
mult = 1


# =============================================================================
#  importing data
# =============================================================================

# these are already scaled to the 100 tonneyear exposure
ar38bkgrd = np.loadtxt('Data_files/ds20k-39Ar_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]
ar38bkgrd_p = np.loadtxt('Data_files/ds20k-39Ar_bkgrd-p.dat',delimiter=' ')[firstbin:lastbin,1]
ar38bkgrd_m = np.loadtxt('Data_files/ds20k-39Ar_bkgrd-m.dat',delimiter=' ')[firstbin:lastbin,1]
ar38bkgrd_dplus = ar38bkgrd_p - ar38bkgrd
ar38bkgrd_dminus = ar38bkgrd_m - ar38bkgrd

gammabkgrd = np.loadtxt('Data_files/ds20k-gamma_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]
gammabkgrd_p = np.loadtxt('Data_files/ds20k-gamma_bkgrd-p.dat',delimiter=' ')[firstbin:lastbin,1]
gammabkgrd_m = np.loadtxt('Data_files/ds20k-gamma_bkgrd-m.dat',delimiter=' ')[firstbin:lastbin,1]
gammabkgrd_dplus = gammabkgrd_p - gammabkgrd
gammabkgrd_dminus = gammabkgrd_m - gammabkgrd


# neutrino signal
cennssig = np.loadtxt('Data_files/ds20k-cenns_bkgrd.dat',delimiter=' ')[firstbin:lastbin,1]*mult #!!! making cenns sig bigger
cennssig_p = np.loadtxt('Data_files/ds20k-cenns_bkgrd-p.dat',delimiter=' ')[firstbin:lastbin,1]*mult #!!!
cennssig_m = np.loadtxt('Data_files/ds20k-cenns_bkgrd-m.dat',delimiter=' ')[firstbin:lastbin,1]*mult #!!!
cennssig_dplus = cennssig_p - cennssig
cennssig_dminus = cennssig_m - cennssig

# get bins
bins = (np.loadtxt('Data_files/ds20k-39Ar_bkgrd.dat',delimiter=' '))[firstbin:lastbin,0]


# =============================================================================
# Defining functions
# =============================================================================


def construct_covariance (uncertainty_amplitudes, correlation_matrix) :
    '''
    Creates a covariance matrix using the uncertinties on the measurements and their correlations.
    If an int is used in place of a correlation_matrix, it is assumed that the uncertainties are either
    fully correlated (Corr = 1) or uncorrelated (Corr = 0)
    
    uncertainty_amplitudes: array of uncertainties
    correlation_matrix: if int, 1 (0) if uncertainties (un)correlated
                        else a full correlation matrix should be supplied
    
    Returns the covariance matrix
    '''
    num_measurements = len(uncertainty_amplitudes)
    # uncertaities uncorrelated - matrix of zeros with diagonal values of one (uncertainties are correlated only with themselves)
    if (type(correlation_matrix) is int) and (correlation_matrix == 0) :
        correlation_matrix = np.eye(num_measurements)
    # uncertainties correlated - matrix of ones (all uncertainties fully correlated)
    elif (type(correlation_matrix) is int) and (correlation_matrix == 1) :
        correlation_matrix = np.ones(shape=(num_measurements,num_measurements))
    # creating covariance matrix - Covij = Corrij x erri x errj
    covariance_matrix = np.zeros(shape=(num_measurements, num_measurements))    
    for i, ex in enumerate(uncertainty_amplitudes) :
        for j, ey in enumerate(uncertainty_amplitudes) :
            covariance_matrix[i, j] = correlation_matrix[i, j]*ex*ey
    return covariance_matrix


def get_prediction (c,energyscaleAr_NP) :
    '''
    Get the background and cenns predicton for a given c and energy scale
    energyscaleAr_NP: energy scale nuisance parameter
    '''
    ar = ar38bkgrd; gam = gammabkgrd; cen = cennssig
    # if energy scale +ve use delta plus
    if energyscaleAr_NP>0 : 
        ar = ar38bkgrd + energyscaleAr_NP * ar38bkgrd_dplus
        gam = gammabkgrd + energyscaleAr_NP * gammabkgrd_dplus
        cen = cennssig + energyscaleAr_NP * cennssig_dplus
    # if energy scale -ve use delta minus
    if energyscaleAr_NP<0 : 
        ar = ar38bkgrd + energyscaleAr_NP * ar38bkgrd_dminus
        gam = gammabkgrd + energyscaleAr_NP * gammabkgrd_dminus
        cen = cennssig + energyscaleAr_NP * cennssig_dminus
        
    bkgrd = ar + gam
    prediction = (c * cen) + bkgrd
    return prediction

def get_ds20k_TNLL (measurement, c, energyscaleAr_NP) :
    '''
    Returns -2logL for a given c and energy scale

    '''
    # calculate chi2
    # chi2 = (measured values - predicted values)^T x total covarince matrix^-1 x (measured values - predicted values)
    residual = measurement - get_prediction(c,energyscaleAr_NP)
    chi2 =  np.matmul(residual, np.matmul(total_cov_inverse, residual))
    # calculate likelihood
    # L = 1/(sqrt(2pi)^nummeasurements |Cov|) x exp(-0.5 chi2) x penalty term: gauss(NP)
    norm_factor = 1. / (np.power(np.sqrt(2*np.pi), len(measurement))*total_cov_det)
    NP_constraint = stats.norm.pdf(energyscaleAr_NP) # penalty term
    L = norm_factor * np.exp(-0.5*chi2) * NP_constraint
    # calculate -2LogL
    NLL = -2. * np.log(L)
    return NLL


def get_loglikelihood_from_fit_params (measurement, params) :
    '''
    Returns -2logL for given values of c and energy scale from params
    
    '''
    c, energyscaleAr_NP = params[0],params[1]
        
    return get_ds20k_TNLL (measurement, c, energyscaleAr_NP)


def get_argmin_TNLL (measurement, c=None, energyscaleAr_NP=None, minos=None, precision=None) :
    '''
    Uses minuit to get an optimised -2logL function
    if c or energyscaleAr_NP are given, they are fixed to that value
    if they aren't provided they float and we maximise the likelihood
    
    Returns optimised -2LogL and the minuit object
    '''
    fix_c, fix_energyscaleAr_NP = False, False
    if type(c) == type(None) : c = 0. # c not provided
    else                     : fix_c = True # c provided so is fixed
    if type(energyscaleAr_NP) == type(None) : energyscaleAr_NP = 0. # energyscale not provided
    else                                    : fix_energyscaleAr_NP = True # energyscale is provided so is fixed
    # defining lambda function which returns -2logL for given params
    fit_func = lambda params : get_loglikelihood_from_fit_params(measurement,params)
    # creating minuit object - parameters (given func: accepts single parameter which is a np array, 
    # starting point for minimisation, parameter names, inital step sizes, fixing values T/F, fcn increment 1SD above min)
    minuit_object = Minuit.from_array_func(fit_func, (c, energyscaleAr_NP), error = (1., 1.), 
                                           limit=None, fix = (fix_c, fix_energyscaleAr_NP),
                                           name = ("c", "energyscaleAr_NP"), errordef  = 1.)
    # minimisation via migrad algorithm
    minuit_object.migrad(precision=precision)
    # MINOS computes the asymetric confidence intervals
    if type(minos) != type(None) :
        minuit_object.minos(minos)
    
    return minuit_object.fmin.fval, minuit_object


def get_Wilks_CL_from_PLR (PLR, dof=1) :
    '''
    Get CL using Wilks' theorem
    The Profile Likelihood Ratio should be distributed like a chi2 distribution - interpolate PLR array to get CL values
    returns array of 1-CL for a given array of PLR values
    '''
    x = np.linspace(0, 20, 1001)
    y = stats.chi2.cdf(x, df=dof)
    return 1. - np.interp(PLR, x, y)

def get_limit_ds20k(alpha=0.1, plot=False) :
    '''
    Returns the upper limit for given alpha
     
    alpha: CL, default 0.1 for 90% limits
    plot: bool, make plots default false 
    '''
    try:
        # getting unconditional likelihood maximum - we find values of c and energyscale NP that minimise TNLL
        TNLL_best_obs, fit_obj = get_argmin_TNLL(measured_values,minos="c")
        c_errs = fit_obj.merrors['c']
        print("Best fit parameters:\n", fit_obj.params)
        
        if (fit_obj.values["c"] < 0) : 
            print('best fit was neg')
            TNLL_best_obs, fit_obj = get_argmin_TNLL(measured_values,c=0)

        c_best = fit_obj.values["c"]
        # multiply by 4 to ensure our range of c includes the 90% CL
        c_scan_values = np.linspace(c_best+4*c_errs.lower, c_best+4*c_errs.upper, 501)
        PLR_scan_tot  = np.array([get_argmin_TNLL(measured_values, c=c)[0] for c in c_scan_values]) - TNLL_best_obs
        CL_scan_tot   = get_Wilks_CL_from_PLR(PLR_scan_tot)
        
        x_lo, x_hi, y_lo, y_hi = [], [], [], []
        # iterating over the 'paired' values of c and CL - adding values to lo until last value in lo is larger than the new one
        for c, CL in zip(c_scan_values, CL_scan_tot) :
            if (len(y_lo) > 0) and (CL < y_lo[-1]) :
                x_hi.append(c)
                y_hi.append(CL)
            else:
                x_lo.append(c)
                y_lo.append(CL)
        # get value of c at CL = aplha for both arrays e.g when y = alpha what is x
        lower_limit = np.interp(alpha, y_lo, x_lo) 
        upper_limit = np.interp(alpha, y_hi[::-1], x_hi[::-1]) # [::-1] reverses array
        
        
        if plot==True:
            # plotting PLR vs C
            PLRfig=plt.figure(0)
            plt.plot(c_scan_values,PLR_scan_tot)
            plt.xlabel("$c$", fontsize=25, labelpad=35)
            plt.ylabel(r"$-2$ $\log$  $\frac{\mathcal{L}\left(c\right)}{\mathcal{L}_{max}}$", fontsize=25)
            plt.show()
            
            #plotting the limits
            CLfig=plt.figure(1)
            plt.axvline(lower_limit, ymax = 0.95-alpha, linestyle='dashed', color='grey')
            plt.axvline(upper_limit, ymax = 0.95-alpha, linestyle='dashed', color='grey')
            plt.axhline(1-alpha, linestyle='dashed', color='grey')
            plt.plot(c_scan_values,1-CL_scan_tot)
            plt.xlabel("c")
            plt.ylabel("1-CL")
            plt.show()
            
        
    except (RuntimeError):
        upper_limit=1e0

    return upper_limit



# =============================================================================
# SM signal, uncertaities and covariance matrix
# =============================================================================

# combining backgrounds - predicting we do not see cenns
prediction = (ar38bkgrd + gammabkgrd)

# as this is a projection, we put the 'observed' = to the predicted backgrounds
measured_values = prediction


# include a 15% correlated uncertainty on the background
ar39_uncertainty = per_err*ar38bkgrd
gamma_uncertainty = per_err*gammabkgrd

# statistical uncertainty on observed
data_stat_uncertainties = measured_values**0.5

# creating covariance matrices
ar39_covariance  = construct_covariance(ar39_uncertainty, 1) # 1 or 0 depending on if uncertainty is correlated or not
gamma_covariance = construct_covariance(gamma_uncertainty, 1)
data_stat_covariance = construct_covariance(data_stat_uncertainties, 0)

# now get the total covariance matrix and inverse/det for use in likelihood calc
total_cov  = ar39_covariance + gamma_covariance + data_stat_covariance
total_cov_inverse = np.linalg.inv(total_cov)
total_cov_det = np.linalg.det(total_cov)


# =============================================================================
# Calulating the limit
# =============================================================================

print("Upper limit: ", get_limit_ds20k(plot=True))

# =============================================================================
# Plotting spectrum for different values of c
# =============================================================================

f=plt.figure(2,figsize=(10,8))

plt.errorbar(bins, measured_values, yerr=data_stat_uncertainties, fmt='o', capsize=3, color='k',label='Background',linewidth=2)
plt.plot(bins, get_prediction(0, 0), linewidth=2,label='c = 0, EnergyScale_NP = 0 ')
plt.plot(bins, get_prediction(0.2, 0), linewidth=2,label='c = 0.2, EnergyScale_NP = 0 ')
plt.plot(bins, get_prediction(1, 0), linewidth=2, color='firebrick',label='c = 1, EnergyScale_NP = 0 ')

plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Number of events',fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)

plt.yscale('log')
#plt.ylim(1e4,5e6)

#plt.xlim(4,50)
plt.xlim(4,14)

plt.legend(fontsize=18,frameon=False,loc=2)
plt.tight_layout()

# =============================================================================
# Plotting background + signal histograms
# =============================================================================

fig=plt.figure(3,figsize=(10,8))
plt.errorbar(bins,prediction,yerr=(prediction)**0.5, fmt='o', capsize=3, color='k',label='Predicted observed',linewidth=2)
plt.bar(bins,gammabkgrd,width=1,color='darkorchid',log=True,alpha=0.7,label='SiPM gamma background') #/widths
plt.bar(bins,ar38bkgrd,width=1,color='deepskyblue',log=True,bottom=gammabkgrd,alpha=0.7,label='39Ar background')
plt.bar(bins,cennssig, width =1, color='firebrick',log=True,bottom=gammabkgrd+ar38bkgrd,alpha=0.7,label=r'CE$\nu$NS spectrum')
plt.plot(bins,0.2*cennssig, linewidth=3, color='firebrick',label=r'CE$\nu$NS spectrum c=0.2')

plt.axvline(x=4,linestyle='dashed',color='k',linewidth=2,label='4e$^-$ analysis threshold')


plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel('Number of events',fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)

plt.yscale('log')
plt.ylim(3e2,5e7)

plt.xlim(firstbin, lastbin)

plt.legend(fontsize=18,frameon=False,loc=9)
plt.tight_layout()


# =============================================================================
# Plotting S/S+B
# =============================================================================

fig=plt.figure(4,figsize=(10,8))
plt.errorbar(bins,prediction/prediction, fmt='o', capsize=3, color='k',label='Background only',linewidth=2)
plt.bar(bins,(cennssig+prediction)/prediction, width =1, color='firebrick',label=r"CE$\nu$NS", alpha=0.8)

plt.axhline(y=1,linestyle='dashed',color='black',linewidth=2)

plt.xlabel('Number of electrons',fontsize=26)
plt.ylabel(r"$\frac{S+B}{B}$",fontsize=26) # per kg day
plt.tick_params(axis='both', which='major', labelsize=20)

plt.xlim(firstbin, lastbin)
plt.ylim(0.95,1.05)

plt.legend(fontsize=18,frameon=False,loc=1)
plt.tight_layout()


