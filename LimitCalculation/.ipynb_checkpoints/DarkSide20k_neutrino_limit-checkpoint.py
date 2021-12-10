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
ds20k_firstbin = 4
ds20k_lastbin = 50

# 15 uncertainty on backgrounds
per_err = 0.15


# =============================================================================
#  importing data
# =============================================================================

# these are already scaled to the 100 tonneyear exposure
ds20k_ar39bkgrd = np.loadtxt('ds20k-39Ar_bkgrd.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_ar39bkgrd_p = np.loadtxt('ds20k-39Ar_bkgrd-p.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_ar39bkgrd_m = np.loadtxt('ds20k-39Ar_bkgrd-m.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_ar39bkgrd_dplus = ds20k_ar39bkgrd_p - ds20k_ar39bkgrd
ds20k_ar39bkgrd_dminus = ds20k_ar39bkgrd_m - ds20k_ar39bkgrd

ds20k_gammabkgrd = np.loadtxt('ds20k-gamma_bkgrd.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_gammabkgrd_p = np.loadtxt('ds20k-gamma_bkgrd-p.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_gammabkgrd_m = np.loadtxt('ds20k-gamma_bkgrd-m.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_gammabkgrd_dplus = ds20k_gammabkgrd_p - ds20k_gammabkgrd
ds20k_gammabkgrd_dminus = ds20k_gammabkgrd_m - ds20k_gammabkgrd

# neutrino signal
ds20k_cennssig = np.loadtxt('ds20k-cenns_bkgrd.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_cennssig_p = np.loadtxt('ds20k-cenns_bkgrd-p.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_cennssig_m = np.loadtxt('ds20k-cenns_bkgrd-m.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_cennssig_dplus = ds20k_cennssig_p - ds20k_cennssig
ds20k_cennssig_dminus = ds20k_cennssig_m - ds20k_cennssig


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


def get_ds20k_prediction (c,energyscaleAr_NP) :
    '''
    Get the background and cenns predicton for a given c and energy scale
    energyscaleAr_NP: energy scale nuisance parameter
    '''
    ar = ds20k_ar39bkgrd; gam = ds20k_gammabkgrd; cen = ds20k_cennssig
    # if energy scale +ve use delta plus
    if energyscaleAr_NP>0 : 
        ar = ds20k_ar39bkgrd + energyscaleAr_NP * ds20k_ar39bkgrd_dplus
        gam = ds20k_gammabkgrd + energyscaleAr_NP * ds20k_gammabkgrd_dplus
        cen = ds20k_cennssig + energyscaleAr_NP * ds20k_cennssig_dplus
    # if energy scale -ve use delta minus
    if energyscaleAr_NP<0 : 
        ar = ds20k_ar39bkgrd + energyscaleAr_NP * ds20k_ar39bkgrd_dminus
        gam = ds20k_gammabkgrd + energyscaleAr_NP * ds20k_gammabkgrd_dminus
        cen = ds20k_cennssig + energyscaleAr_NP * ds20k_cennssig_dminus
        
    bkgrd_shifted = ar + gam
    
    return c * cen + bkgrd_shifted


def get_ds20k_chi2 (measurement, c,energyscaleAr_NP,cov_inverse) :
    '''
    Get the chi2 for given c and energy scale
    '''
    # chi2 = (measured values - predicted values)^T x total covarince matrix^-1 x (measured values - predicted values)
    residual = measurement - get_ds20k_prediction(c,energyscaleAr_NP)
    return np.matmul(residual, np.matmul(cov_inverse, residual))


def get_ds20k_likelihood (measurement, c, energyscaleAr_NP) :
    '''
    Get the likelihood for given c and energy scale

    '''
    # L = 1/(sqrt(2pi)^nummeasurements |Cov|) x exp(-0.5 chi2) x penalty term: gauss(NP)
    norm_factor = 1. / np.power(np.sqrt(2*np.pi), len(measurement))*ds20k_data_plus_SM_cov_det
    NP_constraint = stats.norm.pdf(energyscaleAr_NP) # penalty term
    return norm_factor * np.exp(-0.5*get_ds20k_chi2(measurement, c, energyscaleAr_NP,ds20k_data_plus_SM_cov_inverse)) * NP_constraint


def get_ds20k_TNLL (measurement, c, energyscaleAr_NP) :
    '''
    Returns -2logL for a given c and energy scale

    '''
    return -2. * np.log(get_ds20k_likelihood (measurement, c, energyscaleAr_NP))


def get_loglikelihood_from_fit_params (ds20k_measurement, params) :
    '''
    Returns -2logL for given values of c and enrgy scale from params
    
    '''
    c, energyscaleAr_NP = params[0],params[1]
        
    return get_ds20k_TNLL (ds20k_measurement, c, energyscaleAr_NP)


def get_argmin_TNLL (ds20k_measurement, c=None, energyscaleAr_NP=None, minos=None, precision=None) :
    '''
    Uses minuit to get an optimised -2logL function
    if c or energyscaleAr_NP are given, they are fixed to that value
    if they aren't provided they float and we maximise the likelihood
    
    Returns fval: the value of the cost func at the minimum and the optimised minuit object
    '''
    fix_c, fix_energyscaleAr_NP = False, False
    if type(c) == type(None) : c = 0. # c not provided
    else : fix_c = True # c provided so is fixed
    if type(energyscaleAr_NP) == type(None) : energyscaleAr_NP = 0. # energyscale not provided
    else : fix_energyscaleAr_NP = True # energyscale is provided so is fixed
    # defining lambda function which returns -2logL for given params
    fit_func = lambda params : get_loglikelihood_from_fit_params(ds20k_measurement,params)
    # creating minuit object - parameters (given func: accepts single parameter which is a np array, 
    # starting point for minimisation, parameter names, inital step sizes, fixing values T/F, fcn increment 1SD above min)
    minuit_object = Minuit.from_array_func(fit_func, (c, energyscaleAr_NP), name = ("c", "energyscaleAr_NP"), 
                                           error = (1., 1.), fix = (fix_c, fix_energyscaleAr_NP), errordef  = 1.)
    # minimisation via migrad algorithm
    minuit_object.migrad(precision=precision)
    # MINOS computes the asymetric confidence intervals
    if type(minos) != type(None) :
        minuit_object.minos(minos)
    return minuit_object.get_fmin().fval, minuit_object


def get_Wilks_CL_from_PLR (PLR, dof=1) :
    '''
    Get CL using Wilks' theorem
    The Profile Likelihood Ratio should be distributed like a chi2 distribution
    returns array of CL for a given array of PLR values
    '''
    x = np.linspace(0, 20, 1001)
    y = stats.chi2.cdf(x, df=dof)
    return 1. - np.interp(PLR, x, y)


def get_limits (c_scan_values, CL_scan_tot, alpha, plot = False) :
    '''
    Get limits at a CL of alpha from CL of a PLR scan over values of c
    
    c_scan_values: array of values of c
    CL_scan_tot: array of CL from applying Wilks'theorem to a given PLR scan over the c_scan_values
    alpha: CL
    plot: bool, plot the limits default false 

    Returns the lower and upper limits
    '''
    x_lo, x_hi, y_lo, y_hi = [], [], [], []
    # iterating over the 'paired' values of c and CL - adding values to lo until last value in lo is larger than the new one
    for c, CL in zip(c_scan_values, CL_scan_tot) :
        if (len(y_lo) > 0) and (CL < y_lo[-1]) :
            x_hi.append(c)
            y_hi.append(CL)
        else:
            x_lo.append(c)
            y_lo.append(CL)
    # get value of c at CL = aplha
    lower_limit = np.interp(alpha, y_lo, x_lo) 
    upper_limit = np.interp(alpha, y_hi[::-1], x_hi[::-1])
    
    if plot == True:
        #plotting the limits
        l=np.empty(10); l.fill(lower_limit)
        h=np.empty(10); h.fill(upper_limit)
        concl =np.empty(len(c_scan_values)); concl.fill(1-alpha)
        clarr = np.linspace(0,.9,10)
        plt.plot(l,clarr,'--', color='grey')
        plt.plot(h,clarr,'--', color='grey')
        plt.plot(c_scan_values,concl,'--', color='grey')
        plt.plot(c_scan_values,1-CL_scan_tot)
        plt.xlabel("c")
        plt.ylabel("1-CL")
        plt.show()
    
    return lower_limit, upper_limit


def get_limit_ds20k() :
    '''
    Returns the upper limit for a cenns signal

    '''
    try:
        # 
        TNLL_best_obs, fit_obj = get_argmin_TNLL(ds20k_measured_values,minos="c")
        c_errs = fit_obj.get_merrors()['c']
        
        if (fit_obj.values["c"] < 0) : 
            print('best fit was neg')
            TNLL_best_obs, fit_obj = get_argmin_TNLL(ds20k_measured_values,c=0)

        c_best = fit_obj.values["c"]
        c_scan_values = np.linspace(c_best+4*c_errs.lower, c_best+4*c_errs.upper, 501)
        PLR_scan_tot  = np.array([get_argmin_TNLL(ds20k_measured_values, c=c)[0] for c in c_scan_values]) - TNLL_best_obs
        CL_scan_tot   = get_Wilks_CL_from_PLR(PLR_scan_tot)
        lower_limit, upper_limit = get_limits(c_scan_values, CL_scan_tot, 0.1, True)

        up=upper_limit
    except (RuntimeError):
        up=1e0

    return up



# =============================================================================
# SM signal, uncertaities and covariance matrix
# =============================================================================

# combining backgrounds - predicting we do not see cenns
ds20k_prediction = (ds20k_ar39bkgrd + ds20k_gammabkgrd)

# as this is a projection, we put the 'observed' = to the predicted backgrounds
ds20k_measured_values = ds20k_prediction


# include a 15% correlated uncertainty
ds20k_ar39_uncertainty = per_err*ds20k_ar39bkgrd
ds20k_gamma_uncertainty = per_err*ds20k_gammabkgrd

# statistical uncertainty on observed
ds20k_data_stat_uncertainties = ds20k_measured_values**0.5

# creating covariance matrices
ds20k_ar39_covariance  = construct_covariance(ds20k_ar39_uncertainty, 1) # 1 or 0 depending on if uncertainty is correlated or not
ds20k_gamma_covariance = construct_covariance(ds20k_gamma_uncertainty, 1)
ds20k_data_stat_covariance = construct_covariance(ds20k_data_stat_uncertainties, 0)

# now get the total covariance matrix and inverse/det for use in likelihood calc
ds20k_bkgrd_tot_cov  = ds20k_ar39_covariance + ds20k_gamma_covariance + ds20k_data_stat_covariance
ds20k_data_plus_SM_cov_inverse = np.linalg.inv(ds20k_bkgrd_tot_cov)
ds20k_data_plus_SM_cov_det    = np.linalg.det(ds20k_bkgrd_tot_cov)


# =============================================================================
# Calulating the limit
# =============================================================================

print("Upper limit: ", get_limit_ds20k())


