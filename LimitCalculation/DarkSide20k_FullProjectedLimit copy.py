#!/usr/bin/env python
# coding: utf-8

# **Binned profile likelihood ratio analysis for DS20k using a asymptotic gaussian likelihood**

### imports ###
import numpy as np
import uproot
from scipy import stats
from matplotlib import pyplot as plt
from iminuit import Minuit
import matplotlib
from matplotlib import rc

### defining variables ###
# use a 100tonneyear exposure
ds20k_exposure = 365*1000*100 # exposure in kgday
ds20k_exposure_tonneyear = 100 # exposure in tyr - 5yr run with fiducial mass of 20t

# the data is in number of electrons. use a 4e- threshold
ds20k_firstbin = 4
ds20k_lastbin = 50



# =============================================================================
#  importing the darkside backgrounds
# =============================================================================

# these are already scaled to the 100 tonneyear exposure
ds20k_ar39bkgrd = np.loadtxt('ds20k-39Ar_bkgrd.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1] # getting 2nd col of binned data
ds20k_ar39bkgrd_p = np.loadtxt('ds20k-39Ar_bkgrd-p.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_ar39bkgrd_m = np.loadtxt('ds20k-39Ar_bkgrd-m.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_ar39bkgrd_dplus = ds20k_ar39bkgrd_p - ds20k_ar39bkgrd
ds20k_ar39bkgrd_dminus = ds20k_ar39bkgrd_m - ds20k_ar39bkgrd

ds20k_gammabkgrd = np.loadtxt('ds20k-gamma_bkgrd.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_gammabkgrd_p = np.loadtxt('ds20k-gamma_bkgrd-p.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_gammabkgrd_m = np.loadtxt('ds20k-gamma_bkgrd-m.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_gammabkgrd_dplus = ds20k_gammabkgrd_p - ds20k_gammabkgrd
ds20k_gammabkgrd_dminus = ds20k_gammabkgrd_m - ds20k_gammabkgrd

ds20k_cennsbkgrd = np.loadtxt('ds20k-cenns_bkgrd.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_cennsbkgrd_p = np.loadtxt('ds20k-cenns_bkgrd-p.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_cennsbkgrd_m = np.loadtxt('ds20k-cenns_bkgrd-m.dat',delimiter=' ')[ds20k_firstbin:ds20k_lastbin,1]
ds20k_cennsbkgrd_dplus = ds20k_cennsbkgrd_p - ds20k_cennsbkgrd
ds20k_cennsbkgrd_dminus = ds20k_cennsbkgrd_m - ds20k_cennsbkgrd


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


def get_dm_spectra(op,specnum, exposure, firstbin, lastbin): # op = operator - I'm only interested in 1 and specnum - integer relating to dm mass
    '''
    Loads in dark matter spectra from root files and corrects the binning.
    op: int, EFT operator
    specnum: int, refers to a specfic dm mass
    exposure: experiment exposure in kgday
    firstbin: int, value of the first bin we want to include
    lastbin: int, value of the last bin we want to include
    
    Returns the DM spectra and the +/- energy shifted spectra difference
    '''

    filename = './DS50/Ar_O%s/spectra_ds20k_Ar_O%s_%i.root'%(op,op,specnum)
        
    dmdata = uproot.open(filename)
    hist = dmdata["hSum%i"%specnum]
    dmevents_perc = np.array(hist.values(),dtype='float64')*exposure

    # also have the spectra that correspond to an energy shift of +/-1 sigma
    hist_p = dmdata["hSum_p%i"%specnum]
    dmevents_perc_p = np.array(hist_p.values(),dtype='float64')*exposure

    hist_m = dmdata["hSum_m%i"%specnum]
    dmevents_perc_m = np.array(hist_m.values(),dtype='float64')*exposure

    BSM_prediction_per_c = [] #dmevents_perc
    BSM_prediction_per_c_p = [] #dmevents_perc
    BSM_prediction_per_c_m = [] #dmevents_perc
    
    # the data in these files in binned in half e-, so we need to combine adjacent bins to get 1e- binning
    for i in range(int(len(dmevents_perc)/2)) :
        BSM_prediction_per_c.append((dmevents_perc[i*2]+dmevents_perc[i*2+1])) # joint adjacent bins as they are 0.5 e wide ( to make them 1e- wide)
        BSM_prediction_per_c_p.append((dmevents_perc_p[i*2]+dmevents_perc_p[i*2+1])) 
        BSM_prediction_per_c_m.append((dmevents_perc_m[i*2]+dmevents_perc_m[i*2+1])) 

    # restriciting data to the bins we are interested in
    BSM_prediction_per_c=np.array(BSM_prediction_per_c)[firstbin:lastbin]
    BSM_prediction_per_c_p=np.array(BSM_prediction_per_c_p)[firstbin:lastbin]
    BSM_prediction_per_c_m=np.array(BSM_prediction_per_c_m)[firstbin:lastbin]
    
    # for use in likelihood NP 
    dplus = BSM_prediction_per_c_p-BSM_prediction_per_c
    dminus = BSM_prediction_per_c_m-BSM_prediction_per_c
    
    return BSM_prediction_per_c,dplus,dminus


def get_ds20k_dm_prediction (dm_array,energyscaleAr_NP) :
    '''
    Getting the DM spectrum, shifted for a specific energy scale

    dm_array: array of dm spectrum and the +/- energy shifted spectrum difference
    energyscaleAr_NP: energy scale nuisance parameter
    
    returns the shifted dm spectrum
    '''
    n, dplus, dminus = dm_array
    
    dm_shifted = n
    
    if energyscaleAr_NP>0 : dm_shifted = n + energyscaleAr_NP * dplus
    if energyscaleAr_NP<0 : dm_shifted = n + energyscaleAr_NP * dminus

    return dm_shifted


def get_ds20k_sm_prediction (energyscaleAr_NP) :
    '''
    Getting SM spectrum
    
    energyscaleAr_NP: energy scale nuisance parameter
    
    returns shifted SM spectrum & individual shifted backgrounds
    '''
    ar = ds20k_ar39bkgrd; gam = ds20k_gammabkgrd; cen = ds20k_cennsbkgrd
    if energyscaleAr_NP>0 : 
        ar = ds20k_ar39bkgrd + energyscaleAr_NP * ds20k_ar39bkgrd_dplus
        gam = ds20k_gammabkgrd + energyscaleAr_NP * ds20k_gammabkgrd_dplus
        cen = ds20k_cennsbkgrd + energyscaleAr_NP * ds20k_cennsbkgrd_dplus
        
    if energyscaleAr_NP<0 : 
        ar = ds20k_ar39bkgrd + energyscaleAr_NP * ds20k_ar39bkgrd_dminus
        gam = ds20k_gammabkgrd + energyscaleAr_NP * ds20k_gammabkgrd_dminus
        cen = ds20k_cennsbkgrd + energyscaleAr_NP * ds20k_cennsbkgrd_dminus
        
    sm_shifted = ar + gam + cen
    return sm_shifted, ar, gam, cen


def get_ds20k_prediction (c,dm_array,energyscaleAr_NP) :
    '''
    Get the SM + DM predicton for a given c and energy scale
    energyscaleAr_NP: energy scale nuisance parameter
    '''
    return c*get_ds20k_dm_prediction(dm_array,energyscaleAr_NP) + get_ds20k_sm_prediction(energyscaleAr_NP)[0]


def get_ds20k_chi2 (measurement, c, dm_array,energyscaleAr_NP,cov_inverse) :
    '''
    Get the chi2 for given c and energy scale
    '''
    # chi2 = (measured values - predicted values)^T x total covarince matrix^-1 x (measured values - predicted values)
    residual = measurement - get_ds20k_prediction(c,dm_array,energyscaleAr_NP)
    return np.matmul(residual, np.matmul(cov_inverse, residual))


def get_ds20k_likelihood (measurement, c, dm_array,energyscaleAr_NP) :
    '''
    Get the likelihood for given c and energy scale

    '''
    # L = 1/(sqrt(2pi)^nummeasurements |Cov|) x exp(-0.5 chi2) x penalty term: gauss(NP)
    norm_factor = 1. / np.power(np.sqrt(2*np.pi), len(measurement))*ds20k_data_plus_SM_cov_det
    NP_constraint = stats.norm.pdf(energyscaleAr_NP) # penalty term
    return norm_factor * np.exp(-0.5*get_ds20k_chi2(measurement, c, dm_array,energyscaleAr_NP,ds20k_data_plus_SM_cov_inverse)) * NP_constraint


def get_ds20k_TNLL (measurement, c, dm_array,energyscaleAr_NP) :
    '''
    Returns -2logL for a given c and energy scale

    '''
    return -2. * np.log(get_ds20k_likelihood (measurement, c, dm_array,energyscaleAr_NP))


def get_loglikelihood_from_fit_params (ds20k_measurement, dm_array, params) :
    '''
    Returns -2logL for given values of c and enrgy scale from params
    
    '''
    c, energyscaleAr_NP = params[0],params[1]
        
    return get_ds20k_TNLL (ds20k_measurement, c, dm_array,energyscaleAr_NP)


def get_argmin_TNLL (ds20k_measurement, dm_array, c=None, energyscaleAr_NP=None, minos=None, precision=None) :
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
    fit_func = lambda params : get_loglikelihood_from_fit_params(ds20k_measurement, dm_array ,params)
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


def get_limit_ds20k(operator,specnum,cross_section) :
    '''
    Returns the upper limit for a specfific DM spectrum i.e. a specific mass

    '''
    # get dm spectum
    ds20k_bsm = get_dm_spectra(operator,specnum,ds20k_exposure,ds20k_firstbin,ds20k_lastbin)
    
    try:
        # 
        TNLL_best_obs, fit_obj = get_argmin_TNLL(ds20k_measured_values,ds20k_bsm,minos="c")
        c_errs = fit_obj.get_merrors()['c']
        
        if (fit_obj.values["c"] < 0) : 
            print('best fit was neg')
            TNLL_best_obs, fit_obj = get_argmin_TNLL(ds20k_measured_values,ds20k_bsm,c=0)

        c_best = fit_obj.values["c"]
        c_scan_values = np.linspace(c_best+4*c_errs.lower, c_best+4*c_errs.upper, 501)
        PLR_scan_tot  = np.array([get_argmin_TNLL(ds20k_measured_values,ds20k_bsm, c=c)[0] for c in c_scan_values]) - TNLL_best_obs
        CL_scan_tot   = get_Wilks_CL_from_PLR(PLR_scan_tot)
        lower_limit, upper_limit = get_limits(c_scan_values, CL_scan_tot, 0.1)

        up=upper_limit*cross_section
    except (RuntimeError):
        up=1e0

    return up




def get_full_ds20k_limit(operator,cross_section) :
    '''
    Calculates the upper limit for an array of DM masses.
    Saves the values in a file

    '''
    limit = []
    # looping over all masses and getting the upper limit for each
    for i in range(0,34) :
        limit.append(get_limit_ds20k(operator, i,cross_section))
        
    if operator=='1':
        masses = np.concatenate([np.logspace(-1.2,-1.05,4) , np.logspace(-1,1,30)])
    else : 
        masses = np.logspace(-1.2,1,34) 
    # saving masses and limit to datafile
    np.savetxt('ds20k_O%slimit_energyscale.dat'%operator, np.array([masses, limit]).T, delimiter=', ')
    return np.array(limit)
        




# =============================================================================
# SM signal, uncertaities and covariance matrix
# =============================================================================

# combining backgrounds for SM signal
ds20k_SM_prediction = (ds20k_ar39bkgrd + ds20k_gammabkgrd + ds20k_cennsbkgrd)

# as this is a projection, we put the 'observed' = to the predicted backgrounds
# eg. we assume no dark matter is observed, we just observed the predicted SM rate
ds20k_measured_values = ds20k_SM_prediction


# include a 15% correlated uncertainty on each background component
ds20k_ar39_uncertainty = 0.15*ds20k_ar39bkgrd
ds20k_gamma_uncertainty = 0.15*ds20k_gammabkgrd
ds20k_cenns_uncertainty = 0.15*ds20k_cennsbkgrd

# statistical uncertainty on observed
ds20k_data_stat_uncertainties = ds20k_measured_values**0.5

# creating covariance matrices
ds20k_ar39_covariance  = construct_covariance(ds20k_ar39_uncertainty, 1) # 1 or 0 depending on if uncertainty is correlated or not
ds20k_gamma_covariance = construct_covariance(ds20k_gamma_uncertainty, 1)
ds20k_cenns_covariance  = construct_covariance(ds20k_cenns_uncertainty, 1)
ds20k_data_stat_covariance = construct_covariance(ds20k_data_stat_uncertainties, 0)

# now get the total covariance matrix and inverse/det for use in likelihood calc
ds20k_data_plus_SM_covariance  = ds20k_ar39_covariance + ds20k_gamma_covariance + ds20k_cenns_covariance + ds20k_data_stat_covariance
ds20k_data_plus_SM_cov_inverse = np.linalg.inv(ds20k_data_plus_SM_covariance)
ds20k_data_plus_SM_cov_det    = np.linalg.det(ds20k_data_plus_SM_covariance)


# =============================================================================
# Calulating the limits
# =============================================================================

# cross section for this operator
o1_cross_section = 1e-40 # data in root files produced at this value

# total limit for all masses

get_full_ds20k_limit('1',o1_cross_section)



# =============================================================================
# Plotting
# =============================================================================


fig = plt.figure(figsize=(9,5))
ax = plt.subplot(111)
rc('text', usetex=False)

data = np.loadtxt('ds20k_O1limit_energyscale.dat',delimiter=',')
ax.loglog(data[:,0],data[:,1],'--',label='DS20k $\mathcal{O}_{1}$ projection',linewidth=2.5,color='tab:red')



ax.set_ylim(5e-47,5e-36)
ax.set_xlim(0.06,8)
ax.set_ylabel('Dark Matter-Nucleon $\sigma_{SI}$ [cm$^{2}$]',fontsize=18)
ax.set_xlabel('M$_{\chi}$ [GeV/c$^2$]',fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=15)

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.xaxis.set_ticks([0.1,0.5,1,3,5])

ax.legend(fontsize=15,frameon=False) #loc='center left', bbox_to_anchor=(1, 0.5)

plt.tight_layout()
fig.savefig('Limits-LZ-DS20k-proj-o1.png',transparent=False,facecolor='w')
fig.savefig('Limits-LZ-DS20k-proj-o1.pdf',transparent=False,facecolor='w')




