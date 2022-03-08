############
# cenns.py
# Code has been edited and commented out to remove uncessary parts
# Only interested in CEvNS NR - no ER or WIMP
# sampling more
############

import numpy as np
import numericalunits as nu
import ROOT
import sys, os.path
#import wimp_rate # commented out as I don't have this
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

operator = "O1"

class MakeWimpMigdalSpectra:
  def __init__(self, CS=1e-40, step=0.001, livetime=515.21, mass=46.7 * 0.97, entries=2000, maxene=2, atomic_masss=39.948):
    self.CS        = CS
    self.maxene    = maxene
    self.step      = step
    self.livetime  = livetime*nu.day/nu.year
    self.mass      = mass*nu.kg/(1000 * nu.kg)
    self.entries   = entries
    self.entries2d = entries*10

    #self.wimp      = wimp_rate.WimpNR(atomic_masss)
    # self.det       = response.DetectorResponse()
    self.ene       = np.linspace(0,maxene,int(maxene/step)+1)

    if atomic_masss < 100 :
    #  self.dfmaps  = pd.read_csv("~/Downloads/%sData_Ar_Xe/ArFullData%s.dat"%(operator,operator), sep='\s+', header=None, names=['M','ER','NR','Rate'])
     # self.dfnorm  = pd.read_csv("~/Downloads/%sData_Ar_Xe/Normalisations_Ar%s.dat"%(operator,operator), sep='\s+')
       self.s2bindf     = self.build_map ("data/ds20k_s2_binning_info.csv")
       self.s2_bin_widths = (self.s2bindf['end_ne'] - self.s2bindf['start_ne']).values
       self.s2_bin_edges = np.concatenate([self.s2bindf['start_ne'].values, [self.s2bindf['end_ne'].iloc[-1]]])
       self.s2_bin_centers     = self.s2bindf['linear_center_ne'].values
     # self.responseER = self.build_map ("../detector_response/ds20k_ne_response_er.csv")
       self.responseNR = self.build_map ("data/remake_ds20k_ne_response_nr.csv")
       self.NRlowbound = 0.1  #self.responseNR['energy_bin_start_kev'].iloc[0]
       self.ERlowbound = 0.025 #self.responseER['energy_bin_start_kev'].iloc[0]
# commenting out as atomic mass is <100
#    else :
#      self.dfmaps  = pd.read_csv("~/Downloads/%sData_Ar_Xe/XeFullData%s.dat"%(operator,operator),sep='\s+', header=None, names=['M','ER','NR','Rate'])
#      self.dfnorm  = pd.read_csv("~/Downloads/%sData_Ar_Xe/Normalisations_Xe%s.dat"%(operator,operator), sep='\s+')
#      self.s2bindf     = self.build_map ("../input_spectra/xe1t_s2only_data_release/s2_binning_info.csv")
#      self.s2_bin_widths = (self.s2bindf['end_pe'] - self.s2bindf['start_pe']).values
#      self.s2_bin_edges = np.concatenate([self.s2bindf['start_pe'].values, [self.s2bindf['end_pe'].iloc[-1]]])
#      self.s2_bin_centers     = self.s2bindf['linear_center_pe'].values
#      self.responseER = self.build_map ("../input_spectra/xe1t_s2only_data_release/s2_response_er.csv")
#      self.responseNR = self.build_map ("../input_spectra/xe1t_s2only_data_release/s2_response_nr.csv")
#      self.NRlowbound = 10 #self.responseNR['energy_bin_start_kev'].iloc[0]
#      self.ERlowbound = 0.186 #self.responseER['energy_bin_start_kev'].iloc[0]

#    self.table_masses_GeV                      = self.dfnorm["mass"].unique()
    #self.table_integrals_per_kg_day_1Eneg40cm2 = self.dfnorm["integrated_rate"]  # this is for NR+Migdal at 1E-40 cm2
    #self.table_er_fraction                     = self.dfnorm["ratio"]

    #self.E_EM_max        = self.dfnorm["E_EM_max"]
    #self.E_R_max         = self.dfnorm["E_R_max"]
    #self.E_R_intervals   = self.dfnorm["E_R_interval"]

    #self.E_R_steps       = self.dfnorm["E_R_max"]
    #self.E_EM_steps      = self.dfnorm["E_EM_steps"]
    #self.E_EM_intervals  = self.dfnorm["E_EM_interval"]

    #self.energy_bin_starts_kev_er = self.responseER['energy_bin_start_kev'].values
    #self.energy_bin_end_kev_er    = self.responseER['energy_bin_end_kev'].values
    #self.energy_bins              = np.concatenate([ self.energy_bin_starts_kev_er, [self.energy_bin_end_kev_er[-1]]])

    # self.fout = None

  def build_map (self,name) :
    mymap     = pd.read_csv(name,sep=",")
    return mymap

  def get_sample(self,enes,spectrum):
    '''
    sample from distribution - but only above the threshold
    To save on computation we want to only sample energies which above NRlowbound
    '''
    # find the index of first value in the list above threshold
    index = 0
    while enes[index] < self.NRlowbound:
        index += 1
    # remove data from spec below threshold - and renormalise so it adds up to 1 still
    enes_above_threshold = enes[index:]
    spec_above_threshold = spectrum[index:]
    # renormalise spec
    normaliser = 1 / sum(spec_above_threshold)
    norm_spec_above_threshold = normaliser * spec_above_threshold
    # now sample energies
    ene = np.random.choice(enes_above_threshold,self.entries,p=norm_spec_above_threshold)
    return ene

  def get_cenns(self)  :
    '''
    function which unpacks the CEvNS spectrum and normalises the spectrum by the sum
    '''
    # changed this for my data
    enes , datasp = np.loadtxt('data/rateTOTAr_old_spec_for_comparison.txt', delimiter='\t', unpack=True)
    # dividing by the sum??
    print(sum(datasp))
    datasp /= np.sum(datasp)
    print(sum(datasp))
    return enes, datasp

  def get_wimp_spectrum_keV(self, Mw):
    spectrum  = self.wimp.rate_wimp_nr(erec=self.ene,  mw=Mw, sigma_nucleon = self.CS)*self.mass*self.livetime*self.step
    self.norm  = np.sum(spectrum)
    return    spectrum/self.norm

  def convert_one_energy (self, eNR, response, s2_bin_centers) :
    #find the energy bin
    if eNR > 5. : return 0 
    #index = response['energy_bin_end_kev'].last_valid_index()
    else :
      index = response[(response.energy_bin_start_kev>eNR)]['energy_bin_start_kev'].index[0]
    #extract the s2 distribution for the given energy and normalise
    s2slice  = np.array ( response.iloc[index].tolist()[3:] ) # probabilty dist for the selected energy bin
    if len(s2slice)==0 : return 0
    weight  =  s2slice.sum()
    s2slice /= weight # normalising by dividing by sum (should be 1?? so shouldn't make a diff)
    #sample one value according to the S2 distribution
    s2bin    = np.random.choice( len(s2slice), p=s2slice, size=None)
    #linearly scale inside the bin width
    ne = s2_bin_centers[ s2bin ] * eNR / response['energy_kev'][index]
    if np.random.random() > weight : ne = 0 # confused by this
    return ne

  def convert (self, energies):
    nlist = []
    # loop through the array of energies that was created
    for e in energies :
      nNR = 0
      nER = 0
      eNR = e[1]
      eER = e[0] # always zero so there will be no ER contribution
      if eNR < self.NRlowbound : nNR = 0 # below threshold
      else :
        nNR = self.convert_one_energy (eNR, self.responseNR, self.s2_bin_centers )
      # ER
      if eER < self.ERlowbound  : nER = 0
      else :
        nER = self.convert_one_energy (eER, self.responseER, self.s2_bin_centers )

      ss = nER + nNR
      nlist.append(ss)
    return np.array(nlist)

  def get_wimp_migdal_spectrum_ne(self,  table_index=0):
    '''
    getting the NR spec - main function
    '''

    # extract a sample of energies for the pure NR case
    enes, sp  = self.get_cenns()
    energies0 = self.get_sample( enes, sp ) # sampling from the CEvNS spec 1,000,000 times
    energies  = np.array([ [0,x] for x in energies0]) # making an array of the sampled energies with 0

    #convert the NR energies in detected spectra, for the case of pure NR
    #coult use a convolution, but want to keep the ev by event structure (even if much slower)
    measured_energies = self.convert (energies)

    """ recostrcted spectra """
    hNR    = ROOT.TH1D("hNR",   "hNR",    len(self.s2_bin_edges)-1, self.s2_bin_edges)
    hEne   = ROOT.TH1D("hEne", "hEne",    1000, 0, 1)
    hNR.Sumw2() # sum the squares of the weights

    # filling histogram
    for x in measured_energies :
      if (x>0) : hNR.Fill(x)

    for x in energies0 :
      hEne.Fill(x)

    # 922.472  is the number of events expected per t.yr to produce a signal > 25 eV
    hNR.Scale  ( 922.472 / self.entries ) # shouldn't need this ideally

    c = ROOT.TCanvas ('c','c',780,700)
    hNR.Draw()
    c.Update()
    c.Draw()

    return hNR, hEne



if __name__ == '__main__':


  Mw        = 50           # GeV/c2
  CS        = 1e-40        # cm^2
  maxene    = 2           # keV      #need to control these two
  step      = 0.001         # keV  #need to control these two
  livetime  = 1           # days
  mass      = 1 # kg
  atomic_masss = 131.2
  atomic_masss = 39.9

  # livetime  = livetime*nu.day/nu.year    # transform in years
  # mass      = mass*nu.kg/(1000 * nu.kg)  # transform in tons

  import matplotlib.pyplot as plt
  #wimpMigdal    = MakeWimpMigdalSpectra(CS=CS, atomic_masss=131.2)
  wimpMigdal    = MakeWimpMigdalSpectra(CS=CS, atomic_masss=atomic_masss, step=step, livetime=livetime, mass=mass, maxene=maxene)

  hNR, hEne= wimpMigdal.get_wimp_migdal_spectrum_ne()

  fout  = ROOT.TFile("data/spectra_cenns.root", "recreate")
  hNR.Write("hNR_per_tonne_yr") ;
  hEne.Write("hEne_per_tonne_yr") ;
