#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 11:16:13 2021

@author: EB

Testing out data files
"""

import numpy as np
import uproot

ds20k_exposure = 365*1000*100
ds20k_exposure_tonneyear = 100 

# the data is in number of electrons. use a 4e- threshold
ds20k_firstbin = 4
ds20k_lastbin = 50

dm_array = [1,2,3]

n, dplus, dminus = dm_array

print(n,dplus,dminus)