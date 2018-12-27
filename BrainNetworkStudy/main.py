# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 23:42:01 2018

@author: Mattia
"""

import pyedflib
import mylib 
import pandas as pd
import mvar
import numpy as np
import matplotlib.pyplot as plt

#%% Load data and store into a dataframe

f = pyedflib.EdfReader('../data/S001R01.edf')
# mylib.print_some_stuff(f)

data = mylib.edfToDataFrame(f)
f._close()

#%%

data = data.T.values
N, n = data.shape
p = 3
A_est, sigma = mvar.mvar_fit(data, p)
sigma = np.diag(sigma)  # DTF + PDC support diagonal noise

# compute DTF
D, freqs = mvar.DTF(A_est, sigma)
#mvar.plot_all(freqs, D, 'DTF') it takes a lot to run this

# compute PDC
P, freqs = mvar.PDC(A_est, sigma)
#plot_all(freqs, P, 'PDC') it takes a lot to run this

freqs
