# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 23:42:01 2018

@author: Mattia
"""

import pyedflib
import mylib 
import mvar
import mygraph
import numpy as np
from igraph import Graph

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

# compute PDC
P, freqs = mvar.PDC(A_est, sigma)

# take take a graph relative to a specific frequency

freq = 0.111328 # this has to be chosen better

D = D[np.where(freqs == mylib.find_nearest(freqs, freq)),:,:].reshape(64, 64)

G = Graph.Weighted_Adjacency(D.tolist())

G = mygraph.applyTreshhold(G, 0.2)


