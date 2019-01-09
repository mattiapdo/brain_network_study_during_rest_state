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
import pandas as pd
from igraph import Graph
import re
import igraph

#%% Load data and store into a dataframe
path = './data/S001R01.edf'
print("Loading data from", path)
f = pyedflib.EdfReader(path)
# mylib.print_some_stuff(f)

data = mylib.edfToDataFrame(f)
f._close()
T = 1./f.getSampleFrequencies()[1]

#%%
# 1.1

measures = data.T.values
N, n = measures.shape
p = 3
A_est, sigma = mvar.mvar_fit(measures, p)
sigma = np.diag(sigma)  # The noise for each time series
# compute DTF
print("Computing DFT")
D, freqs = mvar.DTF(A_est, sigma, T)
# compute PDC
print("Computing PDC")
# P, freqs = mvar.PDC(A_est, sigma)
#%%
# take take a graph relative to a specific frequency
freq = 10. # this is in Hz and roughly corresponds to Alpha (mu) ERD
#freq = 28   # this is in Hz and roughly corresponds to Beta ERS
print("Building network to analyse interactions at", freq, "Hz")
# select the correspondent weights
D = D[np.where(freqs == mylib.find_nearest(freqs, freq)),:,:].reshape(64, 64)
G = Graph.Weighted_Adjacency(D.tolist(), mode = 0) # mode=0 is for directed / mode=1 is for indirected graph

# get channel names, cleaning replacing any dot with a void charachter
G.vs["label"] = list(map(lambda x: re.sub('\.', '', x), data.columns.values))
# %%
threshold = 0.05
print("Applying a threshold on network density of", round(threshold*100), "%")
G = mygraph.applyTreshold(G, threshold)

#%%
G = mygraph.add_brain_layout(G)

#%%
print("Making a topological representation")
visual_style = {}
visual_style["vertex_size"] = 20
visual_style["vertex_color"] = "white"
visual_style["vertex_label"] = G.vs["label"]
visual_style["edge_width"] = [1 + 2 * int(weight) for weight in G.es["weight"]]
visual_style["layout"] = G.vs["coordinates"]
igraph.plot(G, **visual_style)

#%%
# 2.1
print("Macroscopic Network Analysis")
# global clustering coefficient
# https://igraph.org/python/doc/igraph.Graph-class.html#transitivity_avglocal_undirected
print("global clustering coefficient = ", round(G.transitivity_avglocal_undirected(), 3))

# average path length
# https://igraph.org/python/doc/igraph.GraphBase-class.html#average_path_length
print("global average path lenght = ", round(G.average_path_length(directed = False, unconn=True), 3))

print("Microscopic Network Analysis")
local_ind = pd.DataFrame.from_dict(      {
                "channel" : G.vs["label"] ,
                "degree" : G.strength(mode = 3),
                "out-degree" : G.strength(mode = 1),
                "in-degree" : G.strength(mode = 2)
                }
        )

local_ind.set_index("channel", inplace = True)
print("Top 10 channels by degree\n", local_ind.sort_values(by = "degree", ascending = False)[1:10]["degree"])
print("Top 10 channels by OUT degree\n", local_ind.sort_values(by = "out-degree", ascending = False)[1:10]["out-degree"])
print("Top 10 channels by IN degree\n", local_ind.sort_values(by = "in-degree", ascending = False)[1:10]["in-degree"])



# 2.7
local_ind = pd.DataFrame.from_dict(      {
                "channel" : G.vs["label"] ,
                "strength" : G.strength(mode = 3,  weights = "weight"),
                "out-strength" : G.strength(mode = 1,  weights = "weight"),
                "in-strength" : G.strength(mode = 2,  weights = "weight")
                }
        )

local_ind.set_index("channel", inplace = True)
print("Top 10 channels by generalized degree\n", local_ind.sort_values(by = "strength", ascending = False)[1:10]["strength"])
print("Top 10 channels by generalized OUT degree\n", local_ind.sort_values(by = "out-strength", ascending = False)[1:10]["out-strength"])
print("Top 10 channels by generalized IN degree\n", local_ind.sort_values(by = "in-strength", ascending = False)[1:10]["in-strength"])

#%%
# 2.4 Study the behaviour of global graph indices in function of network density
print("Studying the behaviour of global graph indices in function of network density...")
densities = [0.05, 0.1, 0.2, 0.3, 0.5] # edit these with the ones provided by the prof

clustering_coeffs = []
average_path_lengths =  []
# there is a specific reason why we start from higher densities and then we
# decrease it: the function applyThreshold takes a graph and simply picks less 
# important edges (to which correspond low weigths) erasing them from the network:
# to do that it has to find the minimum above N* x N* weigths. 
# This is quite computationally expensive in our simple implementation:
# to speed up a bit our analysis we put the function in the condition of having
# at each step a smaller amount N*, avoiding to replicate the same minimizations
# multiple times
G = Graph.Weighted_Adjacency(D.tolist(), mode = 0)
for density in densities[::-1]:
    G = mygraph.applyTreshold(G, density)
    clustering_coeffs.append(G.transitivity_avglocal_undirected())
    average_path_lengths.append(G.average_path_length(directed = False, unconn=True))
    
#%%

print("Ploting the result of the analysis")
mylib.plot_analysis(densities, clustering_coeffs, average_path_lengths)

#%% Motif analysis

path = './data/inputForMA.txt'
print("Writing input file for Motif Analysis to", path)
mylib.write_inputFileForMotifAnalysis(G, file = path)
print("\nDone")
