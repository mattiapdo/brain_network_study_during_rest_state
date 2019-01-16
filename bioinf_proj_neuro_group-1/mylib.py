# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 00:42:40 2018

@author: Mattia
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyedflib
import mylib
import mvar
import mygraph
from igraph import Graph
import re
import igraph
import os

def edfToDataFrame(f):
    # number of rows
    n = f.getNSamples()[0]
    # number of columns
    m = f.signals_in_file
    # init a matrix
    sigbufs = np.zeros([n, m])
    # fill the data
    for i in np.arange(m):
        sigbufs[:,i] = f.readSignal(i)
    # create pandas dataframe
    data = pd.DataFrame(sigbufs)
    # rename columns
    data.columns= f.getSignalLabels()

    return data

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def receive_data(path, freq = 10, threshold = 0.1, DTF = True):
    """
    This function builds a brain network given in input:
    - path of the EDF file
    - freq: the frequency we want to conduct analysis on
    - threshold: a threshold on the graph's density
    - DFT: boolean: if True estimates the Direct Tranfer Function (DFT)
                    if False estimates the Partial Directed Coherence (PDC)
    """
    #%% Load data and store into a dataframe

    print("Loading data from", path, "...")
    f = pyedflib.EdfReader(path)
    print("...done")

    data = mylib.edfToDataFrame(f)
    f._close()
    T = 1./f.getSampleFrequencies()[1]
    measures = data.T.values
    N, n = measures.shape
    p = 3
    A_est, sigma = mvar.mvar_fit(measures, p)
    sigma = np.diag(sigma)  # The noise for each time series
    if DTF:
        # compute DTF
        print("Computing DFT...")
        Adj, freqs = mvar.DTF(A_est, sigma, T)
        print("... done")
    else:
        # compute PDC
        print("Computing PDC")
        Adj, freqs = mvar.PDC(A_est, sigma, T)
        print("... done")

    # take take a graph relative to a specific frequency expressed in Hz
    print("Building network to analyse interactions at", freq, "Hz")
    # select the correspondent weights
    Adj = Adj[np.where(freqs == mylib.find_nearest(freqs, freq)),:,:].reshape(64, 64)
    G = Graph.Weighted_Adjacency(Adj.tolist(), mode = 0) # mode=0 is for directed / mode=1 is for indirected graph

    # get channel names, cleaning replacing any dot with a void charachter
    G.vs["label"] = list(map(lambda x: re.sub('\.', '', x), data.columns.values))
    G = mygraph.add_brain_layout(G)

    print("Applying a threshold on network density of", round(threshold*100), "% ...")
    G = mygraph.applyTreshold(G, threshold)
    print("... done")

    G.es["arrow_size"] = [0.5 for edge in G.es]

    return G

def graph_indices(G, weights = False):

    if weights:
        print("Macroscopic Network Analysis")
        # global clustering coefficient
        # https://igraph.org/python/doc/igraph.Graph-class.html#transitivity_avglocal_undirected
        print("global clustering coefficient = ", round(G.transitivity_avglocal_undirected(weights = "weight"), 3))

        # average path length
        # https://igraph.org/python/doc/igraph.GraphBase-class.html#average_path_length
        print("global average path lenght = ", round(G.average_path_length(directed = True), 3))
        print("\t the latter is calculated regardless from the weights")

        print("Microscopic Network Analysis")
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

    else:
        print("Macroscopic Network Analysis")
        # global clustering coefficient
        # https://igraph.org/python/doc/igraph.Graph-class.html#transitivity_avglocal_undirected
        print("global clustering coefficient = ", round(G.transitivity_avglocal_undirected(), 3))

        # average path length
        # https://igraph.org/python/doc/igraph.GraphBase-class.html#average_path_length
        print("global average path lenght = ", round(G.average_path_length(directed = True), 3))

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

def plot_analysis(densities, clustering_coeffs, average_path_lengths):
    f, axarr = plt.subplots(2, sharex=True, figsize = (9,9))
    axarr[0].plot(densities, clustering_coeffs[::-1], color = "green")
    axarr[0].set_title('Global Indices vs Graph Density')
    axarr[0].set_xlabel("density")
    axarr[0].set_ylabel("clustering coefficient")
    axarr[0].grid(linestyle = ":")
    axarr[1].plot(densities, average_path_lengths[::-1], color = "red")
    axarr[1].set_xlabel("density")
    axarr[1].set_ylabel("average path length")
    axarr[1].grid(linestyle = ":")
    plt.show()

def write_inputFileForMotifAnalysis(G, file):
    with open(file, 'a') as file:
        for source in range(G.vcount()):
            targets = G.get_adjlist()[source]
            for target in targets:
                if target != source:
                    file.write(str(source) + "  " + str(target) +"  " + str(1)+"\n")
    file.close()
    return


def swi(L, Ll, Lr, C, Cl, Cr):
    return (L-Ll)/(Lr-Ll)*(C-Cr)/(Cl-Cr)
