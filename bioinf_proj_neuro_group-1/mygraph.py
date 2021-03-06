# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

#%%
def randargmin(b,**kw):
  """ a random tie-breaking argmax"""
  return np.argmin(np.random.random(b.shape) * (b==b.min()), **kw)

#%%
def applyTreshold(G, threshold):
    while(G.density() > threshold):
        d = np.matrix(G.get_adjacency(attribute = "weight")._get_data())
        masked = np.ma.masked_where(d == 0, d)
        #d = np.matrix([[0,1,0],[1,2,3],[3,2,0]])
        arg_min = np.matrix(np.argwhere(masked == np.amin(masked)))
        nrows = arg_min.shape[0]
        select = np.random.randint(nrows)
        i, j = arg_min[select,][0,0], arg_min[select,][0,1]
        G.delete_edges([(i,j)])

    return G


#%%
def add_brain_layout(g):
    chan_loc = pd.read_csv("./data/channels_topology.csv", sep=';', decimal=',')
    pos = {row[0]: (row[1],row[2]) for row in chan_loc.values}
    g.vs["coordinates"] = [pos[i["label"]] for i in g.vs]
    return(g)

#%%
