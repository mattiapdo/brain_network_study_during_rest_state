# -*- coding: utf-8 -*-

import numpy as np

def randargmin(b,**kw):
  """ a random tie-breaking argmax"""
  return np.argmin(np.random.random(b.shape) * (b==b.min()), **kw)

def applyTreshold(G, threshold):
    while(G.density() > threshold):
        """adj = np.matrix(G.get_adjacency(attribute = "weight")._get_data())
        masked = np.ma.masked_where(adj == 0, adj)
        G.delete_edges([np.unravel_index(randargmin(masked), masked.shape)])"""
        d = np.matrix(G.get_adjacency(attribute = "weight")._get_data())
        masked = np.ma.masked_where(d == 0, d)
        #d = np.matrix([[0,1,0],[1,2,3],[3,2,0]])
        arg_min = np.matrix(np.argwhere(masked == np.amin(masked)))
        nrows = arg_min.shape[0]
        select = np.random.randint(nrows)
        i, j = arg_min[select,][0,0], arg_min[select,][0,1]
        G.delete_edges([(i,j)])
            
    return G



"""d = np.matrix([[0,1,0],[1,2,3],[3,2,0]])
arg_min = np.matrix(np.argwhere(d == np.amin(d)))
nrows = arg_min.shape[0]
select = np.random.randint(nrows)
i, j = arg_min[select,][0,0], arg_min[select,][0,1]"""
