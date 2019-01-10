# README

### Usage

`python main.py`

### Small World Networks

Duncan Watts and Steven Strogatz found that

[1]: https://en.wikipedia.org/wiki/Small-world_network#cite_note-3	"Watts, Duncan J.; Strogatz, Steven H. (June 1998). &quot;Collective dynamics of &#39;small-world&#39; networks&quot;. Nature. 393 (6684): 440–442"

graphs can be classified according to two important global indexes:

1. clustering coefficient
2. average shortest path

We say a network to be a **small world network** when

- for any given node, neighbors of that node are likely to be neighbors of each other (read: clustering coefficient close to 1)

- the average distance L between two nodes (average shortest path) scales as the logarithm of the number of nodes in the network 
  $$
  L \sim \log N
  $$


#### Small World Network Index

https://en.wikipedia.org/wiki/Small-world_network

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3604768/

https://arxiv.org/pdf/cond-mat/9903108.pdf







The bctpy is the python translation of the https://github.com/aestrivex/bctpy

```
bct.randmio_dir(R, iter)¶
    This function randomizes a directed network, while preserving the in- and out-degree distributions. In weighted networks, the function preserves the out-strength but not the in-strength distributions.
    Parameters:	
    W : NxN np.ndarray
        directed binary/weighted connection matrix
    iter : int
        rewiring parameter. Each edge is rewired approximately iter times.
    Returns:	
    R : NxN np.ndarray
        randomized network
    eff : int
        number of actual rewirings carried out
```

```
bct.latmio_dir(R, iter, D=None)¶
    This function “latticizes” a directed network, while preserving the in- and out-degree distributions. In weighted networks, the function preserves the out-strength but not the in-strength distributions.
    Parameters:	
    R : NxN np.ndarray
        directed binary/weighted connection matrix
    iter : int
        rewiring parameter. Each edge is rewired approximately iter times.
    D : np.ndarray | None
        distance-to-diagonal matrix. Defaults to the actual distance matrix if not specified.
    Returns:	
    Rlatt : NxN np.ndarray
        latticized network in original node ordering
    Rrp : NxN np.ndarray
        latticized network in node ordering used for latticization
    ind_rp : Nx1 np.ndarray
        node ordering used for latticization
    eff : int
        number of actual rewirings carried out
```



### Motif Analysis

To perform motif analysis we make use of the software `mfinder` version 1.2.

It is available at http://www.weizmann.ac.il/mcb/UriAlon/research/network-motifs.

It has not a graphical interface, and we must use the propt, running a command like:

`mfinder1.2 ./data/inputForMA.txt -s 3 -r 1000 `

Where..

- `-s` is the desired subgraph size

- `-r` is the number of random networks to be generated

