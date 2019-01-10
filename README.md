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

The `bctpy` is the `Python` translation of the Brain Connectivity Toolbox (natively in `MATLAB`) https://github.com/aestrivex/bctpy

- download the repository
- run  `python setup.py install` 

We make use of two methods:

`bct.randmio_dir(R, iter)`:  this function randomizes a directed network, while preserving the in- and out-degree distributions. In weighted networks, the function preserves the out-strength but not the in-strength distributions.

`bct.latmio_dir(R, iter, D=None)`: this function “latticizes” a directed network, while preserving the in- and out-degree distributions. In weighted networks, the function preserves the out-strength but not the in-strength distributions.



### Motif Analysis

To perform motif analysis we make use of the software `mfinder` version 1.2.

It is available at http://www.weizmann.ac.il/mcb/UriAlon/research/network-motifs.

It has not a graphical interface, and we must use the propt, running a command like:

`mfinder1.2 ./data/inputForMA.txt -s 3 -r 1000 `

Where..

- `-s` is the desired subgraph size

- `-r` is the number of random networks to be generated

