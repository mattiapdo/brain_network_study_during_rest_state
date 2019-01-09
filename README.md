# README

### Usage

`python main.py`



### Motif Analysis

To perform motif analysis we make use of the software `mfinder` version 1.2.

It is available at http://www.weizmann.ac.il/mcb/UriAlon/research/network-motifs.

It has not a graphical interface, and we must use the propt, running a command like:

`mfinder1.2 ./data/inputForMA.txt -s 3 -r 1000 `

Where..

- `-s` is the desired subgraph size

- `-r` is the number of random networks to be generated

