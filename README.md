*Senior Thesis*
**Kyle Campbell**
***Implementing Li and Panigrahi***
Comparing min-cut algorithms - Hao-Orlin and Nagamochi-Ibaraki

Examples of how to construct graphs using the lemon library and run the Li Panigrahi Min-Cut algorithm on them
can be found in the files test_isolating_cut.cc and test_graphs.

Preliminary results comparing runtimes of min cut algorithms are here. As this implementation continues to be improved,
these results will be updated.

Graph is a Bicycle Wheel with 1024 Nodes.

Nagamochi Ibaraki exec. time (ms): 0.19403 Nodes:1024
Hao Orlin exec. time (ms): 115.456 Nodes:1024
Li Panigrahi exec. time (ms): 52.6606 Nodes:1024 Subset Size12
Li Panigrahi exec. time (ms): 100.17 Nodes:1024 Subset Size100
Li Panigrahi exec. time (ms): 189.48 Nodes:1024 Subset Size200
Li Panigrahi exec. time (ms): 429.971 Nodes:1024 Subset Size500
Li Panigrahi exec. time (ms): 598.734 Nodes:1024 Subset Size700
Li Panigrahi exec. time (ms): 697.931 Nodes:1024 Subset Size800
Mincut value: NI | HO | LP 512 | 512 | 512

Graph is a Bicycle Wheel with 10000 nodes.

Nagamochi Ibaraki exec. time (ms): 1.45628 Nodes:10000
Hao Orlin exec. time (ms): 71859.7 Nodes:10000
Li Panigrahi exec. time (ms): 19686 Nodes:10000 Subset Size10
Li Panigrahi exec. time (ms): 3226.73 Nodes:10000 Subset Size50
Li Panigrahi exec. time (ms): 3460.86 Nodes:10000 Subset Size100
Li Panigrahi exec. time (ms): 3121.31 Nodes:10000 Subset Size250
Li Panigrahi exec. time (ms): 4845.52 Nodes:10000 Subset Size500
Li Panigrahi exec. time (ms): 8766.47 Nodes:10000 Subset Size1000
Li Panigrahi exec. time (ms): 17179.3 Nodes:10000 Subset Size2000
Li Panigrahi exec. time (ms): 33588.4 Nodes:10000 Subset Size4000
Li Panigrahi exec. time (ms): 42376 Nodes:10000 Subset Size5000
Li Panigrahi exec. time (ms): 63081 Nodes:10000 Subset Size7500
Mincut value: NI | HO | LP 5000 | 5000 | 512

Graph is a random graph with edge probability 0.5 and 1024 nodes.

Nagamochi Ibaraki exec. time (ms): 161.027 Nodes:1024
Hao Orlin exec. time (ms): 1618.41 Nodes:1024
Li Panigrahi exec. time (ms): 1565.4 Nodes:1024 Subset Size12
Li Panigrahi exec. time (ms): 11719.4 Nodes:1024 Subset Size100
Li Panigrahi exec. time (ms): 25186.9 Nodes:1024 Subset Size200
Li Panigrahi exec. time (ms): 60927.7 Nodes:1024 Subset Size500
Li Panigrahi exec. time (ms): 82212.3 Nodes:1024 Subset Size700
Li Panigrahi exec. time (ms): 95130.3 Nodes:1024 Subset Size800
NI | HO | LP 457 | 457 | 457

Graph is a random graph with edge probability 0.25 and 10000 nodes.
Nagamochi Ibaraki exec. time (ms): 187185 Nodes:10000
Hao Orlin exec. time (ms): 2.23312e+06 Nodes:10000
Li Panigrahi exec. time (ms): 344083 Nodes:1024 Subset Size12

TODO:
 * Recover the nodes on either side of MinCut.
 * Mimic Graph Contraction without overhead of copying
    - Do this with Disjoint Set data structure.
 * Put contraction processes in parallel.
 * Write Li Panigrahi randomized algorithm to have something that is a single call to compute mincut.
 * Don't use Preflow to compute maxflows
