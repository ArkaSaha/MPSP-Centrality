# MPSP Centrality

This repository contains the code for computing MPSPs and Betweenness Centrality in uncertain graphs, as described in the VLDB 2021 paper titled [_Shortest Paths and Centrality in Uncertain Networks_](http://www.vldb.org/pvldb/vol14/p1188-khan.pdf) by _Arkaprava Saha, Ruben Brokkelkamp, Yllka Velaj, Arijit Khan, Francesco Bonchi_.

## Generating Synthetic Graphs and Queries

We generate synthetic graphs according to the Erdos-Renyi (ER) and Barabasi-Albert (BA) models, with various values of number of nodes and edge-node ratio. Also, we generate node pairs as queries for every graph. One such graph and its queries can be generated with the following command:
```
python generate_random_graphs.py <graph-type> <number-of-nodes> <edge-node-ratio>
```
where `graph-type` can be either ER or BA. This command generates a synthetic graph in the location `<graph-type>/<graph-type>_<number-of-nodes>_<number-of-edges>.graph` and the corresponding queries in the location `<graph-type>/<graph-type>_<number-of-nodes>_<number-of-edges>.queries`.

## Usage

### MPSPs

Our proposed method, with hyperparameters `k` (number of MPSPs), `m` (number of Dijkstra+MC runs in Phase 1) and `N` (number of MC samples in Phase 2), can be run using the following command:
```
./mpsp <path-to-graph> <path-to-queries> <path-to-output> <k> <m> <N>
```

The running time and quality results can be respectively obtained using the commands below:
```
grep Average\ Total\ Time <path-to-output>
grep Average\ SP\ Probability <path-to-output>
```

### Centrality

We compare the top-`k` most central nodes according to 4 betweenness centrality methods. Methods 1, 2 and 3 can be run using the following command:
```
./mpsp <path-to-graph> <path-to-output> <k>
```
The code for method 4 can be found [here](https://github.com/XNetLab/ProbGraphBetwn).

The output is first the node ids and betweenness values of the top-`k` produced by method 1, then the top-`k` produced by method 3 and the top-k for method 2.

Method 1 can also be implemented in parallel, i.e., using multiple threads. This version can be run using the following command:
```
./mpsp <path-to-graph> <path-to-output> <k> <number-of-threads>
```

