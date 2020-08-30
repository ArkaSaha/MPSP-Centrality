# MPSP Centrality

This repository contains the code used for the experiments described in the paper titled

_Shortest Paths and Centrality in Uncertain Networks_
by
_Arkaprava Saha, Ruben Brokkelkamp, Yllka Velaj, Arijit Khan & Francesco Bonchi_

## Requirements
The experiments have been run on an Ubuntu 18.04.4 LTS computer with a 3.7 GHz Xeon processor and 16 GB RAM. All programs are written in `C++17` and make use of [the Boost C++ libraries](https://www.boost.org/).

Generating data is done with `Python 3` using the modules
```
[networkx](https://networkx.github.io/)
[python-igraph](https://igraph.org/python/)
```

## Input Data
See [input specifications](data/README.md), [generating test data](data/Synthetic/README.md) and [road data](data/Real/Road/README.md).

## Usage

### Shortest Paths

#### Our Method
Our proposed method, with hyperparameters k (number of MPSPs), m (number of Dijkstra+MC runs in Phase 1) and N (number of MC samples in Phase 2), can be run using the following command:
```
./mpsp <path-to-graph> <path-to-queries> <path-to-output> <k> <m> <N>
```

#### Baseline Method
To run the baseline algorithm, we also need the output of our method in order to set the candidate generation time threshold (the time taken by our method multiplied by a user-provided factor) for each query. Hence we need to run it using the following command:
```
./wise <path-to-graph> <path-to-queries> <path-to-output-mpsp> <path-to-output> <time-multiplying-factor> <k>
```

#### Output
For every method, the running time and quality results can be respectively obtained using the commands below:
```
grep Average\ Total\ Time <path-to-output>
grep Average\ SP\ Probability <path-to-output>
```

### Effect of Each Phase
The effect of each phase on the performance (Section 5.4 of our paper) can be demonstrated using the command below:
```
./dasfaa <path-to-graph> <path-to-queries> <path-to-output>
```
To display the quality of the results obtained, run the following commands:
```
grep Luby-Karp\ vs\ Majority <path-to-output>
grep Luby-Karp\ vs\ Horvitz-Thompson <path-to-output>
```
These commands show the quality of the path returned by Luby-Karp (LK) compared to Majority and Horvitz-Thompson (HT) respectively. Specifically, each of them displays, for every query category, 3 space-separated numbers denoting (in order) the fraction of queries for which LK returns better, the same, or worse paths compared to the other method (Majority or HT).

### Centrality
As mentioned in Section 5.8 of our paper, we compare the top-k most central nodes according to 4 betweenness centrality methods. Methods 0, 1 and 2 can be run using the following command:
```
./mpsp <path-to-graph> <path-to-output> <k>
```
The code for method 3 can be found [here](https://github.com/XNetLab/ProbGraphBetwn).

The output is first the node ids and betweenness values of the top k produced by method 0, then the the top k produced by method 2 and if the graph has at most 1000 nodes also the top k for method 1.

