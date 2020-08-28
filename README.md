# MPSP Centrality

This repository contains the code used for the experiments described in the paper titled

_Shortest Paths and Centrality in Uncertain Networks_
by
_Arkaprava Saha, Ruben Brokkelkamp, Yllka Velaj, Arijit Khan & Francesco Bonchi_

## Requirements
The experiments are written in `c++17` and make use of [the Boost C++ libraries](https://www.boost.org/)

Generating data is done with `python 3` using the modules
```
[networkx](https://networkx.github.io/)
[python-igraph](https://igraph.org/python/)
```

## Input data
See [input specifications](data/README.md) and [generate test data](data/Synthetic/README.md)

## Compiling
Compiling the `c++` code can be done with
```
make all
```

## Usage
### Shortest paths
The mpsp experiments can be done using
```
./mpsp <path-to-graph> <path-to-queries> <path-to-output>
```

To run the baseline algorithm we also need the output from the mpsp experiments
```
./wise <path-to-graph> <path-to-queries> <path-to-output-mpsp> <path-to-output> <time-multiplying-factor>
```

### Betweenness
```
./mpsp <path-to-graph> <path-to-output> <k>
```
