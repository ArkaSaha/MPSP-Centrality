# Generating Synthetic Graphs and Queries

As mentioned in Section 5.2 of our paper, we generate synthetic graphs according to the Erdos-Renyi (ER) and Barabasi-Albert (BA) models, with various values of number of nodes and edge-node ratio. Also, as mentioned in Section 5.1, we generate node pairs as queries for every graph. One such graph and its queries can be generated with the following command:
```
python generate_random_graphs.py <graph-type> <number-of-nodes> <edge-node-ratio>
```
where `graph-type` can be either ER or BA. This command generates a synthetic graph in the location `<graph-type>/<graph-type>_<number-of-nodes>_<number-of-edges>.graph` and the corresponding queries in the location `<graph-type>/<graph-type>_<number-of-nodes>_<number-of-edges>.graph`.

The format of the generated graph and query files is mentioned [here](../README.md).
