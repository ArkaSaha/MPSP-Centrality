# Downloading Road Graphs and Queries

As mentioned in Section 5.3 of our paper, we generate graphs of the road networks of the cities Brno, Porto, Rome and San Francisco. Also, as mentioned in Section 5.1, we generate node pairs as queries for each graph. The graph files are named `<city>.graph` and the query files are named `<city>.queries`. Owing to space limitations, the graph files are compressed. They can be decompressed using the command below:
```
bzip2 -d <city>.graph.bz2
```

The format of the generated graph and query files is mentioned [here](../README.md).
