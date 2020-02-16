## Graph files

The graph files first contain a line with two (space separated) integer n and m.
Here n is the number of vertices and m the number of edges. Vertices are numbered 0 to n-1.
Then follow m lines of the form

u v l p

representing an edge from node u to v of length l and probability p.

u, v and l are integers while p is a double

### Query files

Query files consist of 500 lines of two (space separated) integers u and v.

The first 100 pairs are random nodes u and v.

The next 100 pairs are nodes which are 2 hops apart

The next 100 pairs are nodes which are 4 hops apart

The next 100 pairs are nodes which are 6 hops apart

The next 100 pairs are nodes which are 8 hops apart

