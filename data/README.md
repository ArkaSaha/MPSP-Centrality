## Graph files

The graph files first contain a line with two (space separated) integer n and m.
Here n is the number of vertices and m the number of edges. Vertices are numbered 0 to n-1.
Then follow m lines of the form

u v l p

representing an edge from node u to v of length l and probability p.

u, v and l are integers while p is a double

## Query files

The query files first contain a line with an integer H.
Then the following repeats for H times.
A line with two integer h and n. h stands for how many hops the query nodes are apart and n is the number of queries.
Then n lines follow with two integers u and v, indicating a query for finding a path from node u to v.

