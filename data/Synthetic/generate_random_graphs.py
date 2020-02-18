#!/usr/bin/env python3

import networkx as nx
import random
import os

import igraph

# create directories
for dir_name in ["ER", "BA", "SF", "BP"]:
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


SEED = 12345
random.seed(SEED)
MAX_EDGE_LENGTH = 1000


def generate_queries(g, filename):
    q = open(filename, "w")
    n = g.number_of_nodes()
    m = g.number_of_edges() 

    # 0 represents the completely random pairs, 2, 4, 6 and 8 will contain pairs that are 2, 4, 6, 8 hops apart
    hops_dists = [0, 2, 4, 6, 8]
    queries = {h : set() for h in hops_dists}

    queries_per_category = 100
    queries_so_far = 0
    i = 0
    
    while queries_so_far < len(queries) * queries_per_category:
        i += 1
        if i % 10000 == 0:
            print("{} - {}, ".format(i, queries_so_far), end="", flush=True)
        s = random.randrange(n)
        t = random.randrange(n)
        if s == t:
            continue
        try:
            hops = nx.shortest_path_length(g, source=s, target=t)
            if hops != 0  and hops in queries  \
                          and len(queries[hops]) < queries_per_category \
                          and (s,t) not in queries[hops]:
                queries[hops].add((s, t))
                queries_so_far += 1

            if len(queries[0]) < queries_per_category and (s,t) not in queries[0]:
                queries[0].add((s, t))
                queries_so_far += 1

        except nx.NetworkXNoPath:
            continue
    print()

    # Write the queries to file.
    # The first 100 are the random queries, then 100 two hops queries, ..., lastly 100 eight hop queries
    for hops in hops_dists:
        for s, t in queries[hops]:
            q.write("{} {}\n".format(s, t))

    q.close()



# Another option to (faster) generate queries, but here the pairs won't be picked uniformly at random
# First a start vertex is picked uniformly at random and then the end vertex is picked uniformly at random from
# the vertices at a certain hop distance
def generate_queries_skewed(g, filename):
    q = open(filename, "w")
    n = g.number_of_nodes()
    m = g.number_of_edges() 


    queries = {0: set(), 2:set(), 4:set(), 6:set(), 8:set()}

    queries_per_category = 100
    queries_so_far = 0
    i = 0

    # generate queries at certain hop distance
    for h in [2, 4, 6, 8]:
        queries_so_far = 0
        while queries_so_far < queries_per_category:
            s = random.randrange(n)

            nodes_at_dist_at_most_h = nx.single_source_shortest_path_length(g, source=s, cutoff=h)

            t_options = [k for k in nodes_at_dist_at_most_h if nodes_at_dist_at_most_h[k] == h]

            if len(t_options) == 0: 
                continue

            t = random.choice(t_options)
            if (s,t) not in queries[h]:
                queries[h].add((s,t))
                queries_so_far += 1

    # generate completely random queries
    queries_so_far = 0
    while queries_so_far < queries_per_category:
        s = random.randrange(n)
        t = random.randrange(n)
        if s == t:
            continue
        try:
            if len(queries[0]) < queries_per_category and (s,t) not in queries[0]:
                queries[0].add((s, t))
                queries_so_far += 1

        except nx.NetworkXNoPath:
            continue
    for hops in queries:
        for s, t in queries[hops]:
            q.write("{} {}\n".format(s, t))

    q.close()




def ER(graph_sizes):
    print("Generating ER graphs")
    # Random graph with n vertices and m = 2*n edges, store in ER folder
    # https://networkx.github.io/documentation/stable/reference/generated/networkx.generators.random_graphs.gnm_random_graph.html
    for n in graph_sizes:
        m = 2*n

        print("Generating ER graph with {} nodes and {} edges".format(n, m))
        g = nx.gnm_random_graph(n, m, directed=True)

        f = open("ER/ER_{}_{}.graph".format(n, m), "w")
        f.write("{} {}\n".format(n, m))
        for u, v in g.edges:
            f.write("{} {} {} {}\n".format(u, v, random.randint(1, MAX_EDGE_LENGTH), random.random()))
        f.close()

        generate_queries(g, "ER/ER_{}_{}.queries".format(n, m))



def BA(graph_sizes):
    print("Generating BA graphs")
    # https://igraph.org/python/doc/igraph.GraphBase-class.html#Barabasi
    for n in graph_sizes:
        BA_m = 20 

        print("Generating BA graph with {} nodes".format(n))
        igraph_g = igraph.Graph.Barabasi(n, BA_m , directed=True)
        g = nx.DiGraph(igraph_g.get_edgelist())

        m = g.number_of_edges()

        print("Graph has {} edges".format(m))

        f = open("BA/BA_{}_{}.graph".format(n, m), "w")
        f.write("{} {}\n".format(n, m))
        for u, v in g.edges:
            f.write("{} {} {} {}\n".format(u, v, random.randint(1, MAX_EDGE_LENGTH), random.random()))
        f.close()

        generate_queries(g, "BA/BA_{}_{}.queries".format(n, m))

def BP(graph_sizes):
    print("Generating BP graphs")
    for n in graph_sizes:
        m = 2*n 

        print("Generating BP graph with {} nodes and {} edges".format(n, m))

        g = nx.algorithms.bipartite.generators.gnmk_random_graph(n//2, n//2, m, directed=True)

        m = g.number_of_edges()

        # also add the reverse edges
        g.add_edges_from([(v,u) for (u,v) in g.edges])

        m = g.number_of_edges()

        print("Graph has {} edges".format(m))
        f = open("BP/BP_{}_{}.graph".format(n, m), "w")
        f.write("{} {}\n".format(n, m))
        for u, v in g.edges:
            f.write("{} {} {} {}\n".format(u, v, random.randint(1, MAX_EDGE_LENGTH), random.random()))
            f.write("{} {} {} {}\n".format(v, u, random.randint(1, MAX_EDGE_LENGTH), random.random()))
        f.close()

        generate_queries(g, "BP/BP_{}_{}.queries".format(n, m))

def SF(graph_sizes):
    print("Generating SF graphs")
    # https://igraph.org/python/doc/igraph.GraphBase-class.html
    for n in graph_sizes:

        print("Generating SF graph with {} nodes".format(n))
        g = nx.scale_free_graph(n)
        m = g.number_of_edges()

        print("Graph has {} edges".format(m))

        f = open("SF/SF_{}_{}.graph".format(n, m), "w")
        f.write("{} {}\n".format(n, m))
        for u, v, _ in g.edges:
            f.write("{} {} {} {}\n".format(u, v, random.randint(1, MAX_EDGE_LENGTH), random.random()))
        f.close()

        generate_queries(g, "SF/SF_{}_{}.queries".format(n, m))


graph_sizes = [10000, 20000, 50000, 100000, 500000, 1000000, 5000000, 10000000]

ER(graph_sizes)
#BA(graph_sizes)
#BP(graph_sizes)
#SF(graph_sizes)

    

