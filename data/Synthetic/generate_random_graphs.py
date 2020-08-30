#!/usr/bin/env python3

import networkx as nx
import random
import os
import sys
# from tqdm import tqdm

# import igraph

# create directories
# for dir_name in ["ER", "BA", "SF", "BP"]:
#     if not os.path.isdir(dir_name):
#         os.mkdir(dir_name)


SEED = 12345
random.seed(SEED)
MAX_EDGE_LENGTH = 1000


def generate_queries_random(g, filename):
    """ Generate queries by picking s and t at random and checking if they satisfy the requirements """

    print("Generating queries for {}".format(filename))
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


def generate_queries(g, filename):
    """ Generates queries by first finding all possible pairs s and t at a certain hop distance and then randomly
    selecting a subset of those """

    print("Generating queries for {}".format(filename))
    q = open(filename, "w")
    n = g.number_of_nodes()
    m = g.number_of_edges() 


    candidates = {i : [] for i in [2, 4, 6, 8]}

    queries_per_category = 100

    # generate queries at certain hop distance
    for s in tqdm(range(n)):
        nodes_at_dist_at_most = nx.single_source_shortest_path_length(g, source=s, cutoff=max(candidates))
        
        for t in nodes_at_dist_at_most:
            if nodes_at_dist_at_most[t] in candidates:
                candidates[nodes_at_dist_at_most[t]].append((s, t))

    queries =  {}
    for h in candidates:
        if len(candidates[h]) < queries_per_category:
            print("ERROR: Not enough candidate pairs for hop distance {}".format(h))
        queries[h] = random.sample(candidates[h], queries_per_category)


    queries[0] = []
    # generate completely random queries
    queries_so_far = 0
    while queries_so_far < queries_per_category:
        s = random.randrange(n)
        t = random.randrange(n)
        if s == t:
            continue
        try:
            hops = nx.shortest_path_length(g, source=s, target=t)
            if len(queries[0]) < queries_per_category and (s,t) not in queries[0]:
                queries[0].append((s, t))
                queries_so_far += 1

        except nx.NetworkXNoPath:
            continue
    for hops in sorted(queries):
        for s, t in queries[hops]:
            q.write("{} {}\n".format(s, t))

    q.close()


# Another option to (faster) generate queries, but here the pairs won't be picked uniformly at random
# First a start vertex is picked uniformly at random and then the end vertex is picked uniformly at random from
# the vertices at a certain hop distance
def generate_queries_skewed(g, filename):
    print("Generating queries for {}".format(filename))
    q = open(filename, "w")
    n = g.number_of_nodes()
    m = g.number_of_edges() 


    queries = {0:set(), 2:set(), 3:set(), 4:set(), 5:set(), 6:set()}

    queries_per_category = 100
    queries_so_far = 0
    i = 0

    # generate queries at certain hop distance
    for h in [2, 3, 4, 5, 6]:
        queries_so_far = 0
        nodes = set()
        backups = set()
        while queries_so_far < queries_per_category and len(nodes) < n:
            s = random.randrange(n)
            nodes.add(s)

            nodes_at_dist_at_most_h = nx.single_source_shortest_path_length(g, source=s, cutoff=h)

            t_options = [k for k in nodes_at_dist_at_most_h if nodes_at_dist_at_most_h[k] == h]

            if len(t_options) == 0: 
                continue

            t = random.choice(t_options)
            if (s,t) not in queries[h]:
                queries[h].add((s,t))
                queries_so_far += 1
            backups.update([(s,v) for v in t_options if v != t])

        while queries_so_far < queries_per_category and backups:
            s, t = random.choice(backups)
            backups.discard((s,t))
            queries[h].add((s,t))

        print("Done for {} hops".format(h))

    # generate completely random queries
    queries_so_far = 0
    node_pairs = set()
    while queries_so_far < queries_per_category and len(node_pairs) < n*(n-1)/2:
        s = random.randrange(n)
        t = random.randrange(n)
        if s == t:
            continue
        node_pairs.add((s,t))
        try:
            hops = nx.shortest_path_length(g, source=s, target=t)
            if not any([((s,t) in queries[h]) for h in queries]):
                queries[0].add((s, t))
                queries_so_far += 1
        except nx.NetworkXNoPath:
            continue
    print("Done for 0 hops")

    q.write("{}\n".format(len(queries)))
    for hops in sorted(queries):
        q.write("{} {}\n".format(hops, len(queries[hops])))
        for s, t in queries[hops]:
            q.write("{} {}\n".format(s, t))

    q.close()




def ER(graph_sizes):
    print("Generating ER graphs")
    # Random graph with n vertices and m = 2*n edges, store in ER folder
    # https://networkx.github.io/documentation/stable/reference/generated/networkx.generators.random_graphs.gnm_random_graph.html
    for n in graph_sizes:
        m = 10*n

        print("Generating ER graph with {} nodes and {} edges".format(n, m))
        g = nx.gnm_random_graph(n, m, directed=True)

        f = open("ER/ER_{}_{}.graph".format(n, m), "w")
        f.write("{} {}\n".format(n, m))
        for u, v in g.edges:
            l = random.randint(1, MAX_EDGE_LENGTH)
            f.write("{} {} {} {}\n".format(u, v, l, random.random()))
            # f.write("{} {} {} {}\n".format(u, v, l, l / MAX_EDGE_LENGTH))
        f.close()

        generate_queries_skewed(g, "ER/ER_{}_{}.queries".format(n, m))

        # test_hop_distance("ER/ER_{}_{}".format(n, m))


def BA(graph_sizes,density):
    print("Generating BA graphs")
    # https://igraph.org/python/doc/igraph.GraphBase-class.html#Barabasi
    for n in graph_sizes:
        BA_m = density

        print("Generating BA graph with {} nodes".format(n))
        igraph_g = igraph.Graph.Barabasi(n, BA_m , directed=True)
        g = nx.DiGraph(igraph_g.get_edgelist())

        #g_ba = nx.barabasi_albert_graph(n, 5)
        #g = nx.DiGraph()
        #g.add_edges_from(g_ba.edges())
        ## also add the reverse edges
        #g.add_edges_from([(v, u) for u,v in g_ba.edges()])

        m = g.number_of_edges()

        print("Graph has {} edges".format(m))

        f = open("BA/BA_{}_{}.graph".format(n, m), "w")
        f.write("{} {}\n".format(n, m))
        for u, v in g.edges:
            l = random.randint(1, MAX_EDGE_LENGTH)
            p = random.random()
            # p = l / MAX_EDGE_LENGTH
            f.write("{} {} {} {}\n".format(u, v, l, p))
        f.close()

        generate_queries_skewed(g, "BA_{}_{}.queries".format(n, m))


def convert_edge_entry(x):
    return [int(x[0]), int(x[1]), int(x[2]), float(x[3])]

def test_hop_distance(graph_name):
    print("Testing {}".format(graph_name))
    f_graph = open(graph_name + ".graph", "r")
    f_queries = open(graph_name + ".queries", "r")
    n, m = [int(x) for x in f_graph.readline().split()]

    g = nx.DiGraph()
    g.add_nodes_from(range(n))
    for _ in range(m):
        u, v, l, p = convert_edge_entry(f_graph.readline().split())
        if u < 0 or u >= n:
            print("u is out of range in edge : {} {} {} {}".format(u, v, l, p))
        if v < 0 or v >= n:
            print("v is out of range in edge : {} {} {} {}".format(u, v, l, p))
        if p < 0 or p > 1:
            print("p is out of range in edge : {} {} {} {}".format(u, v, l, p))
        g.add_edge(u, v)

    # first 100 queries can be random, just check if they are within range
    for hops in [0, 2, 4, 6, 8]:
        for _ in range(100):
            s, t = [int(x) for x in f_queries.readline().split()]
            if s < 0 or s >= n:
                print("s is out of range in query : {} {}".format(s, t))
            if t < 0 or t >= n:
                print("t is out of range in query : {} {}".format(s, t))
            try:
                hops_dist = nx.shortest_path_length(g, source=s, target=t)
                if hops != 0 and hops_dist != hops:
                    print("Nodes in query {} {} are {} hops apart while they should be {} hops apart".format(s, t,
                        hops_dist, hops))

            except nx.NetworkXNoPath:
                print("No path between nodes in query : {} {}".format(s, t))
                continue
    print("Done testing".format(graph_name))



assert len(sys.argv) == 4
assert sys.argv[1] in ['ER','BA']
# graph_sizes = [10000, 20000, 50000, 100000, 500000, 1000000, 5000000, 10000000]

if sys.argv[1] == 'ER':
    ER([sys.argv[2]],sys.argv[3])
else:
    BA([sys.argv[2]],sys.argv[3])
