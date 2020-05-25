#!/usr/bin/env python3

import networkx as nx
import random
import os
from tqdm import tqdm
import sys

import igraph


SEED = 12345
random.seed(SEED)
MAX_EDGE_LENGTH = 1000


# Another option to (faster) generate queries, but here the pairs won't be picked uniformly at random
# First a start vertex is picked uniformly at random and then the end vertex is picked uniformly at random from
# the vertices at a certain hop distance
def generate_queries_skewed(g, filename, H = [0, 2, 3, 4, 5, 6]):
    q = open(filename, "w")
    n = g.number_of_nodes()
    m = g.number_of_edges() 


    queries = {h : set() for h in H}

    queries_per_category = 100
    queries_so_far = 0
    i = 0

    # generate queries at certain hop distance
    for h in tqdm([x for x in H if x != 0]):
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

    if 0 in H:
        # generate completely random queries
        queries_so_far = 0
        while queries_so_far < queries_per_category:
            s = random.randrange(n)
            t = random.randrange(n)
            if s == t:
                continue
            try:
                hops = nx.shortest_path_length(g, source=s, target=t)
                if (s,t) not in queries[0]:
                    queries[0].add((s, t))
                    queries_so_far += 1
            except nx.NetworkXNoPath:
                continue

    q.write("{}\n".format(len(queries)))
    for hops in queries:
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
            f.write("{} {} {} {}\n".format(u, v, random.randint(1, MAX_EDGE_LENGTH), random.random()))
        f.close()

        generate_queries_skewed(g, "ER/ER_{}_{}.queries".format(n, m))

        test_hop_distance("ER/ER_{}_{}".format(n, m))

def ER_custom(n, m):
    print("Generating ER graphs")
    # Random graph with n vertices and m = 2*n edges, store in ER folder
    # https://networkx.github.io/documentation/stable/reference/generated/networkx.generators.random_graphs.gnm_random_graph.html
    print("Generating ER graph with {} nodes and {} edges".format(n, m))
    g = nx.gnm_random_graph(n, m, directed=True)

    f = open("ER/ER_{}_{}.graph".format(n, m), "w")
    f.write("{} {}\n".format(n, m))
    for u, v in g.edges:
        f.write("{} {} {} {}\n".format(u, v, random.randint(1, MAX_EDGE_LENGTH), random.random()))
    f.close()

    generate_queries_skewed(g, "ER/ER_{}_{}.queries".format(n, m))

    test_hop_distance("ER/ER_{}_{}".format(n, m))


def BA(graph_sizes):
    print("Generating BA graphs")
    # https://igraph.org/python/doc/igraph.GraphBase-class.html#Barabasi
    for n in graph_sizes:
        BA_m = 6

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
            f.write("{} {} {} {}\n".format(u, v, l, p))
        f.close()

        generate_queries(g, "BA/BA_{}_{}.queries".format(n, m))

        test_hop_distance("BA/BA_{}_{}".format(n, m))

def BP(graph_sizes):
    print("Generating BP graphs")
    for n in graph_sizes:
        m = 2*n

        g = nx.DiGraph()
        g.add_nodes_from(range(n))

        nodes_partition_1 = n//5

        edges = 0
        while edges < m:
            u = random.randrange(nodes_partition_1)
            v = random.randrange(nodes_partition_1, n)

            coin = random.random()
            if coin < 0.5 and not g.has_edge(u, v):
                g.add_edge(u, v)
                edges += 1
            elif not g.has_edge(v,u):
                g.add_edge(v, u)
                edges += 1

        print("Generated BP graph with {} nodes and {} edges".format(n, m))

        f = open("BP/BP_{}_{}.graph".format(n, m), "w")
        f.write("{} {}\n".format(n, m))
        for u, v in g.edges:
            l = random.randint(1, MAX_EDGE_LENGTH)
            p = random.random()
            f.write("{} {} {} {}\n".format(u, v, l, p))
        f.close()

        generate_queries_random(g, "BP/BP_{}_{}.queries".format(n, m))

        test_hop_distance("BP/BP_{}_{}".format(n, m))

#def SF(graph_sizes):
#    print("Generating SF graphs")
#    # https://igraph.org/python/doc/igraph.GraphBase-class.html
#    for n in graph_sizes:
#
#        print("Generating SF graph with {} nodes".format(n))
#        g = nx.scale_free_graph(n)
#        m = g.number_of_edges()
#
#        print("Graph has {} edges".format(m))
#
#        f = open("SF/SF_{}_{}.graph".format(n, m), "w")
#        f.write("{} {}\n".format(n, m))
#        for u, v, _ in g.edges:
#            f.write("{} {} {} {}\n".format(u, v, random.randint(1, MAX_EDGE_LENGTH), random.random()))
#        f.close()
#
#        generate_queries(g, "SF/SF_{}_{}.queries".format(n, m))

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

def read_graph(graph_file):
    f_graph = open(graph_file, "r")
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
    f_graph.close()

    return g


def add_queries(graph_name):
    print("Adding queries {}".format(graph_name))
    f_queries = open(graph_name + ".queries", "r")

    g = read_graph(graph_name + ".graph")
    n, m = g.number_of_nodes(), g.number_of_edge()

    queries = {0:set(), 2:set(), 4:set(), 6:set()}
    # first 100 queries can be random, just check if they are within range
    for hops in [2, 4, 6]:
        N = int(f_queries.readline().strip())
        for _ in range(N):
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
                queries[hops].add((s,t))

            except nx.NetworkXNoPath:
                print("No path between nodes in query : {} {}".format(s, t))
                continue

    f_queries.close()

    # generate completely random queries
    queries_per_category = 100
    queries_so_far = 0
    while queries_so_far < queries_per_category:
        s = random.randrange(n)
        t = random.randrange(n)
        if s == t:
            continue
        try:
            hops = nx.shortest_path_length(g, source=s, target=t)
            if not any((s,t) in queries[hops] for hops in queries):
                queries[0].add((s, t))
                queries_so_far += 1
        except nx.NetworkXNoPath:
            continue
    q = open(graph_name + ".queries", "w")
    for hops in queries:
        q.write("{}\n".format(len(queries[hops])))
        for s, t in queries[hops]:
            q.write("{} {}\n".format(s, t))

    q.close()
    print("Done adding queries".format(graph_name))

g = read_graph("{}.graph".format(sys.argv[1]))
generate_queries_skewed(g, "{}.queries".format(sys.argv[1]))


