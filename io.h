#pragma once
#include "topk.h" 

#include <fstream>
#include <filesystem>

Graph read_from_stdin(){
    // Reads a graph from stdin, file containing graph should look the following
    // First a line with two space separated integers n and m
    // Here n is the number of nodes and m is the number of edges
    // Then follow m lines with of the form 'u v l p' 
    // u, v, l are integers and p is a double
    // representing an edge from node u to node v of length l with probability p
    int n, m;
    cin >> n >> m;
    Graph G = Graph({n, m});
    int u, v, l; double p;
    for(int i=0; i<m; i++){
        cin >> u >> v >> l >> p;
        assert(u >= 0 and u < n);
        assert(v >= 0 and v < n);
        assert(p >= 0.0 and p <= 1.0);

        G.adj[u].push_back(Edge({u, v, l, p, i}));
    }

    G.update_incoming_index2edge();

    return G;
}

Graph read_graph_from_file(string filename)
{
	ifstream graph_in;
	graph_in.open(filename);
    int n, m;
    graph_in >> n >> m;
    Graph G = Graph({n, m});
    int u, v, l; double p;
    for(int i=0; i<m; i++){
        graph_in >> u >> v >> l >> p;
        assert(u >= 0 and u < n);
        assert(v >= 0 and v < n);
        assert(p >= 0.0 and p <= 1.0);

        G.adj[u].push_back(Edge({u, v, l, p, i}));
    }

    G.update_incoming_index2edge();
    graph_in.close();
    return G;
}
