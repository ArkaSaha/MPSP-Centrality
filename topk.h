#pragma once
#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<climits>
#include<list>
#include<map>
#include<vector>
#include<set>
#include<utility>
#include<random>
#include<ctime>

#include<boost/heap/fibonacci_heap.hpp>

#include "statistics.h"


using namespace std;
using namespace boost::heap;


#define ull unsigned long long


struct Edge
{
    int u, v;
    int l;
    double p;
    int index;
    bool available = true;

    bool operator<(const Edge& rhs) const{
        return this->index < rhs.index;
    }

    bool operator==(const Edge& rhs) const{
        return this->index == rhs.index;
    }
    bool operator!=(const Edge& rhs) const{
        return this->index != rhs.index;
    }
};

struct Graph
{
    /*
     * NB: When adj changes the pointers in incoming and index2edge are invalidated and should be recalculated
     */
	int n, m;
    vector< vector<Edge> > adj;
    vector< vector<Edge *> > incoming;
    vector<Edge *> index2edge;

    Graph(int n, int m, vector<vector<Edge> > adj, vector<vector<Edge *> > incoming, vector<Edge *> index2edge) : n(n), m(m), adj(adj), incoming(incoming), index2edge(index2edge) {} ;
    Graph(int n, int m) : Graph(n, m, vector<vector<Edge> >(n, vector<Edge>()), vector<vector<Edge *> >(n, vector<Edge*>()), vector<Edge *>(m)) {} ;
    Graph(){};

    void update_incoming_index2edge(){
        for(int i=0; i<this->n;i++){
            this->incoming[i] = vector<Edge *>();
            for(uint j=0; j<this->adj[i].size(); j++){
                this->incoming[this->adj[i][j].v].push_back(&this->adj[i][j]);
                this->index2edge[this->adj[i][j].index] = &this->adj[i][j];
            } 
        }
    }
};


struct Path
{
    vector<Edge> edges;
    double LB, UB; 
    Path(vector<Edge> edges) : edges(edges){
        LB = 0.0;
        UB = 1.0;
    }

    ull len(){
        ull res = 0;
        for(Edge e: edges){
            res += e.l;
        }
        return res;
    }

    void print() const{
        if(edges.size() > 0){
            cout << edges[0].u;
            for(Edge e: edges){
                cout << " -> " << e.v;
            }
            cout << endl;
        }
        else{
            cout << "Empty Path" << endl;
        }
    }

    double probability() const{
        if(edges.size() == 0) return 0.0;
        double prob = 1.0;
        for(Edge e: edges){
            prob *= e.p;
        }
        return prob;
    }

    ull len() const{
        if(edges.size() == 0) return 0;
        ull len = 0;
        for(Edge e: edges){
            len += e.l;
        }
        return len;
    }

    bool operator==(const Path& rhs) const{
        // Paths are equal if all the edges are the same
        if(this->edges.size() != rhs.edges.size()) return false;
        for(uint i=0; i<this->edges.size(); i++){
            if(this->edges[i] != rhs.edges[i]) return false;
        }
        return true; 
    }

    bool operator!=(const Path& rhs) const{
        return !((*this) == rhs); 
    }

    bool subpath_of(const Path & p2) const{
        if(this->edges.size() > p2.edges.size()) return false;
        for(uint i = 0; i<this->edges.size(); i++){
            if(this->edges[i] != p2.edges[i]) return false;
        }
        return true;
    }

};


// Function declaration
Path dijkstra(const Graph & g, int s, int t);
vector<Path> yen(Graph &g, int s, int t, int k, Statistics & stats, ostream & ofs, double THRESHOLD_NR_OF_SECONDS);
vector<Path> yen(Graph &g, Path p, double THRESHOLD_NR_OF_SECONDS);

pair< vector< pair<Path,double> >, vector<int> > topk(Graph &g, int s, int t, int k, Statistics & stats, ostream & ofs, double THRESHOLD_NR_OF_SECONDS);

double exact_probability(Graph &g, Path p);

double Luby_Karp(const vector<Path> & paths, int n, ull N);

double Luby_Karp(Graph &g, Path p, ull N);





