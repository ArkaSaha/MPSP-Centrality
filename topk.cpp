/*
 * author      : Ruben Brokkelkamp
 * description : Implements the algorithms as described in
 *               Top-K Possible Shortest Path Query over a Large Uncertain Graph
 *               Lei Zou, Peng Peng, and Dongyan Zhao
 *               WISE 2011
 *               
 *               In the comments we refer to this paper by [WISE 2011]
 */


/*
 * Things to consider:
 *                    - What if the probability goes below double precision ~ 10^{-15} ?
 *                          
 */


#include<iostream>
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

using namespace std;
using namespace boost::heap;

#define ull unsigned long long

//mt19937 mersenne{static_cast<mt19937::result_type>(time(nullptr))};}
mt19937 mersenne{static_cast<mt19937::result_type>(12345)};

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

    Graph(int n, int m, vector<vector<Edge> > adj, vector<vector<Edge *>> incoming, vector<Edge *> index2edge) : n(n), m(m), adj(adj), incoming(incoming), index2edge(index2edge) {} ;
    Graph(int n, int m) : Graph(n, m, vector<vector<Edge>>(n, vector<Edge>()), vector<vector<Edge *>>(n, vector<Edge*>()), vector<Edge *>(m)) {} ;
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

    int len(){
        int res = 0;
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
        }
        else{
            cout << "Empty Path";
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
};

set<Edge> p_minus_q(Path p, Path q){
    set<Edge> res = set<Edge>(p.edges.begin(), p.edges.end());
    for(Edge e: q.edges){
        res.erase(e);
    }
    return res;
}

double LB1(vector<Path> &paths, int n){
    // Given a sorted vector of the top n shortest paths
    // computes LB1 of the nth path as in Theorem 1 in [WISE 2011]
    double LB_prob = 1.0;
    for(Edge e: paths[n].edges){
        LB_prob *= e.p;
    }
    
    for(int i=0; i<n; i++){
        set<Edge> edges_Pi = p_minus_q(paths[i], paths[n]);

        double prod_of_probs = 1.0;
        for(Edge e: edges_Pi){
            prod_of_probs *= e.p;
        }
        LB_prob *= (1.0 - prod_of_probs);
    }
    cout << "LB1 of path " << n << " : " << LB_prob << endl;
    return LB_prob;
}

Edge edge_with_most_paths(const map<Edge, set<int>> &edge_to_path){
    Edge max_element;
    size_t most_paths = 0;
    for(auto elt : edge_to_path){
        if(elt.second.size() > most_paths){
            most_paths = elt.second.size();
            max_element = elt.first;
        }
    }

    return max_element;
}

void print_M(const vector<set<Edge>> &M){
    cout << "Size of M : " << M.size() << endl;
    for(auto row: M){
        for(Edge elt: row){
            cout << elt.u << " - " << elt.v << " | ";
        }
        cout << endl;
    }
    cout << "----" << endl;
}

map<Edge, set<int>> build_edge_to_path(const vector<set<Edge>> &M, const set<int> &remaining_indices){
    // Takes an Edge-Path matrix M and turns into into a map giving for each edge in which paths they occur
    // Only taking the edges in remaining into account
    map<Edge, set<int>> edge_to_path;
    for(int i: remaining_indices){
        for(Edge e : M[i]){
            edge_to_path[e].insert(i);
        }
    }

    return edge_to_path;
}

set<Edge> find_set_cover(const vector<set<Edge>> &M){
    // Given an edge path matrix M, computes greedily a set-cover
    
    set<int> remaining_indices = set<int>(); // Will contain the indices of the E(Pi-Pn) which are already covered by our set cover
    for(uint i=0; i<M.size(); i++){
        remaining_indices.insert(i);
    }

    set<Edge> set_cover = set<Edge>();
    while(!remaining_indices.empty()){
        auto edge_to_path = build_edge_to_path(M, remaining_indices);
        Edge max_elt = edge_with_most_paths(edge_to_path); // finds the edge covering most E(Pi-Pn)

        set_cover.insert(max_elt);
        for(int i: edge_to_path[max_elt]){
            remaining_indices.erase(i);
        }
    }
    return set_cover;
}

set<Edge> A_minus_B(set<Edge> A, set<Edge> B){
    // Removes elements in B from A
    for(Edge elt: B){
        A.erase(elt);
    }
    return A;
}

vector<set<Edge>> build_M(const vector<Path> &paths, int n, const set<Edge> &exclude){
    // For computing LB2, computes M as in section 5.1  in [WISE 2011]
    // A 'row' in M contains the edges in E(Pi - Pn)
    vector<set<Edge> > M = vector<set<Edge>>(n, set<Edge>());
    for(int i=0; i<n; i++){
        set<Edge> pi_minus_pn = p_minus_q(paths[i], paths[n]);

        M[i] = A_minus_B(pi_minus_pn, exclude);
    }
    //print_M(M);
    return M;
}

bool M_covers_every_path(const vector<set<Edge>> &M){
    for(auto row: M){
        if(row.size() == 0) return false;
    }
    return true;
}


vector<set<Edge>> find_pairwise_independend_set_covers(const vector<Path> &paths, int n){
    // For computing LB2, finds pairwise independent set covers as desired by Theorem 2 in [WISE 2011]
    set<Edge> exclude = set<Edge>();
    vector<set<Edge>> M = build_M(paths, n, exclude);
    vector<set<Edge>> set_covers = vector<set<Edge>>();

    while(M_covers_every_path(M)){ // while there are set covers to be found
        set<Edge> cur_set_cover = find_set_cover(M);
        set_covers.push_back(cur_set_cover);

        // exclude the edges in the set cover, so that we don't pick them again
        exclude.insert(cur_set_cover.begin(), cur_set_cover.end()); 
        M = build_M(paths, n, exclude);
    }

    return set_covers;

}


double LB2(const vector<Path> &paths, int n){
    // Given a sorted vector of the top n shortest paths
    // computes LB2 of the nth path as in Theorem 2 in [WISE 2011]
    cout << "Computing LB2  of path " << n << " : ";
    paths[n].print(); cout << endl;

    if(n == 0) return 0.0;

    vector<set<Edge>> set_covers = find_pairwise_independend_set_covers(paths, n);

    double prob_of_path_pn = 1.0;
    for(Edge e: paths[n].edges){
        prob_of_path_pn *= e.p;
    }

    double PUES = 0;
    for(uint j=0; j<set_covers.size(); j++){
        double prob_none_of_Sj = 1.0;
        for(Edge e: set_covers[j]){
            prob_none_of_Sj *= (1 - e.p);
        }
        PUES = PUES + (1-PUES) * prob_none_of_Sj;

    }
    return prob_of_path_pn * PUES;
}


void update_LB(vector<Path> &paths, int n){
    // Given a sorted vector of the top n shortest paths, computes and updates the LB value of the nth path
    double LB1_prob = LB1(paths, n);
    double LB2_prob = LB2(paths, n);
    cout << "LB1 : " << LB1_prob << endl;
    cout << "LB2 : " << LB2_prob << endl;
    if(LB2_prob > LB1_prob){
        cout << "LB2 was better!" << endl;
    }
    paths[n].LB = max(LB1_prob, LB2_prob);
}



double UB(const vector<Path> &paths, int n){
    // Given a sorted vector of the top n shortest paths
    // computes UB of the nth path as in Theorem 3 in [WISE 2011]
    double UB_prob = 1.0;
    for(int i=0; i<n; i++){
        UB_prob -= paths[i].LB;
    }
    cout << "UB " << n << " : " << UB_prob << endl;
    return UB_prob;
}

void update_UB(vector<Path> &paths, int n){
    // Given a sorted vector of the top n shortest paths, computes and updates the UB value of the nth path
    paths[n].UB = UB(paths, n);
}

Path dijkstra(const Graph &g, int s, int t){
    // Computes the shortest path between nodes s and t in graph g using Dijkstra's shortest path algorithm
    
	struct node
	{
        int u;
        ull dist;
		node(const int& u, ull dist) : u(u), dist(dist) {}
	};
	struct compare_node
	{
		bool operator()(const node& n1, const node& n2) const
		{
			return n1.dist > n2.dist;
		}
	};
	fibonacci_heap< node, compare<compare_node> > q = fibonacci_heap< node, compare<compare_node> >();

    vector<ull> D = vector<ull>(g.n, ULLONG_MAX); // initialize all distances to infinity
    vector<Edge> prev = vector<Edge>(g.n); // vector to remember the edge on the shortest path 
    typedef fibonacci_heap< node, compare<compare_node> >::handle_type handle_t;
    vector<handle_t> handles = vector<handle_t>(g.n); // handles for priority queue
    vector<bool> visited = vector<bool>(g.n, false); // vector to remember if we have already created a handle

    // initialize node s
    handles[s] = q.push(node(s, 0));
    D[s] = 0;
    visited[s] = true;

    while(!q.empty()){
        node curnode = q.top(); // Take the closest unvisited node from the priority queue
        q.pop();
        if(curnode.u == t){ // if we have arrived at t this means we found the shortest path
            break;
        }
        for(Edge e: g.adj[curnode.u]){
            if(!e.available) continue;
            if(!visited[e.v]){
                // if we have never seen this node, create a handle and set distance and prev
                D[e.v] = D[e.u] + e.l;
                handles[e.v] = q.push(node(e.v, D[e.v]));
                prev[e.v] = e;
                visited[e.v] = true;
            }
            else if(D[e.u] + e.l < D[e.v]){
                // we have already seen this node, but found a shorter path, update handle so that it moves 
                // up in the priority queue and update distance and prev
                D[e.v] = D[e.u] + e.l;
                q.update(handles[e.v], node(e.v, D[e.v]));
                prev[e.v] = e;
            }
        }
    }

    if(D[t] == ULLONG_MAX){
        cout << "There is no path between s = " << s << " and t = " << t << endl;
        // disconnected
        return Path({});
    }
    else{
        // start at t and reverse walk over the shortest path to recreate it
        vector<Edge> edges = {prev[t]};
        int startnode = prev[t].u;
        while(startnode != s){
            edges.push_back(prev[startnode]);
            startnode = prev[startnode].u;
        }
        reverse(edges.begin(), edges.end());
        return Path({edges});
    }
}


bool subpath_of(const Path &p1, const Path &p2){
    // Test if p1 is the same as the beginning of path p2
    for(uint i=0; i<p1.edges.size(); i++){
        if(i > p2.edges.size()) return false;
        if(p1.edges[i] != p2.edges[i]) return false;
    }
    return true;
}


vector<Path> classic_yen(Graph &g, int s, int t, int K){
    // Computes the top k shortest paths using Yen's algorithm
    cout << "\n\n\nYEN\n\n\n" << endl;
    Path p1 = dijkstra(g, s, t);

    vector<Path> A = {p1};
    vector<Path> B = vector<Path>();
    for(int k=1; k<K; k++){
        cout << "\nk : " << k << endl << endl;
        for(uint i=0; i<A.back().edges.size(); i++){

            // we are going find a new path, diverting from the old path from the spurNode
            int spurNode = A[k-1].edges[i].u;
            Path rootPath = Path({A[k-1].edges.begin(), A[k-1].edges.begin() + i});

            // edges we want to avoid
            set<int> edges_to_delete = set<int>();
            for(uint j=0; j<A.size();j++){
                if(subpath_of(rootPath, A[j])){
                    edges_to_delete.insert(A[j].edges[i].index);
                }
            }

            // nodes we want to avoid
            set<int> nodes_to_delete = {};
            for(Edge e: rootPath.edges){
                if(e.u != spurNode){
                    nodes_to_delete.insert(e.u);
                }
            }

            // create a new graph without all the unwanted nodes and edges
            // TODO: do not create a copy but use the old graph (maybe by adding a boolean to edges whether they are in
            // use or not)
            Graph g2;
            g2.n = g.n; g2.m = 0;
            g2.adj = vector<vector<Edge>>(g2.n, vector<Edge>());
            for(int j=0; j<g2.n; j++){
                if(nodes_to_delete.find(j) == nodes_to_delete.end()){
                    for(Edge e: g.adj[j]){
                        if(edges_to_delete.find(e.index) == edges_to_delete.end()){
                            g2.m++;
                            g2.adj[j].push_back(e);
                        }
                    }
                }
            }

            // compute the shortest path in the new graph from the spurnode to the terminal node
            Path spurPath = dijkstra(g2, spurNode, t);

            if(spurPath.edges.size() > 0){
                // if we found a path, concatenate it to the rootpath and add it to the set B
                rootPath.edges.insert(rootPath.edges.end(), spurPath.edges.begin(), spurPath.edges.end());
                if(find(B.begin(), B.end(), rootPath) == B.end()){
                    // only if this path was not in B yet
                    B.push_back(rootPath);
                }
            }
        }

        if(B.size() == 0){
            // there are not enough paths
            break;
        }

        // Get shortest path in B
        int indexmin = 0;
        int lenmin = B[0].len();
        for(uint j=1; j<B.size(); j++){
            int curlen = B[j].len();
            if(curlen < lenmin){
                lenmin = curlen;
                indexmin = j;
            }
        }

        // add the shortest path from B to A as the kth shortst path
        A.push_back(B[indexmin]);
        B.erase(B.begin() + indexmin);
    }

    cout << "Done with yen" << endl;
    for(uint i=0; i<A.size(); i++){
        cout << "Path " << i << " : ";
        A[i].print(); cout << endl;
    }
    return A;

}

double kth_largest(const set<double> &elts, int k){
    // find the kth largest element in a set of doubles
    auto elt = elts.rbegin();
    for(int i=1; i<min((int)elts.size(), k); i++){
        elt++;
    }
    return *elt;
}

vector<Path> yen(Graph &g, int s, int t, int k){
    // Computes the top k_prime shortest paths using Yen's algorithm where k_prime depends on a stop condition
    // instead of being a fixed parameter
    //
    // The stop condition is that the 
    // UB(prob(k_prime'th shortest path)) < k'th largest LB(prob(path Pn is shortest path))
    // where Pn are the top n shortest paths
    cout << "\n\n\nYEN\n\n\n" << endl;
    Path p1 = dijkstra(g, s, t);

    vector<Path> A = {p1};
    vector<Path> B = vector<Path>();

    update_LB(A, 0);
    update_UB(A, 0);    

    set<double> LBs = {A[0].LB};

    int k_prime = 0;
    while(k_prime < k || A[k_prime].UB >= kth_largest(LBs, k)){
        cout << "\nk_prime : " << k_prime << endl << endl;
        for(uint i=0; i<A.back().edges.size(); i++){
            // we are going find a new path, diverting from the old path from the spurNode
            int spurNode = A.back().edges[i].u;
            Path rootPath = Path({A.back().edges.begin(), A.back().edges.begin() + i});


            // nodes we want to avoid
            vector<int> nodes_to_delete = vector<int>();
            for(Edge e: rootPath.edges){
                if(e.u != spurNode){
                    nodes_to_delete.push_back(e.u);

                    for(Edge *e2:  g.incoming[e.u]){
                        e2->available = false;
                    }
                    for(Edge e2: g.adj[e.u]){
                        e2.available = false;
                    }
                }
            }

            // edges we want to avoid
            vector<int> edges_to_delete = vector<int>();
            for(uint j=0; j<A.size();j++){
                if(subpath_of(rootPath, A[j])){
                    edges_to_delete.push_back(A[j].edges[i].index);

                    g.index2edge[A[j].edges[i].index]->available = false;
                }
            }

            // compute the shortest path in the new graph from the spurnode to the terminal node
            Path spurPath = dijkstra(g, spurNode, t);

            if(spurPath.edges.size() > 0){
                // if we found a path, concatenate it to the rootpath and add it to the set B
                rootPath.edges.insert(rootPath.edges.end(), spurPath.edges.begin(), spurPath.edges.end());
                if(find(B.begin(), B.end(), rootPath) == B.end()){
                    // only if this path was not in B yet
                    B.push_back(rootPath);
                }
            }

            // restore graph
            for(int elt: nodes_to_delete){
                for(Edge *e : g.incoming[elt]){
                    e->available = true;
                }
                for(Edge e : g.adj[elt]){
                    e.available = true;
                }
            }
            for(int index: edges_to_delete){
                g.index2edge[index]->available = true;
            }


        }

        if(B.size() == 0){
            // there are not enough paths
            break;
        }

        // Get shortest path in B
        int indexmin = 0;
        int lenmin = B[0].len();
        for(uint j=1; j<B.size(); j++){
            int curlen = B[j].len();
            if(curlen < lenmin){
                lenmin = curlen;
                indexmin = j;
            }
        }

        // add the shortest path from B to A as the kth shortst path
        A.push_back(B[indexmin]);
        B.erase(B.begin() + indexmin);

        k_prime++;

        cout << "New path : "; A[k_prime].print(); cout << endl;
        // compute the LB and UB for the new path
        update_LB(A, k_prime);
        update_UB(A, k_prime);
        LBs.insert(A[k_prime].LB);
        cout << "Current UB     : " << A[k_prime].UB << endl;
        cout << "kth largest LB : " << kth_largest(LBs, k) << endl;
    }
    if(A[k_prime].UB < kth_largest(LBs, k)){
        A.pop_back();
    }
    cout << "Done with yen, found " << A.size() << " paths" << endl;
    for(uint i=0; i<A.size(); i++){
        cout << "Path " << i << " : ";
        A[i].print(); cout << endl;
    }
    return A;
}


double Luby_Karp(const vector<Path> &paths, int n, ull N){
    // Implementation of Luby Karp as in section 7 of [WISE 2011]
    
    if(n == 0) 
        return paths[0].probability();

    ull cnt = 0;
    vector<double> PEPi_Pn = vector<double>(n, 1.0);
    vector<set<Edge>> Pi_Pn = vector<set<Edge>>(n, set<Edge>());
    double S = 0.0;
    for(int i=0; i < n; i++){
        Pi_Pn[i] = p_minus_q(paths[i], paths[n]); // the edges in Pi but not in Pn;
        for(Edge e : Pi_Pn[i]){
            PEPi_Pn[i] *= e.p;
        }
        S += PEPi_Pn[i];
    }
    discrete_distribution<int> dist(PEPi_Pn.begin(), PEPi_Pn.end());
    uniform_real_distribution<double> coin(0.0, 1.0); 
    for(ull i=0; i<N; i++){
        int index = dist(mersenne); // pick a random index according to distribution dist
        map<Edge, bool> truth_assignment = map<Edge, bool>(); // we will generate the truth assignment on the fly

        for(Edge e : Pi_Pn[index]){
            truth_assignment[e] = true;
        }

        bool count_this_iteration = true;
        for(int j = 0; j < index; j++){
            bool path_j_exists = true;
            for(Edge e: Pi_Pn[j]){
                if(truth_assignment.count(e) == 0){
                    // we have not encountered this edge, flip a coin to decide whether it exists or not
                    double coin_flip = coin(mersenne);
                    truth_assignment[e] = (coin_flip < e.p) ? true : false;
                }

                if(!truth_assignment[e]){
                    path_j_exists = false;
                    break;
                }
            }
            if(path_j_exists){
                count_this_iteration = false;
                break;
            }
        } 
        if(count_this_iteration){
            cnt++;
        }
    }

    /*
    cout << "cnt : " << cnt << endl;
    cout << "N   : " << N << endl;
    cout << "S   : " << S << endl;
    */

    double p_tilde = ((double)cnt)/((double)N) * S;

    return (1-p_tilde) * paths[n].probability();
}

vector<pair<Path,double>> topk(Graph &g, int s, int t, int k){
    /* Computes the topk most probably shortest paths. Consists of two steps
     * 1) Use Yen's algorithm with a modified stopping rule to find a set of candidates.
     * 2) Use Luby-Karp Mote Carlo sampling to compute probabilities 
     */

    vector<Path> candidates = yen(g, s, t, k);

    cout << "Found " << candidates.size() << " candidates" << endl;

    // for every path, use Luby-Karp to estimate the probability of it being the shortest path
    vector<pair<double, int>> LK_probabilities = vector<pair<double, int>>(candidates.size()); 
    for(uint i=0; i<candidates.size(); i++){
        double LK = Luby_Karp(candidates, i, 1000000);
        LK_probabilities[i] = {LK, i};
    }

    // sort the paths from high to low probability
    sort(LK_probabilities.begin(), LK_probabilities.end(), 
            [](const auto &lhs, const auto &rhs){
                return lhs.first > rhs.first;
            }
        );

    // take the top k
    vector<pair<Path, double>> topk_mpsp= vector<pair<Path, double>>();
    for(int i=0; i<k; i++){
        int index = LK_probabilities[i].second;
        topk_mpsp.push_back({candidates[index], LK_probabilities[i].first});
    }
    return topk_mpsp;
}

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

        G.adj[u].push_back(Edge({u, v, l, p, i}));
    }

    G.update_incoming_index2edge();

    return G;
}


bool testLB1_UB(){
    vector<Path> paths = vector<Path>();
    paths.push_back({{{0, 1, 3, 0.3, 0}, {1, 3, 1, 0.4, 3}}});
    paths.push_back({{{0, 2, 1, 0.1, 1}, {2, 1, 1, 0.2, 2}, {1, 3, 1, 0.4, 3}}}); 
    paths.push_back({{{0, 2, 1, 0.1, 1}, {2, 3, 2, 0.9, 5}}}); 

    double res = LB1(paths, 1);
    double true_ans =(0.1 * 0.2 * 0.4 * 0.7) ;
    if(res != true_ans){
        cerr << "ERROR, failed test #1" << endl;
        cerr << "Output    : " << res << endl;
        cerr << "Should be : " << true_ans << endl;
        return false;
    }
    else{
        cerr << "Passed test #1" << endl;
    }

    res = LB1(paths, 2);
    true_ans = 0.1*0.9 * ((1-(0.2*0.4)) * (1-(0.3*0.4))) ;
    if(abs(res - true_ans) > 0.000000001){
        cerr << "ERROR, failed test #2" << endl;
        cerr << "Output    : " << res << endl;
        cerr << "Should be : " << true_ans << endl;
        return false;
    }
    else{
        cerr << "Passed test #2" << endl;
    }

    update_LB(paths, 0);
    update_LB(paths, 1);
    update_LB(paths, 2);
    update_UB(paths, 0);
    update_UB(paths, 1);
    update_UB(paths, 2);

    cout << paths[0].LB << endl;
    cout << paths[0].UB << endl;
    cout << paths[1].LB << endl;
    cout << paths[1].UB << endl;
    cout << paths[2].LB << endl;
    cout << paths[2].UB << endl;

    // TODO: check that these values are correct


    return true;
}


int main(){
    Graph G = read_from_stdin();

    Path p = dijkstra(G, 0, 3);

    for(int i=0; i<G.m; i++){
        cout << "edge " << i << " : " << G.index2edge[i]->u << " - " << G.index2edge[i]->v << endl;
    }

    //testLB1_UB();

    vector<Path> paths = yen(G, 0, 5, 4);

    cout << "checkpoint1" << endl;

    for(uint i=0; i<paths.size(); i++){
        double prob = Luby_Karp(paths, i, 100000);
        cout << i << "   : " << prob << endl;
    }


    cout << "random stuff" << endl;


	clock_t t1 = clock();
    // some random graph
    Graph g_random(200000, 500000);

    int MAX_EDGE_LENGTH = 10;
    uniform_int_distribution<int> new_vertex(0, g_random.n-1); 
    uniform_real_distribution<double> new_probability(0.9, 1.0); 
    uniform_int_distribution<int> new_length(1, MAX_EDGE_LENGTH); 

    set<pair<int, int>> generated_edges = set<pair<int, int>>();

    for(int i=0; i<g_random.m; i++){
        int u = new_vertex(mersenne);
        int v = new_vertex(mersenne);

        if(u == v || generated_edges.find({u, v}) != generated_edges.end()){
            // skip this edge
            i--;
            continue;
        } 

        g_random.adj[u].push_back({u, v, new_length(mersenne), new_probability(mersenne), i});
    }
    g_random.update_incoming_index2edge();

    clock_t t2 = clock();
    cout << "Done generating in " << 1000.0 * (t2-t1) / CLOCKS_PER_SEC << " miliseconds" << endl;
    
    //vector<Path> paths11 = classic_yen(&g_random, 0, g_random.n - 1, 10);
    //vector<Path> paths11 = yen(g_random, 0, g_random.n - 1, 3);
    vector<pair<Path,double>> res = topk(g_random, 0, g_random.n-1, 1);

    clock_t t3 = clock();
    cout << "Done in " << 1000.0 * (t3-t2) / CLOCKS_PER_SEC << " miliseconds" << endl;

    for(uint i= 0; i<res.size(); i++){
        cout << "Path " << i+1 << " : probability : " << res[i].second << endl;
        cout <<             "       : length      : " << res[i].first.len() << endl;
        cout <<             "       : edges       : "; res[i].first.print(); cout << endl;
    }

	return 0;
}



