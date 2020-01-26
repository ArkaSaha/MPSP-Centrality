/*
 * author      : Ruben Brokkelkamp
 * description : Implements the algorithms as described in
 *               Top-K Possible Shortest Path Query over a Large Uncertain Graph
 *               Lei Zou, Peng Peng, and Dongyan Zhao
 *               WISE 2011
 *               
 *               In the comments we refer to this paper by [WISE 2011]
 */


#include<iostream>
#include<cmath>
#include<cstdlib>
#include<climits>
#include<list>
#include<map>
#include<vector>
#include<set>
#include<boost/heap/fibonacci_heap.hpp>

using namespace std;
using namespace boost::heap;

#define ull unsigned long long

struct Edge
{
    int u, v;
    int l;
    double p;
    int index;

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
	int n, m;
    vector< vector<Edge> > adj;
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
        for(auto e: edges){
            res += e.l;
        }
        return res;
    }

    void print(){
        if(edges.size() > 0){
            cout << edges[0].u;
            for(auto e: edges){
                cout << " -> " << e.v;
            }
            cout << endl;
        }
        else{
            cout << "Empty Path" << endl;
        }
    }

    bool operator==(const Path& rhs) const{
        // Paths are equal if all the edges are the same
        if(this->edges.size() != rhs.edges.size()) return false;
        for(int i=0; i<this->edges.size(); i++){
            if(this->edges[i] != rhs.edges[i]) return false;
        }
        return true; 
    }

    bool operator!=(const Path& rhs) const{
        return !((*this) == rhs); 
    }
};


double LB1(vector<Path> & paths, int n){
    // Given a sorted vector of the top n shortest paths
    // computes LB1 of the nth path as in Theorem 1 in [WISE 2011]
    double LB_prob = 1.0;
    for(auto e: paths[n].edges){
        LB_prob *= e.p;
    }
    
    for(int i=0; i<n; i++){
        set<Edge> edges_Pi = set<Edge>(paths[i].edges.begin(), paths[i].edges.end());
        for(auto e: paths[n].edges){
            edges_Pi.erase(e);
        }

        double prod_of_probs = 1.0;
        for(auto e: edges_Pi){
            prod_of_probs *= e.p;
        }
        LB_prob *= (1.0 - prod_of_probs);
    }
    cout << "LB1 of path " << n << " : " << LB_prob << endl;
    return LB_prob;
}


void update_LB(vector<Path> & paths, int n){
    // Given a sorted vector of the top n shortest paths, computes and updates the LB value of the nth path
    double LB_prob = LB1(paths, n);
    paths[n].LB = max(paths[n].LB, LB_prob);
}


double UB(vector<Path> &paths, int n){
    // Given a sorted vector of the top n shortest paths
    // computes UB of the nth path as in Theorem 3 in [WISE 2011]
    double UB_prob = 1.0;
    for(int i=0; i<n; i++){
        UB_prob -= paths[i].LB;
    }
    cout << "UB of path " << n << " : " << UB_prob << endl;
    return UB_prob;
}

void update_UB(vector<Path> &paths, int n){
    // Given a sorted vector of the top n shortest paths, computes and updates the UB value of the nth path
    double UB_prob = UB(paths, n);
    paths[n].UB = min(paths[n].UB, UB_prob);
}

Path dijkstra(Graph* g, int s, int t){
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

    vector<ull> D = vector<ull>(g->n, ULLONG_MAX); // initialize all distances to infinity
    vector<Edge> prev = vector<Edge>(g->n); // vector to remember the edge on the shortest path 
    typedef fibonacci_heap< node, compare<compare_node> >::handle_type handle_t;
    vector<handle_t> handles = vector<handle_t>(g->n); // handles for priority queue
    vector<bool> visited = vector<bool>(g->n, false); // vector to remember if we have already created a handle

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
        for(Edge e: g->adj[curnode.u]){
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


bool subpath_of(Path p1, Path p2){
    // Test if p1 is the same as the beginning of path p2
    for(int i=0; i<p1.edges.size(); i++){
        if(i > p2.edges.size()) return false;
        if(p1.edges[i] != p2.edges[i]) return false;
    }
    return true;
}


vector<Path> classic_yen(Graph *g, int s, int t, int K){
    // Computes the top k shortest paths using Yen's algorithm
    cout << "\n\n\nYEN\n\n\n" << endl;
    Path p1 = dijkstra(g, s, t);

    vector<Path> A = {p1};
    vector<Path> B = vector<Path>();
    for(int k=1; k<K; k++){
        cout << "\nk : " << k << endl << endl;
        for(int i=0; i<A.back().edges.size(); i++){

            // we are going find a new path, diverting from the old path from the spurNode
            int spurNode = A[k-1].edges[i].u;
            Path rootPath = Path({A[k-1].edges.begin(), A[k-1].edges.begin() + i});

            // edges we want to avoid
            set<int> edges_to_delete = set<int>();
            for(int j=0; j<A.size();j++){
                if(subpath_of(rootPath, A[j])){
                    edges_to_delete.insert(A[j].edges[i].index);
                }
            }

            // nodes we want to avoid
            set<int> nodes_to_delete = {};
            for(auto e: rootPath.edges){
                if(e.u != spurNode){
                    nodes_to_delete.insert(e.u);
                }
            }

            // create a new graph without all the unwanted nodes and edges
            Graph g2;
            g2.n = g->n; g2.m = 0;
            g2.adj = vector<vector<Edge>>(g2.n, vector<Edge>());
            for(int j=0; j<g2.n; j++){
                if(nodes_to_delete.find(j) == nodes_to_delete.end()){
                    for(auto e: g->adj[j]){
                        if(edges_to_delete.find(e.index) == edges_to_delete.end()){
                            g2.m++;
                            g2.adj[j].push_back(e);
                        }
                    }
                }
            }

            // compute the shortest path in the new graph from the spurnode to the terminal node
            Path spurPath = dijkstra(&g2, spurNode, t);

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
        for(int j=1; j<B.size(); j++){
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
    for(int i=0; i<A.size(); i++){
        cout << "Path " << i << " : ";
        A[i].print();
    }
    return A;

}

double kth_largest(set<double> elts, int k){
    // find the kth largest element in a set of doubles
    auto elt = elts.rbegin();
    for(int i=1; i<min((int)elts.size(), k); i++){
        elt++;
    }
    return *elt;
}

vector<Path> yen(Graph *g, int s, int t, int k){
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
        for(int i=0; i<A.back().edges.size(); i++){
            // we are going find a new path, diverting from the old path from the spurNode
            int spurNode = A.back().edges[i].u;
            Path rootPath = Path({A.back().edges.begin(), A.back().edges.begin() + i});

            // edges we want to avoid
            set<int> edges_to_delete = set<int>();
            for(int j=0; j<A.size();j++){
                if(subpath_of(rootPath, A[j])){
                    edges_to_delete.insert(A[j].edges[i].index);
                }
            }

            // nodes we want to avoid
            set<int> nodes_to_delete = {};
            for(auto e: rootPath.edges){
                if(e.u != spurNode){
                    nodes_to_delete.insert(e.u);
                }
            }

            // create a new graph without all the unwanted nodes and edges
            Graph g2;
            g2.n = g->n; g2.m = 0;
            g2.adj = vector<vector<Edge>>(g2.n, vector<Edge>());
            for(int j=0; j<g2.n; j++){
                if(nodes_to_delete.find(j) == nodes_to_delete.end()){
                    for(auto e: g->adj[j]){
                        if(edges_to_delete.find(e.index) == edges_to_delete.end()){
                            g2.m++;
                            g2.adj[j].push_back(e);
                        }
                    }
                }
            }

            Path spurPath = dijkstra(&g2, spurNode, t);

            // compute the shortest path in the new graph from the spurnode to the terminal node
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
        for(int j=1; j<B.size(); j++){
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

        // compute the LB and UB for the new path
        update_LB(A, k_prime);
        update_UB(A, k_prime);
        LBs.insert(A[k_prime].LB);
    }

    cout << "Done with yen" << endl;
    for(int i=0; i<A.size(); i++){
        cout << "Path " << i << " : ";
        A[i].print();
    }
    return A;

}

Graph read_from_stdin(){
    int n, m;
    cin >> n >> m;
    Graph G = Graph({n, m, vector<vector<Edge>>(n, vector<Edge>())});
    int u, v, l; double p;
    for(int i=0; i<m; i++){
        cin >> u >> v >> l >> p;
        G.adj[u].push_back(Edge({u, v, l, p, i}));
    }
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
    cout << G.n << " " << G.m << endl;
    for(auto a: G.adj){
        for(auto elt: a){
            cout << elt.p << " ";
        }
        cout << endl;
    }

    Path p = dijkstra(&G, 0, 3);
    for(auto e: p.edges){
        cout << e.u << " " << e.v << endl;
    }

    testLB1_UB();


    dijkstra(&G, 0, 5);
    dijkstra(&G, 1, 5);
    dijkstra(&G, 0, 4);
    dijkstra(&G, 2, 3);
    dijkstra(&G, 1, 4);
    dijkstra(&G, 3, 5);



    cout << "YEN" << endl;
    classic_yen(&G, 0, 5, 3);
    classic_yen(&G, 0, 5, 10);

    yen(&G, 0, 5, 2);


	return 0;
}



