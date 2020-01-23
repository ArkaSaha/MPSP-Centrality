#include<iostream>
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

};


double LB1(vector<Path> & paths, int n){
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
    cout << "LB1" << endl;
    cout << LB_prob << endl;
    return LB_prob;
}

void update_LB1(vector<Path> & paths, int n){
    double LB_prob = LB1(paths, n);
    paths[n].LB = max(paths[n].LB, LB_prob);
}


double UB(vector<Path> &paths, int n){
    double UB_prob = 1.0;
    for(int i=0; i<n; i++){
        UB_prob -= paths[i].LB;
    }
    return UB_prob;
}

void update_UB(vector<Path> &paths, int n){
    double UB_prob = UB(paths, n);
    paths[n].UB = min(paths[n].UB, UB_prob);
}

Path dijkstra(Graph* g, int s, int t){

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

    vector<ull> D = vector<ull>(g->n, ULLONG_MAX);
    vector<Edge> prev = vector<Edge>(g->n);
    typedef fibonacci_heap< node, compare<compare_node> >::handle_type handle_t;
    vector<handle_t> handles = vector<handle_t>(g->n);
    vector<bool> visited = vector<bool>(g->n, false);

    handles[s] = q.push(node({s, 0}));
    D[s] = 0;
    visited[s] = true;
    while(!q.empty()){
        node n = q.top();
        cout << "Currently at node " << n.u << endl;
        q.pop();
        if(n.u == t){
            break;
        }
        for(Edge e: g->adj[n.u]){
            cout << "Trying edge " << e.u << " - " << e.v << endl;
            if(!visited[e.v]){
                D[e.v] = D[e.u] + e.l;
                handles[e.v] = q.push(node({e.v, D[e.v]}));
                prev[e.v] = e;
            }
            else if(D[e.u] + e.l < D[e.v]){
                D[e.v] = D[e.u] + e.l;
                q.update(handles[e.v], node({e.v, D[e.v]}));
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

    update_LB1(paths, 0);
    update_LB1(paths, 1);
    update_LB1(paths, 2);
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


	return 0;
}



