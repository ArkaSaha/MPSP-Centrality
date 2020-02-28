# include <iostream>
# include <fstream>
# include <cstdlib>
# include <climits>
# include <list>
# include <queue>
# include <map>
# include <set>
# include <unordered_set>
# include <vector>
# include <algorithm>
# include <random>
# include <ctime>
# include <boost/heap/fibonacci_heap.hpp>
# include "topk.h"

using namespace std;
using namespace boost::heap;

typedef tuple<int,long,double> adj_ent;
typedef tuple<int,int,long,double> edge;

struct AdjGraph
{
	int n, m;
	vector<adj_ent>* adj;
};

tuple< list<edge>,long,double,double,bool > prob_dijkstra(AdjGraph* g, int s, int t, bool prune, double lb, double& elapsed)
{
	struct node
	{
		int vertex;
		long distance;
		node(const int& v, long d) : vertex(v), distance(d) {}
	};
	struct compare_node
	{
		bool operator()(const node& n1, const node& n2) const
		{
			return n1.distance > n2.distance;
		}
	};
	fibonacci_heap< node, compare<compare_node> > heap = fibonacci_heap< node, compare<compare_node> >();
	typedef fibonacci_heap< node, compare<compare_node> >::handle_type handle_t;
	handle_t* handles = new handle_t[g->n];
	bool* exist = new bool[g->n];
	bool* visited = new bool[g->n];
	adj_ent* prev = new adj_ent[g->n];
	double* prob = new double[g->n];
	unordered_set<int> lu = unordered_set<int>();
	for (int i = 0; i < g->n; i++)
	{
		exist[i] = false;
		visited[i] = false;
	}
	handles[s] = heap.push(node(s,0));
	exist[s] = true;
	prob[s] = 1;
	prob[t] = 1;
	long min_dist = 0;
	double prod = 1;
	list<edge> p = list<edge>();
	lu.insert(s);
	clock_t begin = clock();
	while (! (heap.empty() || (prune && lu.empty())))
	{
		node n = heap.top();
		int u = n.vertex;
		long d = n.distance;
		heap.pop();
		exist[u] = false;
		visited[u] = true;
		if (u == t)
		{
			min_dist = d;
			break;
		}
		bool leaf = true;
		for (adj_ent e : g->adj[u])
		{
			int v = get<0>(e);
			if (! visited[v])
			{
				double r = (double)rand() / RAND_MAX, pr = get<2>(e);
				if (r < pr)
				{
					leaf = false;
					long w = get<1>(e);
					long alt = d + w;
					if (! exist[v])
					{
						prev[v] = make_tuple(u,w,pr);
						prob[v] = prob[u] * pr;
						if (prune && prob[v] >= lb)
							lu.insert(v);
						handles[v] = heap.push(node(v,alt));
						exist[v] = true;
					}
					else if (alt < (*handles[v]).distance)
					{
						prev[v] = make_tuple(u,w,pr);
						prob[v] = prob[u] * pr;
						if (prune && prob[v] >= lb)
							lu.insert(v);
						heap.update(handles[v],node(v,alt));
					}
				}
				else
					prod *= (1 - pr);
			}
		}
		if (prune && !leaf)
			lu.erase(u);
	}
	clock_t end = clock();
	if (min_dist != 0)
	{
		if (prob[t] >= lb)
		{
			int v = t;
			while (v != s)
			{
				adj_ent res = prev[v];
				p.push_front(make_tuple(get<0>(res),v,get<1>(res),get<2>(res)));
				v = get<0>(res);
			}
		}
	}
	elapsed += ((double)(end - begin) / CLOCKS_PER_SEC);
	double pr = prob[t];
	delete [] handles;
	delete [] exist;
	delete [] visited;
	delete [] prev;
	delete [] prob;
	return make_tuple(p, min_dist, prod * pr, pr, prune && lu.empty());
}

double approx_prob(vector< list<edge> > cp, list<edge> sp, double exist, double& elapsed)
{
	int C = 0, N = 10000, n = cp.size();
	list<edge>* diff = new list<edge>[n];
	vector<double> pr = vector<double>(n);
	double S = 0;
	clock_t begin = clock();
	for (int i = 0; i < n; i++)
	{
		list<edge> p = cp[i];
		double prob = 1;
		list<edge> l = list<edge>();
		list<edge>::iterator it = set_difference(p.begin(),p.end(),sp.begin(),sp.end(),l.begin());
		l.resize(distance(l.begin(),it));
		for (edge e : l)
			prob *= get<3>(e);
		S += prob;
		diff[i] = l;
		pr[i] = prob;
	}
	clock_t m1 = clock();
	random_device rd;
	mt19937 gen(rd());
	discrete_distribution<> d(pr.begin(), pr.end());
	clock_t m2 = clock();
	for (int k = 1; k <= N; k++)
	{
		int i = d(gen);
		bool f1 = false;
		for (int j = 0; j < i; j++)
		{
			bool f2 = true;
			for (edge e : diff[j])
			{
				double r = (double)rand() / RAND_MAX;
				if (r > get<2>(e))
				{
					f2 = false;
					break;
				}
			}
			if (f2)
			{
				f1 = true;
				break;
			}
		}
		if (f1)
			C++;
	}
	clock_t end = clock();
	elapsed += ((double)(end - m2 + m1 - begin) / CLOCKS_PER_SEC);
	delete [] diff;
	return (1 - C * S / N) * exist;
}

tuple<list<edge>,int,int,bool,int,int> mpsp(AdjGraph* g, int s, int t, int m, double& candidate_time, double& prob_time)
{
	double lb_max = 0, p_max = 0;
	int f_max = 1, n_s = 0, n_d = 0, n_r = 1000, n_p = 0;
	bool match = false, prune = false;
	map< long,vector< tuple<list<edge>,double,double,int> > > paths = map< long,vector< tuple<list<edge>,double,double,int> > >();
	for (int i = 1; i <= n_r; i++)
	{
		if (n_s >= 10 && f_max >= 0.8 * n_s)
		{
			n_r = i;
			break;
		}
		list<edge> p;
		long w;
		double lb, ub;
		bool pruned;
		if (!prune && n_d >= 10)
			prune = true;
		tie(p,w,lb,ub,pruned) = prob_dijkstra(g,s,t,prune,lb_max,candidate_time);
		if (p.empty())
		{
			n_p += pruned;
			continue;
		}
		n_s++;
		if (lb_max < lb)
			lb_max = lb;
		map< long,vector< tuple<list<edge>,double,double,int> > >::iterator it = paths.find(w);
		if (it == paths.end())
		{
			vector< tuple<list<edge>,double,double,int> > ss = vector< tuple<list<edge>,double,double,int> >();
			ss.push_back(make_tuple(p,lb,ub,1));
			paths[w] = ss;
			n_d++;
		}
		else
		{
			vector< tuple<list<edge>,double,double,int> >& vv = it->second;
			bool f = true;
			for (tuple<list<edge>,double,double,int>& tt : vv)
			{
				if (get<0>(tt) == p)
				{
					double& l = get<1>(tt);
					if (lb > l)
						l = lb;
					get<3>(tt)++;
					if (get<3>(tt) > f_max)
						f_max = get<3>(tt);
					f = false;
					break;
				}
			}
			if (f)
			{
				vv.push_back(make_tuple(p,lb,ub,1));
				n_d++;
			}
		}
	}
	f_max = 1;
	list<edge> pp = list<edge>();
	vector< list<edge> > cp = vector< list<edge> >();
	map< int,set< list<edge> > > fp = map< int,set< list<edge> > >();
	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		vector< tuple<list<edge>,double,double,int> > vv = x->second;
		for (tuple<list<edge>,double,double,int> tt : vv)
		{
			list<edge> p = get<0>(tt);
			double lb = get<1>(tt);
			double ub = get<2>(tt);
			int freq = get<3>(tt);
			if (prune && ub < lb_max)
				n_p += freq;
			else
			{
				double prob = approx_prob(cp,p,ub,prob_time);
				cp.push_back(p);
				if (prune && (prob < lb || prob > ub))
					n_p += freq;
				else
				{
					if (prob > p_max)
					{
						p_max = prob;
						pp = p;
					}
				}
				if (fp.find(freq) == fp.end())
					fp[freq] = set< list<edge> >();
				fp[freq].insert(pp);
				if (freq > f_max)
					f_max = freq;
			}
		}
	}
	if (fp[f_max].find(pp) != fp[f_max].end())
		match = true;
	return make_tuple(pp,cp.size(),f_max,match,n_r,n_p);
}


vector<double> betweenness(AdjGraph & g){

    double _t1, _t2;

    vector<double> B = vector<double>(g.n, 0);

    for(int s=0; s<g.n; s++){
        for(int t=0; t<g.n; t++){
            if(s == t) continue;
            auto cur_mpsp = mpsp(&g, s, t, 1000, _t1, _t2);

            if(get<1>(cur_mpsp) > 0){
                cout << "path between " << s << " and " << t << endl;
                // this means that s and t are connected
                // raise the betweenness of every inner node of the path by 1
                list<edge> p = get<0>(cur_mpsp);
                cout << "inner nodes : ";
                for(auto it = next(p.begin()); it != p.end(); it++){
                    cout << get<0>(*it) << " ";
                    B[get<0>(*it)]++;
                }
                cout << endl;
            }
        }
    }

    // normalize betweenness by size of graph
    for(uint i=0; i<B.size(); i++){
        B[i] /= ((g.n-1) * (g.n-2));
    }

    return B;
}


int bfs(AdjGraph* g, int s, int d)
{
	bool* visited = new bool[g->n];
	for (int i=0; i<g->n; i++)
		visited[i] = false;
	queue< tuple<int,int> > q = queue< tuple<int,int> >();
	q.push(make_tuple(s,0));
	while (! q.empty())
	{
		int t, h;
		tie(t,h) = q.front();
		q.pop();
		visited[t] = true;
		if (h == d)
		{
			delete [] visited;
			return t;
		}
		for (adj_ent e : g->adj[t])
		{
			int v = get<0>(e);
			if (! visited[v])
				q.push(make_tuple(v,h+1));
		}
	}
	delete [] visited;
	return -1;
}

AdjGraph generate_er(int n, int m, char* file)
{
	ofstream graph;
	graph.open(file);
	AdjGraph g;
	g.n = n;
	g.m = m;
	graph << g.n << " " << g.m << endl;
	g.adj = new vector<adj_ent>[g.n];
	for (int i=0; i<g.n; i++)
		g.adj[i] = vector<adj_ent>();
	for (int i=0; i<g.m; i++)
	{
		int u = (int)((double)rand() / RAND_MAX * g.n), v = (int)((double)rand() / RAND_MAX * g.n);
		if (u == v || u == g.n || v == g.n)
		{
			i--;
			continue;
		}
		bool f = true;
		for (adj_ent e : g.adj[u])
		{
			if (get<0>(e) == v)
			{
				f = false;
				break;
			}
		}
		if (f && v != g.n)
		{
			long w = rand() % 1000 + 1;
			double p = (double)rand() / RAND_MAX;
			g.adj[u].push_back(make_tuple(v, w, p));
			graph << u << " " << v << " " << w << " " << p << endl;
		}
		else
			i--;
	}
	graph.close();
	return g;
}

AdjGraph read_graph(char* file)
{
	ifstream graph;
	graph.open(file);
	AdjGraph g;
	graph >> g.n >> g.m;
	g.adj = new vector<adj_ent>[g.n];
	for (int i=0; i<g.n; i++)
		g.adj[i] = vector<adj_ent>();
	for (int i = 0; i < g.m; i++)
	{
		int u, v;
		long w;
		double p;
		graph >> u >> v >> w >> p;
		g.adj[u].push_back(make_tuple(v, w, p));
	}
	graph.close();
	return g;
}

vector< tuple<int,int> > generate_queries(AdjGraph* g, int n, int d, char* file)
{
	ofstream queries;
	queries.open(file);
	vector< tuple<int,int> > q = vector< tuple<int,int> >();
	for (int k = 1; k <= d; k++)
	{
		for (int i = 1; i <= n; i++)
		{
			int s = (int)((double)rand() / RAND_MAX * (g->n - 1));
			while (g->adj[s].empty())
				s = (int)((double)rand() / RAND_MAX * (g->n - 1));
			int t = bfs(g, s, k * 2);
			if (t == -1)
			{
				i--;
				continue;
			}
			queries << s << " " << t << endl;
			q.push_back(make_tuple(s,t));
		}
	}
	queries.close();
	return q;
}

vector< tuple<int,int> > read_queries(int n, int d, char* file)
{
	ifstream queries;
	queries.open(file);
	vector< tuple<int,int> > q = vector< tuple<int,int> >();
	for (int i = 1; i <= (d + 1) * n; i++)
	{
		int s, t;
		queries >> s >> t;
		q.push_back(make_tuple(s,t));
	}
	queries.close();
	return q;
}

void experiment_betweenness(char* path_to_graph, char* path_to_output){
    cout << "checkpoint1" << endl;
    AdjGraph g = read_graph(path_to_graph);
    cout << "checkpoint2" << endl;

	ofstream output;
	output.open(path_to_output);

    cout << "checkpoint" << endl;
    auto B = betweenness(g);

    cout << "checkpoint" << endl;
    // sort B from highest to lowest betweenness
    vector<pair<double, int>> sorted_B = vector<pair<double, int>>(B.size());
    for(int i=0; i<B.size(); i++){
        sorted_B[i] = {B[i], i};
    }
    sort(sorted_B.rbegin(), sorted_B.rend());

    output << sorted_B.size() << endl;
    for(auto elt: sorted_B){
        output << elt.second << " " << elt.first << endl;
    }

    output.close();
}

void experiment(char* path_to_graph, char* path_to_queries, char* path_to_output){

	int n = 100, d = 4;
	AdjGraph g = read_graph(path_to_graph);
	Graph G = Graph({g.n, g.m});
	int num = 0;
	for (int i = 0; i < g.n; i++)
	{
		for (adj_ent t : g.adj[i])
		{
			int j = get<0>(t);
			long l = get<1>(t);
			double p = get<2>(t);
			G.adj[i].push_back(Edge({i, j, (int)l, p, num}));
			num++;
		}
	}
	G.update_incoming_index2edge();
	// vector< tuple<int,int> > q = generate_queries(&g,n,d,argv[2]);
	vector< tuple<int,int> > q = read_queries(n,d,path_to_queries);
	ofstream output;
	output.open(path_to_output);
	for (int k = d; k >= 0; k--)
	{
		// cerr << "k = " << k << endl;
		output << "Number of hops = " << k * 2 << endl << endl;
		long a_w = 0;
		double t_c = 0, t_p = 0, a_p = 0, a_r = 0, a_pr = 0;
		int n_m = 0;
		for (int i = 1; i <= n; i++)
		{
			// cerr << "i = " << i << endl;
			int r, s, t;
			tie(s,t) = q.back();
			q.pop_back();
			output << s << "\t" << t << endl;
			long wt = 0;
			double candidate_time = 0, prob_time = 0, pr = 1, prob = 0;
			list<edge> p;
			int c, f, pruned;
			bool m;
			tie(p,c,f,m,r,pruned) = mpsp(&g, s, t, 1000, candidate_time, prob_time);
			if (! p.empty())
			{
				Path path = Path({});
				for (edge e : p)
				{
					Edge ee;
					ee.u = get<0>(e);
					ee.v = get<1>(e);
					ee.l = (int)get<2>(e);
					ee.p = get<3>(e);
					for (Edge eee : G.adj[ee.u])
					{
						if (eee.v == ee.v)
						{
							ee.index = eee.index;
							break;
						}
					}
					path.edges.push_back(ee);
				}
				prob = Luby_Karp(G, path, 10000);
			}
			for (edge e : p)
			{
				output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
				wt += get<2>(e);
				pr *= get<3>(e);
			}
			a_w += wt;
			a_p += pr;
			a_r += r;
			a_pr += pruned;
			t_c += candidate_time;
			t_p += prob_time;
			n_m += m;
			output << "Number of Dijkstra Runs : " << r << endl;
			output << "Number of Samples Pruned : " << pruned << endl;
			output << "Number of Distinct Non-Pruned Candidate Paths : " << c << endl;
			output << "Frequency of Modal Candidate Path : " << f << endl;
			output << "MPSP matches Modal Candidate Path : " << m << endl;
			output << "Length of MPSP : " << wt << endl;
			output << "Probability of MPSP : " << pr << endl;
			output << "Probability of MPSP being the Shortest Path: " << prob << endl;
			output << "Candidate Generation Time : " << candidate_time << " seconds" << endl;
			output << "Probability Computation Time : " << prob_time << " seconds" << endl;
			output << endl;
		}
		output << "Number of Path Matches for " << k * 2 << " hops : " << n_m << endl;
		output << "Average Length of MPSP for " << k * 2 << " hops : " << a_w / 100 << endl;
		output << "Average Probability of MPSP for " << k * 2 << " hops : " << a_p / 100 << endl;
		output << "Average Number of Dijkstra Runs for " << k * 2 << " hops : " << a_r / 100 << endl;
		output << "Average Number of Samples Pruned for " << k * 2 << " hops : " << a_pr / 100 << endl;
		output << "Average Candidate Generation Time for " << k * 2 << " hops : " << t_c / 100 << " seconds" << endl;
		output << "Average Probability Computation Time for " << k * 2 << " hops : " << t_p / 100 << " seconds" << endl;
		output << "Average Total Time for " << k * 2 << " hops : " << (t_c + t_p) / 100 << " seconds" << endl;
		output << endl << endl;
	}
	output.close();
	delete [] g.adj;

}

int main(int argc, char* argv[])
{

    cerr << argc << endl;
	if (argc < 3 or argc > 4)
	{
        cerr << "For mpsp experiment" << endl;
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-queries> <path-to-output>" << endl;
        cerr << "For betweenness experiment" << endl;
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-output>" << endl;
		return 1;
	}
    else if(argc == 3){
        cout << "Doing betweenness experiment" << endl;
        experiment_betweenness(argv[1], argv[2]);
        
    }
    else if(argc == 4){
        experiment(argv[1], argv[2], argv[3]);
    }
	//srand(time(NULL));
	// Graph g = generate_er(40000,80000,argv[1]);

	return 0;
}
