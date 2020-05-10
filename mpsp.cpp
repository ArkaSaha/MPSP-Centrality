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
# include <cmath>
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

double time_difference(timespec begin, timespec end)
{
	return (end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec) / pow(10,9);
}

double dijkstra(AdjGraph* g, int s, int t)
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
	double* prob = new double[g->n];
	for (int i = 0; i < g->n; i++)
	{
		exist[i] = false;
		visited[i] = false;
	}
	handles[s] = heap.push(node(s,0));
	exist[s] = true;
	prob[s] = 1;
	prob[t] = 0;
	list<edge> p = list<edge>();
	while (! heap.empty())
	{
		node n = heap.top();
		int u = n.vertex;
		long d = n.distance;
		heap.pop();
		exist[u] = false;
		visited[u] = true;
		if (u == t)
			break;
		for (adj_ent e : g->adj[u])
		{
			int v = get<0>(e);
			if (! visited[v])
			{
				double pr = get<2>(e);
				long w = get<1>(e);
				long alt = d + w;
				if (! exist[v])
				{
					prob[v] = prob[u] * pr;
					handles[v] = heap.push(node(v,alt));
					exist[v] = true;
				}
				else if (alt < (*handles[v]).distance)
				{
					prob[v] = prob[u] * pr;
					heap.update(handles[v],node(v,alt));
				}
			}
		}
	}
	double pr = prob[t];
	delete [] handles;
	delete [] exist;
	delete [] visited;
	delete [] prob;
	return pr;
}

tuple< list<edge>,long,double,double > prob_dijkstra(AdjGraph* g, int s, int t, double& elapsed_prune, double& elapsed_noprune)
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
	timespec begin, end;
	clock_gettime(CLOCK_MONOTONIC,&begin);
	while (! heap.empty())
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
		for (adj_ent e : g->adj[u])
		{
			int v = get<0>(e);
			if (! visited[v])
			{
				double r = (double)rand() / RAND_MAX, pr = get<2>(e);
				if (r < pr)
				{
					long w = get<1>(e);
					long alt = d + w;
					if (! exist[v])
					{
						prev[v] = make_tuple(u,w,pr);
						prob[v] = prob[u] * pr;
						handles[v] = heap.push(node(v,alt));
						exist[v] = true;
					}
					else if (alt < (*handles[v]).distance)
					{
						prev[v] = make_tuple(u,w,pr);
						prob[v] = prob[u] * pr;
						heap.update(handles[v],node(v,alt));
					}
				}
				else
					prod *= (1 - pr);
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double time = time_difference(begin,end);
	elapsed_noprune += time;
	elapsed_prune += time;
	if (min_dist != 0)
	{
		int v = t;
		while (v != s)
		{
			adj_ent res = prev[v];
			p.push_front(make_tuple(get<0>(res),v,get<1>(res),get<2>(res)));
			v = get<0>(res);
		}
	}
	double pr = prob[t];
	delete [] handles;
	delete [] exist;
	delete [] visited;
	delete [] prev;
	delete [] prob;
	return make_tuple(p, min_dist, prod * pr, pr);
}

double approx_prob(vector< list<edge> > cp, list<edge> sp, double exist, double& elapsed)
{
	int C = 0, N = 1000, n = cp.size();
	list<edge>* diff = new list<edge>[n];
	vector<double> pr = vector<double>(n);
	double S = 0;
	timespec begin, end, m1, m2;
	clock_gettime(CLOCK_MONOTONIC,&begin);
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
	clock_gettime(CLOCK_MONOTONIC,&m1);
	random_device rd;
	mt19937 gen(rd());
	discrete_distribution<> d(pr.begin(), pr.end());
	clock_gettime(CLOCK_MONOTONIC,&m2);
	for (int k = 1; k <= N; k++)
	{
		map<edge,bool> sampled = map<edge,bool>();
		int i = d(gen);
		bool f1 = false;
		for (int j = 0; j < i; j++)
		{
			bool f2 = true;
			for (edge e : diff[j])
			{
				bool s = false;
				map<edge,bool>::iterator it = sampled.find(e);
				if (it == sampled.end())
				{
					double r = (double)rand() / RAND_MAX;
					s = sampled[e] = (r < get<2>(e));
				}
				else
					s = it->second;
				if (!s)
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
	clock_gettime(CLOCK_MONOTONIC,&end);
	elapsed += (time_difference(begin,end) - time_difference(m1,m2));
	delete [] diff;
	return (1 - C * S / N) * exist;
}

tuple<list<edge>,list<edge>,int,int,int,double,double,bool,bool,int,int> mpsp(AdjGraph* g, int s, int t, int m, double& candidate_time_prune, double& candidate_time_noprune, double& prob_time_prune, double& prob_time_noprune)
{
	double lb_max = 0, p_max_p = 0, p_max_np = 0;
	int f_max = 1, n_s = 0, n_d = 0, n_r = 20, n_p = 0;
	bool match_p = false, match_np = false;
	map< long,vector< tuple<list<edge>,double,double,int> > > paths = map< long,vector< tuple<list<edge>,double,double,int> > >();
	for (int i = 1; i <= n_r; i++)
	{
		if (n_s >= 10 && f_max >= 0.5 * n_s)
		{
			n_r = i;
			break;
		}
		list<edge> p;
		long w;
		double lb, ub;
		tie(p,w,lb,ub) = prob_dijkstra(g,s,t,candidate_time_prune,candidate_time_noprune);
		if (p.empty())
			continue;
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
	list<edge> pp_p, pp_np = list<edge>();
	vector< list<edge> > cp_p = vector< list<edge> >(), cp_np = vector< list<edge> >();
	map< int,set< list<edge> > > fp = map< int,set< list<edge> > >();
	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		vector< tuple<list<edge>,double,double,int> > vv = x->second;
		for (tuple<list<edge>,double,double,int> tt : vv)
		{
			list<edge> p = get<0>(tt);
			// double lb = get<1>(tt);
			double ub = get<2>(tt);
			int freq = get<3>(tt);
			double prob_np = approx_prob(cp_np,p,ub,prob_time_noprune);
			cp_np.push_back(p);
			if (prob_np > p_max_np)
			{
				p_max_np = prob_np;
				pp_np = p;
			}
			if (ub < lb_max)
				n_p++;
			else
			{
				double prob_p = approx_prob(cp_p,p,ub,prob_time_prune);
				cp_p.push_back(p);
				// if (prob_p < lb || prob_p > ub)
				// 	n_p += freq;
				if (prob_p > p_max_p)
				{
					p_max_p = prob_p;
					pp_p = p;
				}
			}
			if (fp.find(freq) == fp.end())
				fp[freq] = set< list<edge> >();
			fp[freq].insert(p);
			if (freq > f_max)
				f_max = freq;
		}
	}
	if (fp[f_max].find(pp_p) != fp[f_max].end())
		match_p = true;
	if (fp[f_max].find(pp_np) != fp[f_max].end())
		match_np = true;
	return make_tuple(pp_p,pp_np,cp_p.size(),cp_np.size(),f_max,p_max_p,p_max_np,match_p,match_np,n_r,n_p);
}


vector<double> betweenness(AdjGraph & g)
{
    double _t1, _t2;

    vector<double> B = vector<double>(g.n, 0);

//     clock_t start = clock();

//     for(int s=0; s<g.n; s++)
//     {
//         clock_t t1 = clock();
//         cout << (s);
//         for(int t=0; t<g.n; t++)
//         {
//             if(s == t) continue;
//             auto cur_mpsp = mpsp(&g, s, t, 1000, _t1, _t2);

//             if(get<1>(cur_mpsp) > 0)
//             {
//                 // this means that s and t are connected
//                 // raise the betweenness of every inner node of the path by 1
//                 list<edge> p = get<0>(cur_mpsp);
//                 for(auto it = next(p.begin()); it != p.end(); it++){
//                     B[get<0>(*it)]++;
//                 }
//             }
//         }
//         clock_t t2 = clock();
//         cout << " : " << (double(t2-t1))/CLOCKS_PER_SEC << " seconds" << endl;
//         cout << "remaining : " << (g.n-(s+1)) * (double(t2 - start))/((s+1) * CLOCKS_PER_SEC) << endl;
//     }

//     // normalize betweenness by size of graph
//     for(uint i=0; i<B.size(); i++)
//     {
//         B[i] /= ((g.n-1) * (g.n-2));
//     }

    return B;
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

void experiment_betweenness(char* path_to_graph, char* path_to_output)
{
    clock_t t1 = clock();
    AdjGraph g = read_graph(path_to_graph);

	ofstream output;
	output.open(path_to_output);

    auto B = betweenness(g);

    // sort B from highest to lowest betweenness
    vector<pair<double, size_t>> sorted_B = vector<pair<double, size_t>>(B.size());
    for (size_t i=0; i<B.size(); i++)
    {
        sorted_B[i] = {B[i], i};
    }
    sort(sorted_B.rbegin(), sorted_B.rend());

    output << sorted_B.size() << endl;
    for (auto elt: sorted_B)
    {
        output << fixed << elt.second << " " << elt.first << endl;
    }

    output.close();
    clock_t t2 = clock();
    cout << "Betweenness calculation took " << (double(t2 - t1))/CLOCKS_PER_SEC << " seconds" << endl;

}

void experiment(char* path_to_graph, char* path_to_queries, char* path_to_output)
{
	int num = 0, d;
	AdjGraph g = read_graph(path_to_graph);
	Graph G = Graph({g.n, g.m});
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
	ifstream queries;
	queries.open(path_to_queries);
	ofstream output;
	output.open(path_to_output);
	queries >> d;
	for (int k = 1; k <= d; k++)
	{
		// cerr << "k = " << k << endl;
		int h, n;
		queries >> h >> n;
		output << "Number of hops = " << h << endl;
		output << "Number of queries = " << n << endl << endl;
		long a_w_p = 0, a_w_np = 0;
		double t_c_p = 0, t_c_np = 0, t_p_p = 0, t_p_np = 0, a_c_p = 0, a_c_np = 0, a_p_p = 0, a_p_np = 0, a_r = 0, a_pr = 0;
		int n_m_p = 0, n_m_np = 0, num = 0;
		for (int i = 1; i <= n; i++)
		{
			// cerr << "i = " << i << endl;
			int r, s, t;
			queries >> s >> t;
			output << s << "\t" << t << endl;
			long wt_p = 0, wt_np = 0;
			double candidate_time_prune = 0, candidate_time_noprune = 0, prob_time_prune = 0, prob_time_noprune = 0, pr_p = 1, pr_np = 1, prob_p = 0, prob_np = 0;
			list<edge> p_p, p_np;
			int c_p, c_np, f, pruned;
			bool m_p, m_np;
			tie(p_p,p_np,c_p,c_np,f,prob_p,prob_np,m_p,m_np,r,pruned) = mpsp(&g, s, t, 1000, candidate_time_prune, candidate_time_noprune, prob_time_prune, prob_time_noprune);
			if (! p_p.empty())
			{
				num++;
				output << "With Pruning" << endl;
                for (edge e : p_p)
				{
					output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
					wt_p += get<2>(e);
					pr_p *= get<3>(e);
				}
				a_c_p += c_p;
				a_w_p += wt_p;
				a_p_p += pr_p;
				a_r += r;
				a_pr += pruned;
				t_c_p += candidate_time_prune;
				t_p_p += prob_time_prune;
				n_m_p += m_p;
				Path path = Path({});
				for (edge e : p_p)
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
				prob_p = Luby_Karp(G, path, 1000);
			}
			else
				prob_p = dijkstra(&g,s,t);
			if (! p_np.empty())
			{
				output << "Without Pruning" << endl;
                for (edge e : p_np)
				{
					output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
					wt_np += get<2>(e);
					pr_np *= get<3>(e);
				}
				a_c_np += c_np;
				a_w_np += wt_np;
				a_p_np += pr_np;
				t_c_np += candidate_time_noprune;
				t_p_np += prob_time_noprune;
				n_m_np += m_np;
				if (p_np != p_p)
				{
					Path path = Path({});
					for (edge e : p_np)
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
					prob_np = Luby_Karp(G, path, 1000);
				}
				else
					prob_np = prob_p;
			}
			else
				prob_np = prob_p;
			output << "Number of Dijkstra Runs : " << r << endl;
			output << "Number of Paths Pruned : " << pruned << endl;
			output << "Number of Distinct Candidate Paths with Pruning : " << c_p << endl;
			output << "Number of Distinct Candidate Paths without Pruning : " << c_np << endl;
			output << "Frequency of Modal Candidate Path : " << f << endl;
			output << "MPSP matches Modal Candidate Path with Pruning : " << m_p << endl;
			output << "MPSP matches Modal Candidate Path without Pruning : " << m_np << endl;
			output << "Length of MPSP with Pruning : " << wt_p << endl;
			output << "Length of MPSP without Pruning : " << wt_np << endl;
			output << "Probability of MPSP with Pruning : " << pr_p << endl;
			output << "Probability of MPSP without Pruning : " << pr_np << endl;
			output << "Probability of MPSP being the Shortest Path with Pruning : " << prob_p << endl;
			output << "Probability of MPSP being the Shortest Path without Pruning : " << prob_np << endl;
			output << "Candidate Generation Time with Pruning : " << candidate_time_prune << " seconds" << endl;
			output << "Candidate Generation Time without Pruning : " << candidate_time_noprune << " seconds" << endl;
			output << "Probability Computation Time with Pruning : " << prob_time_prune << " seconds" << endl;
			output << "Probability Computation Time without Pruning : " << prob_time_noprune << " seconds" << endl;
			output << endl;
		}
		output << "Number of Non-Empty Paths for " << h << " hops : " << num << endl;
		output << "Number of Path Matches with Pruning for " << h << " hops : " << n_m_p << endl;
		output << "Number of Path Matches without Pruning for " << h << " hops : " << n_m_np << endl;
		output << "Average Length of MPSP with Pruning for " << h << " hops : " << a_w_p / num << endl;
		output << "Average Length of MPSP without Pruning for " << h << " hops : " << a_w_np / num << endl;
		output << "Average Probability of MPSP with Pruning for " << h << " hops : " << a_p_p / num << endl;
		output << "Average Probability of MPSP without Pruning for " << h << " hops : " << a_p_np / num << endl;
		output << "Average Number of Dijkstra Runs for " << h << " hops : " << a_r / num << endl;
		output << "Average Number of Paths Pruned for " << h << " hops : " << a_pr / num << endl;
		output << "Average Number of Distinct Candidate Paths with Pruning for " << h << " hops : " << a_c_p / num << endl;
		output << "Average Number of Distinct Candidate Paths without Pruning for " << h << " hops : " << a_c_np / num << endl;
		output << "Average Candidate Generation Time with Pruning for " << h << " hops : " << t_c_p / num << " seconds" << endl;
		output << "Average Candidate Generation Time without Pruning for " << h << " hops : " << t_c_np / num << " seconds" << endl;
		output << "Average Probability Computation Time with Pruning for " << h << " hops : " << t_p_p / num << " seconds" << endl;
		output << "Average Probability Computation Time without Pruning for " << h << " hops : " << t_p_np / num << " seconds" << endl;
		output << "Average Total Time with Pruning for " << h << " hops : " << (t_c_p + t_p_p) / num << " seconds" << endl;
		output << "Average Total Time without Pruning for " << h << " hops : " << (t_c_np + t_p_np) / num << " seconds" << endl;
		output << endl << endl;
	}
	queries.close();
	output.close();
	delete [] g.adj;
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc < 3 or argc > 4)
	{
        cerr << "For mpsp experiment" << endl;
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-queries> <path-to-output>" << endl;
        cerr << "For betweenness experiment" << endl;
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-output>" << endl;
		return 1;
	}
    else if(argc == 3)
    {
        cout << "Doing betweenness experiment" << endl;
        experiment_betweenness(argv[1], argv[2]);
    }
    else if(argc == 4)
    {
        experiment(argv[1], argv[2], argv[3]);
    }

	return 0;
}
