# include <iostream>
# include <fstream>
# include <cstdlib>
# include <climits>
# include <list>
# include <map>
# include <set>
# include <unordered_set>
# include <vector>
# include <algorithm>
# include <random>
# include <ctime>
# include <cmath>
# include <boost/heap/fibonacci_heap.hpp>

using namespace std;
using namespace boost::heap;

mt19937 mersenne_mpsp{static_cast<mt19937::result_type>(12345)};

# define NUM_MPSP 1
# define DIJKSTRA_RUNS 20 
# define LUBY_KARP_SAMPLES 1000


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

tuple< vector< list<edge> >, vector<long>, vector<double> > prob_dijkstra(AdjGraph* g, int s, double& elapsed)
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
	vector<long> length = vector<long>(g->n, 0);
	vector<double> prob = vector<double>(g->n ,0);
	for (int i = 0; i < g->n; i++)
	{
		exist[i] = false;
		visited[i] = false;
	}
	handles[s] = heap.push(node(s,0));
	exist[s] = true;
	prob[s] = 1;
	vector< list<edge> > p = vector< list<edge> >(g->n, list<edge>());
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
		length[u] = d;
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
					double ch = prob[u] * pr;
					if (! exist[v])
					{
						prev[v] = make_tuple(u,w,pr);
						prob[v] = ch;
						handles[v] = heap.push(node(v,alt));
						exist[v] = true;
					}
					else
					{
						long dist = (*handles[v]).distance;
						if (alt < dist || (alt == dist && ch > prob[v]))
						{
							prev[v] = make_tuple(u,w,pr);
							prob[v] = ch;
							heap.update(handles[v],node(v,alt));
						}
					}
				}
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double time = time_difference(begin,end);
	elapsed += time;
	for (int t = 0; t < g->n; t++)
	{
		if (! visited[t])
			continue;
		int v = t;
		while (v != s)
		{
			adj_ent res = prev[v];
			p[t].push_front(make_tuple(get<0>(res),v,get<1>(res),get<2>(res)));
			v = get<0>(res);
		}
	}
	delete [] handles;
	delete [] exist;
	delete [] visited;
	delete [] prev;
	return make_tuple(p, length, prob);
}

double approx_prob(vector< pair<list<edge>,double> > cp, list<edge> sp, int N, double exist, double& elapsed)
{
	int C = 0, n = cp.size();
	if (n == 0)
		return exist;
	auto diff = vector< list<edge> >(n);
	vector<double> pr = vector<double>(n);
	double S = 0;
	timespec begin, end, m1, m2;
	clock_gettime(CLOCK_MONOTONIC,&begin);
	for (int i = 0; i < n; i++)
	{
		list<edge> p = cp[i].first;
		double prob = 1;
		list<edge> l = list<edge>();
		for (edge e : p)
		{
			bool keep = true;
			for (edge ee : sp)
			{
				if (ee == e)
				{
					keep = false;
					break;
				}
			}
			if (keep)
				l.push_back(e);
		}
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
		for (edge e : diff[i])
			sampled[e] = true;
		bool f1 = true;
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
				f1 = false;
				break;
			}
		}
		if (f1)
			C++;
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	elapsed += (time_difference(begin,end) - time_difference(m1,m2));
	return (1 - C * S / N) * exist;
}

void mpsp(AdjGraph* g, int s, size_t k, int m, int N, double& candidate_time, double& prob_time)
{
	// double t_pr = 0;
	vector< map< long,vector< tuple<list<edge>,double> > > > paths = vector< map< long,vector< tuple<list<edge>,double> > > >(g->n, map< long,vector< tuple<list<edge>,double> > >());
	for (int i = 1; i <= m; i++)
	{
		auto res = prob_dijkstra(g,s,candidate_time);
		vector< list<edge> >& p = get<0>(res);
		vector<long>& w = get<1>(res);
		vector<double>& ub = get<2>(res);
		for (int t = 0; t < g->n; t++)
		{
			if (p[t].empty())
				continue;
			map< long,vector< tuple<list<edge>,double> > >::iterator it = paths[t].find(w[t]);
			if (it == paths[t].end())
			{
				vector< tuple<list<edge>,double> > ss = vector< tuple<list<edge>,double> >();
				ss.push_back(make_tuple(p[t],ub[t]));
				paths[t][w[t]] = ss;
			}
			else
			{
				vector< tuple<list<edge>,double> >& vv = it->second;
				bool f = true;
				for (tuple<list<edge>,double>& tt : vv)
				{
					if (get<0>(tt) == p[t])
					{
						f = false;
						break;
					}
				}
				if (f)
					vv.push_back(make_tuple(p[t],ub[t]));
			}
		}
	}

	for (int t = 0; t < g->n; t++)
	{
		vector< pair<list<edge>,double> > cp = vector< pair<list<edge>,double> >();

		for (auto x = paths[t].begin(); x != paths[t].end(); x++)
		{
			vector< tuple<list<edge>,double> > vv = x->second;
			for (tuple<list<edge>,double> tt : vv)
			{
				list<edge> p = get<0>(tt);
				double ub = get<1>(tt);
				double prob = approx_prob(cp,p,N,ub,prob_time);
				cp.push_back(make_pair(p,prob));
			}
		}
		struct { bool operator() (pair<list<edge>,double> x, pair<list<edge>,double> y) { return x.second > y.second; } } comp;
		timespec t1, t2;
		clock_gettime(CLOCK_MONOTONIC,&t1);
		sort(cp.begin(),cp.end(),comp);
		clock_gettime(CLOCK_MONOTONIC,&t2);
		prob_time += time_difference(t1,t2);

		// vector< list<edge> > pp = vector< list<edge> >();
		// double prob = 0;
		// size_t len = min(k,cp.size());
		// for (size_t i = 0; i < len; i++)
		// {
		// 	pp.push_back(cp[i].first);
		// 	prob += cp[i].second;
		// }
		// if (len)
		// 	prob /= len;
		// t_pr += prob;
	}

	// return t_pr / g->n;
}

AdjGraph read_graph(char* file, bool source)
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
		if (source)
			g.adj[u].push_back(make_tuple(v, w, p));
		else
			g.adj[v].push_back(make_tuple(u, w, p));
	}
	graph.close();
	return g;
}

void experiment(char* path_to_graph, char* path_to_queries, char* path_to_output, size_t k = NUM_MPSP, int m = DIJKSTRA_RUNS, int N = LUBY_KARP_SAMPLES)
{
	ifstream queries;
	queries.open(path_to_queries);
	ofstream output;
	output.open(path_to_output);
	bool source;
	queries >> source;
	output << "Single " << (source ? "source" : "target") << endl;
	AdjGraph g = read_graph(path_to_graph, source);
	int n;
	queries >> n;
	output << "Number of queries = " << n << endl << endl;
	double t_c = 0, t_p = 0, a_pr = 0;
	for (int i = 1; i <= n; i++)
	{
		int s;
		queries >> s;
		output << s << endl;
		double candidate_time = 0, prob_time = 0;
		mpsp(&g, s, k, m, N, candidate_time, prob_time);
		t_c += candidate_time;
		t_p += prob_time;
		// a_pr += prob;
		// output << "SP Probability : " << prob << endl;
		output << "Candidate Generation Time : " << candidate_time << " seconds" << endl;
		output << "Probability Computation Time : " << prob_time << " seconds" << endl;
		output << endl;
	}
	output << "Average Candidate Generation Time : " << t_c / n << " seconds" << endl;
	output << "Average Probability Computation Time : " << t_p / n << " seconds" << endl;
	output << "Average Total Time : " << (t_c + t_p) / n << " seconds" << endl;
	// output << "Average SP Probability : " << a_pr / n << endl;
	queries.close();
	output.close();
	delete [] g.adj;
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 7)
	{
		cerr << "Usage: ./single <path-to-graph> <path-to-queries> <path-to-output> <k> <m> <N>" << endl;
		return EXIT_FAILURE;
	}
    else
    {
        experiment(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
    }
	return EXIT_SUCCESS;
}
