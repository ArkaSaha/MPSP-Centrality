# include <iostream>
# include <fstream>
# include <cstdlib>
# include <climits>
# include <queue>
# include <map>
# include <set>
# include <list>
# include <vector>
# include <algorithm>
# include <random>
# include <ctime>
# include <boost/heap/fibonacci_heap.hpp>

using namespace std;
using namespace boost::heap;

typedef tuple<int,long,double> adj_ent;
typedef tuple<int,int,long,double> edge;

struct Graph
{
	int n, m;
	vector<adj_ent>* adj;
};

tuple< list<edge>,long,double > prob_dijkstra(Graph* g, int s, int t, double& elapsed)
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
	list<edge> p = list<edge>();
	clock_t begin = clock();
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
			}
		}
	}
	clock_t end = clock();
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
	elapsed += ((double)(end - begin) / CLOCKS_PER_SEC);
	double pr = prob[t];
	delete [] handles;
	delete [] exist;
	delete [] visited;
	delete [] prev;
	delete [] prob;
	return make_tuple(p,min_dist,pr);
}

double approx_prob(Graph* g, vector< list<edge> > cp, list<edge> sp, double exist, double& elapsed)
{
	int C = 0, N = 1000000, n = cp.size();
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

tuple<list<edge>,int,int,bool> mpsp(Graph* g, int s, int t, int m, double& candidate_time, double& prob_time)
{
	double p_max = 0;
	int f_max = 1, n_s = 0;
	bool match = false;
	map< long,vector< tuple<list<edge>,double,int> > > paths = map< long,vector< tuple<list<edge>,double,int> > >();
	for (int i = 1; i <= m; i++)
	{
		if (n_s >= 10 && f_max >= 0.8 * n_s)
			break;
		list<edge> p;
		long w;
		double pr;
		tie(p,w,pr) = prob_dijkstra(g,s,t,candidate_time);
		if (p.empty())
			continue;
		n_s++;
		map< long,vector< tuple<list<edge>,double,int> > >::iterator it = paths.find(w);
		if (it == paths.end())
		{
			vector< tuple<list<edge>,double,int> > ss = vector< tuple<list<edge>,double,int> >();
			ss.push_back(make_tuple(p,pr,1));
			paths[w] = ss;
		}
		else
		{
			vector< tuple<list<edge>,double,int> >& vv = it->second;
			bool f = true;
			for (tuple<list<edge>,double,int>& tt : vv)
			{
				if (get<0>(tt) == p)
				{
					get<2>(tt)++;
					if (get<2>(tt) > f_max)
						f_max = get<2>(tt);
					f = false;
					break;
				}
			}
			if (f)
				vv.push_back(make_tuple(p,pr,1));
		}
	}
	list<edge> pp = list<edge>();
	vector< list<edge> > cp = vector< list<edge> >();
	map< int,set< list<edge> > > fp = map< int,set< list<edge> > >();
	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		vector< tuple<list<edge>,double,int> > vv = x->second;
		for (tuple<list<edge>,double,int> tt : vv)
		{
			list<edge> p = get<0>(tt);
			double ub = get<1>(tt);
			int freq = get<2>(tt);
			double prob = approx_prob(g,cp,p,ub,prob_time);
			cp.push_back(p);
			if (prob > p_max)
			{
				p_max = prob;
				pp = p;
			}
			if (fp.find(freq) == fp.end())
				fp[freq] = set< list<edge> >();
			fp[freq].insert(pp);
			if (freq > f_max)
				f_max = freq;
		}
	}
	if (fp[f_max].find(pp) != fp[f_max].end())
		match = true;
	return make_tuple(pp,cp.size(),f_max,match);
}

int bfs(Graph* g, int s, int d)
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

Graph generate_er(int n, int m, char* file)
{
	ofstream graph;
	graph.open(file);
	Graph g;
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

Graph read_graph(char* file)
{
	ifstream graph;
	graph.open(file);
	Graph g;
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

vector< tuple<int,int> > generate_queries(Graph* g, int n, int d, char* file)
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

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 4)
	{
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-queries> <path-to-output>" << endl;
		return 1;
	}
	int n = 100, d = 4;
	// Graph g = generate_er(40000,80000);
	Graph g = read_graph(argv[1]);
	// vector< tuple<int,int> > q = generate_queries(&g,n,d);
	vector< tuple<int,int> > q = read_queries(n,d,argv[2]);
	ofstream output;
	output.open(argv[3]);
	for (int k = d; k >= 0; k--)
	{
		// cerr << "k = " << k << endl;
		output << "Number of hops = " << k * 2 << endl << endl;
		long a_w = 0;
		double t_c = 0, t_p = 0, a_p = 0;
		int n_m = 0, num = 0;
		for (int i = 1; i <= n; i++)
		{
			// cerr << "i = " << i << endl;
			int s, t;
			tie(s,t) = q.back();
			q.pop_back();
			output << s << "\t" << t << endl;
			long wt = 0;
			double candidate_time = 0, prob_time = 0, pr = 1;
			list<edge> p;
			int c, f;
			bool m;
			tie(p,c,f,m) = mpsp(&g, s, t, 1000, candidate_time, prob_time);
			if (! p.empty())
			{
				num++;
				for (edge e : p)
				{
					output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
					wt += get<2>(e);
					pr *= get<3>(e);
				}
				a_w += wt;
				a_p += pr;
				t_c += candidate_time;
				t_p += prob_time;
				n_m += m;
			}
			output << "Number of Distinct Non-Pruned Candidate Paths : " << c << endl;
			output << "Frequency of Modal Candidate Path : " << f << endl;
			output << "MPSP matches Modal Candidate Path : " << m << endl;
			output << "Length of MPSP : " << wt << endl;
			output << "Probability of MPSP : " << pr << endl;
			output << "Candidate Generation Time : " << candidate_time << " seconds" << endl;
			output << "Probability Computation Time : " << prob_time << " seconds" << endl;
			output << endl;
		}
		output << "Number of Non-Empty Paths for " << k * 2 << " hops : " << num << endl;
		output << "Number of Path Matches for " << k * 2 << " hops : " << n_m << endl;
		output << "Average Length of MPSP for " << k * 2 << " hops : " << a_w / num << endl;
		output << "Average Probability of MPSP for " << k * 2 << " hops : " << a_p / num << endl;
		output << "Average Candidate Generation Time for " << k * 2 << " hops : " << t_c / num << " seconds" << endl;
		output << "Average Probability Computation Time for " << k * 2 << " hops : " << t_p / num << " seconds" << endl;
		output << "Average Total Time for " << k * 2 << " hops : " << (t_c + t_p) / num << " seconds" << endl;
		output << endl << endl;
	}
	delete [] g.adj;
	return 0;
}