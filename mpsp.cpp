# include <iostream>
# include <fstream>
# include <cstdlib>
# include <climits>
# include <list>
# include <queue>
# include <map>
# include <set>
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

tuple< list<edge>,long,double,double > prob_dijkstra(Graph* g, int s, int t, double lb, double& elapsed)
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
	set<int> lu = set<int>();
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
	while (! (heap.empty() || lu.empty()))
	{
		node n = heap.top();
		int u = n.vertex;
		long d = n.distance;
		heap.pop();
		exist[u] = false;
		visited[u] = true;
		lu.erase(u);
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
						if (prob[v] >= lb)
							lu.insert(v);
						handles[v] = heap.push(node(v,alt));
						exist[v] = true;
					}
					else if (alt < (*handles[v]).distance)
					{
						prev[v] = make_tuple(u,w,pr);
						prob[v] = prob[u] * pr;
						if (prob[v] >= lb)
							lu.insert(v);
						heap.update(handles[v],node(v,alt));
					}
				}
				else
					prod *= (1 - pr);
			}
		}
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
		elapsed += ((double)(end - begin) / CLOCKS_PER_SEC);
	}
	else if (lu.empty())
		elapsed += ((double)(end - begin) / CLOCKS_PER_SEC);
	double pr = prob[t];
	delete [] handles;
	delete [] exist;
	delete [] visited;
	delete [] prev;
	delete [] prob;
	return make_tuple(p,min_dist,prod*pr,pr);
}

double approx_prob(Graph* g, vector< list<edge> > cp, list<edge> sp, double exist, double& elapsed)
{
	int C = 0, N = 10000, n = cp.size();
	list<edge>* diff = new list<edge>[n];
	vector<double> pr = vector<double>(n);
	double S = 0;
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
	random_device rd;
	mt19937 gen(rd());
	discrete_distribution<> d(pr.begin(), pr.end());
	clock_t begin = clock();
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
	elapsed += ((double)(end - begin) / CLOCKS_PER_SEC);
	delete [] diff;
	return (1 - C * S / N) * exist;
}

list<edge> mpsp(Graph* g, int s, int t, int m, double& candidate_time, double& prob_time)
{
	double lb_max = 0, p_max = 0;
	map< long,vector< tuple<list<edge>,double,double> > > paths = map< long,vector< tuple<list<edge>,double,double> > >();
	for (int i = 1; i <= m; i++)
	{
		list<edge> p;
		long w;
		double lb, ub;
		double start = candidate_time;
		tie(p,w,lb,ub) = prob_dijkstra(g,s,t,lb_max,candidate_time);
		double end = candidate_time;
		if (p.empty())
		{
			if (end == start)
				i--;
			continue;
		}
		map< long,vector< tuple<list<edge>,double,double> > >::iterator it = paths.find(w);
		if (it == paths.end())
		{
			vector< tuple<list<edge>,double,double> > ss = vector< tuple<list<edge>,double,double> >();
			ss.push_back(make_tuple(p,lb,ub));
			paths[w] = ss;
		}
		else
		{
			vector< tuple<list<edge>,double,double> > vv = it->second;
			bool f = true;
			for (tuple<list<edge>,double,double> tt : vv)
			{
				if (get<0>(tt) == p)
				{
					double& l = get<1>(tt);
					if (lb > l)
						l = lb;
					f = false;
					break;
				}
			}
			if (f)
				vv.push_back(make_tuple(p,lb,ub));
		}
	}
	list<edge> pp = list<edge>();
	vector< list<edge> > cp = vector< list<edge> >();
	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		long w = x->first;
		vector< tuple<list<edge>,double,double> > vv = x->second;
		vector< list<edge> > tmp = vector< list<edge> >();
		for (tuple<list<edge>,double,double> tt : vv)
		{
			list<edge> p = get<0>(tt);
			double lb = get<1>(tt);
			double ub = get<2>(tt);
			if (ub >= lb_max)
			{
				double prob = approx_prob(g,cp,p,ub,prob_time);
				tmp.push_back(p);
				if (prob >= lb && prob <= ub)
				{
					if (prob > p_max)
					{
						p_max = prob;
						pp = p;
					}
				}
			}
		}
		cp.insert(cp.end(),tmp.begin(),tmp.end());
	}
	cout << "Number of Candidate Paths : " << cp.size() << endl;
	return pp;
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

Graph generate_er(int n, int m)
{
	ofstream graph;
	graph.open("graph.txt");
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
			long w = rand();
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

Graph read_graph()
{
	ifstream graph;
	graph.open("graph.txt");
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

int main(int argc, char* argv[])
{
	ofstream queries;
	srand(time(NULL));
	// Graph g = generate_er(40000,80000);
	Graph g = read_graph();
	queries.open("queries.txt");
	for (int k = 1; k <= 4; k++)
	{
		cerr << "k = " << k << endl;
		cout << "Number of hops = " << k * 2 << endl << endl;
		long a_w = 0;
		double t_c = 0, t_p = 0, a_p = 0;
		for (int i = 1; i <= 100; i++)
		{
			cerr << "i = " << i << endl;
			int s = (int)((double)rand() / RAND_MAX * (g.n - 1));
			while (g.adj[s].empty())
				s = (int)((double)rand() / RAND_MAX * (g.n - 1));
			int t = bfs(&g, s, k * 2);
			if (t == -1)
			{
				i--;
				continue;
			}
			queries << s << " " << t << endl;
			cout << s << "\t" << t << endl;
			long wt = 0;
			double candidate_time = 0, prob_time = 0, pr = 1;
			list<edge> p = mpsp(&g, s, t, 1000, candidate_time, prob_time);
			for (edge e : p)
			{
				cout << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
				wt += get<2>(e);
				pr *= get<3>(e);
			}
			a_w += wt;
			a_p += pr;
			t_c += candidate_time;
			t_p += prob_time;
			cout << "Length of MPSP : " << wt << endl;
			cout << "Probability of MPSP : " << pr << endl;
			cout << "Candidate Generation Time : " << candidate_time << " seconds" << endl;
			cout << "Probability Computation Time : " << prob_time << " seconds" << endl;
			cout << endl;
		}
		cout << "Average Length of MPSP for " << k * 2 << " hops : " << a_w / 100 << endl;
		cout << "Average Probability of MPSP for " << k * 2 << " hops : " << a_p / 100 << endl;
		cout << "Average Candidate Generation Time for " << k * 2 << " hops : " << t_c / 100 << " seconds" << endl;
		cout << "Average Probability Computation Time for " << k * 2 << " hops : " << t_p / 100 << " seconds" << endl;
		cout << endl << endl;
	}
	queries.close();
	delete [] g.adj;
	return 0;
}