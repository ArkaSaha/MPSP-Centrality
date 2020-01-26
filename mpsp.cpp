# include <iostream>
# include <cstdlib>
# include <climits>
# include <list>
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
	if (min_dist != 0 && prob[t] >= lb)
	{
		int v = t;
		while (v != s)
		{
			adj_ent res = prev[v];
			p.push_front(make_tuple(get<0>(res),v,get<1>(res),get<2>(res)));
			v = get<0>(res);
		}
		elapsed += ((double)(end - begin) / CLOCKS_PER_SEC);
	}
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
	vector<double> pr = vector<double>();
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
		pr.push_back(prob);
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

list<edge> mpsp(Graph* g, int s, int t, int m)
{
	double lb_max = 0;
	map< long,vector< tuple<list<edge>,double,double> > > paths = map< long,vector< tuple<list<edge>,double,double> > >();
	double candidate_time = 0;
	for (int i = 1; i <= m; i++)
	{
		list<edge> p;
		long w;
		double lb, ub;
		tie(p,w,lb,ub) = prob_dijkstra(g,s,t,lb_max,candidate_time);
		if (p.empty())
		{
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
	cout << "Candidate Generation Time : " << candidate_time << " seconds" << endl;
	list<edge> pp = list<edge>();
	double p_max = 0;
	vector< list<edge> > cp = vector< list<edge> >();
	double prob_time = 0;
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
	cout << "Probability Computation Time : " << prob_time << " seconds" << endl;
	return pp;
}

int main()
{
	srand(time(NULL));
	Graph g;
	g.n = 200000;
	g.m = 500000;
	g.adj = new vector<adj_ent>[g.n];
	for (int i=0; i<g.n; i++)
		g.adj[i] = vector<adj_ent>();
	for (int i=0; i<g.m; i++)
	{
		int v1 = (int)((double)rand() / RAND_MAX * g.n);
		int v2 = (int)((double)rand() / RAND_MAX * g.n);
		if (v1==v2)
		{
			i--;
			continue;
		}
		bool f = true;
		for (adj_ent e : g.adj[v1])
		{
			if (get<0>(e) == v2)
			{
				f = false;
				break;
			}
		}
		if (f)
			g.adj[v1].push_back(make_tuple(v2, rand(), (double)rand() / RAND_MAX));
		else
			i--;
	}
	int s = (int)((double)rand() / RAND_MAX * g.n);
	while (g.adj[s].empty())
		s = (int)((double)rand() / RAND_MAX * g.n);
	int u = get<0>(g.adj[s][(int)((double)rand() / RAND_MAX * g.adj[s].size())]);
	while (g.adj[u].empty())
		u = get<0>(g.adj[s][(int)((double)rand() / RAND_MAX * g.adj[s].size())]);
	int t = get<0>(g.adj[u][(int)((double)rand() / RAND_MAX * g.adj[u].size())]);
	while (t == s)
		t = get<0>(g.adj[u][(int)((double)rand() / RAND_MAX * g.adj[u].size())]);
	cout << s << "\t" << t << endl;
	list<edge> p = mpsp(&g, s, t, 10);
	for (edge e : p)
		cout << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
	delete [] g.adj;
	return 0;
}