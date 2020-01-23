# include <iostream>
# include <cstdlib>
# include <climits>
# include <list>
# include <map>
# include <set>
# include <algorithm>
# include <ctime>
# include <boost/heap/fibonacci_heap.hpp>

using namespace std;
using namespace boost::heap;

typedef tuple<string,long,double> adj_ent;
typedef tuple<string,string,long,double> edge;

struct Graph
{
	int n, m;
	set<string> vertices;
	map< string,list<adj_ent> > adj;
};

tuple< list<edge>,long,double,double > prob_dijkstra(Graph* g, string s, string t, double lb)
{
	struct node
	{
		string vertex;
		long distance;
		node(const string& v, long d) : vertex(v), distance(d) {}
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
	map<string,handle_t> handles = map<string,handle_t>();
	map<string,adj_ent> prev = map<string,adj_ent>();
	map<string,double> prob = map<string,double>();
	set<string> lu = set<string>();
	for (string v : g->vertices)
		handles[v] = heap.push(node(v,LONG_MAX));
	heap.update(handles[s],node(s,0));
	prob[s] = 1;
	prob[t] = 1;
	long min_dist = 0;
	double prod = 1;
	list<edge> p = list<edge>();
	lu.insert(s);
	while (! (heap.empty() || lu.empty()))
	{
		node n = heap.top();
		string u = n.vertex;
		long d = n.distance;
		heap.pop();
		handles.erase(u);
		lu.erase(u);
		if (u == t || d == LONG_MAX)
		{
			min_dist = d;
			break;
		}
		for (adj_ent e : g->adj[u])
		{
			string v = get<0>(e);
			double r = (double)rand() / RAND_MAX;
			if (r < get<2>(e))
			{
				long alt = d + get<1>(e);
				if (alt < (*handles[v]).distance)
				{
					prev[v] = make_tuple(u,get<1>(e),get<2>(e));
					prob[v] = prob[u] * get<2>(e);
					if (prob[v] >= lb)
						lu.insert(v);
					heap.update(handles[v],node(v,alt));
				}
			}
			else
				prod *= (1 - get<2>(e));
		}
	}
	if (min_dist != LONG_MAX && min_dist != 0 && prob[t] >= lb)
	{
		string v = t;
		while (v != s)
		{
			adj_ent res = prev[v];
			p.push_front(make_tuple(get<0>(res),v,get<1>(res),get<2>(res)));
			v = get<1>(res);
		}	
	}
	return make_tuple(p,min_dist,prod*prob[t],prob[t]);
}

struct path
{
	list<edge> edges;
	long weight;
	path(const list<edge>& e, long w) : edges(e), weight(w) {}
};

struct compare_path
{
	bool operator()(const path& p1, const path& p2) const
	{
		return p1.weight < p2.weight;
	}
};

double approx_prob(Graph* g, set<path,compare_path>::iterator first, set<path,compare_path>::iterator last, double exist)
{
	int C = 0, N = 1000000;
	path pp = *last;
	map< list<edge>,list<edge> > diff = map< list<edge>,list<edge> >();
	map< list<edge>,double > pr = map< list<edge>,double >();
	double S = 0;
	for (set<path,compare_path>::iterator x = first; x != last; x++)
	{
		path p = *x;
		pr[p.edges] = 1;
		diff[p.edges] = list<edge>();
		set_difference(p.edges.begin(),p.edges.end(),pp.edges.begin(),pp.edges.end(),diff[p.edges].begin());
		for (edge e : diff[p.edges])
			pr[p.edges] *= get<3>(e);
		S += pr[p.edges];
	}
	for (int k = 1; k <= N; k++)
	{
		set<path,compare_path>::iterator y = first;
		double r = (double)rand() / RAND_MAX;
		for (set<path,compare_path>::iterator x = first; x != last; x++)
		{
			path p = *x;
			if (r < pr[p.edges]/S)
			{
				y = x;
				break;
			}
		}
		bool f1 = false;
		for (set<path,compare_path>::iterator x = first; x != y; x++)
		{
			bool f2 = true;
			path p = *x;
			for (edge e : diff[p.edges])
			{
				r = (double)rand() / RAND_MAX;
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
	return (1 - C * S / N) * exist;
}

list<edge> mpsp(Graph* g, string s, string t, int m)
{
	set<path,compare_path> cp = set<path,compare_path>();
	double lb_max = 0;
	map< list<edge>,double > lbs = map< list<edge>,double >();
	map< list<edge>,double > ubs = map< list<edge>,double >();
	for (int i = 1; i <= m; i++)
	{
		list<edge> p;
		long w;
		double lb, ub;
		tie(p,w,lb,ub) = prob_dijkstra(g,s,t,lb_max);
		if (w == LONG_MAX)
			i--;
		if (p.empty())
			continue;
		if (lb > lb_max)
			lb_max = lb;
		if (lbs.find(p) == lbs.end() || lbs[p] < lb)
			lbs[p] = lb;
		if (ubs.find(p) == ubs.end() || ubs[p] > ub)
			ubs[p] = ub;
		cp.insert(path(p,w));
	}
	map< list<edge>,double > pr = map< list<edge>,double >();
	path pp = path(list<edge>(),0);
	double p_max = 0;
	for (set<path,compare_path>::iterator x = cp.begin(); x != cp.end(); x++)
	{
		path p = *x;
		pr[p.edges] = approx_prob(g,cp.begin(),x,ubs[p.edges]);
		if (pr[p.edges] > p_max)
		{
			p_max = pr[p.edges];
			pp = p;
		}
	}
	return pp.edges;
}

int main()
{
	srand(time(NULL));
	Graph g;
	g.n = 100;
	g.m = 200;
	g.vertices = set<string>();
	g.adj = map< string,list<adj_ent> >();
	for (int i=0; i<g.n; i++)
	{
		g.vertices.insert(to_string(i));
		g.adj[to_string(i)] = list<adj_ent>();
	}
	cout << "vertices" << endl;
	for (int i=0; i<g.m; i++)
	{
		int v1 = (int)((double)rand() / RAND_MAX * g.n);
		int v2 = (int)((double)rand() / RAND_MAX * g.n);
		if (v1==v2)
		{
			i--;
			continue;
		}
		g.adj[to_string(v1)].push_back(make_tuple(to_string(v2), rand(), (double)rand() / RAND_MAX));
	}
	cout << "edges" << endl;
	int s = (int)((double)rand() / RAND_MAX * g.n);
	int t = (int)((double)rand() / RAND_MAX * g.n);
	cout << s << "\t" << t << endl;
	list<edge> p = mpsp(&g, to_string(s), to_string(t), 10000);
	for (edge e : p)
		cout << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
	return 0;
}