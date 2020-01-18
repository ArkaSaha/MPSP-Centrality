# include <iostream>
# include <cstdlib>
# include <climits>
# include <list>
# include <map>
# include <boost/heap/fibonacci_heap.hpp>

using namespace std;
using namespace boost::heap;

typedef tuple<string,long,double> edge;

struct Graph
{
	int n, m;
	list<string> vertices;
	map< string,list< edge > > adj;
};

tuple< list<string>,long > prob_dijkstra(Graph* g, string s, string t)
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
	map<string,string> prev = map<string,string>();
	for (string v : g->vertices)
	{
		prev[v] = "";
		handles[v] = heap.push(node(v,LONG_MAX));
	}
	heap.update(handles[s],node(s,0));
	long min_dist = 0;
	while (! heap.empty())
	{
		node n = heap.top();
		string u = n.vertex;
		long d = n.distance;
		heap.pop();
		handles.erase(u);
		if (u == t || d == LONG_MAX)
		{
			min_dist = d;
			break;
		}
		for (edge e : g->adj[u])
		{
			string v = get<0>(e);
			if (handles.find(v) != handles.end())
			{
				double r = (double)rand()/RAND_MAX;
				if (r < get<2>(e))
				{
					long alt = d + get<1>(e);
					if (alt < (*handles[v]).distance)
					{
						prev[v] = u;
						heap.update(handles[v],node(v,alt));
					}
				}
			}
		}
	}
	list<string> p = list<string>();
	if (min_dist != LONG_MAX)
	{
		p.push_front(t);
		string v = t;
		while (v != s)
		{
			v = prev[v];
			p.push_front(v);
		}	
	}
	return make_tuple(p,min_dist);
}

int main()
{
	return 0;
}