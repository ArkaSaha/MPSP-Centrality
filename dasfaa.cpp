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

tuple< list<edge>,long,double,double > prob_dijkstra(AdjGraph* g, int s, int t, double& elapsed)
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
	double prod_ne = 1, prod_e = 1;
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
					prod_e *= pr;
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
				else
					prod_ne *= (1 - pr);
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double time = time_difference(begin,end);
	elapsed += time;
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
	return make_tuple(p, min_dist, pr, prod_ne * prod_e);
}

double approx_prob(vector< list<edge> > cp, list<edge> sp, double exist, double& elapsed)
{
	int C = 0, N = 1000, n = cp.size();
	auto diff = vector< list<edge> >(n);
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
	return (1 - C * S / N) * exist;
}

tuple<list<edge>,list<edge>,list<edge>,int,double,double,double> mpsp(AdjGraph* g, int s, int t, int m, double& candidate_time, double& prob_time)
{
	double spr_max_m = 0, spr_max_h = 0, p_max_p = 0, p_max_m = 0, p_max_h = 0;
	map< long,vector< tuple< list<edge>,double,vector<double> > > > paths = map< long,vector< tuple< list<edge>,double,vector<double> > > >();

	for (int i = 1; i <= m; i++)
	{
		list<edge> p;
		long w;
		double ub, pw;
		tie(p,w,ub,pw) = prob_dijkstra(g,s,t,candidate_time);
		if (p.empty())
			continue;
		map< long,vector< tuple< list<edge>,double,vector<double> > > >::iterator it = paths.find(w);
		if (it == paths.end())
		{
			vector< tuple< list<edge>,double,vector<double> > > ss = vector< tuple< list<edge>,double,vector<double> > >();
			vector<double> vect = vector<double>(1,pw);
			ss.push_back(make_tuple(p,ub,vect));
			paths[w] = ss;
		}
		else
		{
			vector< tuple< list<edge>,double,vector<double> > >& vv = it->second;
			bool f = true;
			for (tuple< list<edge>,double,vector<double> >& tt : vv)
			{
				if (get<0>(tt) == p)
				{
					get<2>(tt).push_back(pw);
					f = false;
					break;
				}
			}
			if (f)
			{
				vector<double> vect = vector<double>(1,pw);
				vv.push_back(make_tuple(p,ub,vect));
			}
		}
	}

	list<edge> pp_p = list<edge>(), pp_m = list<edge>(), pp_h = list<edge>();
	vector< list<edge> > cp = vector< list<edge> >();

	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		vector< tuple< list<edge>,double,vector<double> > > vv = x->second;
		for (tuple< list<edge>,double,vector<double> > tt : vv)
		{
			double spr_h = 0;
			list<edge> p = get<0>(tt);
			double ub = get<1>(tt);
			for (double pr : get<2>(tt))
				spr_h += (pr / (1 - pow(1 - pr, m)));
			double spr_m = ((double) get<2>(tt).size()) / m;
			double prob_p = approx_prob(cp,p,ub,prob_time);
			cp.push_back(p);
			if (prob_p > p_max_p)
			{
				p_max_p = prob_p;
				pp_p = p;
			}
			if (spr_m > spr_max_m)
			{
				spr_max_m = spr_m;
				pp_m = p;
				p_max_m = prob_p;
			}
			if (spr_h > spr_max_h)
			{
				spr_max_h = spr_h;
				pp_h = p;
				p_max_h = prob_p;
			}
		}
	}

	return make_tuple(pp_p,pp_m,pp_h,cp.size(),p_max_p,p_max_m,p_max_h);
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

void experiment(char* path_to_graph, char* path_to_queries, char* path_to_output)
{
	int d;
	AdjGraph g = read_graph(path_to_graph);
	ifstream queries;
	queries.open(path_to_queries);
	ofstream output;
	output.open(path_to_output);
	queries >> d;
	for (int k = 1; k <= d; k++)
	{
		int h, n;
		queries >> h >> n;
		output << "Number of hops = " << h << endl;
		output << "Number of queries = " << n << endl << endl;
		double t_c = 0, t_p = 0, n_b_m = 0, n_s_m = 0, n_w_m = 0, n_b_h = 0, n_s_h = 0, n_w_h = 0;
		int num = 0;
		for (int i = 1; i <= n; i++)
		{
			int s, t;
			queries >> s >> t;
			output << s << "\t" << t << endl;
			double candidate_time = 0, prob_time = 0, prob_p = 0, prob_m = 0, prob_h = 0;
			list<edge> p_p, p_m, p_h;
			int c_p;
			tie(p_p,p_m,p_h,c_p,prob_p,prob_m,prob_h) = mpsp(&g, s, t, 20, candidate_time, prob_time);
			if (! p_p.empty())
			{
				num++;
				output << "Luby-Karp" << endl;
                for (edge e : p_p)
					output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
				t_c += candidate_time;
				t_p += prob_time;
			}
			else
				prob_p = dijkstra(&g,s,t);
			if (! p_m.empty())
			{
				output << "Majority" << endl;
                for (edge e : p_m)
					output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
			}
			else
				prob_m = prob_p;
			if (! p_h.empty())
			{
				output << "Horvitz-Thompson" << endl;
                for (edge e : p_h)
					output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
			}
			else
				prob_h = prob_p;
			if (!p_p.empty())
			{
				if (prob_p > prob_m)
					n_b_m++;
				else if (prob_p == prob_m)
					n_s_m++;
				else
					n_w_m++;
				if (prob_p > prob_h)
					n_b_h++;
				else if (prob_p == prob_h)
					n_s_h++;
				else
					n_w_h++;
			}
			output << "Number of Distinct Candidate Paths : " << c_p << endl;
			output << "SP Probability (Luby-Karp) : " << prob_p << endl;
			output << "SP Probability (Majority) : " << prob_m << endl;
			output << "SP Probability (Horvitz-Thompson) : " << prob_h << endl;
			output << "Candidate Generation Time : " << candidate_time << " seconds" << endl;
			output << "Probability Computation Time : " << prob_time << " seconds" << endl;
			output << endl;
		}

		output << "Number of Non-Empty Paths for " << h << " hops : " << num << endl;
		if (num)
		{
			output << "Luby-Karp vs Majority for " << h << " hops : " << n_b_m / num << " " << n_s_m / num << " " << n_w_m / num << endl;
			output << "Luby-Karp vs Horvitz-Thompson for " << h << " hops : " << n_b_h / num << " " << n_s_h / num << " " << n_w_h / num << endl;
			output << "Average Candidate Generation Time for " << h << " hops : " << t_c / num << " seconds" << endl;
			output << "Average Probability Computation Time for " << h << " hops : " << t_p / num << " seconds" << endl;
			output << "Average Total Time for " << h << " hops : " << (t_c + t_p) / num << " seconds" << endl;
		}
		output << endl << endl;
	}
	queries.close();
	output.close();
	delete [] g.adj;
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 4)
	{
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-queries> <path-to-output>" << endl;
		return EXIT_FAILURE;
	}
    else
    {
        experiment(argv[1], argv[2], argv[3]);
    }
	return EXIT_SUCCESS;
}
