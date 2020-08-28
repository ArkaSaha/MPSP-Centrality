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
# include "io.h"

mt19937 mersenne_mpsp{static_cast<mt19937::result_type>(12345)};

# define NUM_MPSP 1
# define DIJKSTRA_RUNS 20 
# define LUBY_KARP_SAMPLES 1000


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

tuple<vector< list<edge> >,vector< list<edge> >,int,int,int,int> mpsp(AdjGraph* g, int s, int t, size_t k, int m, int N, double& candidate_time_prune, double& candidate_time_noprune, double& prob_time_prune, double& prob_time_noprune)
{
	double lb_kth_max = 0, p_max_p = 0, p_max_np = 0;
	int f_max = 1, n_s = 0, n_d = 0, n_r = m, n_p = 0;
	bool match_p = false, match_np = false;
	set<double> lbs = set<double>();
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
		lbs.insert(lb);
		size_t pos = 0;
		if (lbs.size() < k)
			lb_kth_max = *(lbs.begin());
		else
		{
			for (auto it = lbs.rbegin(); it != lbs.rend(); it++)
			{
				pos++;
				if (pos == k)
				{
					lb_kth_max = *it;
					break;
				}
			}
		}
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
	// list<edge> pp_p, pp_np = list<edge>();
	vector< pair<list<edge>,double> > cp_p = vector< pair<list<edge>,double> >(), cp_np = vector< pair<list<edge>,double> >();
	// map< int,set< list<edge> > > fp = map< int,set< list<edge> > >();

	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		vector< tuple<list<edge>,double,double,int> > vv = x->second;
		for (tuple<list<edge>,double,double,int> tt : vv)
		{
			list<edge> p = get<0>(tt);
			// double lb = get<1>(tt);
			double ub = get<2>(tt);
			int freq = get<3>(tt);
			double prob_np = approx_prob(cp_np,p,N,ub,prob_time_noprune);
			cp_np.push_back(make_pair(p,prob_np));
			// if (prob_np > p_max_np)
			// {
			// 	p_max_np = prob_np;
			// 	pp_np = p;
			// }
			if (ub < lb_kth_max)
				n_p++;
			else
			{
				double prob_p = approx_prob(cp_p,p,N,ub,prob_time_prune);
				cp_p.push_back(make_pair(p,prob_p));
				// if (prob_p < lb || prob_p > ub)
				// 	n_p += freq;
				// if (prob_p > p_max_p)
				// {
				// 	p_max_p = prob_p;
				// 	pp_p = p;
				// }
			}
			// if (fp.find(freq) == fp.end())
			// 	fp[freq] = set< list<edge> >();
			// fp[freq].insert(p);
			// if (freq > f_max)
			// 	f_max = freq;
		}
	}
	struct { bool operator() (pair<list<edge>,double> x, pair<list<edge>,double> y) { return x.second > y.second; } } comp;
	timespec t1, t2;
	clock_gettime(CLOCK_MONOTONIC,&t1);
	sort(cp_p.begin(),cp_p.end(),comp);
	clock_gettime(CLOCK_MONOTONIC,&t2);
	prob_time_prune += time_difference(t1,t2);
	clock_gettime(CLOCK_MONOTONIC,&t1);
	sort(cp_np.begin(),cp_np.end(),comp);
	clock_gettime(CLOCK_MONOTONIC,&t2);
	prob_time_prune += time_difference(t1,t2);

	vector< list<edge> > pp_p = vector< list<edge> >(), pp_np = vector< list<edge> >();
	for (size_t i = 0; i < min(k,cp_p.size()); i++)
		pp_p.push_back(cp_p[i].first);
	for (size_t i = 0; i < min(k,cp_np.size()); i++)
		pp_np.push_back(cp_np[i].first);

	// if (fp[f_max].find(pp_p) != fp[f_max].end())
	// 	match_p = true;
	// if (fp[f_max].find(pp_np) != fp[f_max].end())
	// 	match_np = true;

	return make_tuple(pp_p,pp_np,cp_p.size(),cp_np.size(),n_r,n_p);
}

tuple<vector< list<edge> >,vector< list<edge> >,int,int,int,int,double,double> mpsp_new(AdjGraph* g, int s, int t, size_t k, int m, int N, double& candidate_time_prune, double& candidate_time_noprune, double& prob_time_prune, double& prob_time_noprune)
{
	double lb_kth_max = 0, p_max_p = 0, p_max_np = 0;
	int f_max = 1, n_s = 0, n_d = 0, n_r = m, n_p = 0;
	bool match_p = false, match_np = false;
	set<double> lbs = set<double>();
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
		lbs.insert(lb);
		size_t pos = 0;
		if (lbs.size() < k)
			lb_kth_max = *(lbs.begin());
		else
		{
			for (auto it = lbs.rbegin(); it != lbs.rend(); it++)
			{
				pos++;
				if (pos == k)
				{
					lb_kth_max = *it;
					break;
				}
			}
		}
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
	// list<edge> pp_p, pp_np = list<edge>();
	vector< pair<list<edge>,double> > cp_p = vector< pair<list<edge>,double> >(), cp_np = vector< pair<list<edge>,double> >();
	// map< int,set< list<edge> > > fp = map< int,set< list<edge> > >();

	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		vector< tuple<list<edge>,double,double,int> > vv = x->second;
		for (tuple<list<edge>,double,double,int> tt : vv)
		{
			list<edge> p = get<0>(tt);
			// double lb = get<1>(tt);
			double ub = get<2>(tt);
			int freq = get<3>(tt);
			double prob_np = approx_prob(cp_np,p,N,ub,prob_time_noprune);
			cp_np.push_back(make_pair(p,prob_np));
			// if (prob_np > p_max_np)
			// {
			// 	p_max_np = prob_np;
			// 	pp_np = p;
			// }
			if (ub < lb_kth_max)
				n_p++;
			else
			{
				double prob_p = approx_prob(cp_p,p,N,ub,prob_time_prune);
				cp_p.push_back(make_pair(p,prob_p));
				// if (prob_p < lb || prob_p > ub)
				// 	n_p += freq;
				// if (prob_p > p_max_p)
				// {
				// 	p_max_p = prob_p;
				// 	pp_p = p;
				// }
			}
			// if (fp.find(freq) == fp.end())
			// 	fp[freq] = set< list<edge> >();
			// fp[freq].insert(p);
			// if (freq > f_max)
			// 	f_max = freq;
		}
	}
	struct { bool operator() (pair<list<edge>,double> x, pair<list<edge>,double> y) { return x.second > y.second; } } comp;
	timespec t1, t2;
	clock_gettime(CLOCK_MONOTONIC,&t1);
	sort(cp_p.begin(),cp_p.end(),comp);
	clock_gettime(CLOCK_MONOTONIC,&t2);
	prob_time_prune += time_difference(t1,t2);
	clock_gettime(CLOCK_MONOTONIC,&t1);
	sort(cp_np.begin(),cp_np.end(),comp);
	clock_gettime(CLOCK_MONOTONIC,&t2);
	prob_time_prune += time_difference(t1,t2);

	vector< list<edge> > pp_p = vector< list<edge> >(), pp_np = vector< list<edge> >();
	double prob_p = 0, prob_np = 0;
	size_t len = min(k,cp_p.size());
	for (size_t i = 0; i < len; i++)
	{
		pp_p.push_back(cp_p[i].first);
		prob_p += cp_p[i].second;
	}
	if (len)
		prob_p /= len;
	for (size_t i = 0; i < len; i++)
	{
		pp_np.push_back(cp_np[i].first);
		prob_np += cp_np[i].second;
	}
	if (len)
		prob_np /= len;

	// if (fp[f_max].find(pp_p) != fp[f_max].end())
	// 	match_p = true;
	// if (fp[f_max].find(pp_np) != fp[f_max].end())
	// 	match_np = true;

	return make_tuple(pp_p,pp_np,cp_p.size(),cp_np.size(),n_r,n_p,prob_p,prob_np);
}

vector<double> betweenness_naive(AdjGraph & g, ofstream& output)
{
    double _t1_p, _t1_np, _t2_p, _t2_np;

    vector<double> B = vector<double>(g.n, 0);

    timespec start;
    clock_gettime(CLOCK_MONOTONIC,&start);

    for(int s=0; s<g.n; s++)
    {
        timespec t1, t2;
        clock_gettime(CLOCK_MONOTONIC,&t1);
        for(int t=0; t<g.n; t++)
        {
            if(s == t) continue;
            //output << s << " " << t << endl;
            auto cur_mpsp = mpsp(&g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1_p, _t1_np, _t2_p, _t2_np);

            if(get<0>(cur_mpsp)[0].size() >= 2) // path consists of at least 2 edges
            {
                list<edge> p = get<0>(cur_mpsp)[0];
                /*
                for (edge e : p)
                	output << get<0>(e) << " " << get<1>(e) << " " << get<2>(e) << " " << get<3>(e) << endl;
                output << get<5>(cur_mpsp) << endl << endl;
                */

                // raise the betweenness of every inner node of the path by 1
                for(auto it = next(p.begin()); it != p.end(); it++){
                    B[get<0>(*it)]++;
                }
            }
        }
        clock_gettime(CLOCK_MONOTONIC,&t2);
        //cout << "last iteration : " << time_difference(t1,t2) << " seconds" << endl;
        //cout << "remaining : " << (g.n-(s+1)) * time_difference(start,t2)/(s+1) << endl;
    }

    // normalize betweenness by size of graph
    for(uint i=0; i<B.size(); i++)
    {
        //B[i] /= ((g.n-1) * (g.n-2));
        B[i] /= ((g.n-1) * (g.n)); // this is the same normalization as Riondato
    }

    return B;
}


vector<double> betweenness_sampling(AdjGraph & g, int samples, ofstream& output)
{
    double _t1, _t2, _t3, _t4;

    double r = samples;

    vector<double> B = vector<double>(g.n, 0);

    random_device rd;
    uniform_int_distribution<int> random_node(0, g.n-1);

    for(int sample = 0; sample < samples; sample++){
      // sample s and t
      int s = random_node(mersenne_mpsp);
      int t = random_node(mersenne_mpsp);

      if(s == t){ 
        sample--; 
        continue;
      }

      list<edge> mpsp_st = get<0>(mpsp(&g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, _t3, _t4))[0];

      if(mpsp_st.size() >= 2)  // the path consists of at least 2 edges
      {
        //output << s << " -> " << t << " : ";

        // raise the betweenness of every inner node of the path by 1/samples
        for(auto it = next(mpsp_st.begin()); it != mpsp_st.end(); it++){
          B[get<0>(*it)] += 1/r;
          //output << (get<0>(*it)) << " ";
        }
        //output << endl;
      }
    }
    //output << endl;

    return B;
}

vector<double> betweenness_hoeffding(AdjGraph &g, double epsilon, double delta, ofstream& output){
    double samples =  log(2*g.n / delta) / (2.0 * epsilon * epsilon);
    output << "samples : " << samples << endl;
    return betweenness_sampling(g, (int) samples, output);
}

vector<double> riondato(AdjGraph &g, double epsilon, double delta, ofstream& output){
  // output << "Riondato" << endl;
  double bound_on_VC = g.n; // Can we find a better bound on VC dimension?
  double c = 0.5;

  double samples = c/(epsilon * epsilon) *(floor(log2(bound_on_VC - 2)) + 1 + log(1/delta));

  // output << "Nr of samples : " << (int)samples  << " (vs " << (g.n * (g.n-1)) << ")"<< endl;

  return betweenness_sampling(g, (int) samples, output);
}

vector<double> betweenness_sampling_deterministic(Graph & g, int samples, ofstream& output)
{
    double _t1, _t2, _t3, _t4;

    double r = samples;

    vector<double> B = vector<double>(g.n, 0);

    random_device rd;
    uniform_int_distribution<int> random_node(0, g.n-1);

    for(int sample = 0; sample < samples; sample++){
      // sample s and t
      int s = random_node(mersenne_mpsp);
      int t = random_node(mersenne_mpsp);

      if(s == t){ 
        sample--; 
        continue;
      }

      Path sp = dijkstra(g, s, t);

      if(sp.edges.size() >= 2)  // the path consists of at least 2 edges
      {
        //output << s << " -> " << t << " : ";

        // raise the betweenness of every inner node of the path by 1/samples
        for(auto it = next(sp.edges.begin()); it != sp.edges.end(); it++){
          B[(*it).u] += 1/r;
          //output << (get<0>(*it)) << " ";
        }
        //output << endl;
      }
    }
    //output << endl;

    return B;
}


vector<double> exp_betweenness_with_riondato(Graph &g, double epsilon, double delta, ofstream & output){
    // output << "Expected betweenness with Riondato deterministic" << endl;

    int world_samples = ceil(2.0/(epsilon * epsilon) * log(2/delta));

    // output << "nr of worlds : " << world_samples << endl;

    random_device rd;
    uniform_real_distribution<double> coin(0.0, 1.0); 

    vector<double> B_total = vector<double>(g.n, 0);
    for(int world = 0; world < world_samples; world++){
        
        // world a world by flipping a coin for every edge
        for(auto &edge: g.index2edge){
            if(coin(mersenne_mpsp) < edge->p){
                edge->available = true;
            }
            else{
                edge->available = false;
            }
        }

        // compute betweenness for this deterministic world
        double bound_on_VC = g.n; // replace with bound on VD 
        double c = 0.5;
        double worlds = c/(epsilon * epsilon) *(floor(log2(bound_on_VC - 2)) + 1 + log(1/delta));
        vector<double> B = betweenness_sampling_deterministic(g, (int) worlds, output);

        for(int i=0; i <g.n; i++){
            B_total[i] += B[i];
        }
    }

    // normalize?
    for(auto & elt: B_total){
        elt /= world_samples;
    }

    return B_total;

}

vector< pair<int, double> > get_topk_from_betweenness(vector<double> B, int k){
  auto index_B = vector<pair<int, double>>(B.size(), pair<int, double>());

  for(uint i=0; i < B.size(); i++) index_B[i] = {i, B[i]};

  sort(index_B.begin(), index_B.end(), [](const pair<int, double> & a, const pair<int, double> & b){
      return a.second > b.second; });   // sort on second element (the betweeness) in decreasing order
  
  if(k > (int)index_B.size()) return index_B;
  
  return vector< pair<int, double> >(index_B.begin(), index_B.begin() + k);
}

vector<int> hedge(AdjGraph &g, double epsilon, int k, ofstream & output){
  // output << "Hedge" << endl;
  double samples = k * log(g.n) / (epsilon * epsilon);

  // output << "Nr of samples : " << samples << endl;

  double _t1, _t2, _t3, _t4;

  vector<double> B = vector<double>(g.n, 0);

  random_device rd;
  uniform_int_distribution<int> random_node(0, g.n-1);

  auto H = vector< vector<int> >();
  auto samples_containing_vertex = vector< vector<int> >(g.n, vector<int>());
  auto contained_in_nr_samples = vector<int>(g.n, 0);


  // sampling part
  for(int sample = 0; sample < samples; sample++){
    // sample s and t
    int s = random_node(mersenne_mpsp);
    int t = random_node(mersenne_mpsp);

    if(s == t){ 
      sample--; 
      continue;
    }

    list<edge> mpsp_st = get<0>(mpsp(&g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, _t3, _t4))[0];

    if(mpsp_st.size() >= 2)  // the path consists of at least 2 edges
    {

      auto h = vector<int>(mpsp_st.size() - 1);
      // raise the betweenness of every inner node of the path by 1/samples
      int index = 0;
      for(auto it = next(mpsp_st.begin()); it != mpsp_st.end(); it++){
        int u = get<0>(*it);
        h[index++] = u;
        contained_in_nr_samples[u]++;
        samples_containing_vertex[u].push_back(H.size());
      }
      H.push_back(h);
    }
  }

  vector<int> topk_nodes = vector<int>(k);

  vector<bool> not_deleted = vector<bool>(H.size(), true);

  // selecting topk node part
  for(int i=0; i<k; i++){
    auto max_elt = max_element(contained_in_nr_samples.begin(), contained_in_nr_samples.end());
    int v = distance(contained_in_nr_samples.begin(), max_elt); 

    topk_nodes[i] = v;

    // updating information about samples h in H
    for(const auto &h : samples_containing_vertex[v]){
      if(not_deleted[h]){
        not_deleted[h] = false;
        for(const auto &elt: H[h]){
          contained_in_nr_samples[elt]--;
        }
      }
    }

    // If we have picked a sample set it to -1, so it will never pop up as max_elt
    contained_in_nr_samples[v]--; 
  }

  return topk_nodes;

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

double jaccard(vector<pair<int, double>> const & list1, vector<pair<int, double>> const & list2, int k){
    auto set1 = set<int>(),  set2 = set<int>();
    k = min(min(k, (int)list1.size()), (int)list2.size());
    for(int i=0; i<k; i++){
        set1.insert(list1[i].first);
        set2.insert(list2[i].first);
    }
    auto is = set<int>(), un = set<int>();
    set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(is, is.begin()));
    set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(un, un.begin()));

    return (double) is.size() / (double)un.size();
}

double intersection_over_k(vector<pair<int, double>> const & list1, vector<pair<int, double>> const & list2, int k){
    auto set1 = set<int>(),  set2 = set<int>();
    k = min(min(k, (int)list1.size()), (int)list2.size());
    for(int i=0; i<k; i++){
        set1.insert(list1[i].first);
        set2.insert(list2[i].first);
    }
    auto is = set<int>(), un = set<int>();
    set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(is, is.begin()));

    return (double) is.size() / (double)k;
}


void finding_epsilon_delta(char* path_to_graph, char* path_to_output){
  AdjGraph g = read_graph(path_to_graph);
  ofstream output;
  output.open(path_to_output);
  //output << scientific;

  timespec t_naive_start, t_naive_end, t_start, t_end; 

  int repetitions = 3;
  vector<int> Ks = {5, 10, 20, 50, 100};

  clock_gettime(CLOCK_MONOTONIC, &t_naive_start);
  auto B_naive = betweenness_naive(g, output);
  auto topk_naive = get_topk_from_betweenness(B_naive, g.n);
  clock_gettime(CLOCK_MONOTONIC, &t_naive_end);


  output << endl;

  auto jac_res = map<int, double>();
  auto overk_res = map<int, double>();
  for(auto k : Ks){
    jac_res[k] = 0;
    overk_res[k] = 0;
  }
  
  clock_gettime(CLOCK_MONOTONIC, &t_start);
  for(int i=0; i < repetitions; i++){
      auto B_naive2 = betweenness_naive(g, output);
      auto topk_naive2 = get_topk_from_betweenness(B_naive2, g.n);

      for(auto k: Ks){
        jac_res[k] += jaccard(topk_naive, topk_naive2, k)/repetitions;
        overk_res[k] += intersection_over_k(topk_naive, topk_naive2, k)/repetitions;
      }
  }
  clock_gettime(CLOCK_MONOTONIC, &t_end);
  output << "computation took on average " << time_difference(t_start, t_end)/repetitions << " seconds" << endl;

  output << "Jaccard - intersection/k" << endl;
  for(auto k: Ks){
      output << k << " : " << jac_res[k] << " - " << overk_res[k] << endl;

  }

    /*
  output << "naive took " << time_difference(t_naive_start, t_naive_end) << " seconds" << endl << endl;
  for(auto k: {5, 10, 20, 50, 100}){
      output << k << " : " << jaccard(topk_naive, topk_naive2, k) << " - " << intersection_over_k(topk_naive, topk_naive2, k) << endl;

  }

*/
  output << endl;

  vector<pair<double, double>> parameters = {{0.05, 0.1}, {0.05, 0.05}, {0.05, 0.03}, {0.05, 0.01}, 
                                             {0.04, 0.1}, {0.04, 0.05}, {0.04, 0.03}, {0.04, 0.01},
                                             {0.03, 0.1}, {0.03, 0.05}, {0.03, 0.03}, {0.03, 0.01},
                                             {0.02, 0.1}, {0.02, 0.05}, {0.02, 0.03}, {0.02, 0.01},
                                             {0.01, 0.1}, {0.01, 0.05}, {0.01, 0.03}, {0.01, 0.01}};


  for(auto const & params : parameters){
      output << "epsilon : " << params.first << endl;
      output << "delta   : " << params.second << endl;

      auto jac_res = map<int, double>();
      auto overk_res = map<int, double>();
      for(auto k : Ks){
        jac_res[k] = 0;
        overk_res[k] = 0;
      }
      
      clock_gettime(CLOCK_MONOTONIC, &t_start);
      for(int i=0; i < repetitions; i++){
          auto B = betweenness_hoeffding(g, params.first, params.second, output);
          auto topk = get_topk_from_betweenness(B, g.n);

          for(auto k: Ks){
            jac_res[k] += jaccard(topk_naive, topk, k)/repetitions;
            overk_res[k] += intersection_over_k(topk_naive, topk, k)/repetitions;
          }
      }
      clock_gettime(CLOCK_MONOTONIC, &t_end);
      output << "computation took on average " << time_difference(t_start, t_end)/repetitions << " seconds" << endl;

      output << "Jaccard - intersection/k" << endl;
      for(auto k: Ks){
          output << k << " : " << jac_res[k] << " - " << overk_res[k] << endl;

      }

      output << endl <<  endl;
  }


    /*

  output << "naive took " << time_difference(t_naive_start, t_naive_end) << " seconds" << endl << endl;
  for(const auto &elt: topk_naive){
    output << elt.first << " " << elt.second << endl;
  }
  output << endl;

  output << "hoeffding (eps = 0.05, delta = 0.1) took " << time_difference(t_naive_start, t_naive_end) << " seconds" << endl << endl;
  for(const auto &elt: topk_005_01){
    output << elt.first << " " << elt.second << endl;
  }
  output << endl;
  */


}

void experiment_betweenness(char* path_to_graph, char* path_to_output, int k)
{
  AdjGraph g = read_graph(path_to_graph);

  Graph g2 = read_graph_from_file(path_to_graph);

  ofstream output;
  output.open(path_to_output);
  output << scientific;

  double epsilon = 0.03;
  double delta = 0.01;


  timespec t_riondato_start, t_riondato_end, t_riondato_det_start, t_riondato_det_end, t_hedge_start, t_hedge_end, t_naive_start, t_naive_end, t_h_start, t_h_end;

    /*
  clock_gettime(CLOCK_MONOTONIC,&t_riondato_start);
  auto B_riondato = riondato(g, epsilon, delta, output);
  auto topk_riondato = get_topk_from_betweenness(B_riondato, min(k, g.n));
  clock_gettime(CLOCK_MONOTONIC,&t_riondato_end);

  // output << "Riondato (epsilon = 0.05) took " << time_difference(t_riondato_start, t_riondato_end) << " seconds" << endl << endl;
  for(const auto &elt: topk_riondato){
    output << elt.first << " ";
  }
  output << endl;

    */

  clock_gettime(CLOCK_MONOTONIC,&t_h_start);
  auto B_hoeffding = betweenness_hoeffding(g, epsilon, delta,  output);
  auto topk_hoeffding = get_topk_from_betweenness(B_hoeffding, min(k, g.n));
  clock_gettime(CLOCK_MONOTONIC,&t_h_end);

  output << "Hoefdding (epsilon = " << epsilon << ", delta = " << delta << ") took " << time_difference(t_h_start, t_h_end) << " seconds" << endl << endl;
  for(const auto &elt: topk_hoeffding){
    output << elt.first << " " << elt.second << endl;
  }
  output << endl;

  clock_gettime(CLOCK_MONOTONIC,&t_riondato_det_start);
  auto B_riondato_det = exp_betweenness_with_riondato(g2, epsilon, delta, output);
  auto topk_riondato_det = get_topk_from_betweenness(B_riondato_det, min(k,g.n));
  clock_gettime(CLOCK_MONOTONIC,&t_riondato_det_end);

  output << "Expected betweenness with Riondato took " << time_difference(t_riondato_det_start, t_riondato_det_end) << " seconds" << endl << endl;
  for(const auto &elt: topk_riondato_det){
    output << elt.first << " " << elt.second << endl;
  }
  output << endl;


  if(g.n <= 1000){
      clock_gettime(CLOCK_MONOTONIC,&t_naive_start);
      auto B_naive = betweenness_naive(g, output);
      auto topk_naive = get_topk_from_betweenness(B_naive, min(k, g.n));
      clock_gettime(CLOCK_MONOTONIC,&t_naive_end);

      output << "Naive took " << time_difference(t_naive_start, t_naive_end) << " seconds" << endl << endl;
      for(const auto &elt: topk_naive){
          output << elt.first << " " << elt.second << endl;
      }
      output << endl;
  }

  output.close();

}

void experiment(char* path_to_graph, char* path_to_queries, char* path_to_output, size_t k = NUM_MPSP, int m = DIJKSTRA_RUNS, int N = LUBY_KARP_SAMPLES)
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

	for (int j = 1; j <= d; j++)
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
			vector< list<edge> > cp_p, cp_np;
			int c_p, c_np, f = 1, pruned;
			bool m_p = false, m_np = false;
			// tie(cp_p,cp_np,c_p,c_np,r,pruned) = mpsp(&g, s, t, k, m, N, candidate_time_prune, candidate_time_noprune, prob_time_prune, prob_time_noprune);
			tie(cp_p,cp_np,c_p,c_np,r,pruned,prob_p,prob_np) = mpsp_new(&g, s, t, k, m, N, candidate_time_prune, candidate_time_noprune, prob_time_prune, prob_time_noprune);
			for (auto p_p : cp_p)
			{
				if (! p_p.empty())
				{
					output << "With Pruning" << endl;
					num++;
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
					// Path path = Path({});
					// for (edge e : p_p)
					// {
					// 	Edge ee;
					// 	ee.u = get<0>(e);
					// 	ee.v = get<1>(e);
					// 	ee.l = (int)get<2>(e);
					// 	ee.p = get<3>(e);
					// 	for (Edge eee : G.adj[ee.u])
					// 	{
					// 		if (eee.v == ee.v)
					// 		{
					// 			ee.index = eee.index;
					// 			break;
					// 		}
					// 	}
					// 	path.edges.push_back(ee);
					// }
					// double prob = Luby_Karp(G, path, N);
					// prob_p += prob;
					// if (find(cp_np.begin(),cp_np.end(),p_p) != cp_np.end())
					// 	prob_np += prob;
				}
				else
					prob_p += dijkstra(&g,s,t);
			}
			for (auto p_np : cp_np)
			{
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
					// if (find(cp_p.begin(),cp_p.end(),p_np) == cp_p.end())
					// {
					// 	Path path = Path({});
					// 	for (edge e : p_np)
					// 	{
					// 		Edge ee;
					// 		ee.u = get<0>(e);
					// 		ee.v = get<1>(e);
					// 		ee.l = (int)get<2>(e);
					// 		ee.p = get<3>(e);
					// 		for (Edge eee : G.adj[ee.u])
					// 		{
					// 			if (eee.v == ee.v)
					// 			{
					// 				ee.index = eee.index;
					// 				break;
					// 			}
					// 		}
					// 		path.edges.push_back(ee);
					// 	}
					// 	prob_np += Luby_Karp(G, path, N);
					// }
				}
			}
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
		if (num)
		{
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
	if (argc != 4 && argc != 7)
	{
        cerr << "For mpsp experiment" << endl;
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-queries> <path-to-output> <k> <m> <N>" << endl;
        cerr << "For betweenness experiment" << endl;
		cerr << "Usage: ./mpsp <path-to-graph> <path-to-output> <k>" << endl;
		return EXIT_FAILURE;
	}
    // else if (string(argv[3]).find_first_not_of("0123456789") == std::string::npos)
    else if (argc == 4)
    {
        cout << "Doing betweenness experiment" << endl;
        experiment_betweenness(argv[1], argv[2], atoi(argv[3]));
        //finding_epsilon_delta(argv[1], argv[2]);
        
    }
    else
    {
        experiment(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
    }

	return EXIT_SUCCESS;
}
