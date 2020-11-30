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
# include <pthread.h>

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

tuple< list<edge>,long,double > prob_dijkstra(AdjGraph* g, int s, int t, double& elapsed, mt19937& gen)
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
	timespec begin, end;
	clock_gettime(CLOCK_MONOTONIC,&begin);
    uniform_real_distribution<double> coin(0.0, 1.0);
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
				double r = coin(gen);
                double pr = get<2>(e);
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
	return make_tuple(p, min_dist, pr);
}

double approx_prob(vector< pair<list<edge>,double> > cp, list<edge> sp, int N, double exist, double& elapsed, mt19937& gen)
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
	discrete_distribution<> d(pr.begin(), pr.end());
    uniform_real_distribution<double> coin(0.0, 1.0);
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
                    double r = coin(gen);
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

tuple<vector< list<edge> >,int,double> mpsp(AdjGraph* g, int s, int t, size_t k, int m, int N, double& candidate_time, double& prob_time, mt19937& gen)
{
	map< long,vector< tuple<list<edge>,double> > > paths = map< long,vector< tuple<list<edge>,double> > >();

	for (int i = 1; i <= m; i++)
	{
		list<edge> p;
		long w;
		double ub;
		tie(p,w,ub) = prob_dijkstra(g,s,t,candidate_time, gen);
		if (p.empty())
			continue;
		map< long,vector< tuple<list<edge>,double> > >::iterator it = paths.find(w);
		if (it == paths.end())
		{
			vector< tuple<list<edge>,double> > ss = vector< tuple<list<edge>,double> >();
			ss.push_back(make_tuple(p,ub));
			paths[w] = ss;
		}
		else
		{
			vector< tuple<list<edge>,double> >& vv = it->second;
			bool f = true;
			for (tuple<list<edge>,double>& tt : vv)
			{
				if (get<0>(tt) == p)
				{
					f = false;
					break;
				}
			}
			if (f)
				vv.push_back(make_tuple(p,ub));
		}
	}

	vector< pair<list<edge>,double> > cp = vector< pair<list<edge>,double> >();

	for (auto x = paths.begin(); x != paths.end(); x++)
	{
		vector< tuple<list<edge>,double> > vv = x->second;
		for (tuple<list<edge>,double> tt : vv)
		{
			list<edge> p = get<0>(tt);
			double ub = get<1>(tt);
			double prob = approx_prob(cp,p,N,ub,prob_time, gen);
			cp.push_back(make_pair(p,prob));
		}
	}
	struct { bool operator() (pair<list<edge>,double> x, pair<list<edge>,double> y) { return x.second > y.second; } } comp;
	timespec t1, t2;
	clock_gettime(CLOCK_MONOTONIC,&t1);
	sort(cp.begin(),cp.end(),comp);
	clock_gettime(CLOCK_MONOTONIC,&t2);
	prob_time += time_difference(t1,t2);

	vector< list<edge> > pp = vector< list<edge> >();
	double prob = 0;
	size_t len = min(k,cp.size());
	for (size_t i = 0; i < len; i++)
	{
		pp.push_back(cp[i].first);
		prob += cp[i].second;
	}
	if (len)
		prob /= len;

	return make_tuple(pp,cp.size(),prob);
}

vector<double> betweenness_naive(AdjGraph * g)
{
    double _t1, _t2;

    vector<double> B = vector<double>(g->n, 0);

    timespec start;
    clock_gettime(CLOCK_MONOTONIC,&start);

    random_device rd;
    mt19937 gen(rd());

    for(int s=0; s<g->n; s++)
    {
        timespec t1, t2;
        clock_gettime(CLOCK_MONOTONIC,&t1);
        for(int t=0; t<g->n; t++)
        {
            if(s == t) continue;
            auto cur_mpsp = mpsp(g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, gen);

            if(get<0>(cur_mpsp).size() > 0 && get<0>(cur_mpsp)[0].size() >= 2) // path consists of at least 2 edges
            {
                list<edge> p = get<0>(cur_mpsp)[0];
                // raise the betweenness of every inner node of the path by 1
                for(auto it = next(p.begin()); it != p.end(); it++){
                    B[get<0>(*it)]++;
                }
            }
        }
        clock_gettime(CLOCK_MONOTONIC,&t2);
    }

    // normalize betweenness by size of graph
    for(uint i=0; i<B.size(); i++)
    {
        B[i] /= ((g->n-1) * (g->n));
    }

    return B;
}

struct thread_data {
    int start, end;
    int samples;
    AdjGraph *g;
    vector<double> B;
};

void *betweenness_naive_pthread_helper(void *arg){
    thread_data *tdata = (thread_data *) arg;
    double _t1, _t2;
    random_device rd;
    mt19937 gen(rd());

    for(int s = tdata->start; s < tdata->end; s++){
        for(int t=0; t<tdata->g->n; t++)
        {
            if(s == t) continue;
            //auto cur_mpsp = mpsp(tdata->g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, gen);
            auto cur_mpsp = mpsp(tdata->g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, gen);

            if(get<0>(cur_mpsp).size() > 0 && get<0>(cur_mpsp)[0].size() >= 2) // path consists of at least 2 edges
            {
                list<edge> p = get<0>(cur_mpsp)[0];
                // raise the betweenness of every inner node of the path by 1
                for(auto it = next(p.begin()); it != p.end(); it++){
                    tdata->B[get<0>(*it)]++;
                }
            }
        }
    }
    pthread_exit(NULL);
}

vector<double> betweenness_naive_pthread(AdjGraph * g, int nr_threads)
{
    vector<double> B = vector<double>(g->n, 0);

    auto tdata = vector<thread_data>(nr_threads);
    int nodes_per_block = g->n/nr_threads;
    if(g->n % nr_threads != 0) nodes_per_block++;
    auto tid = vector<pthread_t>(nr_threads);
    int start = 0;
    for(int i=0; i<nr_threads; i++){
        tdata[i].g = g;
        tdata[i].B = vector<double>(g->n, 0);
        tdata[i].start = start;
        tdata[i].end = min(start + nodes_per_block, g->n);
        start += nodes_per_block;
        //cerr << "starting thread with s = " << tdata[i].start << " -> " << tdata[i].end << endl;
        int rc = pthread_create(&tid[i], NULL, betweenness_naive_pthread_helper, (void *)(&tdata[i]));
        if(rc){
            cerr << "ERROR : pthread_create, rc : " + rc << endl;
        }
    }

    for(int i=0; i<nr_threads;i++){
        pthread_join(tid[i], NULL);
    }

    for(int i=0;i <nr_threads;i++){
        for(uint j=0; j<B.size(); j++){
            B[j] += tdata[i].B[j];
        }
    }


    // normalize betweenness by size of graph
    for(uint i=0; i<B.size(); i++)
    {
        B[i] /= ((g->n-1) * (g->n));
    }

    return B;
}


vector<double> betweenness_sampling(AdjGraph * g, int samples)
{
    double _t1, _t2;

    double r = samples;

    vector<double> B = vector<double>(g->n, 0);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> random_node(0, g->n-1);

    for(int sample = 0; sample < samples; sample++){
      // sample s and t
      int s = random_node(gen);
      int t = random_node(gen);

      if(s == t){
        sample--;
        continue;
      }

      auto cur_mpsp = mpsp(g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, gen);

      if(get<0>(cur_mpsp).size() > 0 && get<0>(cur_mpsp)[0].size() >= 2) // path consists of at least 2 edges
      {
          list<edge> p = get<0>(cur_mpsp)[0];
          // raise the betweenness of every inner node of the path by 1/r
          for(auto it = next(p.begin()); it != p.end(); it++){
              B[get<0>(*it)] += 1/r;
          }
      }
    }

    return B;
}

vector<double> betweenness_hoeffding(AdjGraph *g, double epsilon, double delta, ofstream& output){
    double samples =  log(2*g->n / delta) / (2.0 * epsilon * epsilon);
    output << "samples : " << samples << endl;
    return betweenness_sampling(g, (int) samples);
}


void *betweenness_sampling_pthread_helper(void *arg)
{
    thread_data *tdata = (thread_data *) arg;
    timespec t1, t2;
    double _t1, _t2;

    double total_time = 0;
    pthread_t pid = pthread_self();

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> random_node(0, tdata->g->n-1);

    for(int sample = 1; sample <= tdata->samples; sample++){
      clock_gettime(CLOCK_MONOTONIC,&t1);
      // sample s and t
      int s = random_node(gen);
      int t = random_node(gen);

      if(s == t){
        sample--;
        continue;
      }

      auto cur_mpsp = mpsp(tdata->g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, gen);

      if(get<0>(cur_mpsp).size() > 0 && get<0>(cur_mpsp)[0].size() >= 2) // path consists of at least 2 edges
      {
          list<edge> p = get<0>(cur_mpsp)[0];
          // raise the betweenness of every inner node of the path by 1
          for(auto it = next(p.begin()); it != p.end(); it++){
              tdata->B[get<0>(*it)]++;
          }
      }
      clock_gettime(CLOCK_MONOTONIC, &t2);
      total_time += time_difference(t1, t2);
      double avg_time = total_time/sample;
      double remaining_time = ((tdata->samples) - sample) * avg_time;
      cerr << pid << " : (total =) " << total_time << " / (sample =) " << sample << " = " << avg_time;
      cerr << " || remaining = " << remaining_time << endl;

      
    }
    pthread_exit(NULL);
}

vector<double> betweenness_hoeffding_pthread(AdjGraph *g, double epsilon, double delta, ofstream& output, int nr_threads){
    int samples = (int) ( log(2*g->n / delta) / (2.0 * epsilon * epsilon));
    int samples_per_thread = samples/nr_threads;
    if(samples %nr_threads != 0) samples_per_thread++;
    output << "samples : " << samples << endl;
    output << "samples/thread : " << samples_per_thread << endl;

    auto tdata = vector<thread_data>(nr_threads);
    auto tid = vector<pthread_t>(nr_threads);
    int start = 0;
    for(int i=0; i<nr_threads; i++){
        tdata[i].g = g;
        tdata[i].samples = samples_per_thread;
        tdata[i].B = vector<double>(g->n, 0);
        //cerr << "starting thread #" << i << endl;
        int rc = pthread_create(&tid[i], NULL, betweenness_sampling_pthread_helper, (void *)(&tdata[i]));
        if(rc){
            cerr << "ERROR : pthread_create, rc : " + rc << endl;
        }
        else{
            cerr << "Created thread #" << i << endl;
        }
    }

    for(int i=0; i<nr_threads;i++){
        pthread_join(tid[i], NULL);
    }

    auto B = vector<double>(g->n, 0);

    for(int i=0;i <nr_threads;i++){
        for(uint j=0; j<B.size(); j++){
            B[j] += tdata[i].B[j];
        }
    }

    double r = samples_per_thread * nr_threads;

    // normalize by nr of samples 
    for(uint i=0; i<B.size(); i++)
    {
        B[i] /= r;
    }

    return B;
}

vector<double> riondato(AdjGraph *g, double epsilon, double delta){
  double bound_on_VC = g->n; // Can we find a better bound on VC dimension?
  double c = 0.5;

  double samples = c/(epsilon * epsilon) *(floor(log2(bound_on_VC - 2)) + 1 + log(1/delta));

  return betweenness_sampling(g, (int) samples);
}

vector<double> betweenness_sampling_deterministic(Graph & g, int samples)
{
    double r = samples;

    vector<double> B = vector<double>(g.n, 0);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> random_node(0, g.n-1);

    for(int sample = 0; sample < samples; sample++){
      // sample s and t
      int s = random_node(gen);
      int t = random_node(gen);

      if(s == t){
        sample--;
        continue;
      }

      Path sp = dijkstra(g, s, t);

      if(sp.edges.size() >= 2)  // the path consists of at least 2 edges
      {

        // raise the betweenness of every inner node of the path by 1/samples
        for(auto it = next(sp.edges.begin()); it != sp.edges.end(); it++){
          B[(*it).u] += 1/r;
        }
      }
    }

    return B;
}


vector<double> exp_betweenness_with_riondato(Graph &g, double epsilon, double delta){

    int world_samples = ceil(2.0/(epsilon * epsilon) * log(2/delta));

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> coin(0.0, 1.0);

    vector<double> B_total = vector<double>(g.n, 0);
    for(int world = 0; world < world_samples; world++){

        // sample a world by flipping a coin for every edge
        for(auto &edge: g.index2edge){
            if(coin(gen) < edge->p){
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
        vector<double> B = betweenness_sampling_deterministic(g, (int) worlds);

        for(int i=0; i <g.n; i++){
            B_total[i] += B[i];
        }
    }

    // normalize
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

vector<int> hedge(AdjGraph &g, double epsilon, int k){
  double samples = k * log(g.n) / (epsilon * epsilon);

  double _t1, _t2;

  vector<double> B = vector<double>(g.n, 0);

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> random_node(0, g.n-1);

  auto H = vector< vector<int> >();
  auto samples_containing_vertex = vector< vector<int> >(g.n, vector<int>());
  auto contained_in_nr_samples = vector<int>(g.n, 0);


  // sampling part
  for(int sample = 0; sample < samples; sample++){
    // sample s and t
    int s = random_node(gen);
    int t = random_node(gen);

    if(s == t){
      sample--;
      continue;
    }

    auto cur_mpsp = mpsp(&g, s, t, NUM_MPSP, DIJKSTRA_RUNS, LUBY_KARP_SAMPLES, _t1, _t2, gen);

    if(get<0>(cur_mpsp).size() > 0 && get<0>(cur_mpsp)[0].size() >= 2) // path consists of at least 2 edges
    {
      auto mpsp_st = get<0>(cur_mpsp)[0];
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


void experiment_betweenness(char* path_to_graph, char* path_to_output, int k)
{
  AdjGraph g = read_graph(path_to_graph);

  //Graph g2 = read_graph_from_file(path_to_graph);

  ofstream output;
  output.open(path_to_output);
  output << scientific;

  double epsilon = 0.05;
  double delta = 0.1;

  timespec t_riondato_det_start, t_riondato_det_end, t_naive_start, t_naive_end, t_h_start, t_h_end;

  /*

  clock_gettime(CLOCK_MONOTONIC,&t_riondato_det_start);
  auto B_riondato_det = exp_betweenness_with_riondato(g2, epsilon, delta);
  auto topk_riondato_det = get_topk_from_betweenness(B_riondato_det, min(k,g.n));
  clock_gettime(CLOCK_MONOTONIC,&t_riondato_det_end);

  output << "Expected betweenness with Riondato took " << time_difference(t_riondato_det_start, t_riondato_det_end) << " seconds" << endl << endl;
  for(const auto &elt: topk_riondato_det){
    output << elt.first << " " << elt.second << endl;
  }
  output << endl;
  */


/*

  clock_gettime(CLOCK_MONOTONIC,&t_h_start);
  auto B_hoeffding = betweenness_hoeffding(&g, epsilon, delta,  output);
  auto topk_hoeffding = get_topk_from_betweenness(B_hoeffding, min(k, g.n));
  clock_gettime(CLOCK_MONOTONIC,&t_h_end);
  double serial = time_difference(t_h_start, t_h_end);

  output << "Hoefdding (epsilon = " << epsilon << ", delta = " << delta << ") took " << serial << " seconds" << endl;
  for(const auto &elt: topk_hoeffding){
      output << elt.first << " " << elt.second << endl;
  }
  output << endl;
  */

  // PARALLEL
  for(int nr_threads=8; nr_threads <= 8; nr_threads *= 2){
      output << "Parallel Sampling #threads = " << nr_threads << endl;
      clock_gettime(CLOCK_MONOTONIC,&t_naive_start);
      //auto B_naive_p = betweenness_naive_pthread(&g, nr_threads);
      auto B_naive_p = betweenness_hoeffding_pthread(&g, epsilon, delta, output, nr_threads);
      auto topk_naive_p = get_topk_from_betweenness(B_naive_p, min(k, g.n));
      clock_gettime(CLOCK_MONOTONIC,&t_naive_end);
      double parallel = time_difference(t_naive_start, t_naive_end);

      output << "Time elapsed : " << parallel << " seconds" << endl;
      //output << "Speedup : " << serial / parallel << endl;
      for(const auto &elt: topk_naive_p){
          output << elt.first << " " << elt.second << endl;
      }
      output << endl;
  }


/*
  if(g.n <= 10){
      clock_gettime(CLOCK_MONOTONIC,&t_naive_start);
      auto B_naive = betweenness_naive(&g);
      auto topk_naive = get_topk_from_betweenness(B_naive, min(k, g.n));
      clock_gettime(CLOCK_MONOTONIC,&t_naive_end);

      output << "Naive took " << time_difference(t_naive_start, t_naive_end) << " seconds" << endl << endl;
      for(const auto &elt: topk_naive){
          output << elt.first << " " << elt.second << endl;
      }
      output << endl;
  }
  */

  output.close();

}

void experiment(char* path_to_graph, char* path_to_queries, char* path_to_output, size_t k = NUM_MPSP, int m = DIJKSTRA_RUNS, int N = LUBY_KARP_SAMPLES)
{
	int d;
	AdjGraph g = read_graph(path_to_graph);
	ifstream queries;
	queries.open(path_to_queries);
	ofstream output;
	output.open(path_to_output);
	queries >> d;

    random_device rd;
    mt19937 gen(rd());
	for (int j = 1; j <= d; j++)
	{
		int h, n;
		queries >> h >> n;
		output << "Number of hops = " << h << endl;
		output << "Number of queries = " << n << endl << endl;
		double t_c = 0, t_p = 0, a_pr = 0;
		int num = 0;
		for (int i = 1; i <= n; i++)
		{
			int s, t;
			queries >> s >> t;
			output << s << "\t" << t << endl;
			double candidate_time = 0, prob_time = 0, prob = 0;
			vector< list<edge> > cp;
			int nc;
			tie(cp,nc,prob) = mpsp(&g, s, t, k, m, N, candidate_time, prob_time, gen);
			bool f = false;
			for (auto p : cp)
			{
				if (! p.empty())
				{
					f = true;
					output << "Path" << endl;
					for (edge e : p)
						output << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\t" << get<3>(e) << endl;
				}
				else
					prob += dijkstra(&g,s,t);
			}
			if (f)
			{
				num++;
				t_c += candidate_time;
				t_p += prob_time;
				a_pr += prob;
			}
			output << "Number of Distinct Candidate Paths : " << nc << endl;
			output << "SP Probability : " << prob << endl;
			output << "Candidate Generation Time : " << candidate_time << " seconds" << endl;
			output << "Probability Computation Time : " << prob_time << " seconds" << endl;
			output << endl;
		}
		output << "Number of Non-Empty Paths for " << h << " hops : " << num << endl;
		if (num)
		{
			output << "Average Candidate Generation Time for " << h << " hops : " << t_c / num << " seconds" << endl;
			output << "Average Probability Computation Time for " << h << " hops : " << t_p / num << " seconds" << endl;
			output << "Average Total Time for " << h << " hops : " << (t_c + t_p) / num << " seconds" << endl;
			output << "Average SP Probability for " << h << " hops : " << a_pr / num << endl;
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
    else if (argc == 4)
    {
        cout << "Doing betweenness experiment" << endl;
        experiment_betweenness(argv[1], argv[2], atoi(argv[3]));

    }
    else
    {
        experiment(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
    }
	return EXIT_SUCCESS;
}
