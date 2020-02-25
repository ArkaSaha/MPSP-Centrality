#include "topk.h"
#include "io.h"

void do_queries(string graph_file, string query_file, string output_file){
    if(!filesystem::exists(graph_file)){
        cerr << "ERROR: Graph file " << graph_file << " does not exist" <<endl;
        return; 
    }
    if(!filesystem::exists(query_file)){
        cerr << "ERROR: Query file " << query_file << " does not exist" <<endl;
        return; 
    }

    Graph G = read_graph_from_file(graph_file);
    
    ifstream queries;
    queries.open(query_file);

    ofstream ofs;
    ofs.open(output_file);

    cerr << graph_file << endl;

    for (int k = 0; k <= 4; k++){
		cerr << "k = " << k << endl;
		ofs << "Number of hops = " << k * 2 << endl << endl;
		ull a_w = 0;
        Statistics stats;
		double t_c = 0, t_p = 0, a_p = 0;
        int s, t;
        int candidate_generation_timeout = 0;
		for (int i = 1; i <= 100; i++){
            cerr << i << " ";
			queries >> s >> t;
			ofs << s << "\t" << t << endl;
			auto topk_outcome = topk(G, s, t, 1, stats, ofs);
			a_w += topk_outcome[0].first.len();
			a_p += topk_outcome[0].second;
			t_c += stats.candidate_generation;
			t_p += stats.probability_computation;
            candidate_generation_timeout += stats.candidate_generation_timeout; // NB: adding a boolean to an int
			ofs << "Length of MPSP : " << topk_outcome[0].first.len() << endl;
			ofs << "Probability of MPSP : " << topk_outcome[0].second << endl;
			ofs <<fixed << "Candidate Generation Time : " << (1.0*stats.candidate_generation)/CLOCKS_PER_SEC << " seconds" << endl;
			ofs << fixed << "Probability Computation Time : " << (1.0*stats.probability_computation)/CLOCKS_PER_SEC << " seconds" << endl;
			ofs << endl;
		}
        cerr << endl;
		ofs << "Average Length of MPSP for " << k*2 << " hops : " << a_w/100 << endl;
		ofs << "Average Probability of MPSP for " << k*2 << " hops : " << a_p/100 << endl;
		ofs << "Average Candidate Generation Time for " << k * 2 << " hops : " << (1.0*t_c)/(CLOCKS_PER_SEC*100) << " seconds" << endl;
		ofs << "Average Probability Computation Time for " << k * 2 << " hops : " << (1.0*t_p)/(100*CLOCKS_PER_SEC) << " seconds" << endl;
		ofs << "Number of candidate generation timeouts : " << candidate_generation_timeout << endl;
		ofs << endl << endl;
	}
	queries.close();
    ofs.close();
}

void do_queries(string graph_name){
    do_queries(graph_name + ".graph", graph_name + ".queries", graph_name+".output");
}

int main(){


    string folder = "data/Synthetic/";
    vector<string> graph_names = {
                                  "ER/ER_10000_20000", 
                                  //"ER/ER_20000_40000", 
                                  //"ER/ER_40000_80000", 
                                  //"ER/ER_80000_160000",
                                  //"BA/BA_10000_199790",
                                  //"BA/BA_20000_399790",
                                  //"BA/BA_40000_799790",
                                  //"BA/BA_80000_1599790",
                                  //"BP/BP_10000_40000",
                                  //"BP/BP_20000_80000",
                                  //"BP/BP_40000_160000",
                                  //"BP/BP_80000_320000",
                                  //"SF/SF_10000_21752",
                                  //"SF/SF_20000_43813",
                                  //"SF/SF_40000_87008",
                                  //"SF/SF_80000_173923"
    };

    for(string graph_name: graph_names){
        do_queries(folder + graph_name);
    }


    
    Graph G = read_graph_from_file("data/Synthetic/ER/ER_10000_20000.graph");
    Statistics stats;
    auto paths = yen(G, 0, 1000, 1, stats, cout);

    int n = 15;
    if(paths.size() > 20){
        paths = vector<Path>(paths.begin(), paths.begin() + (n+1));
    }
    cout << paths.size() << endl;
    for(Path p: paths){
        p.print(); cout << endl;
    }

    cout << "\n**\n" << endl;
    auto paths2 = yen(G, paths[n]);
    for(Path p: paths2){
        p.print(); cout << endl;
    }

    for(int i=0; i<=n; i++){
        assert(paths[i] == paths2[i]);
    }

    cout << "The Path:" << endl;
    paths[n].print();

    double exact = exact_probability(G, paths[n]);
    cout << "Exact     Probability : " << exact  << endl;

    double LK = Luby_Karp(paths2, n , 1000000);
    cout << "Luby Karp Probability : " << LK << endl;

    double LK2 = Luby_Karp(G, paths[n], 1000000);
    cout << "Luby Karp Probability2: " << LK2 << endl;
    


    return 0;
}

