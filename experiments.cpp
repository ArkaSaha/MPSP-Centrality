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

    // cerr << graph_file << endl;

    for (int k = 0; k <= 4; k++){
		// cerr << "k = " << k << endl;
		ofs << "Number of hops = " << k * 2 << endl << endl;
		ull a_w = 0;
        Statistics stats;
		double t_c = 0, t_p = 0, a_p = 0;
        int s, t;
        int candidate_generation_timeout = 0;
		for (int i = 1; i <= 100; i++){
            // cerr << i << " ";
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
        // cerr << endl;]
		ofs << "Average Length of MPSP for " << k*2 << " hops : " << a_w/100 << endl;
		ofs << "Average Probability of MPSP for " << k*2 << " hops : " << a_p/100 << endl;
		ofs << "Average Candidate Generation Time for " << k * 2 << " hops : " << (1.0*t_c)/(CLOCKS_PER_SEC*100) << " seconds" << endl;
		ofs << "Average Probability Computation Time for " << k * 2 << " hops : " << (1.0*t_p)/(100*CLOCKS_PER_SEC) << " seconds" << endl;
        ofs << "Average Total Time for " << k * 2 << " hops : " << (1.0*(t_c+t_p))/(100*CLOCKS_PER_SEC) << " seconds" << endl;
		ofs << "Number of candidate generation timeouts : " << candidate_generation_timeout << endl;
		ofs << endl << endl;
	}
	queries.close();
    ofs.close();
}

void do_queries(string graph_name){
    do_queries(graph_name + ".graph", graph_name + ".queries", graph_name+".output");
}

int main(int argc, char* argv[]){

    do_queries(argv[1],argv[2],argv[3]);

    return 0;
}

