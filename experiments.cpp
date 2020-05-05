#include "topk.h"
#include "io.h"
#include <algorithm>


namespace fs = std::filesystem;

void do_queries(string graph_file, string query_file, string output_file, double THRESHOLD_NR_OF_SECONDS){
	if(!fs::exists(graph_file)){
		cerr << "ERROR: Graph file " << graph_file << " does not exist" <<endl;
		return; 
	}
	if(!fs::exists(query_file)){
		cerr << "ERROR: Query file " << query_file << " does not exist" <<endl;
		return; 
	}

	Graph G = read_graph_from_file(graph_file);

	ifstream queries;
	queries.open(query_file);

	ofstream ofs;
	ofs.open(output_file);

	// cerr << graph_file << endl;

	for (int k = 0; k <= 3; k++){
		// cerr << "k = " << k << endl;
		ofs << "Number of hops = " << k * 2 << endl;
		ull a_w = 0;
		Statistics stats;
		double t_c = 0, t_p = 0, a_p = 0;
		int s, t, n[6], n_q;
		int candidate_generation_timeout = 0;
		for (int i = 0; i < 6; i++)
			n[i] = 0;
		queries >> n_q;
		ofs << "Number of queries = " << n_q << endl << endl;
		for (int i = 1; i <= n_q; i++){
			// cerr << i << " ";
			queries >> s >> t;
			ofs << s << "\t" << t << endl;
			auto topk_result = topk(G, s, t, 1, stats, ofs, THRESHOLD_NR_OF_SECONDS);
			auto topk_outcome = topk_result.first;
			vector<int> ranks = topk_result.second;
			a_w += topk_outcome[0].first.len();
			a_p += topk_outcome[0].second;
			t_c += stats.candidate_generation;
			t_p += stats.probability_computation;
			candidate_generation_timeout += stats.candidate_generation_timeout; // NB: adding a boolean to an int
			n[min(ranks[0],5)]++;
			ofs << "MPSP:" << endl;
			for(auto e: topk_outcome[0].first.edges){
				ofs << e.u << " " << e.v << " " << e.l << " " << e.p << endl;
			}
			ofs << "Length of MPSP : " << topk_outcome[0].first.len() << endl;
			ofs << "Probability of MPSP : " << topk_outcome[0].second << endl;
			ofs << "Probability of being MPSP : " << topk_outcome[0].second << endl;
			ofs << "Rank of MPSP : " << ranks[0] << endl;
			ofs <<fixed << "Candidate Generation Time : " << (1.0*stats.candidate_generation)/CLOCKS_PER_SEC << " seconds" << endl;
			ofs << fixed << "Probability Computation Time : " << (1.0*stats.probability_computation)/CLOCKS_PER_SEC << " seconds" << endl;
			ofs << endl;
		}
		// cerr << endl;]
		ofs << "Average Length of MPSP for " << k * 2 << " hops : " << a_w/n_q << endl;
		ofs << "Average Probability of being MPSP for " << k * 2 << " hops : " << a_p/n_q << endl;
		ofs << "Average Candidate Generation Time for " << k * 2 << " hops : " << (1.0*t_c)/(CLOCKS_PER_SEC*n_q) << " seconds" << endl;
		ofs << "Average Probability Computation Time for " << k * 2 << " hops : " << (1.0*t_p)/(n_q*CLOCKS_PER_SEC) << " seconds" << endl;
		ofs << "Average Total Time for " << k * 2 << " hops : " << (1.0*(t_c+t_p))/(n_q*CLOCKS_PER_SEC) << " seconds" << endl;
		ofs << "Number of candidate generation timeouts for " << k * 2 << " hops : " << candidate_generation_timeout << endl;
		ofs << "Number of returned paths at each rank for " << k * 2 << " hops : ";
		for (int i = 0; i < 6; i++)
			ofs << n[i] << " ";
		ofs << endl;
		ofs << endl << endl;
	}
	queries.close();
	ofs.close();
}

int main(int argc, char* argv[]){

	if (argc == 3){
		string graph_name(argv[1]);
		string graph_type = graph_name.substr(0, 2);
		string graph = "data/Synthetic/" + graph_type + "/" + argv[1] + ".graph";
		string queries = "data/Synthetic/" + graph_type + "/" + argv[1] + ".queries";
		string output = "output/" + graph_type + "/" + argv[1] + "_WISE_" + argv[2] + ".output";
		cout << "Graph  file : " << graph << endl;
		cout << "Query  file : " << queries << endl;
		cout << "Output file : " << output << endl;
		cout << "Seconds      : " << atof(argv[2]) << endl;
		do_queries(graph, queries, output, atof(argv[2]));
	}
	else if (argc == 5){
		do_queries(argv[1],argv[2],argv[3], atof(argv[4]));
	}
	else{
		cerr << "For WISE experiment" << endl;
		cerr << "Usage: ./experiment <graph_name> <nr_of_seconds>" << endl;
		cerr << "Usage: ./experiment <path-to-graph> <path-to-queries> <path-to-output> <nr_of_seconds>" << endl;
		return 1;
	}

	return 0;
}

