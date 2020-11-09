#include "topk.h"
#include "io.h"
#include <algorithm>


void do_queries(string graph_file, string query_file, string threshold_file, string output_file, double factor, int k){

	Graph G = read_graph_from_file(graph_file);

	ifstream queries;
	queries.open(query_file);

	ofstream ofs;
	ofs.open(output_file);

	ifstream res;
	res.open(threshold_file);
	string line;

	int d;
	queries >> d;

	for (int j = 0; j < d; j++){
		Statistics stats;
		double t_c = 0, t_p = 0, a_p = 0;
		int h, n, s, t;
		int candidate_generation_timeout = 0;
		vector<int> r = vector<int>(6);
		queries >> h >> n;
		ofs << "Number of hops = " << h << endl;
		ofs << "Number of queries = " << n << endl << endl;
		for (int i = 0; i < 6; i++)
			r[i] = 0;
		for (int i = 1; i <= n; i++){
			queries >> s >> t;
			ofs << s << "\t" << t << endl;
			do {
				getline(res,line);
			} while (line.substr(0,44) != "Candidate Generation Time without Pruning : ");
			auto topk_result = topk(G, s, t, k, stats, ofs, stod(line.substr(44),nullptr) * factor);
			auto topk_outcome = topk_result.first;
			vector<int> ranks = topk_result.second;
			t_c += stats.candidate_generation;
			t_p += stats.probability_computation;
			candidate_generation_timeout += stats.candidate_generation_timeout; // NB: adding a boolean to an int
			r[min(ranks[0],5)]++;
			double prob = 0;
			int num = min((int)topk_outcome.size(),k);
			for (auto p : topk_outcome)
			{
				ofs << "MPSP:" << endl;
				for(auto e : p.first.edges){
					ofs << e.u << " " << e.v << " " << e.l << " " << e.p << endl;
				}
				prob += p.second;
			}
			a_p += (prob / num);
			ofs << "SP Probability : " << prob / num  << endl;
			ofs << "Rank of MPSP : " << ranks[0] << endl;
			ofs <<fixed << "Candidate Generation Time : " << (1.0*stats.candidate_generation)/CLOCKS_PER_SEC << " seconds" << endl;
			ofs << fixed << "Probability Computation Time : " << (1.0*stats.probability_computation)/CLOCKS_PER_SEC << " seconds" << endl;
			ofs << endl;
		}
		ofs << "Average Candidate Generation Time for " << h << " hops : " << (1.0*t_c)/(CLOCKS_PER_SEC*n) << " seconds" << endl;
		ofs << "Average Probability Computation Time for " << h << " hops : " << (1.0*t_p)/(n*CLOCKS_PER_SEC) << " seconds" << endl;
		ofs << "Average Total Time for " << h << " hops : " << (1.0*(t_c+t_p))/(n*CLOCKS_PER_SEC) << " seconds" << endl;
		ofs << "Average SP Probability for " << h << " hops : " << a_p / n << endl;
		ofs << "Number of candidate generation timeouts for " << h << " hops : " << candidate_generation_timeout << endl;
		ofs << "Number of returned paths at each rank for " << h << " hops : ";
		for (int i = 0; i < 6; i++)
			ofs << r[i] << " ";
		ofs << endl;
		ofs << endl << endl;
	}
	queries.close();
	ofs.close();
	res.close();
}

int main(int argc, char* argv[]){

	if (argc == 7){
		do_queries(argv[1], argv[2], argv[3], argv[4], atof(argv[5]), atoi(argv[6]));
	}
	else{
		cerr << "For WISE experiment" << endl;
		cerr << "Usage: ./experiments <graph_name> <time-multiplying-factor> <k>" << endl;
		cerr << "Usage: ./experiments <path-to-graph> <path-to-queries> <path-to-time-thresholds> <path-to-output> <time-multiplying-factor> <k>" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
