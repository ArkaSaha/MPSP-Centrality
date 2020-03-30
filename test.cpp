
#include "io.h"


Path i2path(Graph & g, vector<int> I){
    vector<Edge> edges = vector<Edge>();
    for(int i: I){
        edges.push_back(*(g.index2edge[i]));
    }
    return Path(edges);
}


int main(){

    Graph G = read_graph_from_file("data/toy/test.graph");
    // File should contain the following graph
    /*
       5 6
       0 1 10 0.1
       1 2 3 0.4
       1 3 2 0.9
       3 2 4 0.9
       2 4 5 0.6
       3 4 10 0.9
    */
    Path p = dijkstra(G, 0, 4);

    Path ans = Path({*G.index2edge[0], *G.index2edge[1], *G.index2edge[4]});
    Path ans2 = i2path(G, {0, 1, 4});

    assert(ans == p);
    assert(ans2 == p);

    Path p1 = i2path(G, {0, 1, 4});
    Path p2 = i2path(G, {0, 2, 3, 4});
    Path p3 = i2path(G, {0, 2, 5});
    vector<Path> paths_ans = {p1, p2, p3};

    vector<double> prob_ans = {0.024, 0.02916, 0.035316};

    double exact_tolerance = 0.000001;
    double LK_tolerance = 0.001;

    vector<Path> paths = yen(G, p3, 1000);
    for(uint i=0; i<paths.size(); i++){
        Path p = paths[i];

        assert(paths[i] == paths_ans[i]);
        double exact = exact_probability(G, p);
        double LK = Luby_Karp(G, p, 100000);

        assert(prob_ans[i] - exact_tolerance <= exact and exact <= prob_ans[i] + exact_tolerance);
        assert(prob_ans[i] - LK_tolerance <= LK and LK <= prob_ans[i] + LK_tolerance);

    }








    






    return 0;
}

