#include <ctime> 

struct Statistics{
    clock_t candidate_generation = 0;
    clock_t probability_computation = 0;
    bool candidate_generation_timeout = false;

};

