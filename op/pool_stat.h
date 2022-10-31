#ifndef POOL_STAT_H
#define POOL_STAT_H

#include "float_16_t.h"

class pool_stat{
    public:
    pool_stat(int pool_size, float16_t* prior);

    ~pool_stat();

    int get_pool_size(){return pool_size;}

    float16_t* get_prior(){return prior;}

    float16_t* generate_prior_probability_map();

    float generate_prior_probability(int state);

    private: 
    int pool_size;
    float16_t *prior;

};
#endif