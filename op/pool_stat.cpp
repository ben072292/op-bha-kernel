#include "pool_stat.h"

using namespace std;
pool_stat::pool_stat(int pool_size, float16_t* prior){
    this->pool_size = pool_size;
    this->prior = prior;
}

pool_stat::~pool_stat(){
    delete[] prior;
}

float16_t* pool_stat::generate_prior_probability_map(){
    int size = (1 << pool_size) / _v_size;
    float16_t* ret = float16_alloc(size);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < _v_size; j++){
            ret[i][j] = generate_prior_probability(i*_v_size+j);
        }
    }
    return ret;
}

float pool_stat::generate_prior_probability(int state){
    float ret = 1.0;
    for(int i = pool_size-1; i >=0; i--){
        if((state & 1) == 1) ret *= (1-prior[i]);
        else ret *= prior[i];
        state = state >> 1;
    }
    return ret;
}
