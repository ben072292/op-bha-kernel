#ifndef LATTICE_MODEL_H
#define LATTICE_MODEL_H

#include "pool_stat.h"
#include <vector>


class lattice_model{
    private:
    int pool_size;
    float16_t* posterior_probability_map;

    public:
    lattice_model(pool_stat *pool_stat, bool gen_prior);
    ~lattice_model();
    int get_pool_size(){return pool_size;}
    float16_t* get_posterior_probability_map(){return posterior_probability_map;}
    void set_posterior_probability_map(float16_t* val){posterior_probability_map = val;}
    void delete_posterior_probability_map(){delete[] posterior_probability_map;}
    void free_posterior_probability_map(){free(posterior_probability_map);}
    int* get_up_set(int state);
    int* generate_power_set_adder(int* add_index, int state);
    int* get_down_set(int state);
    int* generate_power_set_reducer(int* reduce_index, int state);
    float16_t* update_posterior_probability(int state, int response, float** dilution_matrix);
    float** generate_dilution_matrix(float alpha, float h);
    float16_t* lattice_model::compute_dilution_response(float** dilution_matrix, int e, int start_state);
    int find_halving_state(float prob);
    int find_halving_state_openmp_static_scheduler(float prob);
    int find_halving_state_openmp_dynamic_scheduler(float prob);
    float get_up_set_mass(int state);
    int* sophisticated_selection();
    int binomial(int N, int K);

};

#endif