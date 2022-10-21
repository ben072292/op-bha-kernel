#include "lattice_model.h"
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <mpi.h>
#include <thread>
#include <omp.h>
using namespace std;

bool is_done(bool* checked_set, int size){
    for(int i = 0; i < size; i++){
        if(!checked_set[i]) return false;
    }
    return true;
}

void select_n(bool* checked_set, int* table, int size, int select_size, vector<int>& ret){
    int counter = 0;
    for(int i = 0; i < size; i++){
        if(!checked_set[table[i]]) ret[counter++] = table[i];
        if(counter == select_size) break;
    }
}

void fill_checked_set_up_helper(int n, int state, int* add_index, bool* checked_set){
    int pow_set_size = 1 << n;
    int i, j, temp;
    for(i= 0; i < pow_set_size; i++){
        temp = state;
        for(j = 0; j < n; j++){
            if((i&(1<<j)) > 0) temp += add_index[j];
        }
        checked_set[temp] = true;
    }
}

void fill_checked_set_up(int state, int pool_size, bool* checked_set){
    int* add_index = new int[pool_size -__builtin_popcount(state)];
    int counter=0, i, index;
    for(i = 0; i < pool_size; i++){
        index = (1 << i);
        if((state & index) == 0) add_index[counter++] = index;
    }
    fill_checked_set_up_helper(pool_size-__builtin_popcount(state), state, add_index, checked_set);
    delete[] add_index;
}

void fill_checked_set_down_helper(int n, int state, int* sub_index, bool* checked_set){
    int pow_set_size = (1 << n);
    int* ret = new int[pow_set_size];
    int i, j, temp;
    for(i= 0; i < pow_set_size; i++){
        temp = state;
        for(j = 0; j < n; j++){
            if((i&(1<<j)) > 0) temp -= sub_index[j];
        }
        checked_set[temp] = true;
    }
}

void fill_checked_set_down(int state, int pool_size, bool* checked_set){

    int* sub_index = new int[__builtin_popcount(state)];
    int counter = 0, index;
    for(int i = 0; i < pool_size; i++){
        index = (1<<i);
        if((state&index) == index) sub_index[counter++] = index;
    }
    fill_checked_set_down_helper(__builtin_popcount(state), state, sub_index, checked_set);
    delete[] sub_index;
}

int main(int argc, char* argv[]){
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    int pool_size = atoi(argv[1]);
    double* prior = new double[pool_size];
    
    double prior_val = atof(argv[2]);
    for(int i = 0; i < pool_size; i++) prior[i] = prior_val;
    
    // omp_set_num_threads(atoi(argv[3]));
    pool_stat* pool = new pool_stat(pool_size, prior);

    double run_time = 0.0;

    lattice_model* model = new lattice_model(pool, true);
    double** dilution_matrix = model->generate_dilution_matrix(0.99, 0.005);
    
    run_time -= MPI_Wtime();

    int* sophsticated_table;
    bool* checked_set = (bool*)malloc(sizeof(bool) * (1 << pool_size));
    
    bool* temp;
    bool done = false, done_reduce = false;
    double local_min = 2.0;
    double global_min = 2.0;

    sophsticated_table = model->sophisticated_selection();

     MPI_Request checked_set_request, min_request, done_request;
     int checked_set_ready, min_ready, done_ready;

     int offset = 0;
    
    while(!done){
        // auto start = chrono::high_resolution_clock::now();
        // #pragma omp parallel for schedule(dynamic)
        for(int i = 0; ; i++){
            int index = i * world_size + world_rank + offset;
            if(index >= 1 << pool_size) break;

            int state = sophsticated_table[index];
            if(checked_set[state]) continue;
            
            double val = model->get_up_set_mass(state);
            checked_set[state] = true;
            if(val < 0.5){
                fill_checked_set_up(state, model->get_pool_size(), checked_set);
            }
            else if (val > 0.5){
                fill_checked_set_down(state, model->get_pool_size(), checked_set);
            }
            double temp = abs(val - 0.5);
            if(temp < local_min){
                local_min = temp;
            }

            // MPI_Iallreduce(checked_set, checked_set_reduce, 1<<pool_size, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD, &checked_set_request);
            // MPI_Wait(&checked_set_request, MPI_STATUS_IGNORE);
            MPI_Allreduce(MPI_IN_PLACE, checked_set, 1<<pool_size, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

            // if(world_rank == 0) cout << "Temp min: " << global_min << endl;
            done = is_done(checked_set, 1 << pool_size);
            // MPI_Bcast(&done, 1, MPI_INT, world_rank, MPI_COMM_WORLD);
            // MPI_Iallreduce(&done, &done_reduce, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD, &done_request);
            MPI_Allreduce(MPI_IN_PLACE, &done, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

        }
        // auto stop = chrono::high_resolution_clock::now();
        // auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        // cout << "Time: " << duration.count() << endl;

        offset++;
    }
    // Reduce all of the local min into the global min
    // MPI_Ireduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD, &min_request);
    MPI_Reduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    run_time += MPI_Wtime();

    if(world_rank == 0){
        cout << "Min: " << global_min << endl << "Time: " << run_time << "s" << endl;
    }

    // final cleanup
    if(world_rank == 0){
        delete[] sophsticated_table;
        // free(select);
    }
    else model->delete_posterior_probability_map();
    
    delete[] checked_set;

    // Finalize the MPI environment.
    MPI_Finalize();
}