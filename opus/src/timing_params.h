#pragma once
#include "pso.h"


// === hyperparameters ===

typedef struct {
    int dimensions;
    int population_size;
    int time_max;
    int n_trials;
    int plu_matrix_side_length;
} profiling_hyperparameters;

// === run_pso ===

typedef struct {
    blackbox_fun f;
    double inertia;
    double social;
    double cognition;
    double local_refinement_box_size;
    double min_minimizer_distance;
    int dimensions;
    int population_size;
    int time_max; 
    int n_trials;
    double* bounds_low;
    double* bounds_high;
    double* vmin;
    double* vmax;
    double** initial_positions;
} run_pso_parameters;

void generate_run_pso_parameters(run_pso_parameters* params, profiling_hyperparameters* hyperparams);

void run_pso_wrapper(run_pso_parameters* params);

// // === plu_factorization ===

// typedef struct {
//     int N;
//     double* A;
//     plu_factorization* plu_ft;
// } plu_factorize_parameters;

// void generate_plu_factorize_parameters(plu_factorize_parameters* params, profiling_hyperparameters* hyperparams);

// void plu_factorize_wrapper(plu_factorize_parameters* params);

// // === plu_solve ===

// typedef struct {
//     plu_factorize_parameters* plu_fact_params;
//     double* b;
//     double* x;
// } plu_solve_parameters;

// void generate_plu_solve_parameters(plu_solve_parameters* params, profiling_hyperparameters* hyperparams);

// void plu_solve_wrapper(plu_solve_parameters* params);

// === fit_surrogate === 

void generate_fit_surrogate_parameters(struct pso_data_constant_inertia* params, profiling_hyperparameters* hyperparams);

void fit_surrogate_wrapper(struct pso_data_constant_inertia* params);
