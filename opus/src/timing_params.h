#pragma once
#include "pso.h"

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

void generate_run_pso_parameters(run_pso_parameters* params);

void run_pso_wrapper(run_pso_parameters* params);