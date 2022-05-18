#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "math.h"

#include "timing_params.h"

double randdouble_between(double a, double b) {
    // drand48: https://pubs.opengroup.org/onlinepubs/7908799/xsh/drand48.html
    double len = b - a;
    return a + len * drand48();
}

int randint_between(int a, int b) {
    return a + rand() / (RAND_MAX / (b - a + 1) + 1);
}

void fill_linspace(double * const array, size_t dim, double start, double end) {
    double step = (end - start) / dim;
    array[0] = start;
    for(size_t it=1; it<dim; ++it) {
        array[it] = array[it-1] + step;
    }
}

// === run_pso ===

// x_1^2 + ... + x_d^2
double f1(double const * const x, size_t d)
{
    double acc=0;
    for(size_t it=0; it<d; ++it) {
        acc += pow(x[it], 2);
    }
    return acc;
}

// x_1^3 + ... + x_(d-1)^3 + 3(x_2 + ... + x_d) + 3
double f2(double const * const x, size_t d) {
    double acc=0;
    for(size_t it=0; it<d-1; ++it) {
        acc += pow(x[it], 3) + 3*x[it+1] + 3;
    }
    return acc;
}

// x_1^2 * x_2 + ... x_(d-1)^2 * x_d + x_d^2
double f3(double const * const x, size_t d) {
    double acc=0;
    for(size_t it=0; it<d-1; ++it) {
        acc += pow(x[it],2) * x[it+1];
    }
    acc += pow(x[d-1], 2);
    return acc;
}

void generate_run_pso_parameters(run_pso_parameters* params, profiling_hyperparameters* hyperparams) {
    srand(clock());

    // Parameter space.
    // fixed parameters (flops count is expressed as a function of these terms)
    params->dimensions = hyperparams->dimensions;
    params->population_size = hyperparams->population_size;
    params->time_max = hyperparams->time_max;
    params->n_trials = hyperparams->n_trials;

    // randomizable parameters (flops count is NOT expressed as a function of these terms)
    // blackbox_fun functions[3] = {&f1, &f2, &f3};
    // size_t function_index = randint_between(0, 2);
    // params->f = functions[function_index];
    params->f = &f1;

    double inertia_space[5];
    fill_linspace(inertia_space, 5, 0.5, 1.0);
    size_t inertia_index = randint_between(0, 4);
    params->inertia = inertia_space[inertia_index];

    double social_space[10];
    fill_linspace(social_space, 10, 0.5, 1.5);
    size_t social_index = randint_between(0, 9);
    params->social = social_space[social_index];

    double cognition_space[10];
    fill_linspace(cognition_space, 10, 0.5, 1.5);
    size_t cognition_index = randint_between(0, 9);
    params->cognition = cognition_space[cognition_index];

    double local_refinement_box_size_space[10];
    fill_linspace(cognition_space, 10, 2, 4);
    size_t local_refinement_box_size_index = randint_between(0, 9);
    params->local_refinement_box_size = local_refinement_box_size_space[local_refinement_box_size_index];

    double min_minimizer_distance_space[10];
    fill_linspace(min_minimizer_distance_space, 10, 2, 4);
    size_t min_minimizer_distance_index = randint_between(0, 9);
    params->min_minimizer_distance = min_minimizer_distance_space[min_minimizer_distance_index];

    double bounds_low_space[10];
    fill_linspace(bounds_low_space, 10, -20, -10);
    size_t bounds_low_index = randint_between(0, 9);
    double bounds_low0 = bounds_low_space[bounds_low_index];
    params->bounds_low = malloc(sizeof(double) * params->dimensions);
    for(size_t it=0; it<params->dimensions; ++it) {
        params->bounds_low[it] = bounds_low0;
    }

    double bounds_high_space[10];
    fill_linspace(bounds_high_space, 10, 10, 20);
    size_t bounds_high_index = randint_between(0, 9);
    double bounds_high0 = bounds_high_space[bounds_high_index];
    params->bounds_high = malloc(sizeof(double) * params->dimensions);
    for(size_t it=0; it<params->dimensions; ++it) {
        params->bounds_high[it] = bounds_high0;
    }

    double vmin_space[10];
    fill_linspace(vmin_space, 10, -20, -10);
    size_t vmin_space_index = randint_between(0, 9);
    double vmin0 = vmin_space[vmin_space_index];
    params->vmin = malloc(sizeof(double) * params->dimensions);
    for(size_t it=0; it<params->dimensions; ++it) {
        params->vmin[it] = vmin0;
    }

    double vmax_space[10];
    fill_linspace(vmax_space, 10, 10, 20);
    size_t vmax_space_index = randint_between(0, 9);
    double vmax0 = vmax_space[vmax_space_index];
    params->vmax = malloc(sizeof(double) * params->dimensions);
    for(size_t it=0; it<params->dimensions; ++it) {
        params->vmax[it] = vmax0;
    }

    double initial_position_space[20];
    fill_linspace(initial_position_space, 20, -10, 10);
    double* initial_position_data = malloc(sizeof(double) * params->dimensions * params->population_size);
    for(size_t it=0; it<params->dimensions * params->population_size; ++it) {
        size_t initial_position_space_index = randint_between(0, 19);
        initial_position_data[it] = initial_position_space[initial_position_space_index];
    }

    params->initial_positions = malloc(sizeof(double) * params->population_size);
    for(size_t it = 0; it < params->population_size; ++it) {
        params->initial_positions[it] = initial_position_data + params->dimensions * it;
    }
}

void destruct_run_pso_params(run_pso_parameters* params) {
    free(params->bounds_low);
    free(params->bounds_high);

    free(params->vmin);
    free(params->vmax);

    free(params->initial_positions[0]);
    free(params->initial_positions);
}

void run_pso_wrapper(run_pso_parameters* params) {
    run_pso(
        params->f, 
        params->inertia, params->social, params->cognition,
        params->local_refinement_box_size,
        params->min_minimizer_distance,
        params->dimensions,
        params->population_size, params->time_max, params->n_trials,
        params->bounds_low, params->bounds_high,
        params->vmin, params->vmax,
        params->initial_positions   
    );
}

// === plu_factorization ===

#define IDX(MAT, ROW, COL) (MAT)[N * (ROW) + (COL)]
// this is how we generate random invertible matrices fast:
// https://en.wikipedia.org/wiki/Diagonally_dominant_matrix
void generate_random_invertible_matrix(int N, double* A) {
    for(size_t i=0; i<N; ++i) {
        double budget = randdouble_between(10, 30);
        IDX(A, i, i) = budget;
        budget -= 1; // for strictly diagonally dominant
        for(size_t j=0; j<N; ++j) {
            if(j != i) {
                double pick = randdouble_between(0, budget);
                IDX(A, j, i) = pick;
                budget -= pick;
            }
        }
    }
}

void generate_plu_factorize_parameters(plu_factorize_parameters* params, profiling_hyperparameters* hyperparams) {
    // fixed parameters (flops count is expressed as a function of these terms)
    params->N = hyperparams->plu_matrix_side_length;

    // randomized parameters
    params->plu_ft = malloc(sizeof(plu_factorization));
    alloc_plu_factorization(params->N, params->plu_ft);
    params->A = malloc(params->N * params->N * sizeof(double));
    generate_random_invertible_matrix(params->N, params->A);
}

void plu_factorize_wrapper(plu_factorize_parameters* params) {
    plu_factorize(params->N, params->A, params->plu_ft);
}

// === plu_solve ===

void generate_random_vector(int N, double* b) {
    for(size_t i=0; i<N; ++i) {
        b[i] = randdouble_between(1, 10);
    }
}

void generate_plu_solve_parameters(plu_solve_parameters* params, profiling_hyperparameters* hyperparams) {
    params->plu_fact_params = malloc(sizeof(plu_factorize_parameters));
    generate_plu_factorize_parameters(params->plu_fact_params, hyperparams);

    // generate plu result
    plu_factorize_wrapper(params->plu_fact_params);

    params->b = malloc(params->plu_fact_params->N * sizeof(double));
    generate_random_vector(params->plu_fact_params->N, params->b);
    params->x = malloc(params->plu_fact_params->N * sizeof(double));
}

void plu_solve_wrapper(plu_solve_parameters* params) {
    plu_solve(params->plu_fact_params->N, params->plu_fact_params->plu_ft, params->b, params->x);
}

// === fit_surrogate ===

/* Parameters of pso used in fit_surrogate:
    - pso->population_size
    - pso->time 
    - pso->dimensions
    - pso->x ([time][population_size_index][dimensions])
    - pso->x_eval (=f(x), [time][population_size_index])
    - pso->lambda (= NULL because it is realloc-ed there)
    - pso->p (= malloc((pso->dimensions + 1) * sizeof(double)) because it's overwritten) 
*/

void generate_fit_surrogate_parameters(struct pso_data_constant_inertia* params, profiling_hyperparameters* hyperparams) {
    params->population_size = hyperparams->population_size;
    params->time = hyperparams->time_max-1;
    params->dimensions = hyperparams->dimensions;

    params->time_max = hyperparams->time_max;

    params->f = &f1;

    params->x = (double***)malloc(params->time_max * sizeof(double **));
    for (int t = 0 ; t < params->time_max ; t++)
    {
        params->x[t] = (double**)malloc(params->population_size * sizeof(double*));
        for (int i = 0 ; i < params->population_size ; i++)
        {
            params->x[t][i] = (double*)malloc(params->dimensions * sizeof(double));
            generate_random_vector(params->dimensions, params->x[t][i]);
        }
    }

    params->x_eval = (double**)malloc(params->time_max * sizeof(double *));
    for (int t = 0 ; t < params->time_max ; t++)
    {
        params->x_eval[t] = (double*)malloc(params->population_size * sizeof(double));
        for(int i =0 ; i < params->population_size ; i++) {
            params->x_eval[t][i] = params->f(params->x[t][i], params->dimensions);
        }
    }

    params->lambda = NULL;
    params->p = malloc((params->dimensions + 1) * sizeof(double));
}

void fit_surrogate_wrapper(struct pso_data_constant_inertia * params) {
    fit_surrogate(params);
}