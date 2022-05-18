#pragma once

#include <stdlib.h>

typedef double (*blackbox_fun)(double const * const, size_t d);

// for other inertia choices, see https://ieeexplore.ieee.org/document/6089659
struct pso_data_constant_inertia
{
    blackbox_fun f;
    // positions x_i, saved for all times
    // i.e. x[t][i] is the position vector x_i(t)
    double *** x;
    // f(x_i), for all times
    double ** x_eval;

    // velocity v_i
    double ** v;
    // best position recorded for particle i y_i
    double ** y;
    double * y_eval;
    
    double * y_hat;

    // Keep track of local refinement points.
    // max lenght of list = tmax
    // current length n_past_refinement_points
    double ** past_refinement_points;
    double * past_refinement_points_eval; 


    double * bound_low;
    double * bound_high;
    double * vmin;
    double * vmax;

    // parameters of the surrogate
    double * lambda;
    double * p;

    double inertia;
    double social;
    double cognition;

    double local_refinement_box_size;
    double min_minimizer_distance;

    double y_hat_eval;
    
    int dimensions;
    int population_size;
    int n_trials;

    int n_past_refinement_points;

    int time_max;
    int time;
};

void run_pso(
    blackbox_fun f,
    double inertia,
    double social, double cognition,
    double local_refinement_box_size,
    double min_minimizer_distance,
    int dimensions,
    int population_size, int time_max, int n_trials,
    double * bounds_low, double * bounds_high,
    double * vmin, double * vmax,
    double ** initial_positions
);

int fit_surrogate(struct pso_data_constant_inertia * pso);