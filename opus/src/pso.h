#pragma once

#include <stdlib.h>

typedef double (*blackbox_fun)(double const * const, size_t d);

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