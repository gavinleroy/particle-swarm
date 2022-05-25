
#include <stdio.h>

#include "helpers.h"
#include "perf_testers/perf_fit_surrogate.h"
#include "pso.h"

double my_f(double const *const x)
{
  return (x[0] - 2) * (x[0] - 2) + (x[1] - 5) * (x[1] - 5);
}

int main(int argc, char **argv)
{
  /*
   * Get a valid PSO object
   */

  double inertia = 0.7;
  double social = 1., cognition = 1.;
  double local_refinement_box_size = 3.;
  double min_minimizer_distance = 1.;
  int dimensions = 2;
  int population_size = 5;
  int time_max = 100;
  int n_trials = 5;
  double bounds_low[2] = {-10., -10.};
  double bounds_high[2] = {10., 10.};
  double vmin[2] = {-10., -10.};
  double vmax[2] = {10., 10.};
  double initial_positions[10] = {-1, 8, 5, -3, 5, 6, 7, 3, -9, -2};

  struct pso_data_constant_inertia pso;
  pso_constant_inertia_init(
      &pso, &my_f, inertia, social, cognition, local_refinement_box_size,
      min_minimizer_distance, dimensions, population_size, time_max, n_trials,
      bounds_low, bounds_high, vmin, vmax, initial_positions);

  pso_constant_inertia_first_steps(&pso);

  while (pso.time < pso.time_max - 1)
  {
    pso_constant_inertia_loop(&pso);
  }

  perf_test_fit_surrogate(&pso);
}