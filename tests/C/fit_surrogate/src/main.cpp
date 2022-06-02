

#include <cmath>
#include <cstdio>

#include "perf_testers/PerformanceTester.hpp"

extern "C"
{
#include "helpers.h"
#include "latin_hypercube.h"
#include "pso.h"
#include "steps/fit_surrogate.h"
}

#define POPSIZE 20
#define DIMENSION 20
#define SPACE_FILLING_DESIGN_SIZE 25

static double griewank_Nd(double const *const x)
{
  // see https://mathworld.wolfram.com/GriewankFunction.html
  double r = 0;
  double t = 1;

  double d = 1. / 4000;

  for (size_t i = 0; i < DIMENSION; i++)
  {
    double v = x[i];

    r += v * v;
    t *= cos(v / sqrt((double)i + 1));
  }

  return (1. + d * r - t);
}

int main(int argc, char **argv)
{
  /*
   * Get a valid PSO object
   */
  srand(42);

  double inertia = 0.8;
  double social = 0.1, cognition = 0.2;
  double local_refinement_box_size = 5.;
  double min_dist = 0.01;
  int dimensions = DIMENSION;
  int population_size = POPSIZE;
  int time_max = 40;
  int n_trials = 10;
  double bounds_low[DIMENSION] = {0};
  double bounds_high[DIMENSION] = {0};
  double vmin[DIMENSION] = {0};
  double vmax[DIMENSION] = {0};

  for (size_t k = 0; k < DIMENSION; k++)
  {
    bounds_low[k] = -500, bounds_high[k] = 700;
    vmin[k] = -50, vmax[k] = 50;
  }

  double lh[SPACE_FILLING_DESIGN_SIZE * DIMENSION] = {0};
  latin_hypercube(lh, SPACE_FILLING_DESIGN_SIZE, DIMENSION);

  double space_filling_design[SPACE_FILLING_DESIGN_SIZE * DIMENSION] = {0};
  for (size_t i = 0; i < SPACE_FILLING_DESIGN_SIZE; i++)
  {
    for (size_t k = 0; k < DIMENSION; k++)
    {
      double lo = bounds_low[k], hi = bounds_high[k];
      space_filling_design[i * DIMENSION + k] =
          lo + (hi - lo) * lh[i * DIMENSION + k];
    }
  }

  struct pso_data_constant_inertia pso;
  pso_constant_inertia_init(&pso, &griewank_Nd, inertia, social, cognition,
                            local_refinement_box_size, min_dist, dimensions,
                            population_size, time_max, n_trials, bounds_low,
                            bounds_high, vmin, vmax, SPACE_FILLING_DESIGN_SIZE);

  pso_constant_inertia_first_steps(&pso, SPACE_FILLING_DESIGN_SIZE,
                                   space_filling_design);

  while (pso.time < pso.time_max - 1)
  {
    pso_constant_inertia_loop(&pso);
  }

  PerformanceTester<fit_surrogate_fun_t> perf_tester;

  auto arg_restorer = [&]()
  { pso.x_distinct_idx_of_last_batch = pso.x_distinct_s - 10; };

  perf_tester.add_function(fit_surrogate, "fit_surrogate", 1);
  perf_tester.perf_test_all_registered(std::move(arg_restorer), &pso);
}