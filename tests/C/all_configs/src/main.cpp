

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "perf_testers/PerformanceTester.hpp"

extern "C"
{
#include <getopt.h>

#include "helpers.h"
#include "latin_hypercube.h"
#include "pso.h"
#include "steps/fit_surrogate.h"
#include "steps/surrogate_eval.h"
}

#define POPSIZE 9
#define DIMENSION 13
#define SPACE_FILLING_DESIGN_SIZE 16

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

static struct option long_options[] = {
    {"interval", required_argument, 0, 'I'},
    {"measurements", required_argument, 0, 'M'},
    {"print", no_argument, 0, 'P'},
    {"bench-fit-surrogate", required_argument, 0, 'f'},
    {"bench-surrogate-eval", required_argument, 0, 'e'},
    {"write", required_argument, 0, 'w'},
    {0, 0, 0, 0}};

int main(int argc, char **argv)
{
  int nb_measurements = 20;
  int time_between_measures = 5;
  bool do_bench_fit_surrogate = false;
  char *fit_surrogate_name = NULL;

  bool do_bench_surrogate_eval = false;
  char *surrogate_eval_name = NULL;

  bool do_print_outputs = false;

  char *write_to_fname = NULL;

  while (1)
  {
    int option_index = 0;
    int c = getopt_long(argc, argv, "IMPfew", long_options, nullptr);
    if (c == -1)
      break;

    switch (c)
    {
    case 'I':
      time_between_measures = std::stoi(optarg);
      break;
    case 'M':
      nb_measurements = std::stoi(optarg);
      break;
    case 'P':
      do_print_outputs = true;
      break;
    case 'f':
      do_bench_fit_surrogate = true;
      fit_surrogate_name = optarg;
      break;
    case 'e':
      do_bench_surrogate_eval = true;
      surrogate_eval_name = optarg;
      break;
    case 'w':
      write_to_fname = optarg;
      break;
    case '?':
      break;

    default:
      printf("?? getopt returned character code 0%o ??\n", c);
    }
  }

  std::ofstream outfile;
  if (write_to_fname)
  {
    std::cerr << "Writing results to " << write_to_fname << "\n";
    outfile.open(write_to_fname, outfile.app);
  }

  /*
   * Get a valid PSO object
   */
  srand(42);

  int time_max = nb_measurements * (time_between_measures - 1);
  double inertia = 0.8;
  double social = 0.1, cognition = 0.2;
  double local_refinement_box_size = 5.;
  double min_dist = 0.01;
  int dimensions = DIMENSION;
  int population_size = POPSIZE;
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
  std::cout << "PSO init" << std::endl;
  pso_constant_inertia_init(&pso, &griewank_Nd, inertia, social, cognition,
                            local_refinement_box_size, min_dist, dimensions,
                            population_size, time_max, n_trials, bounds_low,
                            bounds_high, vmin, vmax, SPACE_FILLING_DESIGN_SIZE);
  std::cout << "PSO first steps" << std::endl;

  pso_constant_inertia_first_steps(&pso, SPACE_FILLING_DESIGN_SIZE,
                                   space_filling_design);

  std::cout << "Run pso: " << nb_measurements << " x " << time_between_measures
            << "\n";
  for (int k_input = 0; k_input < nb_measurements; k_input++)
  {
    while (pso.time < k_input * time_between_measures)
    {
      pso_constant_inertia_loop(&pso);
    }

    if (do_print_outputs)
    {
      size_t lambda_p_s = pso.x_distinct_s + (pso.dimensions + 1);
      // print_vectord(pso.lambda_p, lambda_p_s, "lamdbda_p");

      double x[DIMENSION] = {0};
      printf("s(0) = %f\n", surrogate_eval(&pso, x));
    }

    size_t n_A = pso.x_distinct_s + pso.dimensions + 1;

    if (do_bench_fit_surrogate)
    {

      PerformanceTester<fit_surrogate_fun_t> perf_tester;

      auto arg_restorer = [&]()
      { pso.x_distinct_idx_of_last_batch = pso.x_distinct_s - 10; };

      std::stringstream descr;
      descr << fit_surrogate_name << "__" << n_A;

      struct perf_metrics metrics = perf_tester.perf_test(
          fit_surrogate, descr.str(), std::move(arg_restorer), &pso);
      std::stringstream ss;
      ss << fit_surrogate_name << "," << n_A << "," << metrics.cycles << ","
         << metrics.flops << "," << metrics.l3_misses << std::endl;
      auto repr = ss.str();
      std::cout << repr;
      outfile << repr;
    }

    if (do_bench_surrogate_eval)
    {
      PerformanceTester<surrogate_eval_fun_t> perf_tester;

      // nothing to restore, surrogate eval doesn't modify the pso
      auto arg_restorer = [&]() {};

      double x[DIMENSION] = {0};

      std::stringstream descr;
      descr << surrogate_eval_name << "__" << n_A;

      struct perf_metrics metrics = perf_tester.perf_test(
          &surrogate_eval, descr.str(), std::move(arg_restorer), &pso, x);

      std::stringstream ss;
      ss << fit_surrogate_name << "," << n_A << "," << metrics.cycles << ","
         << metrics.flops << "," << metrics.l3_misses << std::endl;
      auto repr = ss.str();
      std::cout << repr;
      outfile << repr;
    }
  }

  outfile.close();
}