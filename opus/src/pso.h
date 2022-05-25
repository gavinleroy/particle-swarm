#pragma once

#include <stdbool.h>
#include <sys/types.h>

#define USE_ROUNDING_BLOOM_FILTER 1

typedef double (*blackbox_fun)(double const *const, size_t);

// PSO_X : pso::pso, i:int -> x_i :double*
#define PSO_X(pso, i) ((pso)->x + (i) * (pso)->dimensions)
#define PSO_FX(pso, i) ((pso)->x_eval[i])

#define PSO_XD(pso, i) ((pso)->x_distinct + (i) * (pso)->dimensions)

#define PSO_FXD(pso, i) (pso)->x_distinct_eval[i]

// PSO_V : pso::pso, i:int -> v_i:double*
#define PSO_V(pso, i) ((pso)->v + (i) * (pso)->dimensions)

// PSO_Y : pso::pso, i:int -> y_i:double*
// PSO_FY : pso::pso, i:int -> f(y_i):double
#define PSO_Y(pso, i) ((pso)->y + (i) * (pso)->dimensions)
#define PSO_FY(pso, i) (pso)->y_eval[(t) * (pso)->population_size + (i)]

// PSO_STEP3_RAND : pso::pso, i:int (population) k:int (dimension) ->
// step3_rands_i:double*
#define PSO_STEP3_RAND(pso, i) ((pso)->step3_rands + (i) * (pso)->dimensions)

// for other inertia choices, see https://ieeexplore.ieee.org/document/6089659
struct pso_data_constant_inertia
{
  blackbox_fun f;
  // positions x_i, saved for all times
  // i.e. PSO_X(pso, i) is the current position vector x_i
  double *x;
  // f(x) for current positions
  double *x_eval;

  // PSO_V(pso,i) = v_i
  double *v;

  // PSO_Y(pso,i) := y_i := best position ever recorded for particle i
  double *y;
  // PSO_FY(pso,i) := f(y_i)
  double *y_eval;

  // ŷ = best position ever recorded over all particles
  double *y_hat;

  double *x_trial;
  double *v_trial;

  double *x_trial_best;
  double *v_trial_best;

  // Used in steps 10 and 11 in local refinement
  double *x_local;

  double *bound_low;
  double *bound_high;
  double *vmin;
  double *vmax;

  // list of all previously seen points with evaluation
  double *x_distinct;
  // starting position of last batch
  size_t x_distinct_idx_of_last_batch;
  // total size of x_distinct_s
  size_t x_distinct_s;

  // fonction evaluation at x_distinct[k]
  double *x_distinct_eval;

#ifdef USE_ROUNDING_BLOOM_FILTER
  struct rounding_bloom *bloom;
#endif

  // parameters of the surrogate
  // store the concatenation lambda_0 ... lambda_i || p_0 ... p(d+1)
  // (as this is the format of the output vector of fit_surrogate)
  double *lambda_p;

  // random numbers precomputed
  double *step3_rands; // population_size * dimensions

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

void run_pso(blackbox_fun f, double inertia, double social, double cognition,
             double local_refinement_box_size, double min_minimizer_distance,
             int dimensions, int population_size, int time_max, int n_trials,
             double *bounds_low, double *bounds_high, double *vmin,
             double *vmax, double *initial_positions);

void pso_constant_inertia_init(
    struct pso_data_constant_inertia *pso, blackbox_fun f, double inertia,
    double social, double cognition, double local_refinement_box_size,
    double min_minimizer_distance, int dimensions, int population_size,
    int time_max, int n_trials, double *bounds_low, double *bounds_high,
    double *vmin, double *vmax, double *initial_positions);
void pso_constant_inertia_first_steps(struct pso_data_constant_inertia *pso);
bool pso_constant_inertia_loop(struct pso_data_constant_inertia *pso);
