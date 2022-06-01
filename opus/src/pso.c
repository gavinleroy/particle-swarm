#include "pso.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gaussian_elimination_solver.h"
#include "helpers.h"
#include "local_refinement.h"

#include "steps/steps.h"

#if DISTINCTIVENESS_CHECK_TYPE == 2
#include "rounding_bloom.h"
#endif

#define DEBUG_TRIALS 0
#define DEBUG_SURROGATE 0

#include "logging.h"


struct pso_data_constant_inertia *alloc_pso_data_constant_inertia()
{
  return malloc(sizeof(struct pso_data_constant_inertia));
}

void random_number_generation(struct pso_data_constant_inertia *pso)
{
  // Step 3
  pso->step3_rands =
      malloc(pso->population_size * pso->dimensions * sizeof(double));
  for (int i = 0; i < pso->population_size; i++)
  {
    for (int k = 0; k < pso->dimensions; k++)
    {
      PSO_STEP3_RAND(pso, i)
      [k] = rand_between(pso->bound_low[k], pso->bound_high[k]);
    }
  }

  // Step 6
  pso->step6_rands = aligned_alloc(32, pso->population_size * pso->n_trials *
                                           pso->dimensions * sizeof(double));

  for (int i = 0; i < pso->population_size; i++)
  {
    for (int l = 0; l < pso->n_trials; l++)
    {
      double *row_ptr = PSO_STEP6_RAND(pso, i, l);
      for (int j = 0; j < pso->dimensions; j++)
      {
        double w1 = (double)rand() / RAND_MAX;
        double w2 = (double)rand() / RAND_MAX;
        row_ptr[2 * j] = w1;
        row_ptr[2 * j + 1] = w2;
      }
    }
  }
}

void pso_constant_inertia_init(
    struct pso_data_constant_inertia *pso, blackbox_fun f, double inertia,
    double social, double cognition, double local_refinement_box_size,
    double min_dist, int dimensions, int population_size,
    int time_max, int n_trials, double *bounds_low, double *bounds_high,
    double *vmin, double *vmax, size_t sfd_size)
{
  pso->f = f;
  pso->inertia = inertia;
  pso->dimensions = dimensions;
  pso->social = social, pso->cognition = cognition;
  pso->local_refinement_box_size = local_refinement_box_size;
  pso->min_dist2 = min_dist * min_dist;
  pso->population_size = population_size;
  pso->time_max = time_max, pso->n_trials = n_trials;

  pso->time = 0;

  pso->x =
      (double *)malloc(pso->population_size * pso->dimensions * sizeof(double));
  pso->x_eval = malloc(pso->population_size * sizeof(double));

  pso->v = malloc(pso->population_size * pso->dimensions * sizeof(double));

  pso->y = malloc(pso->population_size * pso->dimensions * sizeof(double));

  pso->y_eval = malloc(pso->population_size * sizeof(double *));

  // yhat will be a pointer in another array
  pso->y_hat = NULL;

  // size of one vector, rounded up to be 32B aligned
  size_t size_of_one_vec_32 = (((pso->dimensions * sizeof(double)) + 31) & -32);

  pso->v_trial = aligned_alloc(32, size_of_one_vec_32);
  pso->x_trial = aligned_alloc(32, size_of_one_vec_32);

  pso->v_trial_best = malloc(pso->dimensions * sizeof(double));
  pso->x_trial_best = malloc(pso->dimensions * sizeof(double));

  pso->x_local = malloc(pso->dimensions * sizeof(double));

  pso->bound_low =
      (double *)aligned_alloc(32, size_of_one_vec_32);
  pso->bound_high =
      (double *)aligned_alloc(32, size_of_one_vec_32);

  pso->vmin = (double *)aligned_alloc(32, size_of_one_vec_32);
  pso->vmax = (double *)aligned_alloc(32, size_of_one_vec_32);

  for (int j = 0; j < dimensions; j++)
  {
    pso->vmin[j] = vmin[j];
    pso->vmax[j] = vmax[j];
  }

  // pso->population_size particles per iterations, and 1 for the local minimizer
  // + initial space filling design
  size_t x_distinct_max_nb = pso->time_max * (pso->population_size + 1) + sfd_size;
  pso->x_distinct =
      malloc(x_distinct_max_nb * pso->dimensions * sizeof(double));
  pso->x_distinct_idx_of_last_batch = 0;
  pso->x_distinct_s = 0;

  pso->x_distinct_eval = malloc(x_distinct_max_nb * sizeof(double));

#if DISTINCTIVENESS_CHECK_TYPE == 0
  // Unconditionnal accept ; nothing to allocate
#elif DISTINCTIVENESS_CHECK_TYPE == 1
  // Naive distance calculation ; nothing to allocate
#elif DISTINCTIVENESS_CHECK_TYPE == 2
  // Bloom filter
  pso->bloom = malloc(sizeof(struct rounding_bloom));
  int bloom_entries = pso->time_max * pso->population_size;
  if (bloom_entries < 1000)
    bloom_entries = 1000;
  double bloom_false_pos_rate = 0.01;
  double bloom_rounding_eps = min_dist;
  rounding_bloom_init(pso->bloom, bloom_entries, bloom_false_pos_rate,
                      bloom_rounding_eps, dimensions, bounds_low);
#endif

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // max: pop_size * time
  //  TODO: add the refinement points + 1 * time
  //        add the space filling design +?
  size_t max_n_phi = x_distinct_max_nb;
  size_t n_P = pso->dimensions + 1;
  prealloc_fit_surrogate(max_n_phi, n_P);

  // alloc maximum possible size: max_n_phi for lambda and d+1 for P
  size_t lambda_p_s = max_n_phi + (pso->dimensions + 1);
  pso->lambda_p = malloc(lambda_p_s * sizeof(double));

  // setup bounds in space
  for (int k = 0; k < pso->dimensions; k++)
    pso->bound_low[k] = bounds_low[k];
  for (int k = 0; k < pso->dimensions; k++)
    pso->bound_high[k] = bounds_high[k];
  // setup bounds on velocity
  for (int k = 0; k < pso->dimensions; k++)
    pso->vmin[k] = vmin[k];
  for (int k = 0; k < pso->dimensions; k++)
    pso->vmax[k] = vmax[k];

  // precomputed random numbers
  random_number_generation(pso);
}

void pso_constant_inertia_first_steps(struct pso_data_constant_inertia *pso, size_t sfd_size, double * space_filling_design)
{
  step1_2(pso, sfd_size, space_filling_design);

  step3_optimized(pso);

  step4_optimized(pso);
}

bool pso_constant_inertia_loop(struct pso_data_constant_inertia *pso)
{
  step5_optimized(pso);

  step6_optimized(pso);

  step7_optimized(pso);

  step8_optimized(pso);

  step9_optimized(pso);

  step10_optimized(pso);

  step11_optimized(pso);

  return (pso->time < pso->time_max - 1);
}

void run_pso(blackbox_fun f, double inertia, double social, double cognition,
             double local_refinement_box_size, double min_minimizer_distance,
             int dimensions, int population_size, int time_max, int n_trials,
             double *bounds_low, double *bounds_high, double *vmin,
             double *vmax, size_t sfd_size, double * space_filling_design)
{
  struct pso_data_constant_inertia pso;
  pso_constant_inertia_init(
      &pso, f, inertia, social, cognition, local_refinement_box_size,
      min_minimizer_distance, dimensions, population_size, time_max, n_trials,
      bounds_low, bounds_high, vmin, vmax, sfd_size);

  pso_constant_inertia_first_steps(&pso, sfd_size, space_filling_design);

  printf("t=%d  ŷ=[", pso.time);
  for (int j = 0; j < dimensions; j++)
  {
    printf("%f", pso.y_hat[j]);
    if (j < dimensions - 1)
      printf(", ");
  }
  printf("]  f(ŷ)=%f\n", pso.y_hat_eval);

  while (pso.time < pso.time_max - 1)
  {
    pso_constant_inertia_loop(&pso);

    printf("t=%d  ŷ=[", pso.time);
    for (int j = 0; j < dimensions; j++)
    {
      printf("%f", pso.y_hat[j]);
      if (j < dimensions - 1)
        printf(", ");
    }
    printf("]  f(ŷ)=%f\n", pso.y_hat_eval);
  }
}
