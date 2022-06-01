#include "step6.h"
#include <immintrin.h>
#include <string.h>

#include "float.h"
#include "stdlib.h"

#include "surrogate_eval.h"

static double clamp(double v, double lo, double hi)
{
  if (v < lo)
  {
    return lo;
  }
  else if (v > hi)
  {
    return hi;
  }
  else
  {
    return v;
  }
}

void step6_base(struct pso_data_constant_inertia *pso)
{
  // Determine new particle positions

  for (int i = 0; i < pso->population_size; i++)
  {

#if DEBUG_TRIALS
    printf("step6:start trials for particle %d\n", i);
#endif

    double x_trial_best_seval = DBL_MAX;
    for (int l = 0; l < pso->n_trials; l++)
    {

      for (int j = 0; j < pso->dimensions; j++)
      {
        // compute v_i(t+1) from v_i(t)

        double w1 = (double)rand() / RAND_MAX;
        double w2 = (double)rand() / RAND_MAX;
        double v = pso->inertia * PSO_V(pso, i)[j] +
                   pso->cognition * w1 * (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j]) +
                   pso->social * w2 * (pso->y_hat[j] - PSO_X(pso, i)[j]);

        pso->v_trial[j] = clamp(v, pso->vmin[j], pso->vmax[j]);

        pso->x_trial[j] = clamp(PSO_X(pso, i)[j] + pso->v_trial[j],
                                pso->bound_low[j], pso->bound_high[j]);
      }

      double x_trial_seval = surrogate_eval(pso, pso->x_trial);

#if DEBUG_TRIALS
      char trial_name[16] = {0};
      snprintf(trial_name, sizeof(trial_name), "x%d_trial%d", i, l);
      print_vectord(pso->x_trial, pso->dimensions, trial_name);
      printf("score of trial %d = %f", l, x_trial_seval);
#endif

      if (x_trial_seval < x_trial_best_seval)
      {
#if DEBUG_TRIALS
        printf(" (new best!)");
#endif
        // keep x_trial as x_trial_best by swapping the two buffers: the new
        // x_trial will get overwritten in the next iteration
        x_trial_best_seval = x_trial_seval;

        double *t;

        t = pso->x_trial;
        pso->x_trial = pso->x_trial_best;
        pso->x_trial_best = t;

        t = pso->v_trial;
        pso->v_trial = pso->v_trial_best;
        pso->v_trial_best = t;
      }
#if DEBUG_TRIALS
      printf("\n");
#endif
    }

    // set next position and update velocity

    for (int j = 0; j < pso->dimensions; j++)
    {
      PSO_X(pso, i)[j] = pso->x_trial_best[j];
      PSO_V(pso, i)[j] = pso->v_trial_best[j];
    }

#if DEBUG_TRIALS
    printf("set x[t+1][i]:\n");
    char new_x_name[16] = {0};
    snprintf(new_x_name, sizeof(new_x_name), "x%d_(t=%d)", i, t + 1);
    print_vectord(PSO_X(pso, i), pso->dimensions, new_x_name);
#endif
  }
}

/*
 * use precomputed RAND variables
 */
void step6_opt1(struct pso_data_constant_inertia *pso)
{
  // Determine new particle positions
  int pop_size = pso->population_size;
  int dim = pso->dimensions;
  int n_trials = pso->n_trials;

  for (int i = 0; i < pop_size; i++)
  {
    double x_trial_best_seval = DBL_MAX;
    for (int l = 0; l < n_trials; l++)
    {

      double *row_ptr = PSO_STEP6_RAND(pso, i, l);

      for (int j = 0; j < dim; j++)
      {
        // compute v_i(t+1) from v_i(t)
        double w1 = row_ptr[2 * j];     // (double)rand() / RAND_MAX;
        double w2 = row_ptr[2 * j + 1]; // (double)rand() / RAND_MAX;

        double v = pso->inertia * PSO_V(pso, i)[j] +
                   pso->cognition * w1 * (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j]) +
                   pso->social * w2 * (pso->y_hat[j] - PSO_X(pso, i)[j]);

        pso->v_trial[j] = clamp(v, pso->vmin[j], pso->vmax[j]);

        pso->x_trial[j] = clamp(PSO_X(pso, i)[j] + pso->v_trial[j],
                                pso->bound_low[j], pso->bound_high[j]);
      }

      double x_trial_seval = surrogate_eval(pso, pso->x_trial);

      if (x_trial_seval < x_trial_best_seval)
      {
        // keep x_trial as x_trial_best by swapping the two buffers: the new
        // x_trial will get overwritten in the next iteration
        x_trial_best_seval = x_trial_seval;

        double *t;

        t = pso->x_trial;
        pso->x_trial = pso->x_trial_best;
        pso->x_trial_best = t;

        t = pso->v_trial;
        pso->v_trial = pso->v_trial_best;
        pso->v_trial_best = t;
      }
    }

    // set next position and update velocity

    for (int j = 0; j < pso->dimensions; j++)
    {
      PSO_X(pso, i)[j] = pso->x_trial_best[j];
      PSO_V(pso, i)[j] = pso->v_trial_best[j];
    }
  }
}

/*
 * inner-loop unroll by 4
 */
void step6_opt2(struct pso_data_constant_inertia *pso)
{
  // Determine new particle positions
  int pop_size = pso->population_size;
  int dim = pso->dimensions;
  int n_trials = pso->n_trials;

  for (int i = 0; i < pop_size; i++)
  {
    double x_trial_best_seval = DBL_MAX;
    for (int l = 0; l < n_trials; l++)
    {

      double *row_ptr = PSO_STEP6_RAND(pso, i, l);

      int j = 0;
      for (; j < dim; j += 4)
      {
        // compute v_i(t+1) from v_i(t)
        int j2 = j * 2;
        double w0_1 = row_ptr[j2];     // (double)rand() / RAND_MAX;
        double w0_2 = row_ptr[j2 + 1]; // (double)rand() / RAND_MAX;
        double w1_1 = row_ptr[j2 + 2];
        double w1_2 = row_ptr[j2 + 3];
        double w2_1 = row_ptr[j2 + 4];
        double w2_2 = row_ptr[j2 + 5];
        double w3_1 = row_ptr[j2 + 6];
        double w3_2 = row_ptr[j2 + 7];

        double v0 =
            pso->inertia * PSO_V(pso, i)[j] +
            pso->cognition * w0_1 * (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j]) +
            pso->social * w0_2 * (pso->y_hat[j] - PSO_X(pso, i)[j]);

        double v1 =
            pso->inertia * PSO_V(pso, i)[j + 1] +
            pso->cognition * w1_1 *
                (PSO_Y(pso, i)[j + 1] - PSO_X(pso, i)[j + 1]) +
            pso->social * w1_2 * (pso->y_hat[j + 1] - PSO_X(pso, i)[j + 1]);

        double v2 =
            pso->inertia * PSO_V(pso, i)[j + 2] +
            pso->cognition * w2_1 *
                (PSO_Y(pso, i)[j + 2] - PSO_X(pso, i)[j + 2]) +
            pso->social * w2_2 * (pso->y_hat[j + 2] - PSO_X(pso, i)[j + 2]);

        double v3 =
            pso->inertia * PSO_V(pso, i)[j + 3] +
            pso->cognition * w3_1 *
                (PSO_Y(pso, i)[j + 3] - PSO_X(pso, i)[j + 3]) +
            pso->social * w3_2 * (pso->y_hat[j + 3] - PSO_X(pso, i)[j + 3]);

        pso->v_trial[j] = clamp(v0, pso->vmin[j], pso->vmax[j]);
        pso->v_trial[j + 1] = clamp(v1, pso->vmin[j + 1], pso->vmax[j + 1]);
        pso->v_trial[j + 2] = clamp(v2, pso->vmin[j + 2], pso->vmax[j + 2]);
        pso->v_trial[j + 3] = clamp(v3, pso->vmin[j + 3], pso->vmax[j + 3]);

        pso->x_trial[j] = clamp(PSO_X(pso, i)[j] + pso->v_trial[j],
                                pso->bound_low[j], pso->bound_high[j]);
        pso->x_trial[j + 1] =
            clamp(PSO_X(pso, i)[j + 1] + pso->v_trial[j + 1],
                  pso->bound_low[j + 1], pso->bound_high[j + 1]);
        pso->x_trial[j + 2] =
            clamp(PSO_X(pso, i)[j + 2] + pso->v_trial[j + 2],
                  pso->bound_low[j + 2], pso->bound_high[j + 2]);
        pso->x_trial[j + 3] =
            clamp(PSO_X(pso, i)[j + 3] + pso->v_trial[j + 3],
                  pso->bound_low[j + 3], pso->bound_high[j + 3]);
      }

      for (; j < dim; j++)
      {
        // compute v_i(t+1) from v_i(t)
        double w1 = row_ptr[2 * j];     // (double)rand() / RAND_MAX;
        double w2 = row_ptr[2 * j + 1]; // (double)rand() / RAND_MAX;

        double v = pso->inertia * PSO_V(pso, i)[j] +
                   pso->cognition * w1 * (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j]) +
                   pso->social * w2 * (pso->y_hat[j] - PSO_X(pso, i)[j]);

        pso->v_trial[j] = clamp(v, pso->vmin[j], pso->vmax[j]);

        pso->x_trial[j] = clamp(PSO_X(pso, i)[j] + pso->v_trial[j],
                                pso->bound_low[j], pso->bound_high[j]);
      }

      double x_trial_seval = surrogate_eval(pso, pso->x_trial);

      if (x_trial_seval < x_trial_best_seval)
      {
        // keep x_trial as x_trial_best by swapping the two buffers: the new
        // x_trial will get overwritten in the next iteration
        x_trial_best_seval = x_trial_seval;

        double *t;

        t = pso->x_trial;
        pso->x_trial = pso->x_trial_best;
        pso->x_trial_best = t;

        t = pso->v_trial;
        pso->v_trial = pso->v_trial_best;
        pso->v_trial_best = t;
      }
    }

    // set next position and update velocity

    for (int j = 0; j < pso->dimensions; j++)
    {
      PSO_X(pso, i)[j] = pso->x_trial_best[j];
      PSO_V(pso, i)[j] = pso->v_trial_best[j];
    }
  }
}

/*
 * add vector intrinsics
 */
void step6_opt3(struct pso_data_constant_inertia *pso)
{
  // Determine new particle positions
  int pop_size = pso->population_size;
  int dim = pso->dimensions;
  int n_trials = pso->n_trials;

  for (int i = 0; i < pop_size; i++)
  {
    double x_trial_best_seval = DBL_MAX;
    for (int l = 0; l < n_trials; l++)
    {

      double *row_ptr = PSO_STEP6_RAND(pso, i, l);

      int j = 0;
      for (; j < dim - 3; j += 4)
      {
        // compute v_i(t+1) from v_i(t)
        int j2 = j * 2;

        __m256d inertia, v, t0;
        __m256d cognition, w1, pso_y, pso_x, t10, t11;
        __m256d social, w2, y_hat, t20, t21;

        // pso->inertia * PSO_V(pso, i)[j]
        inertia = _mm256_set1_pd(pso->inertia);
        v = _mm256_load_pd(PSO_V(pso, i) + j);

        // pso->cognition * w0_1 * (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j])
        cognition = _mm256_set1_pd(pso->cognition);
        w1 = _mm256_load_pd(row_ptr + j2);
        pso_y = _mm256_load_pd(PSO_Y(pso, i) + j);
        pso_x = _mm256_load_pd(PSO_X(pso, i) + j);

        // pso->social * w2 * (pso->y_hat[j] - PSO_X(pso, i)[j])
        social = _mm256_set1_pd(pso->social);
        w2 = _mm256_load_pd(row_ptr + j2 + 4);
        y_hat = _mm256_load_pd(pso->y_hat + j);

        t0 = _mm256_mul_pd(inertia, v); // pso->inertia * PSO_V(pso, i)[j]

        t10 = _mm256_mul_pd(cognition, w1); // pso->cognition * w1
        t11 =
            _mm256_sub_pd(pso_y, pso_x); // PSO_Y(pso, i)[j] - PSO_X(pso, i)[j]
        t10 = _mm256_mul_pd(
            t10,
            t11); // pso->cognition * w1 * (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j])

        v = _mm256_add_pd(
            t0, t10); // pso->inertia * PSO_V(pso, i)[j] + pso->cognition * w1 *
                      // (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j])

        t20 = _mm256_mul_pd(social, w2);   // pso->social * w2
        t21 = _mm256_sub_pd(y_hat, pso_x); // pso->y_hat[j] - PSO_X(pso, i)[j]
        t20 = _mm256_mul_pd(
            t20, t21); // pso->social * w2 * (pso->y_hat[j] - PSO_X(pso, i)[j]);

        v = _mm256_add_pd(v, t20);

        // Clamps
        __m256d vmin, vmax, pso_bound_low, pso_bound_high;
        vmin = _mm256_load_pd(pso->vmin + j);
        vmax = _mm256_load_pd(pso->vmax + j);
        pso_bound_low = _mm256_load_pd(pso->bound_low + j);
        pso_bound_high = _mm256_load_pd(pso->bound_high + j);

        // clamp(v, pso->vmin[j], pso->vmax[j]);
        v = _mm256_max_pd(vmin, v); // lower bound
        v = _mm256_min_pd(vmax, v); // upper bound

        //        clamp(PSO_X(pso, i)[j] + pso->v_trial[j],
        //              pso->bound_low[j], pso->bound_high[j]);
        pso_x = _mm256_add_pd(pso_x, v);
        pso_x = _mm256_max_pd(pso_bound_low, pso_x);  // lower bound
        pso_x = _mm256_min_pd(pso_bound_high, pso_x); // upper bound

        _mm256_store_pd(pso->v_trial + j, v); // pso->v_trial[j] = clamp(...)
        _mm256_storeu_pd(pso->x_trial + j, pso_x);
      }

      for (; j < dim; j++)
      {
        // compute v_i(t+1) from v_i(t)
        double w1 = row_ptr[2 * j];     // (double)rand() / RAND_MAX;
        double w2 = row_ptr[2 * j + 1]; // (double)rand() / RAND_MAX;

        double v = pso->inertia * PSO_V(pso, i)[j] +
                   pso->cognition * w1 * (PSO_Y(pso, i)[j] - PSO_X(pso, i)[j]) +
                   pso->social * w2 * (pso->y_hat[j] - PSO_X(pso, i)[j]);

        pso->v_trial[j] = clamp(v, pso->vmin[j], pso->vmax[j]);

        pso->x_trial[j] = clamp(PSO_X(pso, i)[j] + pso->v_trial[j],
                                pso->bound_low[j], pso->bound_high[j]);
      }

      double x_trial_seval = surrogate_eval(pso, pso->x_trial);

      if (x_trial_seval < x_trial_best_seval)
      {
        // keep x_trial as x_trial_best by swapping the two buffers: the new
        // x_trial will get overwritten in the next iteration
        x_trial_best_seval = x_trial_seval;

        double *t;

        t = pso->x_trial;
        pso->x_trial = pso->x_trial_best;
        pso->x_trial_best = t;

        t = pso->v_trial;
        pso->v_trial = pso->v_trial_best;
        pso->v_trial_best = t;
      }
    }

    // set next position and update velocity
    memcpy(PSO_X(pso, i), pso->x_trial_best, pso->dimensions * sizeof(double));
    memcpy(PSO_V(pso, i), pso->v_trial_best, pso->dimensions * sizeof(double));
  }
}

void step6_optimized(struct pso_data_constant_inertia *pso)
{
//    step6_base(pso);
  //  step6_opt1(pso);
//    step6_opt2(pso);
  step6_opt3(pso);
}
