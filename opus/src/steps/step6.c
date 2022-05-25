#include "step6.h"

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

void step6_optimized(struct pso_data_constant_inertia *pso) { step6_base(pso); }
