#include "step11.h"

#include <math.h>
#include <string.h>

#include "../helpers.h"

static int
is_far_from_previous_evaluations(struct pso_data_constant_inertia const *pso,
                                 double *x, double min_dist)
{
  // - find way to check "if minimizer of surrogate is far from previous points"
  //    + see https://en.wikipedia.org/wiki/Nearest_neighbor_search
  //    + in high dim naive search can be best

  double delta2 = min_dist * min_dist;

  // TODO: switch to bloom filter? Then cannot specify different delta
  //  but is that really necesary?
  //  Check previous distinct positions (particles and local minimizers)
  for (int i = 0; i < pso->x_distinct_s; i++)
  {
    double d2 = dist2(pso->dimensions, x, PSO_XD(pso, i));
    if (d2 < delta2)
    {
      return 0;
    }
  }

  return 1;
}

void step11_base(struct pso_data_constant_inertia *pso)
{
  double *x_local = pso->x_local;
  // Determine if minimizer of surrogate is far from previous points
  if (is_far_from_previous_evaluations(pso, x_local,
                                       pso->min_minimizer_distance))
  {
    double x_local_eval = pso->f(x_local);

    double *x_local_in_xdistinct =
        pso->x_distinct + (pso->x_distinct_s * pso->dimensions);

    // Add new refinement point (and evaluation) to list of distinct evaluation
    // positions
    memcpy(x_local_in_xdistinct, x_local, pso->dimensions * sizeof(double));
    pso->x_distinct_eval[pso->x_distinct_s] = x_local_eval;
    pso->x_distinct_s++;

    // update overall best if applicable
    if (x_local_eval < pso->y_hat_eval)
    {
      pso->y_hat = x_local_in_xdistinct;
      pso->y_hat_eval = x_local_eval;
    }
  }
}

void step11_optimized(struct pso_data_constant_inertia *pso)
{
  step11_base(pso);
}
