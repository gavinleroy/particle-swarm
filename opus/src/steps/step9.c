#include "step9.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if USE_ROUNDING_BLOOM_FILTER
#include "../rounding_bloom.h"
#endif

#include "fit_surrogate.h"

void step9_base(struct pso_data_constant_inertia *pso)
{
  // Refit surrogate with time = t+1

  // first update the set of distinct points
  size_t t = pso->time;

  // save position of new batch of points in x_distinct
  // TODO: pass as an argument, no need to save it!
  pso->new_x_distinct_at_t[t] = pso->x_distinct_s;


  size_t cur_x_distinct_batch_start = pso->new_x_distinct_at_t[pso->time];
#if USE_ROUNDING_BLOOM_FILTER
  for (int i = 0; i < pso->population_size; i++)
  {
    // add and check proximity to previous points
    if (!rounding_bloom_check_add(pso->bloom, pso->dimensions,
                                  PSO_X(pso, i), 1))
    {
      // copy point to x_distinct
      memcpy(PSO_XD(pso, pso->x_distinct_s), PSO_X(pso, i), pso->dimensions * sizeof(double));
      pso->x_distinct_eval[pso->x_distinct_s] = pso->x_eval[i];
            
      pso->x_distinct_s++;
    }
  }
#else
  // naive implementation with distance computation
  fprintf(stderr, "Not implemented.\n");
  exit(1);
#endif

  if (fit_surrogate(pso) < 0)
  {
    fprintf(stderr, "ERROR: Failed to fit surrogate\n");
    exit(1);
  }

}

void step9_optimized(struct pso_data_constant_inertia *pso) { step9_base(pso); }
