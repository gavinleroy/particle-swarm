#include "step9.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if USE_ROUNDING_BLOOM_FILTER
#include "../rounding_bloom.h"
#endif

#include "fit_surrogate.h"

#include "../logging.h"

void step9_base(struct pso_data_constant_inertia *pso)
{
  // Refit surrogate with time = t+1

  // first update the set of distinct points
  size_t t = pso->time;

#if USE_ROUNDING_BLOOM_FILTER
  for (int i = 0; i < pso->population_size; i++)
  {
    // add and check proximity to previous points
    if (!rounding_bloom_check_add(pso->bloom, pso->dimensions, PSO_X(pso, i),
                                  1))
    {
      // copy point to x_distinct
      memcpy(PSO_XD(pso, pso->x_distinct_s), PSO_X(pso, i),
             pso->dimensions * sizeof(double));
      pso->x_distinct_eval[pso->x_distinct_s] = pso->x_eval[i];

      pso->x_distinct_s++;
    }
  }
#else
  // naive implementation with distance computation
  fprintf(stderr, "Not implemented.\n");
  exit(1);
#endif

  TIMING_INIT();
  if (fit_surrogate(pso) < 0)
  {
    fprintf(stderr, "ERROR: Failed to fit surrogate\n");
    exit(1);
  }
  TIMING_STEP("fit_surrogate", FIT_SURROGATE_VERSION, pso->time);
}

void step9_optimized(struct pso_data_constant_inertia *pso) { step9_base(pso); }
