#include "step5.h"

#include "stdio.h"
#include "stdlib.h"
#include <string.h>

#include "../gaussian_elimination_solver.h"
#include "math.h"

#include "fit_surrogate.h"

#if USE_ROUNDING_BLOOM_FILTER
#include "../rounding_bloom.h"
#endif

#include "../logging.h"

void step5_base(struct pso_data_constant_inertia *pso)
{
  // Step 5.
  // Fit surrogate
  // f already evaluated on x[0..t][0..i-1]

  // Build set of distinct points: add latest x positions
  for (size_t i = 0; i < pso->population_size; i++)
  {
#if USE_ROUNDING_BLOOM_FILTER
    // add and check proximity to previous points
    if (!rounding_bloom_check_add(pso->bloom, pso->dimensions, PSO_X(pso, i),
                                  1))
    {
      // copy point and value to x_distinct
      memcpy(PSO_XD(pso, pso->x_distinct_s), PSO_X(pso, i),
             pso->dimensions * sizeof(double));
      pso->x_distinct_eval[pso->x_distinct_s] = pso->x_eval[i];

      pso->x_distinct_s++;
    }
#else
    // naive implementation with distance computation
    fprintf(stderr, "Not implemented.\n");
    exit(1);
#endif
  }

  TIMING_INIT();
  if (fit_surrogate(pso) < 0)
  {
    fprintf(stderr, "ERROR: Failed to fit surrogate\n");
    exit(1);
  }
  TIMING_STEP("fit_surrogate", FIT_SURROGATE_VERSION, pso->time);
}

void step5_optimized(struct pso_data_constant_inertia *pso) { step5_base(pso); }
