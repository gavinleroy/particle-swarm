#include "step9.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../distincts.h"

#include "fit_surrogate.h"

#include "../logging.h"

void step9_base(struct pso_data_constant_inertia *pso)
{
  // Refit surrogate with time = t+1

  // first update the set of distinct points
  size_t t = pso->time;

  for (int i = 0; i < pso->population_size; i++)
  {
    // add and check proximity to previous points
    add_to_distincts_if_distinct(pso, PSO_X(pso, i), pso->x_eval[i]);
  }

  TIMING_INIT();
  if (fit_surrogate(pso) < 0)
  {
    fprintf(stderr, "ERROR: Failed to fit surrogate\n");
    exit(1);
  }
  TIMING_STEP("fit_surrogate", STR(FIT_SURROGATE_VERSION), pso->time);
}

void step9_optimized(struct pso_data_constant_inertia *pso) { step9_base(pso); }
