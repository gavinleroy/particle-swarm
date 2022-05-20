#include "step7.h"

void step7_base(struct pso_data_constant_inertia *pso)
{
  // Evaluate swarm positions
  for (int i = 0; i < pso->population_size; i++)
  {
    PSO_FX(pso, pso->time + 1, i) = pso->f(PSO_X(pso, pso->time + 1, i));
  }
}

void step7_optimized(struct pso_data_constant_inertia *pso) { step7_base(pso); }
