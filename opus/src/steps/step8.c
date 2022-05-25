#include "step8.h"

void step8_base(struct pso_data_constant_inertia *pso)
{
  // Update the best positions per particle and overall
  for (int i = 0; i < pso->population_size; i++)
  {
    double x_eval = PSO_FX(pso, i);

    if (x_eval < pso->y_eval[i])
    {
      // y_i <- x_i
      for (int k = 0; k < pso->dimensions; k++)
      {
        PSO_Y(pso, i)[k] = PSO_X(pso, i)[k];
      }
      pso->y_eval[i] = x_eval;

      // is x_i(t+1) better than ŷ(t) ?
      if (x_eval < pso->y_hat_eval)
      {
        pso->y_hat = PSO_Y(pso, i);
        pso->y_hat_eval = x_eval;
      }
    }
  }

  pso->time++;
}

void step8_optimized(struct pso_data_constant_inertia *pso) { step8_base(pso); }
