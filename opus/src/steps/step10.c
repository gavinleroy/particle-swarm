#include "step10.h"

#include "../local_refinement.h"

#include "surrogate_eval.h"

// surrogate inferface taking a void pointer
// to pass to local_optimizer
double surrogate_eval_void(double const *x, void const *args)
{
  struct pso_data_constant_inertia const *pso =
      (struct pso_data_constant_inertia const *)args;
  return surrogate_eval(pso, x);
}

void step10_base(struct pso_data_constant_inertia *pso)
{
  // Local refinement
  local_optimization(&surrogate_eval_void, pso->dimensions, pso->y_hat,
                     pso->local_refinement_box_size, pso->bound_low,
                     pso->bound_high, pso, pso->x_local);
}

void step10_optimized(struct pso_data_constant_inertia *pso)
{
  step10_base(pso);
}
