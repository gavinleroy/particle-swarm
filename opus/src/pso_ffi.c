#include "pso.h"

double *get_current_particle_positions(struct pso_data_constant_inertia *pso)
{
  return pso->x;
}

double *get_current_particle_velocities(struct pso_data_constant_inertia *pso)
{
  return pso->v;
}