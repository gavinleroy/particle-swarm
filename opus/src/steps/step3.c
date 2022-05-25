#include "step3.h"

#include <immintrin.h>

#include "../helpers.h"

void step3_base(struct pso_data_constant_inertia *pso)
{
  // Step 3. Initialize particle velocities
  for (int i = 0; i < pso->population_size; i++)
  {
    for (int k = 0; k < pso->dimensions; k++)
    {
      double uk = rand_between(pso->bound_low[k], pso->bound_high[k]);
      PSO_V(pso, i)[k] = 1. / 2 * (uk - PSO_X(pso, i)[k]);
    }
  }
}

/*
  + precomputed 1/2
  + precomputed random values
*/
void step3_opt1(struct pso_data_constant_inertia *pso)
{
  double uk;
  // Step 3. Initialize particle velocities
  for (int i = 0; i < pso->population_size; i++)
  {
    for (int k = 0; k < pso->dimensions; k++)
    {
      uk = PSO_STEP3_RAND(pso, i)[k];
      PSO_V(pso, i)[k] = 0.5 * (uk - PSO_X(pso, i)[k]);
    }
  }
}

/*
  precomputed 1/2
  precomputed random values
  + scalar replacements
*/
void step3_opt2(struct pso_data_constant_inertia *pso)
{
  double uk;
  double *pso_step3_rands_i_ptr;
  double *pso_x_0_i_ptr;
  double *pso_v_i_ptr;
  // Step 3. Initialize particle velocities
  for (int i = 0; i < pso->population_size; i++)
  {
    pso_step3_rands_i_ptr = PSO_STEP3_RAND(pso, i);
    pso_x_0_i_ptr = PSO_X(pso, i);
    pso_v_i_ptr = PSO_V(pso, i);
    for (int k = 0; k < pso->dimensions; k++)
    {
      uk = pso_step3_rands_i_ptr[k];
      pso_v_i_ptr[k] = 0.5 * (pso_step3_rands_i_ptr[k] - pso_x_0_i_ptr[k]);
    }
  }
}

/*
  precomputed 1/2
  precomputed random values
  scalar replacements
  + loop unrolling
*/
void step3_opt3(struct pso_data_constant_inertia *pso)
{
  double uk;
  double *pso_step3_rands_i_ptr;
  double *pso_x_0_i_ptr;
  double *pso_v_i_ptr;
  // Step 3. Initialize particle velocities
  for (int i = 0; i < pso->population_size; i++)
  {
    pso_step3_rands_i_ptr = PSO_STEP3_RAND(pso, i);
    pso_x_0_i_ptr = PSO_X(pso, i);
    pso_v_i_ptr = PSO_V(pso, i);
    int k = 0;
    for (; k < pso->dimensions - 3; k += 4)
    {
      uk = pso_step3_rands_i_ptr[k];
      pso_v_i_ptr[k] = 0.5 * (uk - pso_x_0_i_ptr[k]);

      uk = pso_step3_rands_i_ptr[k + 1];
      pso_v_i_ptr[k + 1] = 0.5 * (uk - pso_x_0_i_ptr[k + 1]);

      uk = pso_step3_rands_i_ptr[k + 2];
      pso_v_i_ptr[k + 2] = 0.5 * (uk - pso_x_0_i_ptr[k + 2]);

      uk = pso_step3_rands_i_ptr[k + 3];
      pso_v_i_ptr[k + 3] = 0.5 * (uk - pso_x_0_i_ptr[k + 3]);
    }

    // leftover
    for (; k < pso->dimensions; k++)
    {
      uk = pso_step3_rands_i_ptr[k];
      pso_v_i_ptr[k] = 0.5 * (pso_step3_rands_i_ptr[k] - pso_x_0_i_ptr[k]);
    }
  }
}

/*
  precomputed 1/2
  precomputed random values
  scalar replacements
  loop unrolling
  + vectorization
*/
void step3_opt4(struct pso_data_constant_inertia *pso)
{
  double *pso_step3_rands_i_ptr;
  double *pso_x_0_i_ptr;
  double *pso_v_i_ptr;

  __m256d u_i, pso_x_0_i, sub_i, pso_v_i, pso_full_half;

  pso_full_half = _mm256_set1_pd(0.5);

  // Step 3. Initialize particle velocities
  for (int i = 0; i < pso->population_size; i++)
  {
    pso_step3_rands_i_ptr = PSO_STEP3_RAND(pso, i);
    pso_x_0_i_ptr = PSO_X(pso, i);
    pso_v_i_ptr = PSO_V(pso, i);
    int k = 0;
    for (; k < pso->dimensions - 3; k += 4)
    {
      u_i = _mm256_loadu_pd(pso_step3_rands_i_ptr + k);
      pso_x_0_i = _mm256_loadu_pd(pso_x_0_i_ptr + k);

      sub_i = _mm256_sub_pd(u_i, pso_x_0_i);
      pso_v_i = _mm256_mul_pd(pso_full_half, sub_i);

      _mm256_storeu_pd(pso_v_i_ptr, pso_v_i);
    }

    // leftover
    for (; k < pso->dimensions; k++)
    {
      pso_v_i_ptr[k] = 0.5 * (pso_step3_rands_i_ptr[k] - pso_x_0_i_ptr[k]);
    }
  }
}

void step3_optimized(struct pso_data_constant_inertia *pso) { step3_base(pso); }
