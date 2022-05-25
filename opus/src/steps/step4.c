#include <string.h>

#include "step4.h"

// Step 4. Initialise y, y_eval, and x_eval for each particle
// with distinct position
void step4_base(struct pso_data_constant_inertia *pso)
{
  for (int i = 0; i < pso->population_size; i++)
  {
    for (int k = 0; k < pso->dimensions; k++)
    {
      PSO_Y(pso, i)[k] = PSO_X(pso, i)[k]; // fix "aliasing issue"
      // unroll
      // get pointer just one time, and inline the macro
    }

    double x_eval = pso->f(PSO_X(pso, i));
    pso->y_eval[i] = x_eval;
    PSO_FX(pso, i) = x_eval;
  }

  // find ŷ
  double *y_hat = PSO_Y(pso, 0);
  double y_hat_eval = pso->y_eval[0];
  for (int i = 1; i < pso->population_size;
       i++) // unroll with mult. (4?) accumulators
  {
    if (pso->y_eval[i] < y_hat_eval)
    {
      y_hat = PSO_Y(pso, i);
      y_hat_eval = pso->y_eval[i];
    }
  }
  pso->y_hat = y_hat;
  pso->y_hat_eval = y_hat_eval;
}

/*
 * 1. Replace all Pointer accesses (combat aliasing)
 * 2. Unroll "copying y(0) = x(0)" innermost loop w/ Factor 4.
 * 3. Unroll "finding the minimum" loop w/ Factor 4.
 */
void step4_opt1(struct pso_data_constant_inertia *pso)
{
  int pop_size = pso->population_size;
  int dim = pso->dimensions;
  double *pso_x = pso->x;
  double *pso_y = pso->y;
  double *pso_y_eval = pso->y_eval;
  blackbox_fun pso_f = pso->f;

  for (int i = 0; i < pop_size; i++)
  {
    int k = 0;
    int i_dim = i * dim;
    double *pso_x_i_dim = pso_x + i_dim; // PSO_X(pso, 0, i)
    double *pso_y_i_dim = pso_y + i_dim;

    for (; k < dim - 3; k += 4)
    {
      int k1 = k + 1, k2 = k + 2, k3 = k + 3;
      // serpate loads and stores
      // vec?
      // => memcpy
      double t0 = pso_x_i_dim[k];  // t = PSO_X(pso, 0, i)[k];
      double t1 = pso_x_i_dim[k1]; // t = PSO_X(pso, 0, i)[k];
      double t2 = pso_x_i_dim[k2]; // t = PSO_X(pso, 0, i)[k];
      double t3 = pso_x_i_dim[k3]; // t = PSO_X(pso, 0, i)[k];

      pso_y_i_dim[k] = t0;  // PSO_Y(pso, i)[k]
      pso_y_i_dim[k1] = t1; // PSO_Y(pso, i)[k]
      pso_y_i_dim[k2] = t2; // PSO_Y(pso, i)[k]
      pso_y_i_dim[k3] = t3; // PSO_Y(pso, i)[k]
    }
    for (; k < dim; k++)
    {
      pso_y_i_dim[k] = pso_x_i_dim[k];
    }

    double x_eval = pso_f(pso_x_i_dim);
    pso_y_eval[i] = x_eval; // pso->y_eval[i] = x_eval;
    PSO_FX(pso, i) = x_eval;
  }

  // find ŷ
  double *min0 = pso_y;
  double *min1 = pso_y + dim;
  double *min2 = pso_y + 2 * dim;
  double *min3 = pso_y + 3 * dim;
  double min0_eval = pso_y_eval[0];
  double min1_eval = pso_y_eval[1];
  double min2_eval = pso_y_eval[2];
  double min3_eval = pso_y_eval[3];

  int i0 = 4;
  for (; i0 < pop_size - 3; i0 += 4) // unroll with mult. (4?) accumulators
  {
    int i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3;

    if (pso_y_eval[i0] < min0_eval)
    {
      min0 = pso_y + dim * i0;    // PSO_Y(pso, i0);
      min0_eval = pso_y_eval[i0]; // pso->y_eval[i0];
    }
    if (pso_y_eval[i1] < min1_eval)
    {
      min1 = pso_y + dim * i1;    // PSO_Y(pso, i0+1);
      min1_eval = pso_y_eval[i1]; // pso->y_eval[i0+1];
    }
    if (pso_y_eval[i2] < min2_eval)
    {
      min2 = pso_y + dim * i2;    // PSO_Y(pso, i0+2);
      min2_eval = pso_y_eval[i2]; // pso->y_eval[i0+2];
    }
    if (pso_y_eval[i3] < min3_eval)
    {
      min3 = pso_y + dim * i3;    // PSO_Y(pso, i0+3);
      min3_eval = pso_y_eval[i3]; // pso->y_eval[i0+3];
    }
  }

  for (; i0 < pop_size; i0++) // unroll with mult. (4?) accumulators
  {
    if (pso_y_eval[i0] < min0_eval)
    {
      min0 = pso_y + dim * i0;
      min0_eval = pso_y_eval[i0];
    }
  }

  double *min1_2, *min3_4;
  double min1_2_eval, min3_4_eval;

  if (min0_eval < min1_eval)
  {
    min1_2 = min0;
    min1_2_eval = min0_eval;
  }
  else
  {
    min1_2 = min1;
    min1_2_eval = min1_eval;
  }
  if (min2_eval < min3_eval)
  {
    min3_4 = min2;
    min3_4_eval = min2_eval;
  }
  else
  {
    min3_4 = min3;
    min3_4_eval = min3_eval;
  }

  if (min1_2_eval < min3_4_eval)
  {
    pso->y_hat = min1_2;
    pso->y_hat_eval = min1_2_eval;
  }
  else
  {
    pso->y_hat = min3_4;
    pso->y_hat_eval = min3_4_eval;
  }
}

/*
 * Use memcpy
 */
void step4_opt2(struct pso_data_constant_inertia *pso)
{
  int pop_size = pso->population_size;
  int dim = pso->dimensions;
  double *pso_x = pso->x;
  double *pso_y = pso->y;
  double *pso_y_eval = pso->y_eval;
  blackbox_fun pso_f = pso->f;

  for (int i = 0; i < pop_size; i++)
  {
    int k = 0;
    int i_dim = i * dim;
    double *pso_x_i_dim = pso_x + i_dim; // PSO_X(pso, 0, i)
    double *pso_y_i_dim = pso_y + i_dim;

    memcpy(PSO_Y(pso, i), PSO_X(pso, i), dim * sizeof(double));

    double x_eval = pso_f(pso_x_i_dim);
    pso_y_eval[i] = x_eval; // pso->y_eval[i] = x_eval;
    PSO_FX(pso, i) = x_eval;
  }

  // find ŷ
  double *min0 = pso_y;
  double *min1 = pso_y + dim;
  double *min2 = pso_y + 2 * dim;
  double *min3 = pso_y + 3 * dim;
  double min0_eval = pso_y_eval[0];
  double min1_eval = pso_y_eval[1];
  double min2_eval = pso_y_eval[2];
  double min3_eval = pso_y_eval[3];

  int i0 = 4;
  for (; i0 < pop_size - 3; i0 += 4) // unroll with mult. (4?) accumulators
  {
    int i1 = i0 + 1, i2 = i0 + 2, i3 = i0 + 3;

    if (pso_y_eval[i0] < min0_eval)
    {
      min0 = pso_y + dim * i0;    // PSO_Y(pso, i0);
      min0_eval = pso_y_eval[i0]; // pso->y_eval[i0];
    }
    if (pso_y_eval[i1] < min1_eval)
    {
      min1 = pso_y + dim * i1;    // PSO_Y(pso, i0+1);
      min1_eval = pso_y_eval[i1]; // pso->y_eval[i0+1];
    }
    if (pso_y_eval[i2] < min2_eval)
    {
      min2 = pso_y + dim * i2;    // PSO_Y(pso, i0+2);
      min2_eval = pso_y_eval[i2]; // pso->y_eval[i0+2];
    }
    if (pso_y_eval[i3] < min3_eval)
    {
      min3 = pso_y + dim * i3;    // PSO_Y(pso, i0+3);
      min3_eval = pso_y_eval[i3]; // pso->y_eval[i0+3];
    }
  }

  for (; i0 < pop_size; i0++) // unroll with mult. (4?) accumulators
  {
    if (pso_y_eval[i0] < min0_eval)
    {
      min0 = pso_y + dim * i0;
      min0_eval = pso_y_eval[i0];
    }
  }

  double *min1_2, *min3_4;
  double min1_2_eval, min3_4_eval;

  if (min0_eval < min1_eval)
  {
    min1_2 = min0;
    min1_2_eval = min0_eval;
  }
  else
  {
    min1_2 = min1;
    min1_2_eval = min1_eval;
  }
  if (min2_eval < min3_eval)
  {
    min3_4 = min2;
    min3_4_eval = min2_eval;
  }
  else
  {
    min3_4 = min3;
    min3_4_eval = min3_eval;
  }

  if (min1_2_eval < min3_4_eval)
  {
    pso->y_hat = min1_2;
    pso->y_hat_eval = min1_2_eval;
  }
  else
  {
    pso->y_hat = min3_4;
    pso->y_hat_eval = min3_4_eval;
  }
}

void step4_optimized(struct pso_data_constant_inertia *pso)
{
  //  step4_base(pso);
  step4_opt2(pso);
}
