#include "surrogate_eval.h"

#include "../helpers.h"

#define QUOTE(x) #x
#define STR(x) QUOTE(x)

double surrogate_eval(struct pso_data_constant_inertia const *pso,
                      double const *x)
{
  return surrogate_eval_0(pso, x);
}

double surrogate_eval_0(struct pso_data_constant_inertia const *pso,
                        double const *x)
{
  // TODO: add the past_local_refinements
  // Note: will also require adding them to the bloom filter / distinctiveness
  // check

  // TODO: I think there is something wrong with the surrogate:
  //  I don't have the same values when evaluating in the Jupyter notebook
  //  To investigate!
  //  VO

  double res = 0;

  double *lambda_p = pso->lambda_p;
  // lambda_p is the concatenation (lambda_0 ... lambda_i || p_0 ... p_(d+1))
  double *lambda = lambda_p;
  double *p_coef = lambda_p + pso->x_distinct_s;

  // iterate directly on x_distinct
  for (size_t k = 0; k < pso->x_distinct_s; k++)
  {
    double *u = pso->x_distinct + k * pso->dimensions;
    double d = dist(pso->dimensions, u, x);
    res += lambda[k] * d * d * d;
  }

  for (int j = 0; j < pso->dimensions; j++)
  {
    res += p_coef[j + 1] * x[j];
  }
  res += p_coef[0];

  return res;
}
