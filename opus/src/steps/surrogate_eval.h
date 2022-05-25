#pragma once

#include "../pso.h"

#ifndef SURROGATE_EVAL_VERSION
#define SURROGATE_EVAL_VERSION 0
#endif

double surrogate_eval(struct pso_data_constant_inertia const *pso,
                      double const *x);

double surrogate_eval_0(struct pso_data_constant_inertia const *pso,
                        double const *x);
