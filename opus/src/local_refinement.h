#pragma once

#include <sys/types.h>

typedef double (*local_optimization_function)(double const *const,
                                              void const *const);

void local_optimization(local_optimization_function f, // R^d -> R
                        size_t dimensions,             // R
                        double const *center,          // R^d
                        double xi,                     // R
                        double const *a,               // R^d
                        double const *b,               // R^d
                        void const *const additionnal_f_args, double *x_min);
