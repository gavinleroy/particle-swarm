#pragma once

#include <sys/types.h>

void set_logging_directory(char const * logdir);
void log_surrogate(char const * fname,
    double * lambda,
    double * p,
    double * x,
    size_t t,
    size_t dim,
    size_t popsize
);