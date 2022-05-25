#pragma once

#include <stdint.h>
#include <stdio.h>

void print_matrixd(double const *M, int N, char const *name);
void print_rect_matrixd(double const *A, int n, int m, char const *name);

void print_vectord(double const *M, int N, char const *name);
void print_vectori(int const *M, int N, char const *name);
void print_vectoru64(uint64_t const *M, int N, char const *name);

double rand_between(double a, double b);

double dist2(size_t dim, double const *x, double const *y);
double dist(size_t dim, double const *x, double const *y);
