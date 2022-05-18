#pragma once

#include <stdint.h>

void print_matrixd(double const * M, int N, char const *name);
void print_rect_matrixd(double const * A, int n, int m, char const *name);

void print_vectord(double const * M, int N, char const * name);
void print_vectori(int const * M, int N, char const * name);
void print_vectoru64(uint64_t const * M, int N, char const * name);