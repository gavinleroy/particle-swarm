#pragma once

/*
 * This file defines a Bloom filter that can be used to test the proximity of
 * points in our problem space.
 * It works by binning the points in a grid of edge size epsilon. The
 * coordinates of the bin is a vector of integers that can be tested for
 * inclusion in a bloom filter.
 *
 *      (bin_id)_k := floor(v_k / eps)
 *
 * Using this simple approach would yield errors for points close together but
 * separated by an edge, e.g. in 1D with epsilon=1.0 for 2.95 and 3.05. They
 * would end up in separate bins.
 *
 * Instead, we can look for the presence of points in all the neighboring bins,
 * e.g. for point 2.95, look at bin 1, 2 and 3.
 *
 * As an optimization, we can restrict ourselves to checking only 2^d bins by
 * looking at the position of the point in the dual grid (i.e. the grid shifted
 * by eps)
 *
 *      (dual_bin_id)_k := floor( v_k / eps - 1. / 2)
 *
 * In 1D:
 *
 *    Grid:
 *
 *   -1       0       1       2       3       4
 *    |-------|-------|-------|-------|-------|
 *
 *    Dual grid
 *       -1       0       1       2       3       4
 *    ----|-------|-------|-------|-------|-------|
 *
 *                                   ^
 *                                   |
 *                                   ---- 2.95
 *
 *     For x=(2.95), bin_id=2 and dual_bin_id=3
 *
 * The bloom filter is checked for the presence of points in (dual_bin_id-1) or
 * dual_bin_id-1, if the check comes back negative, bin_id is added to the bloom
 * filter. (A true collision implies d(x,y) <= sqrt(d * (3/2 esp)^2) ) Note:
 * this implies 2^d bloom filter checks ! To decrease the complexity, it's
 * possible to only check only at bin_id (but then close points separated by a
 * grid edge will not be found).
 *
 * To make the bin_id space larger, the points are translated by -lower_bound
 * (lower limit of the bounding box) to be saved as uint64_t.
 */

#include "bloom.h"

struct rounding_bloom
{
  struct bloom bloom;
  double *lower_bound;
  double epsilon;
  int dims;
};

int rounding_bloom_init(struct rounding_bloom *bloom, int entries, double error,
                        double epsilon, int dims, double const *lower_bound);
int rounding_bloom_check_add(struct rounding_bloom *bloom, int dims,
                             double *const x, int add);

void rounding_bloom_print(struct rounding_bloom *bloom);

void rounding_bloom_free(struct rounding_bloom *bloom);
