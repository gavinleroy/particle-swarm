#include "rounding_bloom.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "helpers.h"

#define DEBUG_BLOOM 0

// For FFI calls
struct rounding_bloom * rounding_bloom_alloc()
{
    return malloc(sizeof(struct rounding_bloom));
}

int rounding_bloom_init(struct rounding_bloom * bloom, int entries, double error, double epsilon, int dims, double const * lower_bound)
{
    int ret;
    if( (ret = bloom_init(&bloom->bloom, entries, error) > 0) )
    {
        return ret;
    }

    bloom->dims = dims;
    bloom->lower_bound = malloc(dims * sizeof(double));
    for (int i = 0 ; i < dims ; i++)
    {
        // Keep some space for the dual grid
        bloom->lower_bound[i] = lower_bound[i] - epsilon;
    }

    bloom->epsilon = epsilon;
    return 0;
}



int rounding_bloom_check_add(struct rounding_bloom * bloom, int dims, double * const x, int add)
{
    int ret = 0;
    uint64_t * bin_id = malloc(dims * sizeof(uint64_t));
    uint64_t * dual_bin_id = malloc(dims * sizeof(uint64_t));


    for (int i = 0 ; i < dims ; i++)
    {
        assert(x[i] >= bloom->lower_bound[i]);
        double y = x[i] - bloom->lower_bound[i];
        double k = y / bloom->epsilon;
        assert(k < (double)UINT64_MAX);
        bin_id[i] = floor(k);
        dual_bin_id[i] = floor(k - 0.5);
    }

    uint64_t max_off = 1 << dims;


    uint64_t * neighbor_dual_bin_id = malloc(dims * sizeof(uint64_t));

    #if DEBUG_BLOOM
    print_vectord(x, dims, "x");
    print_vectoru64(bin_id, dims, "bin_id");
    print_vectoru64(dual_bin_id, dims, "dual_bin_id");
    #endif

    for(uint64_t off = 0 ; off < max_off ; off++)
    {
        for(int i = 0 ; i < dims ; i++)
        {
            neighbor_dual_bin_id[i] = dual_bin_id[i] + ((off >> i) & 1);
        }

        if (bloom_check(&bloom->bloom, (char*)neighbor_dual_bin_id, dims * sizeof(uint64_t)))
        {
            ret=1;
            #if DEBUG_BLOOM
            printf("Found collision with\n");
            print_vectoru64(neighbor_dual_bin_id, dims, "neighbor_dual_bin_id");
            #endif
            break;
        }
    }

    if (add && ret == 0)
    {
        #if DEBUG_BLOOM
        printf("No collision. Adding to bloom filter.\n");
        #endif

        bloom_add(&bloom->bloom, (char*)bin_id, sizeof(uint64_t) * dims);
    }

    free(neighbor_dual_bin_id);
    free(dual_bin_id);
    free(bin_id);
    return ret;
}


int rounding_bloom_check(struct rounding_bloom * bloom, int dims, double * const x)
{
    return rounding_bloom_check_add(bloom, dims, x, 0);
}




void rounding_bloom_print(struct rounding_bloom * bloom)
{
  bloom_print(&bloom->bloom);
  printf(" ->eps = %f\n", bloom->epsilon);
}


void rounding_bloom_free(struct rounding_bloom * bloom)
{
  free(bloom->lower_bound);
  bloom_free(&bloom->bloom);
}
