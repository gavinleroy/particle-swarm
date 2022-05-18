#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "pso.h"

#define NR 32
#define CYCLES_REQUIRED 1e8
#define REP 10
#define EPS (1e-3)

// #define DEBUG_TIMING

#include "timing.h"
#include "tsc_x86.h"

double perf_counter(enum profiled_function_t profiled_function, profiling_hyperparameters * hyperparams) {
    #ifdef DEBUG_TIMING
    printf("Starting perf_counter\n");
    #endif
    double cycles = 0.;
    long num_runs = 10;
    double multiplier = 1;
    myInt64 start, end;

    // Parameter building phase
    #ifdef DEBUG_TIMING
    printf("Starting to build parameters\n");
    #endif
    void* params;

    switch (profiled_function)
    {
    case run_pso_f:
        params = malloc(sizeof(run_pso_parameters));
        generate_run_pso_parameters((run_pso_parameters*) params, hyperparams);
        break;

    case plu_factorization_f:
        params = malloc(sizeof(plu_factorize_parameters));
        generate_plu_factorize_parameters((plu_factorize_parameters*) params, hyperparams);
        break;

    case plu_solve_f:
        params = malloc(sizeof(plu_solve_parameters));
        generate_plu_solve_parameters((plu_solve_parameters*) params, hyperparams);
        break;

    case fit_surrogate_f:
        params = malloc(sizeof(struct pso_data_constant_inertia));
        generate_fit_surrogate_parameters((struct pso_data_constant_inertia*) params, hyperparams);
        break;
    
    default:
        break;
    }

#ifdef DEBUG_TIMING
    printf("Starting to warmup\n");
#endif
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        #ifdef DEBUG_TIMING
        printf("New warmup round\n");
        #endif
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            switch (profiled_function)
            {
            case run_pso_f:
                run_pso_wrapper((run_pso_parameters*) params);
                break;

            case plu_factorization_f:
                plu_factorize_wrapper((plu_factorize_parameters*) params);
                break;

            case plu_solve_f:
                plu_solve_wrapper((plu_solve_parameters*) params);
                break;

            case fit_surrogate_f:
                fit_surrogate_wrapper((struct pso_data_constant_inertia*) params);
                break;
            
            default:
                break;
            }
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);

    #ifdef DEBUG_TIMING
    printf("Starting to measure\n");
    #endif
    // Actual performance measurements repeated REP times.
    // The following allows us to find the average number of cycles.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        #ifdef DEBUG_TIMING
        printf("Measurement rep number %d\n", j);
        #endif

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            switch (profiled_function)
            {
            case run_pso_f:
                run_pso_wrapper((run_pso_parameters*) params);
                break;

            case plu_factorization_f:
                plu_factorize_wrapper((plu_factorize_parameters*) params);
                break;

            case plu_solve_f:
                plu_solve_wrapper((plu_solve_parameters*) params);
                break;

            case fit_surrogate_f:
                fit_surrogate_wrapper((struct pso_data_constant_inertia*) params);
                break;
            
            default:
                break;
            }
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

    }

    #ifdef DEBUG_TIMING
    printf("End of perf_counter\n");
    #endif
    total_cycles /= REP;

    cycles = total_cycles;
    return  cycles; //round((100.0 * flops) / cycles) / 100.0;
}