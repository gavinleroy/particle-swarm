#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "pso.h"

#define NR 32
#define CYCLES_REQUIRED 1e8
#define REP 10
#define EPS (1e-3)

#include "timing.h"
#include "tsc_x86.h"

double perf_counter(enum profiled_function_t profiled_function) {
    printf("Starting perf_counter\n");
    double cycles = 0.;
    long num_runs = 10;
    double multiplier = 1;
    myInt64 start, end;

    // Parameter building phase
    printf("Starting to build parameters\n");
    void* params;

    switch (profiled_function)
    {
    case run_pso_f:
        params = malloc(sizeof(run_pso_parameters));
        generate_run_pso_parameters((run_pso_parameters*) params);
        break;
    
    default:
        break;
    }

    printf("Starting to warmup\n");
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        printf("New warmup round\n");
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            switch (profiled_function)
            {
            case run_pso_f:
                run_pso_wrapper((run_pso_parameters*) params);
                break;
            
            default:
                break;
            }
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);

    printf("Starting to measure\n");
    // Actual performance measurements repeated REP times.
    // The following allows us to find the average number of cycles.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        printf("Measurement rep number %d\n", j);

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            switch (profiled_function)
            {
            case run_pso_f:
                run_pso_wrapper((run_pso_parameters*) params);
                break;
            
            default:
                break;
            }
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

    }

    printf("End of perf_counter\n");
    total_cycles /= REP;

    cycles = total_cycles;
    return  cycles; //round((100.0 * flops) / cycles) / 100.0;
}