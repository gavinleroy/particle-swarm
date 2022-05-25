#pragma once

#include <inttypes.h>

#ifndef LOG_TIMING
#define LOG_TIMING 0
#endif

void set_logging_directory(char const *logdir);
void stop_logging(void);

void log_timing(char const *step, int version, uint64_t time, uint64_t cycles);
uint64_t timing_step(char const *step, int version, uint64_t time,
                     uint64_t tsc_start);

#if LOG_TIMING
#include "perf_testers/tsc_x86.h"

#define TIMING_INIT() uint64_t log_timing_tsc_ctr = start_tsc()
#define TIMING_STEP(STEP, VERSION, TIME)                                       \
  log_timing_tsc_ctr = timing_step(STEP, VERSION, TIME, log_timing_tsc_ctr)
#else
#define TIMING_INIT()                                                          \
  do                                                                           \
  {                                                                            \
  } while (0)
#define TIMING_STEP(STEP, VERSION, TIME)                                       \
  do                                                                           \
  {                                                                            \
  } while (0)
#endif