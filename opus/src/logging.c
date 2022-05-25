#include "logging.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <inttypes.h>

#include <stdio.h>
#include <stdlib.h>

#include "perf_testers/tsc_x86.h"

static int logging_dir_fd = -1;
static FILE *timing_fp = NULL;

void set_logging_directory(char const *logdir)
{
  if ((logging_dir_fd = open(logdir, O_RDONLY | O_DIRECTORY)) < 0)
  {
    perror("Error opening logdir");
    exit(1);
  }

  int timing_fd = -1;

  if ((timing_fd = openat(logging_dir_fd, "timing.csv",
                          O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0)
  {
    perror("Error opening file descriptor");
    exit(1);
  }

  if ((timing_fp = fdopen(timing_fd, "w")) == NULL)
  {
    perror("Error opening file pointer");
    exit(1);
  }
}

void log_timing(char const *step, int version, uint64_t time, uint64_t cycles)
{
  if (timing_fp == NULL)
  {
    return;
  }

  fprintf(timing_fp, "%s,%d,%" PRIu64 ",%" PRIu64 "\n", step, version, time,
          cycles);
}

uint64_t timing_step(char const *step, int version, uint64_t time,
                     uint64_t tsc_start)
{
  uint64_t cycles = stop_tsc(tsc_start);
  log_timing(step, version, time, cycles);
  return start_tsc();
}

void stop_logging(void) { fclose(timing_fp); }