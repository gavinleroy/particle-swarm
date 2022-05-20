#include <papi.h>
#include <stdio.h>

#include "handle_error.h"

int main()
{
  int retval;

  /* Do some computation here */
  double a = 1.1;
  double b = 1.2;
  double c = 0;

  retval = PAPI_hl_region_begin("computation");
  if ( retval != PAPI_OK )
    handle_error(1);

  for (int i; i <= 1000000; i++) {
    c += a*b;
  }

  retval = PAPI_hl_region_end("computation");
  if ( retval != PAPI_OK )
    handle_error(1);

  printf("c: %f\n", c);
}