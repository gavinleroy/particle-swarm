extern "C"
{
#include "perf_fit_surrogate.h"
#include "../steps/fit_surrogate.h"
};

#include <iostream>

#include "cstring"

#include "PerformanceTester.hpp"

template <> class ArgumentRestorer<fit_surrogate_fun_t>
{
public:
  ArgumentRestorer<fit_surrogate_fun_t>(struct pso_data_constant_inertia *pso)
  {
  }

  void restore_arguments(struct pso_data_constant_inertia *pso)
  {
    pso->x_distinct_idx_of_last_batch = pso->x_distinct_s - 1;
  }

  ~ArgumentRestorer<fit_surrogate_fun_t>() {}
};

static PerformanceTester<fit_surrogate_fun_t> perf_tester;

extern "C" void add_function_FIT_SURROGATE(fit_surrogate_fun_t f,
                                           char const *name, int flop)
{
  std::cout << "Register " << name << "\n";
  perf_tester.add_function(f, name, flop);
}

extern "C" int perf_test_fit_surrogate(struct pso_data_constant_inertia *pso)
{
  register_functions_FIT_SURROGATE();
  return perf_tester.perf_test_all_registered(pso);
}

// NOTE I bet we can put these in the templated perf framework and remove the
// weirdness of compiling seperately and linking.
extern "C" void register_functions_FIT_SURROGATE()
{
  add_function_FIT_SURROGATE(&fit_surrogate_0, "Fit surrogate Base", 1);

  add_function_FIT_SURROGATE(&fit_surrogate_1, "Fit surrogate without pso->",
                             1);
  add_function_FIT_SURROGATE(&fit_surrogate_2,
                             "Fit surrogate distances precomputed", 1);
  add_function_FIT_SURROGATE(&fit_surrogate_3, "Fit surrogate memcpy phi", 1);

  add_function_FIT_SURROGATE(&fit_surrogate_4, "Fit surrogate early exit", 1);
}
