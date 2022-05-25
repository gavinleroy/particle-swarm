#pragma once

/**
 *      _________   _____________________  ____  ______
 *     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
 *    / /_  / /| | \__ \ / / / /   / / / / / / / __/
 *   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
 *  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
 *
 *  http://www.acl.inf.ethz.ch/teaching/fastcode
 *  How to Write Fast Numerical Code 263-2300 - ETH Zurich
 *  Copyright (C) 2019
 *                   Tyler Smith        (smitht@inf.ethz.ch)
 *                   Alen Stojanov      (astojanov@inf.ethz.ch)
 *                   Gagandeep Singh    (gsingh@inf.ethz.ch)
 *                   Markus Pueschel    (pueschel@inf.ethz.ch)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see http://www.gnu.org/licenses/.
 */

#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <vector>

extern "C"
{
#include "tsc_x86.h"
}

//#define PERF_TESTER_NR 32
#define PERF_TESTER_CYCLES_REQUIRED 1e8
#define PERF_TESTER_REP 50
//#define PERF_TESTER_EPS (1e-3)

// destructuring function pointer template type
// into return type and argument type:
// see https://devblogs.microsoft.com/oldnewthing/20200713-00/?p=103978

// Argument saver/restorer.
// Need to provide a template specialization for your particular
// function under test.
template <typename F> class ArgumentRestorer;

// eg.
/*
template <>
class ArgumentRestorer<my_fun_to_test_ptr_t>
{
private:
  // private storage
  arg1 arg1_0;
  arg2 arg2_0;
public:
  ArgumentRestorer<my_fun_to_test_ptr_t>(arg1, arg2)
  {
    // allocate private storage and save arg1 arg2
  }

  void restore_arguments(arg1, arg2)
  {
    // copy back arg1_0 and arg2_0 into the args
  }

  ~ArgumentRestorer<my_fun_to_test_ptr_t>(arg1, arg2, arg3)
  {
    // free private storage
  }
};
*/

// The PerformanceTester template parametrized by a function
// pointer. To allow explicit type parametrization, needs an
// empty primary template.
template <typename F> class PerformanceTester;

// Partial specialization specialized for function pointers.
template <typename Ret_T, typename... Args_T>
class PerformanceTester<Ret_T (*)(Args_T...)>
{
  // the signature of the functions(s) we want to benchmark
  typedef Ret_T (*fun_T)(Args_T...);

private:
  std::vector<fun_T> userFuncs;
  std::vector<std::string> funcNames;
  std::vector<int> funcFlops;
  int numFuncs = 0;

public:
  /*
   * Registers a user function to be tested by the driver program. Registers a
   * string description of the function as well
   */
  void add_function(fun_T f, std::string name, int flop)
  {
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    funcFlops.push_back(flop);
    numFuncs++;
  };

  /*
   * Checks the given function for validity (XXX does it?) (RE: no, it doesn't).
   * If valid, then computes and reports and returns the number of cycles
   * required per iteration
   */
  double perf_test(fun_T f, std::string desc, int flops, Args_T... args)
  {
    double cycles = 0.;
    long num_runs = 10;
    double multiplier = 1;
    myInt64 start, end;

    // TODO: save initial paramters with arg_saver and arg_restorer
    // if saver == nulll; restorer == null
    ArgumentRestorer<fun_T> arg_restorer(args...);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do
    {
      num_runs = num_runs * multiplier;
      start = start_tsc();
      for (size_t i = 0; i < num_runs; i++)
      {
        f(args...);
      }
      end = stop_tsc(start);

      cycles = (double)end;

      arg_restorer.restore_arguments(args...);

      multiplier = (PERF_TESTER_CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    std::list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < PERF_TESTER_REP; j++)
    {
      for (size_t i = 0; i < num_runs; ++i)
      {
        start = start_tsc();
        f(args...);
        end = stop_tsc(start);
        cycles = ((double)end);
        total_cycles += cycles;
        cyclesList.push_back(cycles); // XXX why have this

        arg_restorer.restore_arguments(args...);
      }
    }
    total_cycles /= (PERF_TESTER_REP * num_runs);
    cycles = total_cycles;

    return cycles;
  }

  int perf_test_all_registered(Args_T... args)
  {
    std::cout << "Starting performance tests.";
    double perf;
    int i;

    if (numFuncs == 0)
    {
      std::cout << std::endl;
      std::cout << "No functions registered - nothing for driver to do"
                << std::endl;
      std::cout << "Register functions by calling register_func(f, name)"
                << std::endl;
      std::cout << "in register_funcs()" << std::endl;
      return -1;
    }

    std::cout << numFuncs << " functions registered." << std::endl;

    for (i = 0; i < numFuncs; i++)
    {
      perf = perf_test(userFuncs[i], funcNames[i], 1, args...);
      std::cout << std::endl << "Running: " << funcNames[i] << std::endl;
      std::cout << perf << " cycles" << std::endl;
    }

    return 0;
  }
};
