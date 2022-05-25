# Particle swarm fast implementation

Article: [Particle swarm with radial basis function surrogates for expensive black-box optimization](https://acl.inf.ethz.ch/teaching/fastcode/2022/project/project-ideas/particle-swarm.pdf)

## Code Formatting 

Running `git config --local core.hooksPath .githooks` will enable the githooks which use Clang Format to cleanup modified files. Please opt-in.

# Correctness / Performance

In order to test the code, for both performance and correctness, we have a system set in place to modularly test each component. To build the necessary shared library do the following:

```bash
$ cd opus
$ DEBUG=1 PSO_SHARED=1 make
```

This will make a shared library `pso.so` (`libpso.dylib` on MacOS) which you should put somewhere on your `LD_LIBRARY` (`DYLD_LIBRARY_PATH` on MacOS) such that Julia can find it. 
In order to run the suite, do the following:

```bash
$ cd /tests
$ julia
julia> ]
(tests) pkg> activate .
<backspace> # this means hit the backspace/delete key to leave Pkg mode
julia> include("src/tests.jl") # this will run the enabled tests, Gavin plans on improving this interface later
```

Occassionally you might see some strange behavior from the testing suite. This usually occurs if you make changes to the shared object and re-run the suite. If this is the case, leave the Julia REPL and reload things. This *reboot* tends to fix things. 



# Great ideas

## fit surrogate: (VO)
1. Cache distance calculations for phi
2. Precompute P ?
3. Blocking + scalar replacement in copy to Ab ; to lambda ; to p
4. Vectorize dist calculations


## pso_constant_inertia_loop (XS)
  0. Move steps to individual functions
  1. blocking, scalar repl, vectorize...
  2. Move mallocs to preallocation phase

## RNG

1. precomputed random number generation!


## Use x_distinct as the main place to store positions

- currently all positions stored in x
  + but really only accessed through x_distinct
  + so non-contiguous accesses !
- keep only current positions in x
- keep f(x) only for x_distinct
  + 

## Q° for Tomaso:

- can we have the dimension as a compile time contstant?
  * rationale: the code is going to run for days, 


# TODO to follow the paper

- the space filling design
- use the local refinements in fit_surrogate
