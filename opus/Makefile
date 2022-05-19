######
# Makefile options
######

# Exemple: build a library in debug mode
#	$ make PSO_SHARED=1 DEBUG=1

# Build using the debug configuration? Default=yes
DEBUG ?= 1

# Build the pso solver as a library? Default=no
PSO_SHARED ?= 0

######
# Compilation flags
######

# Debug flags: no optimisations, debug symbols, extra warnings
# march=native to allow vectorisation
DEBUG_FLAGS:= \
-march=native \
-O0 -g \
-Wall -Wextra -Wpedantic -Wformat=2 -Wswitch-default -Wswitch-enum -Wfloat-equal \
-pedantic-errors -Werror=format-security \
-Werror=vla


# Release flags: optimisations, NDEBUG (disable assertions), disable hardening, disable PIE (only for executable)
RELEASE_FLAGS:=-DNDEBUG -O2 -flto -march=native -fcf-protection=none -fno-stack-clash-protection -fno-stack-protector

# Flags common to debug and release
# Currently starts empty
# COMMON_FLAGS := 

# Change these variables if setting the C/C++ standard is necessary
# CFLAGS += -std=gnu17
# CXXFLAGS += -std=gnu++17

# Libraries to link
LDLIBS+=-lm


######
# C/C++ source files
######

# Object files required for the library and the executable
# (indifferently C or C++, `make` will use the correct rule based on the
# source file extension)
OBJ_COMMON := src/helpers.o src/local_refinement.o src/logging.o \
                 src/plu_factorization.o src/pso.o src/bloom.o src/murmurhash.o \
		 src/rounding_bloom.o src/gaussian_elimination_solver.o \
		 src/timing.o src/timing_params.o

# Object files required for the the library
OBJ_LIB := src/perf_ge_solve.o src/pso_ffi.o
# Object files required for the the executable
OBJ_EXE := src/main.o


######
# Library vs executable compilation options
######

# Detect the system
UNAME_S := $(shell uname -s)


ifeq ($(PSO_SHARED), 1)
	# Option 1. Build a shared library

	# For the shared file we also need to pass
	# a preprocessor definition
	# TEST_PERF=1
	CPPFLAGS += -DTEST_PERF

	ifeq ($(UNAME_S),Darwin)
	target := libpso.dylib
	COMMON_FLAGS += -dynamiclib
	LDFLAGS += -dynamiclib
	else
	target := libpso.so
	COMMON_FLAGS += -fPIC
	LDFLAGS += -shared
	endif

	OBJFILES := $(OBJ_COMMON) $(OBJ_LIB)

else # PSO_SHARED
	# Option 2. Build an executable
	
	# Disable PIE
	COMMON_FLAGS+=-fno-PIE
	LDFLAGS+=-no-pie
	
	# Optionnal sanitizers.
	# Only enable int the executable as loading sanitized libraries
	# from a non-sanitized executable results in an error
	DEBUG_FLAGS+=-fsanitize=undefined -fsanitize=address
	
	target := pso
	OBJFILES := $(OBJ_COMMON) $(OBJ_EXE)

endif


######
# Compilation flags
######


ifeq ($(DEBUG), 0)
	ALL_FLAGS := $(COMMON_FLAGS) $(RELEASE_FLAGS)
else
	ALL_FLAGS := $(COMMON_FLAGS) $(DEBUG_FLAGS)
endif

CXXFLAGS+= $(ALL_FLAGS)
CFLAGS+= $(ALL_FLAGS)
# LTO and sanitizers: the compilation flags are needed at link time
LDFLAGS+= $(ALL_FLAGS)


######
# Build rules
######

# Note that for the compilation we rely on the implicit rules
# https://www.gnu.org/software/make/manual/html_node/Implicit-Rules.html

# Linking rule
# Note that we might be linking C and C++, use CXX as the linker
$(target): $(OBJFILES)
	$(CXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)



.PHONY: clean
clean:
	rm $(OBJ_COMMON) $(OBJ_EXE) $(OBJ_LIB) ||:
	rm pso libpso.so libpso.dylib ||:
