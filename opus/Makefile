CFLAGS+=\
-O0 -ggdb3 \
-march=native \
-Wall -Wextra -Wpedantic -Wformat=2 -Wswitch-default -Wswitch-enum -Wfloat-equal \
-pedantic-errors -Werror=format-security \
-Werror=vla \
-Wp,-D_GLIBCXX_ASSERTIONS \
-fexceptions \
-flto -fPIE

# LDFLAGS+=-Wl,-z,noexecstack,-z,relro,-z,defs,-z,now,-pie
LDLIBS+=-lm

# optionnal sanitizers
CFLAGS += -fsanitize=undefined -fsanitize=address
LDFLAGS += -fsanitize=undefined -fsanitize=address


CFILES := src/helpers.c src/local_refinement.c src/logging.c src/main.c src/plu_factorization.c src/pso.c src/timing.c src/timing_params.c
OBJFILES := $(CFILES:.c=.o)

pgm: $(OBJFILES)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

.PHONY: clean
clean:
	rm $(OBJFILES)
	rm pgm
