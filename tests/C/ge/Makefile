
# Debug flags
CFLAGS+=-O0 -ggdb3 \
-Wall -Wextra -Wpedantic -Wformat=2 -Wswitch-default -Wswitch-enum -Wfloat-equal \
-pedantic-errors -Werror=format-security \
-Werror=vla \
-I../../../opus/src

# Release flags
#CFLAGS += -O2 -flto -march=native

CPPFLAGS+=-DTEST_PERF

LDLIBS+=-lpso -L../../../opus

CFILES := src/main.c
OBJFILES := $(CFILES:.c=.o)

# Optionnal sanitizers
CFLAGS += -fsanitize=undefined -fsanitize=address
LDFLAGS += -fsanitize=undefined -fsanitize=address

test: $(OBJFILES)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


.PHONY: clean
clean:
	rm $(OBJFILES) ||:
	rm test ||: