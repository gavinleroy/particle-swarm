
# Debug flags
CPPFLAGS+=-I../../../opus/src

# Release flags
CFLAGS:=-fno-PIE -DNDEBUG -O2 -flto -march=native -fcf-protection=none -fno-stack-protector
LDFLAGS+=-no-pie $(CFLAGS)


LDLIBS+=-lpso -L../../../opus

CFILES := src/main.c
OBJFILES := $(CFILES:.c=.o)

# Optionnal sanitizers
#CFLAGS += -fsanitize=undefined -fsanitize=address
#LDFLAGS += -fsanitize=undefined -fsanitize=address

test: $(OBJFILES)
	$(CXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


.PHONY: clean
clean:
	rm $(OBJFILES) ||:
	rm test ||: