CFLAGS   = -o2 -Wall -pedantic -Wshadow -Wunreachable-code -Wpointer-arith -Wcast-qual -Wcast-align -Wstrict-prototypes -Wmissing-prototypes -Wformat-security -Wstack-protector -Wconversion -std=c99 -lm

all: bcalc

bcalc: bcalc.c bcalc.h
	gcc bcalc.c $(CFLAGS) -o bcalc
clean:
	rm -rf *~ *.o bcalc *tmpdatafile*
