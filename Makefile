CFLAGS   = -o2 -Wall -lm

all: bcalc

bcalc: bcalc.c bcalc.h
	gcc bcalc.c $(CFLAGS) -o bcalc
clean:
	rm -rf *~ *.o bcalc *tmpdatafile*
