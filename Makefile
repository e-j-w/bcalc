CFLAGS   = -O -o2 -g -Wall -fPIC -lm

all: bcalc

bcalc: bcalc.c bcalc.h
	gcc bcalc.c $(CFLAGS) -o bcalc
clean:
	rm -rf *~ *.o bcalc *tmpdatafile*
