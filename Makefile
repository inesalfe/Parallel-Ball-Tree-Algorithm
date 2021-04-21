CC=gcc
CFLAGS=-O3 -Wall

ballAlg: gen_points.o ballAlg.c
	$(CC) $(CFLAGS) ballAlg.c -o ballAlg gen_points.o -lm -fopenmp

ballAlg_1: gen_points.o ballAlg_1.c
	$(CC) $(CFLAGS) ballAlg_1.c -o ballAlg_1 gen_points.o -lm -fopenmp

gen_points.o: gen_points.c gen_points.h
	$(CC) -c $(CFLAGS) gen_points.c

clean:
	rm -f *.o core a.out proj *~