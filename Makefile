CC=mpicc
CFLAGS=-g 

# ballAlg: gen_points.o ballAlg.c
# 	$(CC) $(CFLAGS) ballAlg.c -o ballAlg gen_points.o -lm -fopenmp

ballAlg: gen_points.o ballAlg.c
	$(CC) $(CFLAGS) ballAlg.c -o ballAlg gen_points.o -lm

gen_points.o: gen_points.c gen_points.h
	$(CC) -c $(CFLAGS) gen_points.c

clean:
	rm -f *.o core a.out proj *~