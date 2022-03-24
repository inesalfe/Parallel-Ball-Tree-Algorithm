CC=mpicc
CFLAGS=-O3 

ballAlg:  ballAlg.c
	$(CC) $(CFLAGS) ballAlg.c -o ballAlg -lm -fopenmp

clean:
	rm -f *.o core a.out proj *~