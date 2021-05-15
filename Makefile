CC=mpicc
CFLAGS=-g 

ballAlg:  ballAlg.c
	$(CC) $(CFLAGS) ballAlg.c -o ballAlg -lm

clean:
	rm -f *.o core a.out proj *~