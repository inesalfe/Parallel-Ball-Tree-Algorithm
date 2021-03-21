CFLAGS=-O3 -Wall
CPP = /usr/local/opt/llvm/bin/clang
CPPFLAGS = -I/usr/local/opt/llvm/include -fopenmp
LDFLAGS = -L/usr/local/opt/llvm/lib

ballAlg: gen_points.o ballAlg.c
	$(CPP) $(CPPFLAGS) $(CFLAGS) ballAlg.c -o ballAlg gen_points.o -lm $(LDFLAGS)

gen_points.o: gen_points.c gen_points.h
	$(CPP) -c $(CPPFLAGS) $(CFLAGS) gen_points.c

clean:
	rm -f *.o core a.out proj *~
