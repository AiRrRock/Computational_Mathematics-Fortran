CC=gfortran-4.9
CFLAGS=-Wall -std=f2008ts -fimplicit-none -Wno-maybe-uninitialized -static-libgfortran -flto

all:
	$(CC) $(CFLAGS) -c src/SOLVE.F -J obj/ -o obj/sol.o
	$(CC) $(CFLAGS) -c src/DECOMP.F -J obj/ -o obj/dec.o
	$(CC) $(CFLAGS) -c src/main.f95 -I obj/ -o obj/test.o
	$(CC) $(CFLAGS) -o bin/lab.darn ./obj/dec.o obj/sol.o obj/test.o
clean:
	rm -f bin/*.*
	rm -f obj/*.*
	rm -f src/*.swp
	rm -f src/*.swm
	rm -f src/*.swn
run:
	cd ./bin; ./lab.darn
	
