# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra -O3

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

all:
	$(CC) testbessel.cc -o testbessel	-Wl,-rpath='$$ORIGIN/../lib/libcon2020' -L ../lib/libcon2020 -lcon2020

	./testbessel
	python3 testbessel.py

clean:
	-rm -v testbessel
