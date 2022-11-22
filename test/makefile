# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra -O3

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

OS=$(shell uname -s)
ifeq ($(OS),Linux) 
	LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib/libcon2020' -L ../lib/libcon2020 -lcon2020
else 
	LDFLAGS=-L ../lib/libcon2020 -lcon2020
endif

all:
	$(CC) testbessel.cc -o testbessel	$(LDFLAGS) 
	gcc -fPIC -lm -O3 testc.c -o testc $(LDFLAGS) 
	$(CC) testc.cc -o testcc $(LDFLAGS) 
	./testc
	./testcc
	./testbessel
	python3 testbessel.py
	python3 lmic.py

clean:
	-rm -v testbessel
