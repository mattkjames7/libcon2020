# CFLAGS for CC
CFLAGS=-std=c++17 -Wextra -O3

#includeflags
IFLAGS=-I../include

# Compiler
CCo=g++ -fPIC -c $(CFLAGS) $(IFLAGS)
CC=g++ -fPIC $(CFLAGS) $(IFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS) $(IFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS) $(IFLAGS)

OS=$(shell uname -s)
ifeq ($(OS),Linux) 
	LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib' -L ../lib -lcon2020 -lgtest -lgtest_main -lm
else 
	LDFLAGS=-L ../lib -lcon2020 -lgtest -lgtest_main -lm
endif


json:
	git clone https://github.com/nlohmann/json.git
	cd json && \
	git checkout v3.12.0

all: json
	$(CC) bessel.cc -o bessel $(LDFLAGS)
	./bessel

	$(CC) testdata.cc tests.cc -o tests $(LDFLAGS)
	./tests

	$(CC) save_test_data.cc -o save_test_data $(LDFLAGS) -Ijson/include -lz
	./save_test_data
	$(CC) testbessel.cc -o testbessel	$(LDFLAGS) 
	gcc $(IFLAGS) -fPIC -lm -O3 testc.c -o testc $(LDFLAGS) 
	$(CC) testc.cc -o testcc $(LDFLAGS) 
	./testc
	./testcc
	./testbessel
	python3 testbessel.py
	python3 lmic.py
	python3 testsmooth.py


save:
	$(CC) testdata.cc savetestdata.cc -o savetestdata $(LDFLAGS) 
	./savetestdata

clean:
	-rm -v testbessel
