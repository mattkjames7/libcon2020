# CFLAGS for CC
CFLAGS=-std=c++23 -Wextra -O3

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

save:
	$(CC) save_test_data.cc -o save_test_data $(LDFLAGS) -Ijson/include -lz
	./save_test_data


all: json
	$(CC) unit/model_output.cc -o model_output_tests $(LDFLAGS) -Ijson/include -lz
	./model_output_tests

	$(CC) unit/bessel.cc -o bessel_tests $(LDFLAGS) -Ijson/include -lz
	./bessel_tests



	$(CC) bessel.cc -o bessel $(LDFLAGS)
	./bessel

	$(CC) testdata.cc tests.cc -o tests $(LDFLAGS)
	./tests


	$(CC) testbessel.cc -o testbessel	$(LDFLAGS) 
	gcc $(IFLAGS) -fPIC -lm -O3 testc.c -o testc $(LDFLAGS) 
	$(CC) testc.cc -o testcc $(LDFLAGS) 
	./testc
	./testcc
	./testbessel
	python3 testbessel.py
	python3 lmic.py
	python3 testsmooth.py


clean:
	-rm -v testbessel
