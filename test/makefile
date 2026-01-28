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
	LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib' -L ../lib -lcon2020 -Lgoogletest/install/lib -lgtest -lgtest_main -lm
else 
	LDFLAGS=-L ../lib -lcon2020 -Lgoogletest/install/lib -lgtest -lgtest_main -lm
endif


json:
	git clone https://github.com/nlohmann/json.git
	cd json && \
	git checkout v3.12.0

gtest:
	git clone https://github.com/google/googletest.git googletest
	cd googletest && \
	git checkout v1.17.0 && \
	mkdir -p install && \
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="install" -DBUILD_GMOCK=OFF -DINSTALL_GTEST=OFF && \
	cmake --build build -j$(nproc) && \
	cmake --install build

save:
	$(CC) save_test_data.cc -o save_test_data $(LDFLAGS) -Ijson/include -Igoogletest/install/include -lz
	./save_test_data


all: json gtest
	$(CC) unit/model_output.cc -o model_output_tests $(LDFLAGS) -Ijson/include -Igoogletest/install/include -lz
	./model_output_tests

	$(CC) unit/bessel.cc -o bessel_tests $(LDFLAGS) -Ijson/include -Igoogletest/install/include -lz
	./bessel_tests

	python3 testbessel.py
	python3 lmic.py
	python3 testsmooth.py


clean:
	-rm -v bessel_tests model_output_tests save_test_data
	-rm -rf json googletest