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
	LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib' -L ../lib -lcon2020 -Lgoogletest/build/lib -lgtest -lgtest_main -lm
else ifeq ($(OS),Darwin)
	LDFLAGS=-L ../lib -lcon2020 -Lgoogletest/build/lib -lgtest -lgtest_main -lm -Wl,-rpath,@loader_path/../lib
	CFLAGS+=-fvisibility=default
else
	LDFLAGS=-L ../lib -lcon2020 -Lgoogletest/build/lib -lgtest -lgtest_main -lm
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
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_GMOCK=OFF -DINSTALL_GTEST=OFF && \
	cmake --build build -j$(nproc)

save:
	$(CC) save_test_data.cc -o save_test_data $(LDFLAGS) -Ijson/include -lz
	./save_test_data


all: json gtest
	$(CC) unit/model_output.cc -o model_output_tests $(LDFLAGS) -Ijson/include -Igoogletest/googletest/include -lz
	./model_output_tests

	$(CC) unit/bessel.cc -o bessel_tests $(LDFLAGS) -Ijson/include -Igoogletest/googletest/include -lz
	./bessel_tests


clean:
	-rm -v bessel_tests model_output_tests save_test_data
	-rm -rf json googletest