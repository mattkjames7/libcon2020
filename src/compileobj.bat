mkdir ..\build
g++ -c -lm -fPIC -std=c++17 -Wextra -O3 con2020.cc -o ..\build\con2020.o 
g++ -c -lm -fPIC -std=c++17 -Wextra -O3 bessel.cc -o ..\build\bessel.o 
g++ -c -lm -fPIC -std=c++17 -Wextra -O3 trap.cc -o ..\build\trap.o 
g++ -c -lm -fPIC -std=c++17 -Wextra -O3 polyeval.cc -o ..\build\polyeval.o 
g++ -c -lm -fPIC -std=c++17 -Wextra -O3 libcon2020.cc -o ..\build\libcon2020.o
g++ -c -lm -fPIC -std=c++17 -Wextra -O3 smoothd.cc -o ..\build\smoothd.o

