mkdir ..\build
g++ -c -lm -fPIC -std=c++17 -Wextra con2020.cc -o ..\build\con2020.o 
g++ -c -lm -fPIC -std=c++17 -Wextra bessel.cc -o ..\build\bessel.o 
g++ -c -lm -fPIC -std=c++17 -Wextra trap.cc -o ..\build\trap.o 
g++ -c -lm -fPIC -std=c++17 -Wextra polyeval.cc -o ..\build\polyeval.o 
g++ -c -lm -fPIC -std=c++17 -Wextra libcon2020.cc -o ..\build\libcon2020.o
g++ -c -lm -fPIC -std=c++17 -Wextra smoothd.cc -o ..\build\smoothd.o

