call compileobj.bat

mkdir ..\lib\libcon2020
g++ -lm -fPIC -std=c++17 -Wextra ..\build\*.o -shared -o ..\lib\libcon2020\libcon2020.dll
