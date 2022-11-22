@echo "Testing con2020 installation"
cd test
g++ testc_installed.cc -o testc_installed.exe -lcon2020
.\testc_installed.exe
rm testc_installed.exe
cd ..