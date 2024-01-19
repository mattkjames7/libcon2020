#include "savetestdata.h"

void saveBesselCephes() {

    int i;
    std::vector<double> x(1001), z0(1001), z1(1001);

    for (i=0;i<1001;i++) {
        x[i] = 0.1*((double) i);
    }

    j0(1001,x.data(),z0.data());
    j1(1001,x.data(),z1.data());

    std::filesystem::path filePath = "cephes-bessel.bin";
    std::ofstream file(filePath,std::ios::binary);

    saveVector(file,x);
    saveVector(file,z0);
    saveVector(file,z1);

    file.close();

}

void saveBesselMaass() {

    int i;
    std::vector<double> x(1001), z0(1001), z1(1001);

    for (i=0;i<1001;i++) {
        x[i] = 0.1*((double) i);
    }

    j0m(1001,x.data(),z0.data());
    j1m(1001,x.data(),z1.data());

    std::filesystem::path filePath = "maass-bessel.bin";
    std::ofstream file(filePath,std::ios::binary);

    saveVector(file,x);
    saveVector(file,z0);
    saveVector(file,z1);

    file.close();

}


int main() {

    saveBesselCephes();
    saveBesselMaass();

    return 0;
}