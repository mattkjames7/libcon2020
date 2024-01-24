#include "tests.h"



posFieldVectors readTestVectorFile(std::filesystem::path filePath) {

    std::ifstream file(filePath,std::ios::binary);

    std::vector<double> p0 = readVector(file);
    std::vector<double> p1 = readVector(file);
    std::vector<double> p2 = readVector(file);

    std::vector<double> b0 = readVector(file);
    std::vector<double> b1 = readVector(file);
    std::vector<double> b2 = readVector(file);

    tupleVectors P = {p0,p1,p2};
    tupleVectors B = {b0,b1,b2};

    return {P,B};
}

tupleVectors getModelFields(
    tupleVectors pos, 
    std::string equationType, 
    bool cartesian
    ) {

    Con2020 model;

    model.SetCartIn(cartesian);
    model.SetCartOut(cartesian);

    model.SetEqType(equationType.c_str());

    std::vector<double> p0 = std::get<0>(pos);
    std::vector<double> p1 = std::get<1>(pos);
    std::vector<double> p2 = std::get<2>(pos);

    int n = p0.size();

    std::vector<double> b0(n);
    std::vector<double> b1(n);
    std::vector<double> b2(n);

    for (int i=0;i<n;i++) {
        con2020.Field(p0[i],p1[i],p2[i],&b0[i],&b1[i],&b2[i]);
    }

    return {b0,b1,b2};

}

void testModel(
    std::filesystem::path filePath,
    std::string equationType,
    bool cartesian 
) {

    int ndot;
    std::ostringstream out;
    out << "Testing model type: \"" << equationType << "\" (";
    if (cartesian) {
        out << "Cartesian";
    } else {
        out << "Polar";
    }
    out << ")";
    ndot = 46 - out.str().size();

    for (int i=0;i<ndot;i++) {
        out << ".";
    }
    std::cout << out.str();

    posFieldVectors posField = readTestVectorFile(filePath);
    tupleVectors pos = std::get<0>(posField);
    tupleVectors field0 = std::get<1>(posField);

    tupleVectors field1 = getModelFields(pos,equationType,cartesian);

    bool result = compareVectors(field0,field1);


    
    if (result) {
        std::cout << "PASS";
    } else {
        std::cout << "FAIL";
    }
    std::cout << std::endl;

}

void compareModelOutputs() {
    testModel("con2020-analytic-rtp.bin","analytic",false);
    testModel("con2020-analytic-xyz.bin","analytic",true);
    testModel("con2020-hybrid-rtp.bin","hybrid",false);
    testModel("con2020-hybrid-xyz.bin","hybrid",true);
    testModel("con2020-integral-rtp.bin","integral",false);
    testModel("con2020-integral-xyz.bin","integral",true);
}


int main() {
    compareModelOutputs();

    return 0;
}