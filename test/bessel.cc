#include "con2020.h"
#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <vector>

bool compareVectors(
    const std::vector<double> &v0,
    const std::vector<double> &v1,
    double epsilon = 1e-12
) {
    if (v0.size() != v1.size()) {
        return false;
    }
    size_t i;
    for (i=0;i<v0.size();i++) {
        if (std::fabs(v0[i] - v1[i]) > epsilon) {
            return false;
        }
    }
    return true;
}

TEST_CASE("Bessel function j0(x) test cases", "[j0(x)]") {

    std::vector<double> x = {0.0,1.0,5.0,10.0,50.0,100.0};
    std::vector<double> expected = {
        1.0,
        0.7651976865579665,
        -0.1775967713143383,
        -0.24593576445134832,
        0.055812327669252086,
        0.01998585030422333
    };
    double result;

    size_t i;
    for (i=0;i<x.size();i++) {
        result = j0(x[i]);
        REQUIRE(result == ApproxEq(expected[i]).epsilon(1e-12));
    }

}

TEST_CASE("Bessel function j0(n,x,j) test cases", "[j0(n,x,j)]") {

    std::vector<double> x = {0.0,1.0,5.0,10.0,50.0,100.0};
    std::vector<double> expected = {
        1.0,
        0.7651976865579665,
        -0.1775967713143383,
        -0.24593576445134832,
        0.055812327669252086,
        0.01998585030422333
    };
    std::vector<double> result(x.size());

    /* this bit will be simpler once I use vectors all around */
    size_t n = x.size();
    double xa[n], j[n];
    size_t i;
    for (i=0;i<n;i++) {
        xa[i] = x[i]; 
    }
    j0(n,xa,j);
    for (i=0;i<n;i++) {
        result[i] = j[i]; 
    }

    for (i=0;i<x.size();i++) {
        REQUIRE(compareVectors(result,expected));
    }

}


TEST_CASE("Bessel function j1(x) test cases", "[j1(x)]") {

    std::vector<double> x = {0.0,1.0,5.0,10.0,50.0,100.0};
    std::vector<double> expected = {
        0.0,
        0.44005058574493355,
        -0.3275791375914653,
        0.04347274616886141,
        -0.09751182812517509,
        -0.0771453520141123,
    };
    double result;

    size_t i;
    for (i=0;i<x.size();i++) {
        result = j1(x[i]);
        REQUIRE(result == ApproxEq(expected[i]).epsilon(1e-12));
    }

}

TEST_CASE("Bessel function j1(n,x,j) test cases", "[j1(n,x,j)]") {

    std::vector<double> x = {0.0,1.0,5.0,10.0,50.0,100.0};
    std::vector<double> expected = {
        0.0,
        0.44005058574493355,
        -0.3275791375914653,
        0.04347274616886141,
        -0.09751182812517509,
        -0.0771453520141123,
    };
    std::vector<double> result(x.size());

    /* this bit will be simpler once I use vectors all around */
    size_t n = x.size();
    double xa[n], j[n];
    size_t i;
    for (i=0;i<n;i++) {
        xa[i] = x[i]; 
    }
    j1(n,xa,j);
    for (i=0;i<n;i++) {
        result[i] = j[i]; 
    }

    for (i=0;i<x.size();i++) {
        REQUIRE(compareVectors(result,expected));
    }

}