#include <con2020.h>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <zlib.h>
#include <fstream>

nlohmann::json read_json_gz(const std::string &filename) {
    gzFile gz = gzopen(filename.c_str(), "rb");
    if (!gz) {
        throw std::runtime_error("Failed to open gzip file");
    }

    std::string json_str;
    char buffer[4096];
    int bytes_read;
    while ((bytes_read = gzread(gz, buffer, sizeof(buffer))) > 0) {
        json_str.append(buffer, bytes_read);
    }
    gzclose(gz);

    return nlohmann::json::parse(json_str);
}


TEST(BesselFunctionTests, TestBesselJ0_lt5) {
    auto test_data = read_json_gz("bessel_data.json.gz");

    auto x_vals0 = test_data["x0"].get<std::vector<double>>();
    auto j0_expected0 = test_data["j0_0"].get<std::vector<double>>();

    std::vector<double> j0_computed0(x_vals0.size());
    con2020::bessel::j0(x_vals0.size(), x_vals0.data(), j0_computed0.data());
    for (size_t i = 0; i < x_vals0.size(); ++i) {
        EXPECT_NEAR(j0_computed0[i], j0_expected0[i], 1e-10) << "at x = " << x_vals0[i];
    }
}


TEST(BesselFunctionTests, TestBesselJ0_gt5) {
    auto test_data = read_json_gz("bessel_data.json.gz");

    auto x_vals1 = test_data["x1"].get<std::vector<double>>();
    auto j0_expected1 = test_data["j0_1"].get<std::vector<double>>();

    std::vector<double> j0_computed1(x_vals1.size());
    con2020::bessel::j0(x_vals1.size(), x_vals1.data(), j0_computed1.data());
    for (size_t i = 0; i < x_vals1.size(); ++i) {
        EXPECT_NEAR(j0_computed1[i], j0_expected1[i], 1e-10) << "at x = " << x_vals1[i];
    }
}


TEST(BesselFunctionTests, TestBesselJ1_lt5) {
    auto test_data = read_json_gz("bessel_data.json.gz");

    auto x_vals0 = test_data["x0"].get<std::vector<double>>();
    auto j1_expected0 = test_data["j1_0"].get<std::vector<double>>();
    std::vector<double> j1_computed0(x_vals0.size());
    con2020::bessel::j1(x_vals0.size(), x_vals0.data(), j1_computed0.data());
    for (size_t i = 0; i < x_vals0.size(); ++i) {
        EXPECT_NEAR(j1_computed0[i], j1_expected0[i], 1e-10) << "at x = " << x_vals0[i];
    }
}


TEST(BesselFunctionTests, TestBesselJ1_gt5) {
    auto test_data = read_json_gz("bessel_data.json.gz");

    auto x_vals1 = test_data["x1"].get<std::vector<double>>();
    auto j1_expected1 = test_data["j1_1"].get<std::vector<double>>();

    std::vector<double> j1_computed1(x_vals1.size());
    con2020::bessel::j1(x_vals1.size(), x_vals1.data(), j1_computed1.data());
    for (size_t i = 0; i < x_vals1.size(); ++i) {
        EXPECT_NEAR(j1_computed1[i], j1_expected1[i], 1e-10) << "at x = " << x_vals1[i];
    }
}
