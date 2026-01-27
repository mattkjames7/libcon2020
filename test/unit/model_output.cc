#include <vector>
#include <tuple>
#include <string>
#include <zlib.h>
#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include "con2020.h"


typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> Vec3D;


Vec3D get_model_vectors(
    const std::string& eqtype,
    const bool is_cartesian,
    std::vector<double>& pos0,  // x or r
    std::vector<double>& pos1,  // y or theta
    std::vector<double>& pos2   // z or phi
) {
    std::size_t N = pos0.size();
    std::vector<double> B0(N);
    std::vector<double> B1(N); 
    std::vector<double> B2(N);

    con2020.SetEqType(eqtype.c_str());
    con2020.SetCartIn(is_cartesian);
    con2020.SetCartOut(is_cartesian);
    con2020.Field(
        N,
        pos0.data(), pos1.data(), pos2.data(),
        B0.data(), B1.data(), B2.data()
    );

    return std::make_tuple(B0, B1, B2);
}


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


// Helper function to convert array format ["key", value, "key", value, ...] to object
nlohmann::json array_to_object(const nlohmann::json& arr) {
    nlohmann::json obj;
    for (size_t i = 0; i < arr.size(); i += 2) {
        std::string key = arr[i].get<std::string>();
        obj[key] = arr[i + 1];
    }
    return obj;
}


TEST(ModelOutputTests, TestCartesianFieldAnalytic) {
    const nlohmann::json test_data = read_json_gz("test_data.json.gz");

    const std::string eqtype = "analytic";
    auto cartesian = array_to_object(test_data[eqtype]["cartesian"]);
    auto x_vals = cartesian["x"].get<std::vector<double>>();
    auto y_vals = cartesian["y"].get<std::vector<double>>();
    auto z_vals = cartesian["z"].get<std::vector<double>>();
    auto Bx_expected = cartesian["Bx"].get<std::vector<double>>();
    auto By_expected = cartesian["By"].get<std::vector<double>>();
    auto Bz_expected = cartesian["Bz"].get<std::vector<double>>();
    const size_t N = x_vals.size();

    std::vector<double> Bx_computed(N);
    std::vector<double> By_computed(N);
    std::vector<double> Bz_computed(N);

    std::tie(Bx_computed, By_computed, Bz_computed) = get_model_vectors(
        eqtype,
        true,
        x_vals,
        y_vals,
        z_vals
    );

    for (size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(Bx_computed[i], Bx_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(By_computed[i], By_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Bz_computed[i], Bz_expected[i], 1e-6) << "at index " << i;
    }
}


TEST(ModelOutputTests, TestSphericalFieldAnalytic) {
    const nlohmann::json test_data = read_json_gz("test_data.json.gz");

    const std::string eqtype = "analytic";
    auto spherical = array_to_object(test_data[eqtype]["spherical"]);
    auto r_vals = spherical["r"].get<std::vector<double>>();
    auto theta_vals = spherical["theta"].get<std::vector<double>>();
    auto phi_vals = spherical["phi"].get<std::vector<double>>();
    auto Br_expected = spherical["Br"].get<std::vector<double>>();
    auto Btheta_expected = spherical["Btheta"].get<std::vector<double>>();
    auto Bphi_expected = spherical["Bphi"].get<std::vector<double>>();
    const size_t N = r_vals.size();

    std::vector<double> Br_computed(N);
    std::vector<double> Btheta_computed(N);
    std::vector<double> Bphi_computed(N);

    std::tie(Br_computed, Btheta_computed, Bphi_computed) = get_model_vectors(
        eqtype,
        false,
        r_vals,
        theta_vals,
        phi_vals
    );

    for (size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(Br_computed[i], Br_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Btheta_computed[i], Btheta_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Bphi_computed[i], Bphi_expected[i], 1e-6) << "at index " << i;
    }
}


TEST(ModelOutputTests, TestCartesianFieldHybrid) {
    const nlohmann::json test_data = read_json_gz("test_data.json.gz");

    const std::string eqtype = "hybrid";
    auto cartesian = array_to_object(test_data[eqtype]["cartesian"]);
    auto x_vals = cartesian["x"].get<std::vector<double>>();
    auto y_vals = cartesian["y"].get<std::vector<double>>();
    auto z_vals = cartesian["z"].get<std::vector<double>>();
    auto Bx_expected = cartesian["Bx"].get<std::vector<double>>();
    auto By_expected = cartesian["By"].get<std::vector<double>>();
    auto Bz_expected = cartesian["Bz"].get<std::vector<double>>();
    const size_t N = x_vals.size();

    std::vector<double> Bx_computed(N);
    std::vector<double> By_computed(N);
    std::vector<double> Bz_computed(N);

    std::tie(Bx_computed, By_computed, Bz_computed) = get_model_vectors(
        eqtype,
        true,
        x_vals,
        y_vals,
        z_vals
    );

    for (size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(Bx_computed[i], Bx_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(By_computed[i], By_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Bz_computed[i], Bz_expected[i], 1e-6) << "at index " << i;
    }
}


TEST(ModelOutputTests, TestSphericalFieldHybrid) {
    const nlohmann::json test_data = read_json_gz("test_data.json.gz");

    const std::string eqtype = "hybrid";
    auto spherical = array_to_object(test_data[eqtype]["spherical"]);
    auto r_vals = spherical["r"].get<std::vector<double>>();
    auto theta_vals = spherical["theta"].get<std::vector<double>>();
    auto phi_vals = spherical["phi"].get<std::vector<double>>();
    auto Br_expected = spherical["Br"].get<std::vector<double>>();
    auto Btheta_expected = spherical["Btheta"].get<std::vector<double>>();
    auto Bphi_expected = spherical["Bphi"].get<std::vector<double>>();
    const size_t N = r_vals.size();

    std::vector<double> Br_computed(N);
    std::vector<double> Btheta_computed(N);
    std::vector<double> Bphi_computed(N);

    std::tie(Br_computed, Btheta_computed, Bphi_computed) = get_model_vectors(
        eqtype,
        false,
        r_vals,
        theta_vals,
        phi_vals
    );

    for (size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(Br_computed[i], Br_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Btheta_computed[i], Btheta_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Bphi_computed[i], Bphi_expected[i], 1e-6) << "at index " << i;
    }
}


TEST(ModelOutputTests, TestCartesianFieldIntegral) {
    const nlohmann::json test_data = read_json_gz("test_data.json.gz");

    const std::string eqtype = "integral";
    auto cartesian = array_to_object(test_data[eqtype]["cartesian"]);
    auto x_vals = cartesian["x"].get<std::vector<double>>();
    auto y_vals = cartesian["y"].get<std::vector<double>>();
    auto z_vals = cartesian["z"].get<std::vector<double>>();
    auto Bx_expected = cartesian["Bx"].get<std::vector<double>>();
    auto By_expected = cartesian["By"].get<std::vector<double>>();
    auto Bz_expected = cartesian["Bz"].get<std::vector<double>>();
    const size_t N = x_vals.size();

    std::vector<double> Bx_computed(N);
    std::vector<double> By_computed(N);
    std::vector<double> Bz_computed(N);

    std::tie(Bx_computed, By_computed, Bz_computed) = get_model_vectors(
        eqtype,
        true,
        x_vals,
        y_vals,
        z_vals
    );

    for (size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(Bx_computed[i], Bx_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(By_computed[i], By_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Bz_computed[i], Bz_expected[i], 1e-6) << "at index " << i;
    }
}


TEST(ModelOutputTests, TestSphericalFieldIntegral) {
    const nlohmann::json test_data = read_json_gz("test_data.json.gz");

    const std::string eqtype = "integral";
    auto spherical = array_to_object(test_data[eqtype]["spherical"]);
    auto r_vals = spherical["r"].get<std::vector<double>>();
    auto theta_vals = spherical["theta"].get<std::vector<double>>();
    auto phi_vals = spherical["phi"].get<std::vector<double>>();
    auto Br_expected = spherical["Br"].get<std::vector<double>>();
    auto Btheta_expected = spherical["Btheta"].get<std::vector<double>>();
    auto Bphi_expected = spherical["Bphi"].get<std::vector<double>>();
    const size_t N = r_vals.size();

    std::vector<double> Br_computed(N);
    std::vector<double> Btheta_computed(N);
    std::vector<double> Bphi_computed(N);

    std::tie(Br_computed, Btheta_computed, Bphi_computed) = get_model_vectors(
        eqtype,
        false,
        r_vals,
        theta_vals,
        phi_vals
    );

    for (size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(Br_computed[i], Br_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Btheta_computed[i], Btheta_expected[i], 1e-6) << "at index " << i;
        EXPECT_NEAR(Bphi_computed[i], Bphi_expected[i], 1e-6) << "at index " << i;
    }
}
