#include <con2020.h>
#include <nlohmann/json.hpp>
#include <array>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <zlib.h>


template <std::size_t N, typename T = double>
std::array<T, N> random_array(T min, T max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(min, max);

    std::array<T, N> arr;
    for (std::size_t i = 0; i < N; ++i) {
        arr[i] = dis(gen);
    }
    return arr;
}


nlohmann::json collect_model_data(std::string eqtype) {
    const std::size_t N = 100;

    Con2020 model;

    // Generate random test data (cartesian)
    auto x_vals = random_array<N>(-60.0, 60.0);
    auto y_vals = random_array<N>(-60.0, 60.0);
    auto z_vals = random_array<N>(-20.0, 20.0);

    // Prepare arrays to hold the magnetic field results
    std::array<double, N> Bx_vals;
    std::array<double, N> By_vals;
    std::array<double, N> Bz_vals;

    // Generate random test data (spherical polar)
    auto r_vals = random_array<N>(0.1, 60.0);
    auto theta_vals = random_array<N>(0.0, M_PI);
    auto phi_vals = random_array<N>(0.0, 2 * M_PI);

    // Arrays to hold spherical polar results
    std::array<double, N> Br_vals;
    std::array<double, N> Btheta_vals;
    std::array<double, N> Bphi_vals;

    // set model type
    model.SetEqType(eqtype.c_str());

    // Compute magnetic field for each set of coordinates
    model.SetCartIn(true);
    model.SetCartOut(true);
    model.Field(
        N,
        x_vals.data(), y_vals.data(), z_vals.data(),
        Bx_vals.data(), By_vals.data(), Bz_vals.data()
    );

    model.SetCartIn(false);
    model.SetCartOut(false);
    model.Field(
        N,
        r_vals.data(), theta_vals.data(), phi_vals.data(),
        Br_vals.data(), Btheta_vals.data(), Bphi_vals.data()
    );

    // Collect data into JSON
    nlohmann::json model_data;
    model_data["cartesian"] = {
        "x", x_vals,
        "y", y_vals,
        "z", z_vals,
        "Bx", Bx_vals,
        "By", By_vals,
        "Bz", Bz_vals
    };
    model_data["spherical"] = {
        "r", r_vals,
        "theta", theta_vals,
        "phi", phi_vals,
        "Br", Br_vals,
        "Btheta", Btheta_vals,
        "Bphi", Bphi_vals
    };
    return model_data;
}


nlohmann::json collect_bessel_data() {
    const std::size_t N = 1000;
    auto x_vals0 = random_array<N>(0.0, 5.0);
    auto x_vals1 = random_array<N>(5.0, 10.0);

    std::array<double, N> j0_vals0;
    std::array<double, N> j0_vals1;
    std::array<double, N> j1_vals0;
    std::array<double, N> j1_vals1;

    j0(N, x_vals0.data(), j0_vals0.data());
    j0(N, x_vals1.data(), j0_vals1.data());
    j1(N, x_vals0.data(), j1_vals0.data());
    j1(N, x_vals1.data(), j1_vals1.data());

    // Collect data into JSON
    nlohmann::json bessel_data;
    bessel_data["x0"] = x_vals0;
    bessel_data["x1"] = x_vals1;
    bessel_data["j0_0"] = j0_vals0;
    bessel_data["j0_1"] = j0_vals1;
    bessel_data["j1_0"] = j1_vals0;
    bessel_data["j1_1"] = j1_vals1;
    return bessel_data;
}


void save_json_gz(const nlohmann::json& data, const std::string& filename) {

    // Serialize JSON to string
    std::string json_str = data.dump(2);

    // Open gzip file
    gzFile gz = gzopen(filename.c_str(), "wb");
    if (!gz) {
        throw std::runtime_error("Failed to open gzip file");
    }

    // Write compressed data
    gzwrite(gz, json_str.data(), json_str.size());
    gzclose(gz);
}


void save_test_data() {
    auto eqtypes = {std::string("hybrid"), std::string("analytic"), std::string("integral")};
    nlohmann::json all_data;
    for (const auto& eqtype : eqtypes) {
        all_data[eqtype] = collect_model_data(eqtype);
    }

    save_json_gz(all_data, "test_data.json.gz");
}


void save_bessel_data() {
    nlohmann::json bessel_data = collect_bessel_data();
    save_json_gz(bessel_data, "bessel_data.json.gz");
}


int main() {
    save_test_data();
    save_bessel_data();
    return 0;
}
