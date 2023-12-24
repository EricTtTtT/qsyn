#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <fstream>
#include <string>
#include <complex>
#include <iostream>

#include "./solovay_kitaev.hpp"
#include "util/dvlab_string.hpp"
#include "util/phase.hpp"

namespace qsyn::sk_decomp {

bool SKD::read_skd_file(std::filesystem::path const& filepath) {
    auto const extension = filepath.extension();

    if (extension == ".tex")
        return read_tex(filepath);
    else if (extension == "") {
        return true;
    } else {
        spdlog::error("File format \"{}\" is not supported!!", extension);
        return false;
    }
}

bool SKD::read_tex(std::filesystem::path const& filepath) {
    std::fstream tex_file{filepath};
    if (!tex_file.is_open()) {
        spdlog::error("Cannot open the TEX file \"{}\"!!", filepath);
        return false;
    }
    // set _input_matrix[2][2]
    _input_matrix.resize(2);
    _input_matrix[0].resize(2);
    _input_matrix[1].resize(2);

    size_t row_id = 0;

    while (!tex_file.eof()) {
        std::string line;
        std::getline(tex_file, line);
        if (line.empty()) continue;
        if (line[0] == '%') continue;
        if (line[0] == '\\') {
            if (line.find("\\begin{bmatrix}") != std::string::npos) {
                spdlog::info("Found the beginning of the circuit");
                continue;
            } else if (line.find("\\end{bmatrix}") != std::string::npos) {
                spdlog::info("Found the end of the circuit");
                break;
            } else {
                spdlog::error("Unknown command \"{}\"!!", line);
                return false;
            }
        }
        spdlog::info("Found a line: {}", line);
        std::replace(line.begin(), line.end(), '&', ' ');
        std::replace(line.begin(), line.end(), '\\', ' ');
        
        std::istringstream iss(line);
        double realPart, imagPart;
        char i, plus;
        
        size_t col_id = 0;
        while (iss >> realPart >> plus >> imagPart >> i) {
            std::complex<double> num(realPart, imagPart);
            _input_matrix[row_id][col_id] = num;
            ++col_id;
        }
        ++row_id;
    }

    spdlog::info("The input matrix is:");
    for (auto const& row : _input_matrix) {
        spdlog::info("{} + {}i, {} + {}i", row[0].real(), row[0].imag(), row[1].real(), row[1].imag());
    }

    return true;
}

void SKD::set_basic_approximations() {
    // TODO: check if the basic approximations are correct
    _basic_approximations["h"] = {0, 0.7071, -0.7071, 0, -0.7071, -0.7071, -1, 0, 0};
    _basic_approximations["t"] = {1, 0, 0, 0, 0.7071, 0.7071, 0, 0, 1};
    _basic_approximations["s"] = {1, 0, 0, 0, 0, 1, 0, 1, 0};
}

void SKD::run() {

}

}  // namespace qsyn::sk_decomp
