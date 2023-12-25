#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <fstream>
#include <string>
#include <complex>
#include <iostream>
#include <numbers>

#include "./solovay_kitaev.hpp"
#include "util/dvlab_string.hpp"
#include "util/phase.hpp"
#include "solovay_kitaev.hpp"

namespace qsyn::sk_decomp {

Matrix operator*(Matrix const& lhs, Matrix const& rhs) {
    Matrix result(lhs.size(), std::vector<Complex>(rhs[0].size(), Complex(0.0, 0.0)));
    for (size_t i = 0; i < lhs.size(); ++i) {
        for (size_t j = 0; j < rhs[0].size(); ++j) {
            for (size_t k = 0; k < lhs[0].size(); ++k) {
                result[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }
    return result;
}

double distance(Matrix const& lhs, Matrix const& rhs) {
    double result = 0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        for (size_t j = 0; j < lhs[0].size(); ++j) {
            result += std::abs(lhs[i][j] - rhs[i][j]);
        }
    }
    return result;
}

Matrix adjoint(Matrix const& matrix) {
    Matrix result(matrix[0].size(), std::vector<Complex>(matrix.size(), Complex(0.0, 0.0)));
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            result[j][i] = std::conj(matrix[i][j]);
        }
    }
    return result;
}

Matrix diagonalize(Matrix const& matrix) {
    Complex a = matrix[0][0], b = matrix[0][1], c = matrix[1][0], d = matrix[1][1];
    Complex trace = a + d;
    Complex det = a * d - b * c;
    
    Complex lambda1 = (trace + std::sqrt(trace * trace - Complex(4.0) * det)) / Complex(2.0);
    Complex lambda2 = (trace - std::sqrt(trace * trace - Complex(4.0) * det)) / Complex(2.0);

    Matrix result(matrix.size(), std::vector<Complex>(matrix[0].size(), Complex(0.0, 0.0)));
    result[0][0] = lambda1;
    result[1][1] = lambda2;
    return result;
}

std::pair<Vector3, double> u_to_bloch(Matrix const& matrix) {
    Complex a = matrix[0][0], b = matrix[0][1], c = matrix[1][0], d = matrix[1][1];
    double angle = std::acos((a + d) / 2.0).real();
    double sin = std::sin(angle);
    if (sin < 1e-10) {
        return {{0.0, 0.0, 1.0}, 2.0 * angle};
    } else {
        double nx = ((b + c) / Complex(2.0 * sin)).real();
        double ny = ((b - c) / Complex(2.0 * sin)).real();
        double nz = ((a - d) / Complex(2.0 * sin)).real();
        return {{nx, ny, nz}, 2.0 * angle};
    }
}

std::pair<Matrix, Matrix> group_comm_decomp(Matrix const& matrix) {
    Vector3 axis;
    double angle;
    std::tie(axis, angle) = u_to_bloch(matrix);
    double phi = 2.0 * std::asin(std::sqrt(std::sqrt(0.5 - 0.5 * std::cos(angle / 2.0))));
    Matrix v = {
        {Complex(std::cos(phi / 2.0), 0.0), Complex(0.0, -std::sin(phi / 2.0))},
        {Complex(0.0, -std::sin(phi / 2.0)), Complex(std::cos(phi / 2.0), 0.0)}
    };
    Matrix w;
    if (axis[2] > 0) {
        w = {
            {Complex(std::cos(2.0 * std::numbers::pi - phi / 2.0), 0.0), Complex(0.0, -std::sin(2.0 * std::numbers::pi - phi / 2.0))},
            {Complex(0.0, -std::sin(2.0 * std::numbers::pi - phi / 2.0)), Complex(std::cos(2.0 * std::numbers::pi - phi / 2.0), 0.0)}
        };
    } else {
        w = {
            {Complex(std::cos(phi / 2.0), 0.0), Complex(0.0, -std::sin(phi / 2.0))},
            {Complex(0.0, -std::sin(phi / 2.0)), Complex(std::cos(phi / 2.0), 0.0)}
        };
    }

    Matrix ud = diagonalize(matrix);
    Matrix vwvdwd = diagonalize(v * w * adjoint(v) * adjoint(w));
    Matrix s = ud * adjoint(vwvdwd);

    Matrix v_hat = s * v * adjoint(s);
    Matrix w_hat = s * w * adjoint(s);
    return {v_hat, w_hat};
}


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

void SKD::set_basis(std::vector<std::string> const& basis) {
    for (std::string b : basis) {
        if (b == "h") {
            _basis_gates["h"] = {
                {Complex(1.0 / std::sqrt(2)), Complex(1.0 / std::sqrt(2))},
                {Complex(1.0 / std::sqrt(2)), Complex(-1.0 / std::sqrt(2))}
            };
        } else if (b == "t") {
            _basis_gates["t"] = {
                {Complex(1.0), Complex(0.0)},
                {Complex(0.0), std::exp(Complex(0, std::numbers::pi / 4))}
            };
        } else if (b == "s") {
            _basis_gates["s"] = {
                {Complex(1.0), Complex(0.0)},
                {Complex(0.0), Complex(0.0, 1.0)}
            };
        } else {
            spdlog::error("Unknown basis gate \"{}\"!!", b);
        }
    }
    spdlog::info("The basic approximations are:");
    for (auto const& [name, matrix] : _basis_gates) {
        spdlog::info("{}: {} + {}i, {} + {}i", name, matrix[0][0].real(), matrix[0][0].imag(), matrix[0][1].real(), matrix[0][1].imag());
        spdlog::info("{}: {} + {}i, {} + {}i", name, matrix[1][0].real(), matrix[1][0].imag(), matrix[1][1].real(), matrix[1][1].imag());
    }
}

void SKD::create_basic_approximations(int depth) {
    _basis_approximations.clear();
    std::unordered_map<std::string, Matrix> current = _basis_gates;
    for (int i = 0; i < depth; ++i) {
        for (auto const& [name, matrix] : current) {
            _basis_approximations[name] = matrix;
        }
        std::unordered_map<std::string, Matrix> next;
        for (auto const& [name, matrix] : current) {
            for (auto const& [name2, matrix2] : _basis_gates) {
                next[name + name2] = matrix * matrix2;
            }
        }
        current = next;
    }

    spdlog::info("The basic approximations are:");
    for (auto const& [name, matrix] : _basis_approximations) {
        spdlog::info("{}: {} + {}i, {} + {}i", name, matrix[0][0].real(), matrix[0][0].imag(), matrix[0][1].real(), matrix[0][1].imag());
        spdlog::info("{}: {} + {}i, {} + {}i", name, matrix[1][0].real(), matrix[1][0].imag(), matrix[1][1].real(), matrix[1][1].imag());
    }
}

std::string SKD::find_closest_approximation(Matrix const& matrix) {
    double min_distance = 1e9;
    std::string min_name;
    for (auto const& [name, approx] : _basis_approximations) {
        double d = distance(matrix, approx);
        if (d < min_distance) {
            min_distance = d;
            min_name = name;
        }
    }
    
    return min_name;
}


Matrix SKD::sk_decomp(Matrix const& u, size_t depth) {
    // """Solovay-Kitaev Algorithm."""
    // if n == 0:
    //     return find_closest_u(gates, u)
    // else:
    //     u_next = sk_algo(u, gates, n - 1)
    //     v, w = gc_decomp(u @ u_next.adjoint())
    //     v_next = sk_algo(v, gates, n - 1)
    //     w_next = sk_algo(w, gates, n - 1)
    //     return v_next @ w_next @ v_next.adjoint() @ w_next.adjoint() @ u_next
    if (depth == 0) {
        std::string approx_name = find_closest_approximation(u);
        return _basis_approximations[approx_name];
    }
    Matrix u_next = sk_decomp(u, depth - 1);
    Matrix v, w;
    std::tie(v, w) = group_comm_decomp(u * adjoint(u_next));
    Matrix v_next = sk_decomp(v, depth - 1);
    Matrix w_next = sk_decomp(w, depth - 1);
    return v_next * w_next * adjoint(v_next) * adjoint(w_next) * u_next;
}

void SKD::run() {
    set_basis({"h", "t", "s"});
    create_basic_approximations(_depth);

    Matrix approx = sk_decomp(_input_matrix, _depth);

    double d = distance(_input_matrix, approx);

    spdlog::info("The approximation is:");
    for (auto const& row : approx) {
        spdlog::info("{} + {}i, {} + {}i", row[0].real(), row[0].imag(), row[1].real(), row[1].imag());
    }
    spdlog::info("The distance is: {}", d);
}

}  // namespace qsyn::sk_decomp
