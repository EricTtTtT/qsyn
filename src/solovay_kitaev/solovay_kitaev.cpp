#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <fstream>
#include <string>
#include <complex>
#include <iostream>
#include <numbers>
#include <queue>

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

Matrix operator-(Matrix const& lhs, Matrix const& rhs) {
    Matrix result(lhs.size(), std::vector<Complex>(lhs[0].size(), Complex(0.0, 0.0)));
    for (size_t i = 0; i < lhs.size(); ++i) {
        for (size_t j = 0; j < lhs[0].size(); ++j) {
            result[i][j] = lhs[i][j] - rhs[i][j];
        }
    }
    return result;
}

Matrix sqrt(Matrix const& matrix) {
    // Calculating the trace and determinant
    Complex trace = matrix[0][0] + matrix[1][1];
    Complex det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    Complex s = std::sqrt(trace * trace / 4.0 - det);
    Complex lambda1 = trace / 2.0 + s;
    Complex lambda2 = trace / 2.0 - s;

    // Compute the square roots of the eigenvalues
    lambda1 = std::sqrt(lambda1);
    lambda2 = std::sqrt(lambda2);

    // If the matrix is diagonal, return the square roots of the diagonal elements
    if (matrix[0][1] == Complex(0, 0) && matrix[1][0] == Complex(0, 0)) {
        return {{lambda1, Complex(0, 0)}, {Complex(0, 0), lambda2}};
    }

    // If the matrix is not diagonal, use a simplification for the square root
    // This part is an approximation and may not be accurate for all matrices
    Complex a = matrix[0][0];
    Complex b = matrix[0][1];
    Complex c = matrix[1][0];
    Complex d = matrix[1][1];

    // Solving the system for the square root matrix elements
    Complex x = std::sqrt((a + lambda1) / 2.0);
    Complex y = b / (2.0 * x);
    Complex z = c / (2.0 * x);
    Complex w = std::sqrt((d + lambda2) / 2.0);

    return {{x, y}, {z, w}};
}

double trace(Matrix const& matrix) {
    double result = 0;
    for (size_t i = 0; i < matrix.size(); ++i) {
        result += matrix[i][i].real();
    }
    return result;
}

double trace_dist(Matrix const& lhs, Matrix const& rhs) {
    // """Compute trace distance between two 2x2 matrices."""
    // return np.real(0.5 * np.trace(np.sqrt((u - v).adjoint() @ (u - v))))
    
    // Matrix diff = lhs - rhs;
    // Matrix adj = adjoint(diff);
    // Matrix product = diff * adj;
    // Matrix sqrt_product = sqrt(product);
    // double result = trace(sqrt_product) / 2.0;
    // return result;

    double distance = 0.0;

    assert(lhs.size() == 2 && rhs.size() == 2);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::complex<double> diff = lhs[i][j] - rhs[i][j];
            distance += std::norm(diff);  // Using std::norm to calculate squared magnitude for complex numbers
        }
    }

    return std::sqrt(distance);
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
        double nx = ((b + c) / Complex(0.0, 2.0 * sin)).real();
        double ny = ((b - c) / Complex(2.0 * sin)).real();
        double nz = ((a - d) / Complex(0.0, 2.0 * sin)).real();
        return {{nx, ny, nz}, 2.0 * angle};
    }
}

bool is_unitary(Matrix const& matrix) {
    // python3: abs(np.linalg.det(U2x2)) - 1 < 1e-4
    Complex x = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
    double det = std::abs(x);
    return det - 1.0 < 1e-4;
}

bool is_single_qubit(Matrix const& matrix) {
    return matrix.size() == 2 && matrix[0].size() == 2;
}

std::pair<Matrix, Matrix> group_comm_decomp(Matrix const& matrix) {
    Vector3 axis;
    double angle;
    std::tie(axis, angle) = u_to_bloch(matrix);
    double phi = 2.0 * std::asin(std::sqrt(std::sqrt(0.5 - 0.5 * std::sqrt(1 - std::sin(angle / 2.0) * std::sin(angle / 2.0)))));
    Matrix v = {
        {Complex(std::cos(phi / 2.0), 0.0), Complex(0.0, -std::sin(phi / 2.0))},
        {Complex(0.0, -std::sin(phi / 2.0)), Complex(std::cos(phi / 2.0), 0.0)}
    };
    Matrix w;
    if (axis[2] > 0) {
        w = {
            {Complex(std::cos((2.0 * kPI - phi) / 2.0), 0.0), Complex(-std::sin((2.0 * kPI - phi) / 2.0), 0.0)},
            {Complex(std::sin((2.0 * kPI - phi) / 2.0), 0.0), Complex(std::cos((2.0 * kPI - phi) / 2.0), 0.0)}
        };
    } else {
        w = {
            {Complex(std::cos(phi / 2.0), 0.0), Complex(-std::sin(phi / 2.0), 0.0)},
            {Complex(std::sin(phi / 2.0), 0.0), Complex(std::cos(phi / 2.0), 0.0)}
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
    if (extension.compare(".tex") == 0) {
        return read_tex(filepath);
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
                continue;
            } else if (line.find("\\end{bmatrix}") != std::string::npos) {
                break;
            } else {
                spdlog::error("Unknown command \"{}\"!!", line);
                return false;
            }
        }
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

    spdlog::debug("The input matrix is:");
    for (auto const& row : _input_matrix) {
        spdlog::debug("{} + {}i, {} + {}i", row[0].real(), row[0].imag(), row[1].real(), row[1].imag());
    }

    return true;
}

void SKD::print_input_matrix() const {
    std::cout << "The input matrix is:" << std::endl;
    for (auto const& row : _input_matrix) {
        std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    }
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
                {Complex(0.0), std::exp(Complex(0, kPI / 4))}
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
    spdlog::debug("The basis gates are:");
    for (auto const& [name, matrix] : _basis_gates) {
        spdlog::debug("{}: {} + {}i, {} + {}i", name, matrix[0][0].real(), matrix[0][0].imag(), matrix[0][1].real(), matrix[0][1].imag());
        spdlog::debug("{}: {} + {}i, {} + {}i", name, matrix[1][0].real(), matrix[1][0].imag(), matrix[1][1].real(), matrix[1][1].imag());
    }
}

void SKD::create_basic_approximations(size_t length) {
    _basis_approximations.clear();
    std::unordered_map<std::string, Matrix> current = _basis_gates;
    for (size_t i = 0; i < length; ++i) {
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

    spdlog::debug("The basic approximations are:");
    for (auto const& [name, matrix] : _basis_approximations) {
        spdlog::debug("{}: {} + {}i, {} + {}i", name, matrix[0][0].real(), matrix[0][0].imag(), matrix[0][1].real(), matrix[0][1].imag());
        spdlog::debug("{}: {} + {}i, {} + {}i", name, matrix[1][0].real(), matrix[1][0].imag(), matrix[1][1].real(), matrix[1][1].imag());
    }
}

std::string SKD::find_closest_approximation(Matrix const& matrix, std::unordered_map<std::string, Matrix> const& approximations, bool print) {
    size_t report_num = 20;
    
    std::priority_queue<std::pair<std::string, double>,
        std::vector<std::pair<std::string, double>>,
        std::function<bool(std::pair<std::string, double> const&, std::pair<std::string, double> const&)>>
        approx_max_heap([](std::pair<std::string, double> const& lhs, std::pair<std::string, double> const& rhs) {
            return lhs.second < rhs.second;
        });

    for (auto const& [name, approx] : approximations) {
        double d = trace_dist(matrix, approx);
        if (approx_max_heap.size() < report_num) {
            approx_max_heap.push({name, d});
        } else if (d < approx_max_heap.top().second) {
            approx_max_heap.pop();
            approx_max_heap.push({name, d});
        }
    }
    
    std::string min_name;
    if (print) {
        std::cout << "Report the closest " << report_num << " approximations:" << std::endl;
    }
    while (!approx_max_heap.empty()) {
        if (print) {
            std::cout << "  " << std::setw(5) << approx_max_heap.top().first << ": distance=" << approx_max_heap.top().second << std::endl;
        }
        min_name = approx_max_heap.top().first;
        approx_max_heap.pop();
    }

    // spdlog::debug("The approximation is: {}", min_name);
    
    return min_name;
}

void SKD::report_basis() {
    std::string approx_name = find_closest_approximation(_input_matrix, _basis_approximations, true);
}

void SKD::report_decomp_result() const {
    double d = trace_dist(_input_matrix, _approx_matrix.matrix);

    // std::cout << "Input matrix:" << std::endl;
    // for (auto const& row : _input_matrix) {
    //     std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    // }
    std::cout << "The approximation is: " << std::endl;
    std::cout << _approx_matrix.name << std::endl;
    for (auto const& row : _approx_matrix.matrix) {
        std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    }
    std::cout << "The distance is: " << d << std::endl;
}

InfoMatrix SKD::sk_decomp(Matrix const& u, size_t depth) {
    InfoMatrix ret;
    if (depth == 0) {
        std::string approx_name = find_closest_approximation(u, _basis_approximations);
        ret.matrix = _basis_approximations[approx_name];
        ret.name = approx_name;
        return ret;
    }
    InfoMatrix u_next = sk_decomp(u, depth - 1);
    Matrix v, w;
    std::tie(v, w) = group_comm_decomp(u * adjoint(u_next.matrix));
    InfoMatrix v_next = sk_decomp(v, depth - 1);
    InfoMatrix w_next = sk_decomp(w, depth - 1);
    ret.matrix = v_next.matrix * w_next.matrix * adjoint(v_next.matrix) * adjoint(w_next.matrix) * u_next.matrix;
    ret.name = " -> (v:" + v_next.name + ", w:" + w_next.name + ", u+: " + u_next.name + ")";
    return ret;
}

void SKD::run() {
    spdlog::info("parameters: depth={}, length={}, param={}", _depth, _length, _param);
    InfoMatrix approx = sk_decomp(_input_matrix, _depth);
    _approx_matrix.matrix = approx.matrix;
    _approx_matrix.name = approx.name;

    double d = trace_dist(_input_matrix, approx.matrix);

    // std::cout << "Input matrix:" << std::endl;
    // for (auto const& row : _input_matrix) {
    //     std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    // }
    std::cout << "The approximation is: " << std::endl;
    for (auto const& row : approx.matrix) {
        std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    }
    std::cout << "The distance is: " << d  << std::endl;
}

}  // namespace qsyn::sk_decomp
