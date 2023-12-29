#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <string>
#include <complex>
#include <iostream>
#include <numbers>
#include <queue>
#include <cassert>

#include "./solovay_kitaev.hpp"
#include "util/dvlab_string.hpp"
#include "util/phase.hpp"
#include "solovay_kitaev.hpp"

namespace qsyn::sk_decomp {

std::string GateSequence::to_string() const {
    std::string ret;
    for (auto const& gate : gates) {
        ret += gate + " ";
    }
    return ret;
}

// Matrix operator*(Matrix const& lhs, Matrix const& rhs) {
//     Matrix result(lhs.size(), std::vector<Complex>(rhs[0].size(), Complex(0.0, 0.0)));
//     for (size_t i = 0; i < lhs.size(); ++i) {
//         for (size_t j = 0; j < rhs[0].size(); ++j) {
//             for (size_t k = 0; k < lhs[0].size(); ++k) {
//                 result[i][j] += lhs[i][k] * rhs[k][j];
//             }
//         }
//     }
//     return result;
// }

// Matrix operator*(Complex const& lhs, Matrix const& rhs) {
//     Matrix result(rhs.size(), std::vector<Complex>(rhs[0].size(), Complex(0.0, 0.0)));
//     for (size_t i = 0; i < rhs.size(); ++i) {
//         for (size_t j = 0; j < rhs[0].size(); ++j) {
//             result[i][j] = lhs * rhs[i][j];
//         }
//     }
//     return result;
// }

// Matrix operator-(Matrix const& lhs, Matrix const& rhs) {
//     Matrix result(lhs.size(), std::vector<Complex>(lhs[0].size(), Complex(0.0, 0.0)));
//     for (size_t i = 0; i < lhs.size(); ++i) {
//         for (size_t j = 0; j < lhs[0].size(); ++j) {
//             result[i][j] = lhs[i][j] - rhs[i][j];
//         }
//     }
//     return result;
// }

// Matrix sqrt(Matrix const& matrix) {
//     // Calculating the trace and determinant
//     Complex trace = matrix[0][0] + matrix[1][1];
//     Complex det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
//     Complex s = std::sqrt(trace * trace / 4.0 - det);
//     Complex lambda1 = trace / 2.0 + s;
//     Complex lambda2 = trace / 2.0 - s;

//     // Compute the square roots of the eigenvalues
//     lambda1 = std::sqrt(lambda1);
//     lambda2 = std::sqrt(lambda2);

//     // If the matrix is diagonal, return the square roots of the diagonal elements
//     if (matrix[0][1] == Complex(0, 0) && matrix[1][0] == Complex(0, 0)) {
//         return {{lambda1, Complex(0, 0)}, {Complex(0, 0), lambda2}};
//     }

//     // If the matrix is not diagonal, use a simplification for the square root
//     // This part is an approximation and may not be accurate for all matrices
//     Complex a = matrix[0][0];
//     Complex b = matrix[0][1];
//     Complex c = matrix[1][0];
//     Complex d = matrix[1][1];

//     // Solving the system for the square root matrix elements
//     Complex x = std::sqrt((a + lambda1) / 2.0);
//     Complex y = b / (2.0 * x);
//     Complex z = c / (2.0 * x);
//     Complex w = std::sqrt((d + lambda2) / 2.0);

//     return {{x, y}, {z, w}};
// }

// double trace(Matrix const& matrix) {
//     double result = 0;
//     for (size_t i = 0; i < matrix.size(); ++i) {
//         result += matrix[i][i].real();
//     }
//     return result;
// }

// double trace_dist(Matrix const& lhs, Matrix const& rhs) {
//     double distance = 0.0;

//     for (size_t i = 0; i < lhs.size(); ++i) {
//         for (size_t j = 0; j < lhs[0].size(); ++j) {
//             std::complex<double> diff = lhs[i][j] - rhs[i][j];
//             distance += std::norm(diff);  // Using std::norm to calculate squared magnitude for complex numbers
//         }
//     }

//     return std::sqrt(distance);
// }

double distance(MatrixX const& lhs, MatrixX const& rhs) {
    double distance = 0.0;
    MatrixX diff = lhs - rhs;
    distance = diff.norm();
    distance = std::sqrt(distance);
    return distance;
}

// Matrix adjoint(Matrix const& matrix) {
//     Matrix result(matrix[0].size(), std::vector<Complex>(matrix.size(), Complex(0.0, 0.0)));
//     for (size_t i = 0; i < matrix.size(); ++i) {
//         for (size_t j = 0; j < matrix[0].size(); ++j) {
//             result[j][i] = std::conj(matrix[i][j]);
//         }
//     }
//     return result;
// }

// Matrix diagonalize(Matrix const& matrix) {
//     Complex a = matrix[0][0], b = matrix[0][1], c = matrix[1][0], d = matrix[1][1];
//     Complex trace = a + d;
//     Complex det = a * d - b * c;
    
//     Complex lambda1 = (trace + std::sqrt(trace * trace - Complex(4.0) * det)) / Complex(2.0);
//     Complex lambda2 = (trace - std::sqrt(trace * trace - Complex(4.0) * det)) / Complex(2.0);

//     Matrix result(matrix.size(), std::vector<Complex>(matrix[0].size(), Complex(0.0, 0.0)));
//     result[0][0] = lambda1;
//     result[1][1] = lambda2;
//     return result;
// }

Matrix2 u2_to_su2(Matrix2 const& matrix) {
    Complex det = matrix.determinant();
    return std::sqrt(Complex(1.0) / det) * matrix;
}

Matrix3 su2_to_so3(Matrix2 const& matrix) {
    double a = matrix(0, 0).real();
    double b = matrix(0, 0).imag();
    double c = -matrix(0, 1).real();
    double d = -matrix(0, 1).imag();
    Matrix3 rotation {
        {a * a - b * b - c * c + d * d, 2 * a * b + 2 * c * d, -2 * a * c + 2 * b * d},
        {-2 * a * b + 2 * c * d, a * a - b * b + c * c - d * d, 2 * a * d + 2 * b * c},
        {2 * a * c + 2 * b * d, 2 * b * c - 2 * a * d, a * a + b * b - c * c - d * d}
    };
    return rotation;
}

// std::pair<Vector3, double> u_to_bloch(Matrix const& matrix) {
//     Complex a = matrix[0][0], b = matrix[0][1], c = matrix[1][0], d = matrix[1][1];
//     double angle = std::acos((a + d) / 2.0).real();
//     double sin = std::sin(angle);
//     if (sin < 1e-10) {
//         return {{0.0, 0.0, 1.0}, 2.0 * angle};
//     } else {
//         double nx = ((b + c) / Complex(0.0, 2.0 * sin)).real();
//         double ny = ((b - c) / Complex(2.0 * sin)).real();
//         double nz = ((a - d) / Complex(0.0, 2.0 * sin)).real();
//         return {{nx, ny, nz}, 2.0 * angle};
//     }
// }

bool is_unitary(Matrix2 const& matrix) {
    // Complex x = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
    Complex x = matrix.determinant();
    double det = std::abs(x);
    return det - 1.0 < 1e-4;
}

bool is_single_qubit(Matrix2 const& matrix) {
    return matrix.size() == 4;
}

void group_comm_decomp(Matrix3 const& matrix, GateSequence& v, GateSequence& w) {
    // Vector3 axis;
    // double angle;
    // std::tie(axis, angle) = u_to_bloch(matrix);
    // double phi = 2.0 * std::asin(std::sqrt(std::sqrt(0.5 - 0.5 * std::sqrt(1 - std::sin(angle / 2.0) * std::sin(angle / 2.0)))));
    // Matrix v = {
    //     {Complex(std::cos(phi / 2.0), 0.0), Complex(0.0, -std::sin(phi / 2.0))},
    //     {Complex(0.0, -std::sin(phi / 2.0)), Complex(std::cos(phi / 2.0), 0.0)}
    // };
    // Matrix w;
    // if (axis[2] > 0) {
    //     w = {
    //         {Complex(std::cos((2.0 * kPI - phi) / 2.0), 0.0), Complex(-std::sin((2.0 * kPI - phi) / 2.0), 0.0)},
    //         {Complex(std::sin((2.0 * kPI - phi) / 2.0), 0.0), Complex(std::cos((2.0 * kPI - phi) / 2.0), 0.0)}
    //     };
    // } else {
    //     w = {
    //         {Complex(std::cos(phi / 2.0), 0.0), Complex(-std::sin(phi / 2.0), 0.0)},
    //         {Complex(std::sin(phi / 2.0), 0.0), Complex(std::cos(phi / 2.0), 0.0)}
    //     };
    // }

    // Matrix ud = diagonalize(matrix);
    // Matrix vwvdwd = diagonalize(v * w * adjoint(v) * adjoint(w));
    // Matrix s = ud * adjoint(vwvdwd);
    // Matrix v_hat = s * v * adjoint(s);
    // Matrix w_hat = s * w * adjoint(s);
    // return {v_hat, w_hat};

    // double phi = _solve_phi(matrix);
    // Matrix vx = _rotate_angle(phi, {1.0, 0.0, 0.0});
    // Matrix vy = _rotate_angle(phi, {0.0, 1.0, 0.0});
    // Matrix commutator = _get_commutator(vx, vy);

    // Vector3 matrix_axis = _get_rotation_axis(matrix);
    // Vector3 commutator_axis = _get_rotation_axis(commutator);

    // Matrix sim_matrix = _get_rotation_between(matrix_axis, commutator_axis);
    // Matrix sim_matrix_dagger = adjoint(sim_matrix);

    // Matrix v = sim_matrix * vx * sim_matrix_dagger;
    // Matrix w = sim_matrix * vy * sim_matrix_dagger;

    // return {v, w};
    std::cout << matrix << v.to_string() << w.to_string() << std::endl;
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
            _input_matrix(row_id, col_id) = num;
            ++col_id;
        }
        ++row_id;
    }

    spdlog::debug("The input matrix is:");
    spdlog::debug("  {} + {}i, {} + {}i", _input_matrix(0, 0).real(), _input_matrix(0, 0).imag(), _input_matrix(0, 1).real(), _input_matrix(0, 1).imag());
    spdlog::debug("  {} + {}i, {} + {}i", _input_matrix(1, 0).real(), _input_matrix(1, 0).imag(), _input_matrix(1, 1).real(), _input_matrix(1, 1).imag());
    return true;
}

void SKD::print_input_matrix() const {
    std::cout << "The input matrix is:" << std::endl;
    // for (auto const& row : _input_matrix) {
    //     std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    // }
    std::cout << _input_matrix << std::endl;
}

void SKD::set_basis(std::vector<std::string> const& basis) {
    for (std::string b : basis) {
        if (b == "h") {
            _basis_gates["h"] = Matrix2 {
                {Complex(1.0 / std::sqrt(2)), Complex(1.0 / std::sqrt(2))},
                {Complex(1.0 / std::sqrt(2)), Complex(-1.0 / std::sqrt(2))}
            };
        } else if (b == "t") {
            _basis_gates["t"] = Matrix2 {
                {Complex(1.0), Complex(0.0)},
                {Complex(0.0), std::exp(Complex(0, kPI / 4))}
            };
        } else if (b == "s") {
            _basis_gates["s"] = Matrix2 {
                {Complex(1.0), Complex(0.0)},
                {Complex(0.0), Complex(0.0, 1.0)}
            };
        } else {
            spdlog::error("Unknown basis gate \"{}\"!!", b);
        }
    }
}

void SKD::create_basic_approximations(size_t length) {
    _basis_approximations.clear();
    std::vector<GateSequence> curr_;
    for (auto const& [name, matrix] : _basis_gates) {
        curr_.emplace_back(GateSequence{{name}, u2_to_su2(matrix), su2_to_so3(u2_to_su2(matrix))});
    }

    for (size_t i = 0; i < length; ++i) {
        for (auto const& gs : curr_) {
            _basis_approximations.emplace_back(gs);
        }
        std::vector<GateSequence> next_;
        for (auto const& gs : curr_) {
            for (auto const& [name2, matrix2] : _basis_gates) {
                std::vector<std::string> gates = gs.gates;
                Matrix2 su2 = gs.matrix2 * matrix2;
                gates.emplace_back(name2);
                next_.emplace_back(GateSequence{
                    gates,
                    su2,
                    su2_to_so3(su2) * gs.matrix3
                });
            }
        }
        curr_ = next_;
    }
}

size_t SKD::find_closest_approximation(GateSequence const& gs, bool print) {
    size_t report_num = 20;
    
    std::priority_queue<std::pair<size_t, double>,
        std::vector<std::pair<size_t, double>>,
        std::function<bool(std::pair<size_t, double> const&, std::pair<size_t, double> const&)>>
        approx_max_heap([](std::pair<size_t, double> const& lhs, std::pair<size_t, double> const& rhs) {
            return lhs.second < rhs.second;
        });

    for (size_t i=0; i<_basis_approximations.size(); ++i) {
        double d = distance(gs.matrix3, _basis_approximations[i].matrix3);
        if (approx_max_heap.size() < report_num) {
            approx_max_heap.push({i, d});
        } else if (d < approx_max_heap.top().second) {
            approx_max_heap.pop();
            approx_max_heap.push({i, d});
        }
    }

    std::string min_name;
    if (print) {
        std::cout << "Report the closest " << report_num << " approximations:" << std::endl;
    }
    std::vector<std::pair<size_t, double>> id_dist_sorted;
    while (!approx_max_heap.empty()) {
        size_t gs_idx = approx_max_heap.top().first;
        std::string gs_name = _basis_approximations[gs_idx].to_string();
        id_dist_sorted.emplace_back(std::make_pair(gs_idx, approx_max_heap.top().second));
        approx_max_heap.pop();
    }
    if (print) {
        for (auto it = id_dist_sorted.rbegin(); it != id_dist_sorted.rend(); ++it) {
            size_t gs_idx = it->first;
            std::string gs_name = _basis_approximations[gs_idx].to_string();
            std::cout << "  gate: " << std::setw(20) << gs_name << ", distance=" << it->second << std::endl;
        }
    }

    return id_dist_sorted.back().first;
}

void SKD::report_basis() {
    // std::string approx_name = find_closest_approximation(_input_matrix, _basis_approximations, true);
}

void SKD::report_decomp_result() const {
    // double d = trace_dist(_input_matrix, _approx_matrix.matrix);

    // // std::cout << "Input matrix:" << std::endl;
    // // for (auto const& row : _input_matrix) {
    // //     std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    // // }
    // std::cout << "The approximation is: " << std::endl;
    // std::cout << _approx_matrix.name << std::endl;
    // for (auto const& row : _approx_matrix.matrix) {
    //     std::cout << "  " << row[0].real() << " + " << row[0].imag() << "i, " << row[1].real() << " + " << row[1].imag() << "i" << std::endl;
    // }
    // std::cout << "The distance is: " << d << std::endl;
}

GateSequence SKD::sk_decomp(GateSequence const& u, size_t depth) {
    if (depth == 0) {
        size_t approx_idx = find_closest_approximation(u);
        return _basis_approximations[approx_idx];
    }

    // u_n1 = self._recurse(sequence, n - 1, check_input=check_input)
    // # print(f"n: {n}, ori_distance: {np.linalg.norm(np.subtract(u_n1.product, sequence.product))}")

    // v_n, w_n = commutator_decompose(
    //     sequence.dot(u_n1.adjoint()).product, check_input=check_input
    // )

    // v_n1 = self._recurse(v_n, n - 1, check_input=check_input)
    // w_n1 = self._recurse(w_n, n - 1, check_input=check_input)
    // approx = v_n1.dot(w_n1).dot(v_n1.adjoint()).dot(w_n1.adjoint()).dot(u_n1)

    GateSequence u_next = sk_decomp(u, depth - 1);
    GateSequence v, w;
    group_comm_decomp(u.matrix3, v, w);
    GateSequence v_next = sk_decomp(v, depth - 1);
    GateSequence w_next = sk_decomp(w, depth - 1);
    // ret.matrix = v_next.matrix * w_next.matrix * adjoint(v_next.matrix) * adjoint(w_next.matrix) * u_next.matrix;
    // ret.name = " -> (v:" + v_next.name + ", w:" + w_next.name + ", u+: " + u_next.name + ")";
    GateSequence ret;
    ret.matrix3 = v_next.matrix3 * w_next.matrix3 * v_next.matrix3.adjoint() * w_next.matrix3.adjoint() * u_next.matrix3;
    ret.matrix2 = v_next.matrix2 * w_next.matrix2 * v_next.matrix2.adjoint() * w_next.matrix2.adjoint() * u_next.matrix2;
    ret.gates = {"v:" + v_next.to_string(), "w:" + w_next.to_string(), "u+:" + u_next.to_string()};
    return ret;
}

void SKD::run() {
    spdlog::info("parameters: depth={}, length={}, param={}", _depth, _length, _param);
    GateSequence input_gs;
    input_gs.gates = {""};
    input_gs.matrix2 = u2_to_su2(_input_matrix);
    input_gs.matrix3 = su2_to_so3(input_gs.matrix2);
    _approx_matrix = sk_decomp(input_gs, _depth);

    double d = distance(_input_matrix, _approx_matrix.matrix2);
    std::cout << "The approximation is: " << std::endl;
    std::cout << _approx_matrix.to_string() << std::endl;
    std::cout << _approx_matrix.matrix2 << std::endl;
    std::cout << "The distance is: " << d << std::endl;
}

}  // namespace qsyn::sk_decomp
