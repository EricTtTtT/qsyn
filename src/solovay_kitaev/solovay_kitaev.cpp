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

#include <fsolve/fsolve.hpp>

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



Matrix3 rotate_angle(double const& angle, Vector3 const& axis) {
    Matrix3 cross_product_matrix {
        {0., -axis(2), axis(1)},
        {axis(2), 0., -axis(0)},
        {-axis(1), axis(0), 0.}
    };
    Matrix3 res = std::cos(angle) * Matrix3::Identity() + std::sin(angle) * cross_product_matrix;
    res += (1 - std::cos(angle)) * axis * axis.adjoint();
    return res;
}

double solve_phi(Matrix3 const& matrix) {
    // TODO
    std::cout << matrix << std::endl;
    return 0;
}

Matrix3 get_rotation_between(Vector3 const& from, Vector3 const& to) {
    Vector3 from_vector = from.normalized();
    Vector3 to_vector = to.normalized();
    double dot = from_vector.dot(to_vector);

    Vector3 cross_v = from_vector.cross(to_vector).normalized();
    Matrix3 cross {
        {0., -cross_v(2), cross_v(1)},
        {cross_v(2), 0., -cross_v(0)},
        {-cross_v(1), cross_v(0), 0.}
    };
    Matrix3 rotation_matrix = Matrix3::Identity() + cross + cross * cross / (1 + dot);
    return rotation_matrix;
}

Matrix3 get_commutator(Matrix3 const& lhs, Matrix3 const& rhs) {
    Matrix3 a_dagger = lhs.adjoint();
    Matrix3 b_dagger = rhs.adjoint();
    return lhs * rhs * a_dagger * b_dagger;
}

Vector3 get_rotation_axis(Matrix3 const& matrix) {
    double trace = matrix.trace();
    double theta = std::acos(0.5 * (trace - 1));
    double x, y, z;
    if (std::sin(theta) > 1e-10) {
        x = 1 / (2 * std::sin(theta)) * (matrix(2, 1) - matrix(1, 2));
        y = 1 / (2 * std::sin(theta)) * (matrix(0, 2) - matrix(2, 0));
        z = 1 / (2 * std::sin(theta)) * (matrix(1, 0) - matrix(0, 1));
    } else {
        x = 1.0;
        y = 0.0;
        z = 0.0;
    }
    return {x, y, z};
}


double distance(Matrix3 const& lhs, Matrix3 const& rhs) {
    double distance = 0.0;
    Matrix3 diff = lhs - rhs;
    distance = diff.norm();
    distance = std::sqrt(distance);
    return distance;
}

double distance(Matrix2 const& lhs, Matrix2 const& rhs) {
    double distance = 0.0;
    Matrix2 diff = lhs - rhs;
    distance = diff.norm();
    distance = std::sqrt(distance);
    return distance;
}

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

Matrix2 so3_to_su2(Matrix3 const& matrix) {
    Matrix3 m_rounded = matrix;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            m_rounded(i, j) = std::round(m_rounded(i, j) * 1e10) / 1e10;
        }
    }

    double theta, psi, phi;
    if (m_rounded(2, 0) != 1 && m_rounded(2, 1) != -1) {
        theta = -std::asin(m_rounded(2, 0));
        psi = std::atan2(m_rounded(2, 1) / std::cos(theta), m_rounded(2, 2) / std::cos(theta));
        phi = std::atan2(m_rounded(1, 0) / std::cos(theta), m_rounded(0, 0) / std::cos(theta));
    } else {
        phi = 0;
        if (m_rounded(2, 0) == 1) {
            theta = kPI / 2;
            psi = phi + std::atan2(m_rounded(0, 1), m_rounded(0, 2));
        } else {
            theta = -kPI / 2;
            psi = -phi + std::atan2(-m_rounded(0, 1), -m_rounded(0, 2));
        }
    }

    Matrix2 uz_phi {
        {std::exp(Complex(0.,-0.5) * phi), 0},
        {0, std::exp(Complex(0.5) * phi)}
    };
    Matrix2 uy_theta {
        {std::cos(theta / 2), std::sin(theta / 2)},
        {-std::sin(theta / 2), std::cos(theta / 2)}
    };
    Matrix2 ux_psi {
        {std::cos(psi / 2), std::sin(psi / 2) * Complex(0., 1.)},
        {std::sin(psi / 2) * Complex(0., 1.), std::cos(psi / 2)}
    };
    return uz_phi * uy_theta * ux_psi;
}

bool is_unitary(Matrix2 const& matrix) {
    // Complex x = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
    Complex x = matrix.determinant();
    double det = std::abs(x);
    return det - 1.0 < 1e-4;
}

bool is_single_qubit(Matrix2 const& matrix) {
    return matrix.size() == 4;
}

void group_comm_decomp(Matrix3 const& matrix3, GateSequence& v, GateSequence& w) {
    double phi = solve_phi(matrix3);
    Matrix3 vx = rotate_angle(phi, {1.0, 0.0, 0.0});
    Matrix3 vy = rotate_angle(phi, {0.0, 1.0, 0.0});
    Matrix3 commutator = get_commutator(vx, vy);

    Vector3 matrix_axis = get_rotation_axis(matrix3);
    Vector3 commutator_axis = get_rotation_axis(commutator);

    Matrix3 sim_matrix = get_rotation_between(matrix_axis, commutator_axis);
    Matrix3 sim_matrix_dagger = sim_matrix.adjoint();

    v.gates = {"Rz(" + std::to_string(phi) + ")"};
    v.matrix2 = so3_to_su2(v.matrix3);
    v.matrix3 = sim_matrix * vx * sim_matrix_dagger;
    w.gates = {"Ry(" + std::to_string(phi) + ")"};
    w.matrix2 = so3_to_su2(w.matrix3);
    w.matrix3 = sim_matrix * vy * sim_matrix_dagger;
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

    _basis_gs.clear();
    for (auto const& [name, matrix2] : _basis_gates) {
        GateSequence gs;
        gs.gates = {name};
        gs.matrix2 = u2_to_su2(matrix2);
        gs.matrix3 = su2_to_so3(gs.matrix2);
        _basis_gs.emplace_back(gs);
    }

    std::vector<GateSequence> curr_ = _basis_gs;
    _basis_approximations.clear();
    for (size_t i = 0; i < length; ++i) {
        for (auto const& gs : curr_) {
            _basis_approximations.emplace_back(gs);
        }
        std::cout << "generating approximations of length " << i + 1 << std::endl;
        std::vector<GateSequence> next_;
        for (auto const& gs_c : curr_) {
            for (GateSequence const& gs_b : _basis_gs) {
                GateSequence gs;
                gs.gates = gs_c.gates;
                gs.gates.emplace_back(gs_b.gates[0]);
                gs.matrix2 = gs_b.matrix2 * gs_c.matrix2;
                gs.matrix3 = gs_b.matrix3 * gs_c.matrix3;
                next_.emplace_back(gs);
            }
        }
        curr_ = next_;
    }

    spdlog::debug("The number of approximations is: {}", _basis_approximations.size());
    for (auto const& gs : _basis_approximations) {
        std::cout << "  " << gs.to_string() << std::endl;
        std::cout << gs.matrix3 << std::endl;
        std::cout << gs.matrix2 << std::endl;
    }
}

size_t SKD::find_closest_approximation(GateSequence const& gs, bool print) {
    size_t report_num = 3;
    
    std::priority_queue<std::pair<size_t, double>,
        std::vector<std::pair<size_t, double>>,
        std::function<bool(std::pair<size_t, double> const&, std::pair<size_t, double> const&)>>
        approx_max_heap([](std::pair<size_t, double> const& lhs, std::pair<size_t, double> const& rhs) {
            return lhs.second < rhs.second;
        });

    for (size_t i=0; i<_basis_approximations.size(); ++i) {
        // double d = distance(gs.matrix2, _basis_approximations[i].matrix2);
        double d = distance(gs.matrix3, _basis_approximations[i].matrix3);
        if (approx_max_heap.size() < report_num) {
            approx_max_heap.push({i, d});
        } else if (d < approx_max_heap.top().second) {
            approx_max_heap.pop();
            approx_max_heap.push({i, d});
        }
    }

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
    std::vector<std::pair<size_t, double>> id_dist;
    for (size_t i=0; i<_basis_approximations.size(); ++i) {
        double d = distance(_input_matrix, _basis_approximations[i].matrix2);
        id_dist.emplace_back(std::make_pair(i, d));
    }
    std::sort(id_dist.begin(), id_dist.end(), [](std::pair<size_t, double> const& lhs, std::pair<size_t, double> const& rhs) {
        return lhs.second < rhs.second;
    });

    size_t max_count = 20;
    for (auto it = id_dist.begin(); it != id_dist.end(); ++it) {
        size_t gs_idx = it->first;
        std::string gs_name = _basis_approximations[gs_idx].to_string();
        std::cout << "  gate: " << std::setw(20) << gs_name << ", distance=" << it->second << std::endl;
        if (max_count-- == 0) break;
    }
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

    GateSequence u_next = sk_decomp(u, depth - 1);
    GateSequence v, w;
    group_comm_decomp(u.matrix3, v, w);
    GateSequence v_next = sk_decomp(v, depth - 1);
    GateSequence w_next = sk_decomp(w, depth - 1);
    GateSequence ret;
    ret.matrix3 = v_next.matrix3 * w_next.matrix3 * v_next.matrix3.adjoint() * w_next.matrix3.adjoint() * u_next.matrix3;
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

    _approx_matrix.matrix2 = so3_to_su2(_approx_matrix.matrix3);

    double d = distance(_input_matrix, _approx_matrix.matrix2);
    std::cout << "The approximation is: " << std::endl;
    std::cout << _approx_matrix.to_string() << std::endl;
    std::cout << _approx_matrix.matrix2 << std::endl;
    std::cout << "The distance is: " << d << std::endl;
}

}  // namespace qsyn::sk_decomp
