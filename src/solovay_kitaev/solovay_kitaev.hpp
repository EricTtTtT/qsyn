#pragma once

#include <cstddef>
#include <filesystem>
#include <ranges>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <complex>

#include "qcir/qcir_gate.hpp"
#include "qcir/gate_type.hpp"
#include "qcir/qcir_qubit.hpp"
#include "qsyn/qsyn_type.hpp"
#include "spdlog/common.h"

#include <Eigen/Dense>

/*
The Solovay Kitaev discrete decomposition algorithm.
*/

namespace qsyn::sk_decomp {
    const double kPI = 3.14159265358979323846;
    const double init_e = 0.14; // precision required of the initial approximations
    const double c_approx = 4 * std::sqrt(2);   // magic constant
    using Complex = std::complex<double>;

    using Matrix3 = Eigen::Matrix3cd;
    using Matrix2 = Eigen::Matrix2cd;
    using MatrixX = Eigen::MatrixXcd;

    struct GateSequence {
        std::vector<std::string> gates;
        Matrix2 matrix2;
        Matrix3 matrix3;

        std::string to_string() const;
        GateSequence operator*(GateSequence const& other) const;
    };

    // Matrix operator*(Matrix const& lhs, Matrix const& rhs);
    // Matrix operator*(Complex const& lhs, Matrix const& rhs);
    // Matrix operator-(Matrix const& lhs, Matrix const& rhs);
    // double trace_dist(Matrix const& lhs, Matrix const& rhs);
    // double trace(Matrix const& matrix);
    // Matrix adjoint(Matrix const& matrix);
    // Matrix diagonalize(Matrix const& matrix);
    // Matrix sqrt(Matrix const& matrix);
    Matrix2 u2_to_su2(Matrix2 const& matrix);
    Matrix3 su2_to_so3(Matrix2 const& matrix);
    bool is_unitary(Matrix2 const& matrix);
    bool is_single_qubit(Matrix2 const& matrix);

    double distance(MatrixX const& lhs, MatrixX const& rhs);
    // double distance2(Matrix2 const& lhs, Matrix2 const& rhs);


    void group_comm_decomp(Matrix3 const& matrix, GateSequence& v, GateSequence& w);
    // std::pair<Vector3, double> u_to_bloch(Matrix const& matrix);

    class SKD {
        // TODO: define the main structure
    public:
        SKD(): _depth(0), _length(0), _param(0.14), _filename(""), _procedures({}) {}
        ~SKD() = default;
        SKD(SKD const& other);
        SKD(SKD&& other) noexcept = default;

        std::vector<std::string> const& get_procedures() const { return _procedures; }
        void add_procedure(std::string const& p) { _procedures.emplace_back(p); }

        void set_filename(std::string f) { _filename = std::move(f); }
        std::string get_filename() const { return _filename; }

        bool read_skd_file(std::filesystem::path const& filepath);
        bool read_tex(std::filesystem::path const& filepath);

        bool is_input_empty() const { return _input_matrix.isZero(); }
        bool is_input_unitary() const { return is_unitary(_input_matrix); }
        bool is_input_single_qubit() const { return is_single_qubit(_input_matrix); }

        void set_depth(size_t num)  { _depth = num;  }
        void set_length(size_t num) { _length = num; }
        void set_param(double num)  { _param = num; }
        size_t get_length() { return _length; }
        size_t get_depth()  { return _depth;  }
        size_t get_param()  { return _param;  }

        void set_basis(std::vector<std::string> const& basis);
        bool is_generated_approximations() const { return !_basis_approximations.empty(); }
        void create_basic_approximations(size_t length);

        size_t find_closest_approximation(GateSequence const& gs, bool print = false);
        
        void print_input_matrix() const;
        void report_basis();
        void report_decomp_result() const;

        GateSequence sk_decomp(GateSequence const& u, size_t depth);

        void run();

    private:
        size_t _depth;  // the maximal recursion depth
        size_t _length; // the length of approximation sequence
        double _param;  // decomposition parameters ε
        std::string _filename;
        std::vector<std::string> _procedures;

        Matrix2 _input_matrix;
        GateSequence _approx_matrix;

        std::unordered_map<std::string, Matrix2> _basis_gates;  
        std::vector<GateSequence> _basis_approximations;
    };
}