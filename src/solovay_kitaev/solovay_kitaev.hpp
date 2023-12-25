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

/*
The Solovay Kitaev discrete decomposition algorithm.
*/

namespace qsyn::sk_decomp {
    const double kPI = 3.14159265358979323846;
    using Complex = std::complex<double>;
    using Matrix = std::vector<std::vector<Complex>>;
    using Vector3 = std::array<double, 3>;

    // present sequence of gate
    using Gate = qsyn::qcir::QCirGate;
    using GateSequence = std::vector<Gate>;

    Matrix operator*(Matrix const& lhs, Matrix const& rhs);
    double distance(Matrix const& lhs, Matrix const& rhs);
    Matrix adjoint(Matrix const& matrix);
    Matrix diagonalize(Matrix const& matrix);
    bool is_unitary(Matrix const& matrix);
    bool is_single_qubit(Matrix const& matrix);
    std::pair<Matrix, Matrix> group_comm_decomp(Matrix const& matrix);
    std::pair<Vector3, double> u_to_bloch(Matrix const& matrix);

    class SKD {
        // TODO: define the main structure
    public:
        SKD(): _depth(0), _length(0), _param(0.0), _filename(""), _procedures({}) {}
        ~SKD() = default;
        SKD(SKD const& other);
        SKD(SKD&& other) noexcept = default;

        std::vector<std::string> const& get_procedures() const { return _procedures; }
        void add_procedure(std::string const& p) { _procedures.emplace_back(p); }

        void set_filename(std::string f) { _filename = std::move(f); }
        std::string get_filename() const { return _filename; }

        bool read_skd_file(std::filesystem::path const& filepath);
        bool read_tex(std::filesystem::path const& filepath);

        bool is_input_empty() const { return _input_matrix.empty(); }
        bool is_input_unitary() const { return is_unitary(_input_matrix); }
        bool is_input_single_qubit() const { return is_single_qubit(_input_matrix); }

        void set_depth(size_t num) { _depth = num; }
        size_t get_depth() { return _depth; }
        void set_length(size_t num) { _length = num; }
        size_t get_length() { return _length; }
        void set_param(size_t num) { _param = num; }
        size_t get_param() { return _param; }

        void set_basis(std::vector<std::string> const& basis);
        bool is_generated_approximations() const { return !_basis_approximations.empty(); }
        void create_basic_approximations(int depth);
        std::string find_closest_approximation(Matrix const& matrix, std::unordered_map<std::string, Matrix> const& approximations, bool print = false);
        
        void report_basis();
        void report_decomp_result() const;

        Matrix sk_decomp(Matrix const& u, size_t depth);

        void run();

    private:
        size_t _depth;  // the maximal recursion depth
        size_t _length; // the length of approximation sequence
        double _param;  // decomposition parameters ε
        std::string _filename;
        std::vector<std::string> _procedures;

        Matrix _input_matrix;

        std::unordered_map<std::string, Matrix> _basis_gates;  
        std::unordered_map<std::string, Matrix> _basis_approximations;
    };




}