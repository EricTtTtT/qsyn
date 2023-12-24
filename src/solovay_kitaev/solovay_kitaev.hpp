#pragma once

#include <cstddef>
#include <filesystem>
#include <ranges>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "qcir/qcir_gate.hpp"
#include "qcir/gate_type.hpp"
#include "qcir/qcir_qubit.hpp"
#include "qsyn/qsyn_type.hpp"
#include "spdlog/common.h"

/*
The Solovay Kitaev discrete decomposition algorithm.
*/

namespace qsyn::sk_decomp {

    // present sequence of gate
    using GateSequence = std::vector<qsyn::qcir::GateType>;

    class SKD {
        // TODO: define the main structure
    public:
        SKD() {}
        ~SKD() = default;
        SKD(SKD const& other);
        SKD(SKD&& other) noexcept = default;

        std::vector<std::string> const& get_procedures() const { return _procedures; }
        void add_procedure(std::string const& p) { _procedures.emplace_back(p); }

        void set_filename(std::string f) { _filename = std::move(f); }
        std::string get_filename() const { return _filename; }

        bool read_skd_file(std::filesystem::path const& filepath);
        bool read_txt(std::filesystem::path const& filepath);

        void set_depth(size_t num) { _depth = num; }
        size_t get_depth() { return _depth; }
        void set_param(size_t num) { _param = num; }
        size_t get_param() { return _param; }



    private:
        size_t _depth;  // the maximal recursion depth
        size_t _param;  // decomposition parameters ε
        std::string _filename;
        std::vector<std::string> _procedures;
    };

}