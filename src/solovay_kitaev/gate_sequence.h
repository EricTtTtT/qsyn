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

namespace qsyn::sk_decomp {

class GateSequence {
public:
    GateSequence();

    // Add a gate to the sequence
    void addGate(const qsyn::qcir::QCirGate& gate);

    // Remove a gate from the sequence
    void removeGate(std::size_t index);

    // Get the sequence of gates
    std::vector<qsyn::qcir::QCirGate> getSequence() const;


private:
    std::vector<qsyn::qcir::QCirGate> gates;
    
};

}  // namespace qsyn::sk_decomp