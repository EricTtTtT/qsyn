#pragma once

#include <cstddef>
#include <vector>

#include "./cli/cli.hpp"
#include "./solovay_kitaev.hpp"
#include "util/data_structure_manager.hpp"

namespace qsyn::sk_decomp {

// Solovay-Kitaev Decomposiotn manager
using SKDMgr = dvlab::utils::DataStructureManager<SKD>;

bool skd_mgr_not_empty(SKDMgr const& skd_mgr);

}  // namespace qsyn::sk_decomp

template <>
inline std::string dvlab::utils::data_structure_info_string(qsyn::sk_decomp::SKD const& t) {
    return fmt::format("{:<19} {}", t.get_filename().substr(0, 19),
                       fmt::join(t.get_procedures(), " ➔ "));
}

template <>
inline std::string dvlab::utils::data_structure_name(qsyn::sk_decomp::SKD const& t) {
    return t.get_filename();
}