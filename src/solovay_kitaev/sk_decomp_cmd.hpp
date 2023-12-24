#pragma once

#include "./sk_decomp_mgr.hpp"
#include "argparse/argparse.hpp"
#include "cli/cli.hpp"

namespace qsyn::sk_decomp {

std::function<bool(size_t const&)> valid_recursion_depth();

bool add_sk_decomp_cmds(dvlab::CommandLineInterface& cli, sk_decomp::SKDMgr& skd_mgr);

}