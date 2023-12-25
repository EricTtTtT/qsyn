#pragma once

#include "./sk_decomp_mgr.hpp"
#include "argparse/argparse.hpp"
#include "cli/cli.hpp"

namespace qsyn::sk_decomp {

bool valid_recursion_depth(size_t const& n, size_t const& e);
bool add_sk_decomp_cmds(dvlab::CommandLineInterface& cli, sk_decomp::SKDMgr& skd_mgr);

}