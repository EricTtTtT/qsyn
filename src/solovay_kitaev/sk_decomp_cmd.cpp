#include "solovay_kitaev/sk_decomp_cmd.hpp"

#include <fmt/ostream.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <filesystem>
#include <ostream>
#include <string>

#include "argparse/arg_parser.hpp"
#include "argparse/arg_type.hpp"
#include "cli/cli.hpp"
#include "solovay_kitaev/solovay_kitaev.hpp"
#include "solovay_kitaev/sk_decomp_mgr.hpp"
#include "util/cin_cout_cerr.hpp"
#include "util/data_structure_manager_common_cmd.hpp"
#include "util/dvlab_string.hpp"
#include "util/phase.hpp"
#include "util/util.hpp"

using namespace dvlab::argparse;
using dvlab::CmdExecResult;
using dvlab::Command;

namespace qsyn::sk_decomp {

std::function<bool(size_t const&)> valid_recursion_depth() {
    return [&](size_t const& n) {
        if (n >= 1) return true;
        /* 
        TODO: If n is too small so that the achievable theoretical difference d(U, S) is greater 
        than ε, issue an error message and forbid the users from executing the following commands.
        */
        /*
        if (qcir_mgr.is_id(id)) return true;
        spdlog::error("QCir {} does not exist!!", id);
        return false;
        */
        return false;
    };
}

bool skd_mgr_not_empty(SKDMgr const& skd_mgr) {
    if (skd_mgr.empty()) {
        spdlog::error("SK decomp. is empty. Please create a SK Decomp. first!!");
        spdlog::info("Use sk_decomp read to add a new matrix from a file.");
        return false;
    }
    return true;
}

dvlab::Command sk_decomp_read_cmd(SKDMgr& skd_mgr) {
    return {"read",
            [](ArgumentParser& parser) {
                parser.description("read in a single-qubit unitary matrix, and report the basic approximation of the gate set");

                parser.add_argument<std::string>("filepath")
                    .constraint(path_readable)
                    .constraint(allowed_extension({".txt", ""}))
                    .help("the filepath to a single-qubit unitary matrix. Supported extension: .txt");
                    
            },
            [&](ArgumentParser const& parser) {
                SKD buffer_skd;
                auto filepath = parser.get<std::string>("filepath");
                if (!buffer_skd.read_skd_file(filepath)) {
                    fmt::println("Error: the format in \"{}\" has something wrong!!", filepath);
                    return CmdExecResult::error;
                }
                if (skd_mgr.empty()) {
                    skd_mgr.add(skd_mgr.get_next_id(), std::make_unique<SKD>(std::move(buffer_skd)));
                } else {
                    skd_mgr.set(std::make_unique<SKD>(std::move(buffer_skd)));
                }
                skd_mgr.get()->set_filename(std::filesystem::path{filepath}.stem());
                return CmdExecResult::done;
            }};
}

dvlab::Command sk_decomp_run_cmd(SKDMgr& skd_mgr) {
    return {"run",
            [](ArgumentParser& parser) {
                parser.description("Recursively decompose U into a series of single-qubit gates");
                
                parser.add_argument<size_t>("-e")
                    .help("the decomposition parameters, approximate the desired quantum gate to an accuracy within ε > 0");

                parser.add_argument<size_t>("-n")
                    .constraint(valid_recursion_depth())
                    .help("the maximal recursion depth");
            },
            [&](ArgumentParser const& parser) {
                if (!skd_mgr_not_empty(skd_mgr)) return CmdExecResult::error;

                skd_mgr.get()->set_param(parser.parsed("-e") ? parser.get<size_t>("-e") : 1);
                skd_mgr.get()->set_depth(parser.parsed("-n") ? parser.get<size_t>("-n") : 1);
                return CmdExecResult::done;
            }};
}

Command sk_decomp_cmd(SKDMgr& skd_mgr) {
    
    auto cmd = dvlab::utils::mgr_root_cmd(skd_mgr);

    cmd.add_subcommand(sk_decomp_read_cmd(skd_mgr));
    cmd.add_subcommand(sk_decomp_run_cmd(skd_mgr));
    return cmd;
}

bool add_sk_decomp_cmds(dvlab::CommandLineInterface& cli, qsyn::sk_decomp::SKDMgr& skd_mgr) {

    if (!cli.add_command(sk_decomp_cmd(skd_mgr))) {
        spdlog::error("Registering \"sk_decomp\" commands fails... exiting");
        return false;
    }
    return true;
}

}  // namespace qsyn::sk_decomp
