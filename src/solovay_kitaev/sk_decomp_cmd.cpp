#include "solovay_kitaev/sk_decomp_cmd.hpp"

#include <fmt/ostream.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <filesystem>
#include <ostream>
#include <string>
#include <cmath>

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

bool valid_recursion_depth(int const& n, double const& e) { 
    // If n is too small so that the achievable theoretical difference d(U, S) is greater 
    // than ε, issue an error message and forbid the users from executing the following commands.
    // double d = std::pow((init_e * c_approx * c_approx), std::pow(1.5, n)) / c_approx / c_approx;
    if (0.5 < e || n >= 0) {
        spdlog::error("e is too large!! please select an e < 0.14");
        return false;
    }
    else {
        return true;
    }
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
                    .constraint(allowed_extension({".tex", ""}))
                    .help("the filepath to a single-qubit unitary matrix. Supported extension: .tex");
                    
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
                if (!skd_mgr.get()->is_input_unitary()) {
                    spdlog::error("The input matrix is not a unitary matrix.");
                    return CmdExecResult::error;
                }
                if (!skd_mgr.get()->is_input_single_qubit()) {
                    spdlog::error("The input matrix is not a single-qubit unitary matrix.");
                    return CmdExecResult::error;
                }

                skd_mgr.get()->set_filename(std::filesystem::path{filepath}.stem());
                return CmdExecResult::done;
            }};
}

dvlab::Command sk_decomp_set_cmd(SKDMgr& skd_mgr) {
    return {"set",
            [](ArgumentParser& parser) {
                parser.description("set the parameters of the decomposition algorithm");

                parser.add_argument<double>("-e")
                    .constraint([](double const& e) { return e > 0; })
                    .help("the decomposition parameters, approximate the desired \
                           quantum gate to an accuracy within ε > 0");

                parser.add_argument<int>("-l")
                    .constraint([](int const& l) { return l > 0; })
                    .help("the length of approximation sequence");

                parser.add_argument<int>("-n")
                    .constraint([](int const& n) { return n != -1; })
                    .help("maximal recursion depth");
            },
            [&](ArgumentParser const& parser) {
                if (!skd_mgr_not_empty(skd_mgr)) {
                    return CmdExecResult::error;
                }
                if (parser.get<double>("-e")) {
                    skd_mgr.get()->set_param(parser.get<double>("-e"));
                } else {
                    spdlog::error("Decomposition parameter (\"-e\") is not defined yet!");
                    return CmdExecResult::error;
                }
                if (parser.get<int>("-l")) {
                    skd_mgr.get()->set_length(parser.get<int>("-l"));
                } else {
                    spdlog::error("The length of approximation sequence (\"-l\") is not defined yet!");
                    return CmdExecResult::error;
                }
                if (parser.get<int>("-n")) {
                    skd_mgr.get()->set_depth(parser.get<int>("-n"));
                } else {
                    spdlog::error("Maximal recursion depth (\"-n\") is not defined yet!");
                    return CmdExecResult::error;
                }
                if (!valid_recursion_depth(parser.get<int>("-n"), parser.get<double>("-e"))) {
                    return CmdExecResult::error;
                }
                return CmdExecResult::done;
            }};
}

dvlab::Command sk_decomp_report_basis_cmd(SKDMgr& skd_mgr) {
    return {"report_basis",
            [](ArgumentParser& parser) {
                parser.description("report the basic approximation of the gate set");
            },
            [&](ArgumentParser const& parser) {
                if (!skd_mgr_not_empty(skd_mgr)) return CmdExecResult::error;
                if (parser.num_parsed_args()) {}
                if (skd_mgr.get()->is_input_empty()) {
                    spdlog::error("The input matrix is empty.");
                    return CmdExecResult::error;
                }
                if (!skd_mgr.get()->is_generated_approximations()) {
                    skd_mgr.get()->set_basis({"h", "t", "s"}); 
                    skd_mgr.get()->create_basic_approximations(skd_mgr.get()->get_length());
                }
                skd_mgr.get()->report_basis();
                return CmdExecResult::done;
            }};
}

dvlab::Command sk_decomp_report_decomp_result_cmd(SKDMgr& skd_mgr) {
    return {"report_decomp_result",
            [](ArgumentParser& parser) {
                parser.description("report the decomposition circuit");
            },
            [&](ArgumentParser const& parser) {
                if (!skd_mgr_not_empty(skd_mgr)) return CmdExecResult::error;
                if (parser.num_parsed_args()) {}
                skd_mgr.get()->report_decomp_result();
                return CmdExecResult::done;
            }};
}

dvlab::Command sk_decomp_run_cmd(SKDMgr& skd_mgr) {
    return {"run",
            [](ArgumentParser& parser) {
                parser.description("Recursively decompose U into a series of single-qubit gates");
            },
            [&](ArgumentParser const& parser) {
                if (!skd_mgr_not_empty(skd_mgr)) return CmdExecResult::error;
                if (parser.num_parsed_args()) {}

                if (!skd_mgr.get()->is_generated_approximations()) {
                    skd_mgr.get()->set_basis({"h", "t", "s"}); 
                    skd_mgr.get()->create_basic_approximations(skd_mgr.get()->get_length());
                }
                skd_mgr.get()->run();
                return CmdExecResult::done;
            }};
}

Command sk_decomp_cmd(SKDMgr& skd_mgr) {
    auto cmd = dvlab::utils::mgr_root_cmd(skd_mgr);

    cmd.add_subcommand(sk_decomp_read_cmd(skd_mgr));
    cmd.add_subcommand(sk_decomp_set_cmd(skd_mgr));
    cmd.add_subcommand(sk_decomp_report_basis_cmd(skd_mgr));
    cmd.add_subcommand(sk_decomp_report_decomp_result_cmd(skd_mgr));
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
