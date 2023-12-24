#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cstddef>
#include <fstream>
#include <string>

#include "./solovay_kitaev.hpp"
#include "util/dvlab_string.hpp"
#include "util/phase.hpp"

namespace qsyn::sk_decomp {

bool SKD::read_skd_file(std::filesystem::path const& filepath) {
    auto const extension = filepath.extension();

    if (extension == ".txt")
        return read_txt(filepath);
    else if (extension == "") {
        return true;
    } else {
        spdlog::error("File format \"{}\" is not supported!!", extension);
        return false;
    }
}

bool SKD::read_txt(std::filesystem::path const& filepath) {
    // TODO: read in the text file to a matrix
    std::fstream txt_file{filepath};
    if (!txt_file.is_open()) {
        spdlog::error("Cannot open the TXT file \"{}\"!!", filepath);
        return false;
    }
    return true;
}

}
