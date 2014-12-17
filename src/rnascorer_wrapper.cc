// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "RNAScoring.h"
#include <string>
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include <gmpxx.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char * argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("structure", po::value<std::string>()->required(), "Structure file")
        ("paramdir,p", po::value<std::string>()->default_value("rna-scoring/data/Turner99"), "Turner99 parameter directory")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
        ("help,h", "Display this help message")
        ;

    po::positional_options_description p;
    p.add("structure", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    if (vm.count("help") or argc == 1) {
        std::cout << desc << std::endl;
        return 1;
    };

    po::notify(vm);

    // Process file-related options
    fs::path struct_file, param_dir;

    struct_file = fs::path(vm["structure"].as<std::string>());
    param_dir = fs::path(vm["paramdir"].as<std::string>());

    // Setup dangle model
    int dangle_model = vm["dangle-model"].as<int>();

    mpq_class result = rnascoring::get_classical_score(struct_file.native(), param_dir.native(), dangle_model);
    printf("Computed energy %s = %5.3f\n", result.get_str(10).c_str(), result.get_d());
    return(0);
}
