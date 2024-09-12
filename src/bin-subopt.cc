// Copyright (c) 2015 Andrew Gainer-Dewar.

#include <iostream>
#include <omp.h>
#include <stdexcept>

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/program_options.hpp"

#include "pmfe_types.h"
#include "subopt.h"
#include "rational.h"

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("verbose,v", po::bool_switch()->default_value(false), "Write verbose debugging output")
        ("outfile,o", po::value<std::string>(), "Output file")
        ("consoleout,C", po::bool_switch()->default_value(false), "Output to console")
        ("delta", po::value<std::string>()->default_value("0"), "Energy delta value")
        ("multiloop-penalty,a", po::value<std::string>(), "Multiloop penalty parameter")
        ("unpaired-penalty,b", po::value<std::string>(), "Unpaired base penalty parameter")
        ("branch-penalty,c", po::value<std::string>(), "Branching helix penalty parameter")
        ("dummy-scaling,d", po::value<std::string>(), "Dummy scaling parameter")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
        ("sorted,s", po::bool_switch(), "Sort results in increasing energy order")
        ("num-threads,t", po::value<int>()->default_value(0), "Number of threads")
        ("transformed-input,I", po::bool_switch()->default_value(false), "Input a, b, c, d is transformed")
        ("transform-output,O", po::bool_switch()->default_value(false), "Transform structure output")
        ("help,h", "Display this help message")
        ;

    po::positional_options_description p;
    p.add("sequence", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    if (vm.count("help") or argc == 1) {
        std::cout << desc << std::endl;
        return 1;
    };

    po::notify(vm);

    // Process thread-related options
    size_t num_threads = (vm["num-threads"].as<int>());
    omp_set_num_threads(num_threads);

    // Process logging-related options
    bool verbose = vm["verbose"].as<bool>();
    if (verbose) {
        boost::log::core::get()->set_filter(
            boost::log::trivial::severity >= boost::log::trivial::debug);
    } else {
        boost::log::core::get()->set_filter
            (boost::log::trivial::severity >= boost::log::trivial::warning);
    }

    // Process file-related options
    fs::path seq_file(vm["sequence"].as<std::string>());
    fs::path out_file;
    if (vm.count("outfile")) {
        out_file = fs::path(vm["outfile"].as<std::string>());
    } else {
        out_file = seq_file;
        out_file.replace_extension(".rnasubopt");
    }

    // Set up the parameters
    pmfe::Rational delta = pmfe::get_rational_from_word(vm["delta"].as<std::string>());

    pmfe::ParameterVector params = pmfe::ParameterVector();

    if (vm.count("multiloop-penalty")) {
        params.multiloop_penalty = pmfe::get_rational_from_word(vm["multiloop-penalty"].as<std::string>());
    };

    if (vm.count("unpaired-penalty")) {
        params.unpaired_penalty = pmfe::get_rational_from_word(vm["unpaired-penalty"].as<std::string>());
    };

    if (vm.count("branch-penalty")) {
        params.branch_penalty = pmfe::get_rational_from_word(vm["branch-penalty"].as<std::string>());
    };

    if (vm.count("dummy-scaling")) {
        params.dummy_scaling = pmfe::get_rational_from_word(vm["dummy-scaling"].as<std::string>());
    };

    if (vm["transformed-input"].as<bool>()) {
        params.untransform_params();
    };

    params.canonicalize();

    bool sorted = vm["sorted"].as<bool>();
    bool transform = vm["transform-output"].as<bool>();

    // Set up dangle model
    pmfe::dangle_mode dangles = pmfe::convert_to_dangle_mode(vm["dangle-model"].as<int>());

    // Get results
    
    std::vector<pmfe::RNAStructureWithScore> structures = suboptimal_structures(seq_file, params, dangles, delta, sorted, transform);

    // Print some status information
    std::cout << "Found " << structures.size() << " suboptimal structures." << std::endl;

    // If input was transformed, transform back
    if (vm["transformed-input"].as<bool>()) {
        params.transform_params();
    };

    std::stringstream outbuffer;
        
    outbuffer << "#\tSuboptimal secondary structures within " << delta.get_d() << " of minimum energy." << std::endl;
    outbuffer << "#\tCoefficients:\t" <<
        "a = " << params.multiloop_penalty << " ≈ " << params.multiloop_penalty.get_d() << ",\t" <<
        "b = " << params.unpaired_penalty << " ≈ " << params.unpaired_penalty.get_d() << ",\t" <<
        "c = " << params.branch_penalty << " ≈ " << params.branch_penalty.get_d() << ",\t" <<
        "d = " << params.dummy_scaling << " ≈ " << params.dummy_scaling.get_d() << "." << std::endl;

    pmfe::RNASequence seq(seq_file);
    outbuffer << "#\t" << seq << "\tM\tU\tB\tw\tEnergy" << std::endl << std::endl;

    for (unsigned int i = 0; i < structures.size(); ++i) {
        outbuffer << i << "\t" << structures[i] << "\t≅ " << structures[i].score.energy.get_d() << std::endl;
    }

    if (vm["consoleout"].as<bool>()) {
        std::cout << outbuffer.str();
    } else {
        // Write out suboptimal structures
        fs::ofstream outfile(out_file);

        if (!outfile.is_open()) {
            std::stringstream error_message;
            error_message << "Output file " << out_file << "is invalid.";
            throw std::invalid_argument(error_message.str());
        }

        outfile << outbuffer.str();
    }
}
