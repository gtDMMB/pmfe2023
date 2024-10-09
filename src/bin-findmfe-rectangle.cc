// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "mfe.h"
#include "pmfe_types.h"
#include "rational.h"

#include <iostream>
#include <omp.h>
#include <string>
#include <unordered_set>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char * argv[]) {
    // Set up argument processing
    po::options_description desc("Options");
    desc.add_options()
        ("sequence", po::value<std::string>()->required(), "Sequence file")
        ("outfile,o", po::value<std::string>(), "Output file")
        ("verbose,v", po::bool_switch()->default_value(false), "Write verbose debugging output")
        ("multiloop-penalty-min,a", po::value<std::string>()->required(), "Multiloop penalty parameter min")
        ("multiloop-penalty-max,A", po::value<std::string>()->required(), "Multiloop penalty parameter max")
        ("unpaired-penalty,b", po::value<std::string>()->default_value("0"), "Unpaired base penalty parameter")
        ("branch-penalty-min,c", po::value<std::string>()->required(), "Branching helix penalty parameter min")
        ("branch-penalty-max,C", po::value<std::string>()->required(), "Branching helix penalty parameter max")
        ("dummy-scaling,d", po::value<std::string>()->default_value("1"), "Dummy scaling parameter")
        ("step-size,s", po::value<std::string>()->default_value("0.1"), "Step size for checking pmfe")
        ("dangle-model,m", po::value<int>()->default_value(1), "Dangle model")
        ("num-threads,t", po::value<int>()->default_value(0), "Number of threads")
        ("transform-input,I", po::bool_switch()->default_value(false), "Input a, b, c, d is transformed")
        ("transform-output,O", po::bool_switch()->default_value(false), "Transform structure output")
        ("parameter-output,P", po::bool_switch()->default_value(false), "Output parameters where each Structure is found")
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

    //Process step-size
    pmfe::Rational step_size = pmfe::get_rational_from_word(vm["step-size"].as<std::string>());
    if (step_size < 0.1) {
        std::cout << "Warning, step-size " << step_size << " may result in long computation." << std::endl;
    }

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

    //Process Outfile
    fs::path out_file;

    if (vm.count("outfile")) {
        out_file = fs::path(vm["outfile"].as<std::string>());
    } else {
        out_file = seq_file;
        out_file.replace_extension(".rnarect");
    }

    fs::ofstream outfile;

    if (vm["parameter-output"].as<bool>()){
        // Write out structure with params
        outfile.open(out_file);
    }


    // Set up the parameters to run findmfe     
    pmfe::Rational a = pmfe::get_rational_from_word(vm["multiloop-penalty-min"].as<std::string>());
    pmfe::Rational A = pmfe::get_rational_from_word(vm["multiloop-penalty-max"].as<std::string>());
    pmfe::Rational b = pmfe::get_rational_from_word(vm["unpaired-penalty"].as<std::string>());
    pmfe::Rational c = pmfe::get_rational_from_word(vm["branch-penalty-min"].as<std::string>());
    pmfe::Rational C = pmfe::get_rational_from_word(vm["branch-penalty-max"].as<std::string>());
    pmfe::Rational d = pmfe::get_rational_from_word(vm["dummy-scaling"].as<std::string>());

    // Setup dangle model
    pmfe::dangle_mode dangles = pmfe::convert_to_dangle_mode(vm["dangle-model"].as<int>());

    //Run findmfe
    pmfe::ParameterVector params = pmfe::ParameterVector();

    // float a_loc = a.to_double();
    // float A_loc = A.to_double();
    // float c_loc = c.to_double();
    // float C_loc = C.to_double();
    pmfe::Rational temp_c = c;
    pmfe::RNAStructureWithScore result = pmfe::mfe(seq_file, params, dangles);
    std::unordered_set<std::string> structure_set;
    pmfe::RNAStructureWithScore prev_result = result;

    while (a <= A) {
        temp_c = c;
        while(temp_c <= C) {    
            params.multiloop_penalty = a;
            params.unpaired_penalty = b;
            params.branch_penalty = temp_c;
            params.dummy_scaling = d;
            
            // Check if input is transformed
            if (vm["transform-input"].as<bool>()) {
                params.untransform_params();
            };

            params.canonicalize();

            result = pmfe::mfe(seq_file, params, dangles);
            result.transformed = vm["transform-output"].as<bool>();
            
            // Print strucure for params if it hasn't already been found
            if (structure_set.emplace(result.string()).second) {
                std::cout << result << std::endl;
            }
            
            if (vm["parameter-output"].as<bool>()){
                // Write out structure with params
                // fs::ofstream outfile(out_file);

                if (!outfile.is_open()) {
                    std::stringstream error_message;
                    error_message << "Output file " << out_file << "is invalid.";
                    throw std::invalid_argument(error_message.str());
                }
                
                outfile << params.multiloop_penalty << ", " << params.unpaired_penalty 
                        << ", " << params.branch_penalty << ", " << result << std::endl;
            }

            temp_c += step_size;
        }
        a += step_size;
    }


    return(0);
}
