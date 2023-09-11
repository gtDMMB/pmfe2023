// Copyright (c) 2015 Andrew Gainer-Dewar.

#ifndef RNA_POLYTOPE_H
#define RNA_POLYTOPE_H

#include "pmfe_types.h"
#include "rna_polytope.h"
#include "BBPolytope.h"

#include <map>
#include "boost/filesystem/fstream.hpp"

#include <CGAL/Gmpq.h>

namespace fs = boost::filesystem;

namespace pmfe {
    typedef CGAL::Gmpq Q; // We'll do all the geometry over Q
    typedef iB4e::BBPolytope<Q> BBP;

    class compare_fp {
    public:
        bool operator() (const BBP::FPoint& left, const BBP::FPoint& right) const {
            // Custom comparitor for the results map in the polytope
            for (int i = 0; i < left.dimension(); ++i) {
                if (left[i] < right[i]) {
                    return true;
                } else if (left[i] > right[i]) {
                    return false;
                }
                // Otherwise, continue
            }

            // If we make it out of the loop, the points are equal
            return false;
        }
    };

    class RNAPolytope: public BBP {
    public:
        ScoreVector classical_scores;
        RNASequence sequence;
        dangle_mode dangles;
        std::map<FPoint, RNAStructureWithScore, compare_fp> structures;
        Rational multiloop_weight;
        bool scale_b_param;

        RNAPolytope(RNASequence sequence, dangle_mode dangles, Rational multiloop_weight);
        RNAPolytope(RNASequence sequence, dangle_mode dangles);

        BBP::FPoint vertex_oracle(BBP::FVector objective);
        void write_to_file(const fs::path poly_file) const;
        BBP::FPoint remove_b_param(BBP::FPoint point, ParameterVector vec);

    protected:
        void hook_preinit();
        void hook_postinit();
        void hook_perloop(size_t confirmed);
        void hook_postloop();
    };
}
#endif
