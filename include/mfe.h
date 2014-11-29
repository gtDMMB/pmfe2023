#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

#include "parametrizer_types.h"
#include <gmpxx.h>

ScoreVector mfe(std::string seq_file, std::string output_file, std::string param_dir, ParameterVector params, int dangle_model = 1);

ScoreVector mfe(std::string seq_file, std::string output_file, std::string param_dir, int dangle_model = 1);

void init_fold(const char* seq, ParameterVector params);
void free_fold(int len);

#endif
