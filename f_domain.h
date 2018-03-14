#pragma once

#include "matrixfns.h"
#include "word.h"
/*
class F_Domain
{
  public:
    F_Domain();

    F_Domain(arma::cx_dvec& base_vector, std::vector<Word> words);
}
*/

void build_f_domain(arma::cx_dvec& base_vector, std::vector<Word> words, CompMat3 H);
