#include "f_domain.h"

void build_f_domain(arma::cx_dvec& base_vector, std::vector<Word> words, CompMat3 H)
{
  std::vector<arma::cx_dvec> vec_images;
  vec_images.reserve(words.size());
  
  for (size_t i = 0; i < words.size(); ++i )
  {
    arma::cx_dvec temp = words[i].get_matrix() * base_vector;
    // sanity check
    std::cout << arma::trans(temp) * H * temp << "\n";
    vec_images.push_back(temp);
  }

}
