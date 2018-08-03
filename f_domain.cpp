#include "f_domain.h"

void build_f_domain(const point& base_vector,
                    const std::vector<Word>& words,
                    const CompMat3& H)
{
  std::vector<arma::cx_dvec> vec_images;
  vec_images.reserve(words.size());
  
  for (size_t i = 0; i < words.size(); ++i)
  {
    const point temp = words[i].get_matrix() * base_vector;
    // sanity check
    const comp_d vec_sign = h_n(temp, H, temp);
//    comp_d vec_sign = (arma::trans(temp) * H * temp).eval()(0, 0);
    if ((vec_sign.real() > 1e-16) || (abs(vec_sign.imag()) > 1e-6))
    {
      std::cout << vec_sign.real() << " + " << vec_sign.imag() << "i ~ ";
      std::cout << words[i].as_string() << "\n";
    }
    //std::cout << "\n";
    vec_images.push_back(temp);
  }
  
}
