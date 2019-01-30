// g++ -o test.out *.cpp -O1 -larmadillo -std=c++11

#include "getwordsupton.h"
#include "getbraidrelns.h"
#include "f_domain.h"

#include "reduction.h"

#define DEBUG(x) std::cout << x << "\n";

Word get_generator(const Generator gen, const std::vector<Word>& generators)
{
  for (size_t i = 0; i < generators.size(); ++i)
  {
    if (gen == generators[i].first_element())
      return generators[i];
  }
  return Word();
}

std::vector<Word> is_power(const std::vector<Word>& seen_words,
                           const unsigned int max_power=100,
                           const bool v=true)
{
  std::vector<Word> add_words;
  for (size_t i = 0; i < seen_words.size(); ++i)
  {
    const Word word_to_check = seen_words[i];
    bool found_match = false;
    for (size_t j = 0; j < add_words.size(); ++j)
    {
      const Word base_word = add_words[j];
      if (word_to_check.is_equal(base_word))
        continue;
      if (word_to_check.word_length() < base_word.word_length())
        continue;
      found_match = false;
      for (int p=1; p <= base_word.get_order(); ++p)
      {
        const Word pow_word = power(base_word, p);
        if (word_to_check.is_equal(pow_word))
        {
          if (v)
            std::cout << word_to_check.as_string() << " = (" << base_word.as_string() << ")^" << p
                      << " | " << word_to_check.get_order() << " ~ " << base_word.get_order() << "\n";
          found_match = true;
          break;
        }
      }
      if (found_match)
        break;
    }
    if (!found_match)
      add_words.push_back(word_to_check);
  }
  std::cout << add_words.size() << "\n";
  return add_words;
}

std::vector<Word> reduce_words(const std::vector<Word>& seen_words,
                               const std::vector<Word>& generators,
                               const bool keep_loxo=false)
{
auto start = std::chrono::high_resolution_clock::now();

  std::vector<Word> reduced_words;
  reduced_words.reserve(seen_words.size());
  for (size_t i = 0; i < seen_words.size(); ++i)
  {
    Word comp_word = seen_words[i];
    if (comp_word.get_isom_class() == IsomClass::Loxodromic && !keep_loxo)
      continue;

    while (comp_word.first_element() == inverse(comp_word.last_element()) &&
           comp_word.word_length() > 0)
    {
      Word P = get_generator(inverse(comp_word.first_element()), generators);
      comp_word = conjugate(comp_word, P);
    }

    size_t rot = comp_word.word_length();
    for (size_t i = 0; i < rot; ++i)
    {
      Word P = get_generator(inverse(comp_word.first_element()), generators);
      comp_word = conjugate(comp_word, P);     
    }

    bool matched = false;
    for (size_t j = 0; j < reduced_words.size(); ++j)
    {
      const Word test_word = reduced_words[j];
      // is it equal to a seen word or the inverse of a seen word
      if (comp_word.is_equal(test_word)) { matched = true; break; }
      if (comp_word.is_equal_inverse(test_word)) { matched = true; break; }

      // is it equal to a the conjugate of a seen word or the conjugate of an inverse of a seen word
      for (size_t k = 0; k < generators.size(); ++k)
      {
        const Word conj_test_word = conjugate(test_word, generators[k]);
        if (comp_word.is_equal(conj_test_word)) { matched = true; break; }
        if (comp_word.is_equal_inverse(conj_test_word)) { matched = true; break; }
      }

      if (matched)
        break;
    }
    if (!matched)
      reduced_words.push_back(comp_word);
  }
  return reduced_words;
}

std::vector<Word> kill_conjugates(const std::vector<Word>& seen_words,
                                  const std::vector<Word>& generators)
{
  std::vector<Word> reduced_words;
  reduced_words.reserve(seen_words.size());
  for (size_t i = 0; i < seen_words.size(); ++i)
  {
    const Word comp_word = seen_words[i];
    Word conj_test_word = seen_words[i];
    bool matched = false;
    for (size_t j = 0; j < reduced_words.size(); ++j)
    {
      Word test_word = reduced_words[j];
      size_t rot_len = conj_test_word.word_length();
      for (size_t k = 0; k < rot_len; ++k)
      {
        const Word last_char = get_generator(conj_test_word.last_element(), generators);          
        conj_test_word = conjugate(conj_test_word, last_char);
        if (conj_test_word.is_equal(test_word)) { matched = true; break; }
        if (conj_test_word.is_equal_inverse(test_word)) { matched = true; break; }
      }
      if (matched)
        break;
    }
    if (!matched)
      reduced_words.push_back(comp_word);
  }
  std::cout << reduced_words.size() << "\n";
  return reduced_words;
}

void print_word_vector(const std::vector<Word>& word_vec,
                       const bool print_loxo=false)
{
  std::cout << "|\tWord\t|\tOrder\t|\n";
  for (size_t i = 0; i < word_vec.size(); ++i)
  {
    if (word_vec[i].get_isom_class() == IsomClass::Loxodromic && !print_loxo)
      continue;

    std::cout << "|\t" << word_vec[i].as_string() << "\t|\t"
              << word_vec[i].str_isom_class() << "\t|";

    std::cout << "\n";
  }
}

void print_c_and_p(const std::vector<Word>& word_vec)
{
  for (size_t i = 0; i < word_vec.size(); ++i)
  {
    if (word_vec[i].get_isom_class() == IsomClass::Loxodromic)
      continue;

    for (int k = 0; k < word_vec[i].get_order(); ++k)
      std::cout << word_vec[i].as_string();

    std::cout << "\n";
  }
}

int main(int argc, char *argv[])
{
  unsigned int upto = 5;
  unsigned int ref_order = 2;
  if (argc > 1)
  {
    unsigned int val;

    // get upto
    std::istringstream iss_upto( argv[1] );
    if (iss_upto >> val)
      upto = val;

    // get reflection order
    std::istringstream iss_order( argv[2] );
    if (iss_order >> val)
      ref_order = val;
  }
  else
  {
    std::cout << "enter a number\n";
    return 0;
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Set up the generators R1, R2, R3, E1, E2, E3
  //
  //
  const double ref_angle = 2.0 * PI / (double)ref_order;
  const comp_d psi(std::cos(ref_angle / 3.0), std::sin(ref_angle / 3.0));
  const comp_d conj_root_psi(std::cos(-ref_angle / 6.0), std::sin(-ref_angle / 6.0));
  const comp_d comp_zero(0.0, 0.0);

  const comp_d rho(1.0 / 2.0, std::sqrt(7.0) / 2.0);
  const comp_d sigma(1.0, 0.0);
  const comp_d tau = sigma;

/*
  const comp_d rho(1.0, 0.0);
  const comp_d sigma(std::sqrt(2), 0.0);
  const comp_d tau = sigma;
*/
/*
  const comp_d rho = 2.0 * std::cos(PI / 5.0) * comp_d(std::cos(2.0 * PI / 5.0), std::sin(2.0 * PI / 5.0));
  const comp_d sigma = comp_d(1.0, 0.0);
  const comp_d tau = sigma;
*/
/*
  const comp_d rho = 0.25 * comp_d(1.0, std::sqrt(3.0)) * comp_d(std::sqrt(5.0), -std::sqrt(3.0));
  const comp_d sigma(1.0 / 2.0, std::sqrt(5.0) / 2.0);
  const comp_d tau = sigma;
*/
/*
  const comp_d rho = comp_d(1.0, 1.0) * comp_d(1.0, -1.0 / std::sqrt(2.0));
  const comp_d sigma = std::sqrt(2 + std::sqrt(2));
  const comp_d tau = sigma;
*/
/*
  const double k = 4.0;
  const comp_d rho = 2.0 * std::cos(PI / k) * comp_d(std::cos(PI / k), std::sin(PI / k));
  const comp_d sigma = 2.0 * std::cos(PI / k);
  const comp_d tau = sigma;
*/
/*
  const comp_d u(std::cos(2.0 * PI / 7.0), std::sin(2.0 * PI / 7.0));
  const comp_d rho = u + u*u + u*u*u*u;
  const comp_d sigma = u*u*u*u*u;
  const comp_d tau = sigma;
*/

  const comp_d alpha = std::sqrt(2.0 - 2.0 * std::real(psi * psi * psi));
  const comp_d beta_1 = comp_d(0.0, -1.0) * conj_root_psi * rho;
  const comp_d beta_2 = comp_d(0.0, -1.0) * conj_root_psi * sigma;
  const comp_d beta_3 = comp_d(0.0, -1.0) * conj_root_psi * tau;
  CompMat3 mat_H(3, 3);
  mat_H << alpha << beta_1 << std::conj(beta_3) << arma::endr
        << std::conj(beta_1) << alpha << beta_2 << arma::endr
        << beta_3 << std::conj(beta_2) << alpha << arma::endr;
  print_e_structure(mat_H, "Hermitian form H");

  CompMat3 mat_S(3, 3);
  mat_S << rho << psi * std::sqrt(1 - 2 * std::real(rho)) << psi * psi * std::sqrt(2 * std::real(rho)) << arma::endr
        << std::conj(psi) << 0.0 << 0.0 << arma::endr
        << 0 << std::conj(psi) * std::sqrt(2 * std::real(rho)) << -1.0 << arma::endr;

  ref_order = 2;
  const comp_d u(std::cos(2.0 * PI / 7.0), std::sin(2.0 * PI / 7.0));

  CompMat3 mat_R1(3, 3);
  //mat_R1 << psi * psi << rho            << -psi * std::conj(tau) << arma::endr
  //       << 0.0       << std::conj(psi) << 0.0                   << arma::endr
  //       << 0.0       << 0.0            << std::conj(psi)        << arma::endr;
  mat_R1 << 1.0 << u + u*u + u*u*u*u << u*u  << arma::endr
         << 0.0 << -1.0              << 0.0  << arma::endr
         << 0.0 << 0.0               << -1.0 << arma::endr;
  std::vector<Generator> vec_R1{Generator::R1};
  Word word_R1(vec_R1, mat_R1, IsomClass::Reflection,
               arma::trace(mat_R1), (int)ref_order);

  std::vector<Generator> vec_E1{Generator::E1};
  Word word_E1(vec_E1, arma::inv(mat_R1), IsomClass::Reflection,
               arma::trace(arma::inv(mat_R1)), (int)ref_order);

  CompMat3 mat_R2(3, 3);
  //mat_R2 << std::conj(psi)        << 0.0       << 0.0            << arma::endr
  //       << -psi * std::conj(rho) << psi * psi << sigma          << arma::endr
  //       << 0.0                   << 0.0       << std::conj(psi) << arma::endr;
  mat_R2 << -1.0                         << 0.0 << 0.0       << arma::endr
         << std::conj(u + u*u + u*u*u*u) << 1.0 << u*u*u*u*u << arma::endr
         << 0.0                          << 0.0 << -1.0      << arma::endr;
  std::vector<Generator> vec_R2{Generator::R2};
  Word word_R2(vec_R2, mat_R2, IsomClass::Reflection,
               arma::trace(mat_R2), (int)ref_order);

  std::vector<Generator> vec_E2{Generator::E2};
  Word word_E2(vec_E2, arma::inv(mat_R2), IsomClass::Reflection,
               arma::trace(arma::inv(mat_R2)), (int)ref_order);

  CompMat3 mat_R3(3, 3);
  //mat_R3 << std::conj(psi) << 0.0                     << 0.0       << arma::endr
  //       << 0.0            << std::conj(psi)          << 0.0       << arma::endr
  //       << tau            << -psi * std::conj(sigma) << psi * psi << arma::endr;
  mat_R3 << -1.0      << 0.0  << 0.0 << arma::endr
         << 0.0       << -1.0 << 0.0 << arma::endr
         << u*u*u*u*u << u*u  << 1.0 << arma::endr;
  std::vector<Generator> vec_R3{Generator::R3};
  Word word_R3(vec_R3, mat_R3, IsomClass::Reflection,
               arma::trace(mat_R3), (int)ref_order);

  std::vector<Generator> vec_E3{Generator::E3};
  Word word_E3(vec_E3, arma::inv(mat_R3), IsomClass::Reflection,
               arma::trace(arma::inv(mat_R3)), (int)ref_order);
  //////////////////////////////////////////////////////////////////////////////
  std::vector<Word> generators;
  generators.reserve(6);
  generators.push_back(word_R1);
  generators.push_back(word_R2);
  generators.push_back(word_R3);
  generators.push_back(word_E1);
  generators.push_back(word_E2);
  generators.push_back(word_E3);

  KnuthBendix kb;

  // add inverse pairs to knuth bendix class
  std::vector<Generator> R1E1 = {Generator::R1, Generator::E1};
  WordVector d_R1E1(R1E1);
  kb.add(Relation(d_R1E1));

  std::vector<Generator> R2E2 = {Generator::R2, Generator::E2};
  WordVector d_R2E2(R2E2);
  kb.add(Relation(d_R2E2));

  std::vector<Generator> R3E3 = {Generator::R3, Generator::E3};
  WordVector d_R3E3(R3E3);
  kb.add(Relation(d_R3E3));

  std::vector<Generator> E1R1 = {Generator::E1, Generator::R1};
  WordVector d_E1R1(E1R1);
  kb.add(Relation(d_E1R1));

  std::vector<Generator> E2R2 = {Generator::E2, Generator::R2};
  WordVector d_E2R2(E2R2);
  kb.add(Relation(d_E2R2));

  std::vector<Generator> E3R3 = {Generator::E3, Generator::R3};
  WordVector d_E3R3(E3R3);
  kb.add(Relation(d_E3R3));

  // std::cout << "Checking for braid relations in the group presentation\n";
  // get_braid_relns(generators, kb);
  // std::cout << "Found " << kb.size() << " braid relations\n";

  // auto start = std::chrono::high_resolution_clock::now();
  // std::cout << "finding all unique words upto length " << upto << "\n";
  // std::vector<Word> unique_words = get_words_upto_n(upto, generators, kb);
  // auto finish = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> elapsed = finish - start;
  // std::cout << "Found " << unique_words.size() << " words. Elapsed time: " << elapsed.count() << " s\n";
  // std::cout << "Looking for relevant relations in group presentation\n";
  // std::vector<Word> relevant_words = reduce_words(unique_words, generators);
  // std::cout << relevant_words.size() << "\n";

  // std::cout << "Removing conjugate words from group presentation\n";
  // relevant_words = kill_conjugates(relevant_words, generators);

  // std::cout << "Removing power words from group presentation\n";
  // relevant_words = is_power(relevant_words, 100, false);
  // print_word_vector(relevant_words, false); //print_c_and_p(relevant_words);

  // std::cout << "Running Knuth-bendix\n";
  // kb.run_algo();
  // std::cout << "Total relations found " << kb.size() <<"\n";
  // kb.print();

  // arma::cx_dvec base_vector(3);
  // base_vector << 1.0 << arma::endr
  //             << 1.0 << arma::endr
  //             << 1.0 << arma::endr;

  ComFunDomain f_dom(word_R1 * word_R2 * word_R3);
  Face start_face(word_R1, word_R2, word_R3);
  std::cout << f_dom.face_count() << std::endl;

  f_dom.add_face(start_face);
  std::cout << f_dom.face_count() << std::endl;

  f_dom.build_f_domain(100);

  return 0;
};
