#include "getwordsupton.h"

std::vector<Word> get_words_upto_n(const unsigned int n, const std::vector<Word>& gens,
                                   KnuthBendix& kb)
{
  size_t i = 0;
  size_t count = 0;
  std::vector<Word> seen_words;
  seen_words.reserve(1000000);
  seen_words.push_back(Word());

  //inverse pairs
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

  while (i < n)
  {
    std::vector<Word> seen_words_n_minus_1;
    seen_words_n_minus_1.reserve(1000000);
    for (size_t j = 0; j < seen_words.size(); ++j)
    {
      const Word test_word = seen_words[j];
      if (test_word.word_length() == i)
        seen_words_n_minus_1.push_back(test_word);
    }

    for (size_t j = 0; j < seen_words_n_minus_1.size(); ++j)
    {
      const Word base_word = seen_words_n_minus_1[j];
      const Generator base_last = base_word.last_element();

      for (size_t k = 0; k < gens.size(); ++k)
      {
        const Word new_gen = gens[k];
        // if the final element is the inverse of the generator we're about to
        // do nothing as they will cancel
        if (inverse(base_last) == new_gen.last_element())
          continue;

        // get a new word by multiplying the generator to the base word
        const Word new_word = base_word * new_gen;
        // is the new word the identity?
        // this is almost certainly not getting hit
        if (new_word.get_isom_class() == IsomClass::Identity)
        {
          kb.add(Relation(&new_word));
          if (VERBOSE) 
            std::cout << new_word.as_string() << " is the identity!\n";
        }
        // is the new word elliptic?
        if (new_word.get_order() > 0)
        {
          const Word pow_word = power(new_word, new_word.get_order());     
          kb.add(Relation(&pow_word));        
          if (VERBOSE)          
            std::cout << new_word.as_string() << " has order " << new_word.get_order() << "\n";
        }
        // check if we've seen this before!
        bool seen = false;
/*        const Word *wd = nullptr; // uninitialized pointer exciting!
        size_t l = 0;
        for (wd = seen_words.data(), l = 0; l < seen_words.size(); ++l, ++wd)
        {
        count += 1;
          if (new_word.is_equal(wd))// || new_word.is_equal_inverse(wd))
          {
            seen = true;
            kb.add(Relation(&new_word, wd));
            if (VERBOSE)
              std::cout << new_word.as_string() << " is equal to "
                        << wd->as_string() << "\n";
            break;
          }
        }*/

        size_t l = 0;
        for (unsigned int w = 0; w < seen_words.size(); ++w)
        {
          count += 1;
          const Word wd = seen_words[w];
          if (new_word.is_equal(wd))// || new_word.is_equal_inverse(wd))
          {
            seen = true;
            kb.add(Relation(&new_word, &wd));
            if (VERBOSE)
              std::cout << new_word.as_string() << " is equal to "
                        << wd.as_string() << "\n";
            break;
          }
        }

        if (!seen)
          seen_words.push_back(new_word);
      }

    }
    int max_pos = (1 - std::pow(6, i+2)) / (-5);
    std::cout << i << " | " << seen_words_n_minus_1.size() << " -> "
              << seen_words.size() << " / " << max_pos
              << " (" << 100.0*seen_words.size() / max_pos << "%)\n";
    i += 1;
  }
  std::cout << "Total unique words seen " << seen_words.size() <<"\n";
  std::cout << "Total relations found " << kb.size() <<"\n";
//  kb.print();
//  kb.run_algo();
//  std::cout << "Total relations found " << kb.size() <<"\n";
//  kb.print();
  std::cout << "Comparisons made " << count << "\n";

  return seen_words;
}
