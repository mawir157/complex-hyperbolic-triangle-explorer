#include "getbraidrelns.h"

Pair::Pair(std::vector<Generator> l, std::vector<Generator> r) : 
    LHS(l)
  , RHS(r) {}

std::vector<Pair> get_braid_relns(const std::vector<Word>& gens, KnuthBendix& kb,
                                  unsigned int max_braid)
{
  std::vector<Pair> pair_vec;
  for (size_t i = 0; i < gens.size(); ++i)
  {
    const Word gen_a = gens[i];
    for (size_t j = i+1; j < gens.size(); ++j)
    {
      const Word gen_b = gens[j];
      Word LHS = gen_a;
      Word RHS = gen_b;
      if (gen_b.is_equal(gen_a) || gen_b.is_equal_inverse(gen_a))
        continue;
      
      bool alt = true;
      unsigned int count = 0;
      while (!LHS.is_equal(RHS) && count <  max_braid)
      {
        if (alt)
        {
          LHS = LHS * gen_b;
          RHS = RHS * gen_a;
        }
        else
        {
          LHS = LHS * gen_a;
          RHS = RHS * gen_b;
        }
        alt = !alt;
        count += 1;
      }
      if ((count >= max_braid) && (VERBOSE))
        std::cout << "No braid relation between " << gen_a.as_string() 
                  << " and " << gen_b.as_string() << "\n";
      else
      {
        if (VERBOSE)
          std::cout << LHS.as_string() << " = " << RHS.as_string() << "\n";
        Pair pr(LHS.get_gen_vec(), RHS.get_gen_vec());
        pair_vec.push_back(pr);
        Relation braid_rel(&LHS, &RHS);
        kb.add(braid_rel);
      }
    }  
  }
  return pair_vec;
}

int get_braid_relation(const Word& a, const Word& b,
                       const unsigned int max_braid,
                       const unsigned int min_braid)
{
  Word lhs = a;
  Word rhs = b;
  int braid = 1;
  bool parity = true;

  while (braid < max_braid)
  {
    if ((lhs.is_equal(rhs)) && (braid >= min_braid))
      return braid;

    if (parity)
    {
      lhs = lhs * b;
      rhs = rhs * a;
      parity = false;
    }
    else
    {
      lhs = lhs * a;
      rhs = rhs * b;
      parity = true;
    }
    ++braid;
  }
  return - 1;
}