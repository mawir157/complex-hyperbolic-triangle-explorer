#pragma once

#include "word.h"
#include "reduction.h"

struct Pair 
{
  std::vector<Generator> LHS;
  std::vector<Generator> RHS;

  Pair(std::vector<Generator> l, std::vector<Generator> r);
};

std::vector<Pair> get_braid_relns(const std::vector<Word>& gens, KnuthBendix& kb,
	                              unsigned int max_braid=10);
