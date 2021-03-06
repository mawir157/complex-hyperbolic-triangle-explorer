#pragma once

#include "word.h"
//#include "reduction.h"
#include "knuth_bendix.h"

std::vector<Word> get_words_upto_n(const unsigned int n, const std::vector<Word>& gens, 
	                               KnuthBendix& kb);
