#pragma once

#include <algorithm>
#include "word.h"

class WordDeq
{
  public:
    // constructors
    WordDeq();
    WordDeq(const std::vector<Generator>& word_vec);
    WordDeq(const Word* wd);

    // overloaded operators
    bool operator<(const WordDeq& wv) const;
    bool operator>(const WordDeq& wv) const;
    bool operator==(const WordDeq& wv) const;
    
    // for use in Knuth-Bendix
    int  overlap(const WordDeq& wd) const;
    int  contains(const WordDeq& w) const;
    bool reduce(const WordDeq& new_rel, size_t old_len, size_t at);

    size_t size() const { return m_word_deq.size(); }
    std::vector<Generator> get_vector() const { return m_word_deq; }
    std::string as_string() const;

  private:
    std::vector<Generator> m_word_deq;
};

class Relation
{
  public:
    // constructor
    Relation();
    Relation(const Word* wd);
    Relation(const WordDeq& wd);
    Relation(const Word* wd1, const Word* wd2);
    Relation(const WordDeq& wd1, const WordDeq& wd2);
    // overloaded operators
    bool operator==(const Relation& r) const;
    bool operator<(const Relation& r) const;

    void print() const;   
      
    bool reduce(const Relation& rd, bool is_left=true);
    bool is_redundant(const Relation& r) const;
    bool is_trivial() const { return lhs == rhs;}
    int  get_overlap(const Relation& rhs) const;
    Relation k_b(const Relation& r, size_t overlap) const;

  // we must have lhs > rhs (hence lhs -> rhs is a Relation)
  private:
    WordDeq lhs;
    WordDeq rhs;
};

class KnuthBendix
{
  public:
    bool   add(Relation r);
    bool   run_algo();
    bool   apply_reductions();
    bool   remove_repeated();
    bool   remove_trivial();
    bool   remove_redundant();
    void   order();
    size_t size() const { return m_relations.size(); }
    void   print() const;

  private:
    std::vector<Relation> m_relations;
};
