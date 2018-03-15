#pragma once

#include <algorithm>
#include "word.h"

class WordVector
{
  public:
    // constructors
    WordVector();
    WordVector(const std::vector<Generator>& word_vec);
    WordVector(const Word* wd);

    // overloaded operators
    bool operator<(const WordVector& wv) const;
    bool operator>(const WordVector& wv) const;
    bool operator==(const WordVector& wv) const;
    
    // for use in Knuth-Bendix
    int  overlap(const WordVector& wd) const;
    int  contains(const WordVector& w) const;
    bool reduce(const WordVector& new_rel, size_t old_len, size_t at);

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
    Relation(const WordVector& wd);
    Relation(const Word* wd1, const Word* wd2);
    Relation(const WordVector& wd1, const WordVector& wd2);
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
    WordVector lhs;
    WordVector rhs;
};
