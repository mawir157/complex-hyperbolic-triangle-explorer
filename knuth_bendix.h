#pragma once

#include "reduction.h"

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