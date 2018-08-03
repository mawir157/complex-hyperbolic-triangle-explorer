#include "knuth_bendix.h"

bool KnuthBendix::add(Relation r)
{
  for (size_t i = 0; i < m_relations.size(); ++i)
    if ((m_relations[i] == r) || 
        r.is_redundant(m_relations[i]) ||
        r.is_trivial())
      return false;

  m_relations.push_back(r);
  return true;
}

void KnuthBendix::print() const
{
  for (size_t i = 0; i < m_relations.size(); ++i)
     m_relations[i].print();
}

bool KnuthBendix::remove_repeated()
{
  bool removed = false;
  bool carry_on = true;
  while (carry_on)
  {
    carry_on = false;
    for (size_t i = 0; i < m_relations.size(); ++i)
    {
      for (size_t j = i + 1; j < m_relations.size(); ++j)
      {
        if (m_relations[i] == m_relations[j])
        {
          carry_on = true;
          removed = true;
          m_relations.erase(m_relations.begin()+j);
          break;
        }
      }
    }
  }
  return removed;
}

bool KnuthBendix::remove_trivial()
{
  bool removed = false;
  bool carry_on = true;
  while (carry_on)
  {
    carry_on = false;
    for (size_t i = 0; i < m_relations.size(); ++i)
    {
      if (m_relations[i].is_trivial())
      {
        carry_on = true;
        removed = true;
        m_relations.erase(m_relations.begin()+i);
        break;
      }
    }
  }
  return removed;
}

bool KnuthBendix::remove_redundant()
{
  bool removed = false;
  bool carry_on = true;  
  while (carry_on)
  { 
    carry_on = false;
    for (size_t i = 0; i < m_relations.size(); ++i)
    {
      for (size_t j = 0; j < m_relations.size(); ++j)
      {
        if (i == j)
          continue;

        if (m_relations[i].is_redundant(m_relations[j]))
        {
          carry_on = true;
          removed = true;
          m_relations.erase(m_relations.begin()+i);
          break;          
        }
      }
      if (carry_on)
        break;
    }
  }
  return removed;
}

bool KnuthBendix::apply_reductions()
{
  bool carry_on = true;  
  bool removed = false;
  while (carry_on)
  { 
    carry_on = false;
    for (size_t i = 0; i < m_relations.size(); ++i)
    {
      for (size_t j = 0; j < m_relations.size(); ++j)
      {
        if (m_relations[i] == m_relations[j])
          continue;

//        carry_on |= m_relations[i].reduce(m_relations[j], true);
        if (m_relations[i].reduce(m_relations[j], false))
        {
          carry_on = true;
          removed = true;
          break;                
        }
      }
      if (carry_on)
        break;
    }
//    remove_redundant();
//    return false;
//    std::cout << size() << "\n";
//    std::cout << "******************************************\n";
  }
  return removed;
}

void KnuthBendix::order()
{
  std::sort(m_relations.begin(), m_relations.end());
}

bool KnuthBendix::run_algo()
{
  bool carry_on = true;
  size_t super_init = m_relations.size();
  while (carry_on)
  {
    carry_on = false;
    carry_on |= apply_reductions();
    //step 1 reduce all word where possible
    carry_on |= remove_redundant();
    // remove repeated relations
    carry_on |= remove_repeated();
    // remove A -> A
    carry_on |= remove_trivial();
    // put the realtions in smallest to largest order
    order();
  }
  std::cout << "done preprocessing Knuth-Bendix\n";


  carry_on = true;
  size_t count = 0;
  while (carry_on || count < 10)
  {
    carry_on = false;
//if (m_relations.size() > 200)
//  return true;
    size_t init_i_loop = m_relations.size();
    for (size_t i = 0; i < init_i_loop; ++i)
    {
//      std::cout << m_relations.size() << " | " << i << "\n";
      size_t init_j_loop = m_relations.size();
      for (size_t j = 0; j < init_j_loop; ++j)
      {
        int overlap = m_relations[i].get_overlap(m_relations[j]);
        if (overlap > 0)
        {
          Relation new_rel = m_relations[i].k_b(m_relations[j], (size_t)overlap);
          carry_on |= add(new_rel); // if (carry_on) { new_rel.print(); }

          if (carry_on)
            break;
        }
      }
      if (carry_on)
        break;
    } 
    carry_on |= apply_reductions();
    // remove relations where the lhs contains a the lhs of another
    carry_on |= remove_redundant();
    // remove repeated relations
    carry_on |= remove_repeated();
    // remove A -> A
    carry_on |= remove_trivial();
    ++count;

    if (carry_on)
    {
      std::cout << m_relations.size() << " (" << count << ") [" 
                << (int)m_relations.size() - (int)super_init << "]\n";
      order();
    }
  }
  return true;
}
