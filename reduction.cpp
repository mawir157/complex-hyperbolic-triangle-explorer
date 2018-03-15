#include "reduction.h"

Relation::Relation()
{
  lhs = WordVector();
  rhs = WordVector();
}

Relation::Relation(const WordVector& w1, const WordVector& w2)
{
  if (w1 == w2)
    return;  //without intialising anything?

  if (w1 < w2)
  {
    lhs = w2; rhs = w1;
  }
  else
  {
    lhs = w1; rhs = w2;
  }  
}

Relation::Relation(const WordVector& w1)
{
  lhs = w1; 
}

Relation::Relation(const Word* wd)
{
  WordVector w1(wd);
  WordVector w2;

  if (w1 == w2)
    return;  //without intialising anything?

  if (w1 < w2)
  {
    lhs = w2; rhs = w1;
  }
  else
  {
    lhs = w1; rhs = w2;
  }
}

Relation::Relation(const Word* wd1, const Word* wd2)
{
  WordVector w1(wd1);
  WordVector w2(wd2);

  if (w1 == w2)
    return;  //without intialising anything?

  if (w1 < w2)
  {
    lhs = w2; rhs = w1;
  }
  else
  {
    lhs = w1; rhs = w2;
  }
}

bool Relation::operator==(const Relation& r) const
{
  return (lhs == r.lhs) && (rhs == r.rhs);
}

bool Relation::operator<(const Relation& r) const
{
  return lhs < r.lhs;
}

void Relation::print() const
{
  std::cout << lhs.as_string() << " -> " << rhs.as_string() << "\n" ;
}

bool Relation::reduce(const Relation& rd, bool is_left)
{
  if (is_left)
  {
    int index = lhs.contains(rd.lhs);
    if (index >= 0)
    {
//      std::cout << "+-+-+-+- left +-+-+-+-+\n";
//      print(); rd.print();
      lhs.reduce(rd.rhs, rd.lhs.size(), (size_t)index);
      if (lhs < rhs)
        std::swap(lhs, rhs);
//      print();
//      std::cout << "+-+-+-+-+-+-+-+-+-+-+-+\n";
      return true;
    }
  }
  else
  {
    int index = rhs.contains(rd.lhs);
    if (index >= 0)
    {
//      std::cout << "+-+-+- right -+-+-+-+-+\n";
//      print(); rd.print();
      rhs.reduce(rd.rhs, rd.lhs.size(), (size_t)index);
      if (lhs < rhs)
        std::swap(lhs, rhs);
//      print(); 
//      std::cout << "+-+-+-+-+-+-+-+-+-+-+-+\n";
      return true;
    }    
  }

  return false;
}

bool Relation::is_redundant(const Relation& r) const
{
  int index = lhs.contains(r.lhs);
  return (index >= 0);
}

int Relation::get_overlap(const Relation& r) const
{
  return lhs.overlap(r.lhs);
}

Relation Relation::k_b(const Relation& r, size_t overlap) const
{
  std::vector<Generator> L_1 = lhs.get_vector();
  std::vector<Generator> R_1 = rhs.get_vector();

  std::vector<Generator> L_2 = r.lhs.get_vector();
  std::vector<Generator> R_2 = r.rhs.get_vector();

  size_t B_length = L_1.size() - overlap;

//  R1+C <----> A+R2
  std::vector<Generator> new_lhs;
  for (size_t i = 0; i < R_1.size(); ++i)
    new_lhs.push_back(R_1[i]);
  for (size_t i = B_length; i < L_2.size(); ++i)
    new_lhs.push_back(L_2[i]);

  std::vector<Generator> new_rhs;
  for (size_t i = 0; i < overlap; ++i)
    new_rhs.push_back(L_1[i]);
  for (size_t i = 0; i < R_2.size(); ++i)
    new_rhs.push_back(R_2[i]);

  Relation new_relation(new_lhs, new_rhs);

  return new_relation;
}

////////////////////////////////////////////////////////////////////////////////
WordVector::WordVector()
{
  std::vector<Generator> id_deq;
  m_word_deq = id_deq;
}

WordVector::WordVector(const std::vector<Generator>& word_vec) :
  m_word_deq(word_vec) {}

WordVector::WordVector(const Word* wd)
{
  m_word_deq =  wd->get_gen_vec(); 
}

bool WordVector::operator<(const WordVector& w) const
{
  if (m_word_deq.size() < w.size())
    return true;

  if (m_word_deq.size() > w.size())
    return false;
  
  // at this point the two words are the same length
  std::vector<Generator> comp_deq = w.get_vector();
  for (size_t i = 0; i < m_word_deq.size(); ++i)
  {
    if (m_word_deq[i] < comp_deq[i])
      return true;
    else if (m_word_deq[i] > comp_deq[i])
      return false;
  }
  // if we get here then the words are equal
  return false;
}


bool WordVector::operator>(const WordVector& w) const
{
  if (m_word_deq.size() > w.size())
    return true;

  if (m_word_deq.size() < w.size())
    return false;
  
  // at this point the two words are the same length
  std::vector<Generator> comp_deq = w.get_vector();
  for (size_t i = 0; i < m_word_deq.size(); ++i)
  {
    if (m_word_deq[i] > comp_deq[i])
      return true;
    else if (m_word_deq[i] < comp_deq[i])
      return false;
  }
  // if we get here then the words are equal
  return false;
}

bool WordVector::operator==(const WordVector& w) const
{
  if (m_word_deq.size() != w.size())
    return false;
  
  // at this point the two words are the same length
  std::vector<Generator> comp_deq = w.get_vector();
  for (size_t i = 0; i < m_word_deq.size(); ++i)
  {
    if (m_word_deq[i] != comp_deq[i])
      return false;
  }
  // if we get here then the words are equal
  return true;
}
// template this out TODO
std::string WordVector::as_string() const
{
  if (m_word_deq.size() == 0) 
    return "e";

  std::string word_str;
  for (size_t i = 0; i < m_word_deq.size(); ++i)
  {
    Generator gen = m_word_deq[i];
    switch (gen)
    {
      case Generator::R1: word_str.append("R1"); break;
      case Generator::R2: word_str.append("R2"); break;
      case Generator::R3: word_str.append("R3"); break;
      case Generator::E1: word_str.append("E1"); break;
      case Generator::E2: word_str.append("E2"); break;
      case Generator::E3: word_str.append("E3"); break;
      case Generator::ID:
                 default: word_str.append("[]");
    }
  }
  return word_str;
}

bool WordVector::reduce(const WordVector& new_rel, size_t old_len, size_t at)
{
  // remove the old word    
  m_word_deq.erase(m_word_deq.begin() + at, 
                   m_word_deq.begin() + at + old_len) ;
  // and insert the new word
  size_t new_len = new_rel.size();
  std::vector<Generator> new_vec = new_rel.get_vector();
  if (new_len > 0)
    m_word_deq.insert(m_word_deq.begin() + at,
                      new_vec.begin(),
                      new_vec.end());  
  return true;
}

int WordVector::contains(const WordVector& w) const
{
  // if w is the identity crash!
  if (w.size() == 0)
    return -1;  

  // if w is longer than this bail out  
  if (w.size() > this->size())
    return -1;

  std::vector<Generator> gen_vec = w.get_vector();
  size_t i_small = 0;
  for (size_t i_large = 0; i_large < this->size(); ++i_large)
  {
    if (gen_vec[i_small] == m_word_deq[i_large])
      ++i_small;
    else
      i_small = 0;

    if (i_small == w.size())
      return i_large - i_small + 1;
  }
  return -3;
}

int WordVector::overlap(const WordVector& w) const
{
  if ((w.size() == 0) | (this->size() == 0))
    return -1;

  std::vector<Generator> left_vec = m_word_deq;
  std::vector<Generator> right_vec = w.get_vector();

  size_t overlap = 0;
  size_t max_overlap = std::min(left_vec.size(), right_vec.size());
  size_t best_overlap = 0;
  size_t best_overlap_index = 0;

  for (overlap = 1; overlap < max_overlap; ++overlap)
  {
    size_t l_start =  left_vec.size() - overlap;
    unsigned int count = 0;
    size_t i;
    for (i = 0; i < overlap; ++i)
    {
        if (left_vec[l_start + i] == right_vec[i])
            count +=1;
        else
            break;
    }
    if (count == overlap)
    {
        best_overlap = overlap;
        best_overlap_index = l_start;
    }
  }
  return best_overlap_index;
}
////////////////////////////////////////////////////////////////////////////////

