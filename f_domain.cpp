#include "f_domain.h"

Face::Face(const Word& w_1, const Word& w_2, const Word& w_3,
           const bool Face, const bool side) :
    m_word_side_a(w_1)
  , m_word_side_b(w_2)
  , m_word_side_c(w_3)
  , m_face_matched(Face)
  , m_side_paired(side)
  , m_braid_order(get_braid_relation(w_2, w_3))
  , m_word_bc(w_2 * w_3)
{}

std::string Face::as_string() const
{
  std::string word_str = "F(";
  word_str.append(m_word_side_a.as_string());
  word_str.append("; ");
  word_str.append(m_word_side_b.as_string());
  word_str.append(", ");
  word_str.append(m_word_side_c.as_string());
  word_str.append(")");
  if (m_face_matched)
    word_str.append("o");
  else
    word_str.append(" ");

  if (m_side_paired)
    word_str.append("o");
  else
    word_str.append(" ");

  return word_str;
}

ComFunDomain::ComFunDomain(const Word& centre) :
    m_123(centre)
  , kill_switch(false)
  {}

bool ComFunDomain::seen_before(const Face& f) const
{
  for (unsigned int i = 0; i < ms_faces.size(); ++i)
  {
    Face g = ms_faces[i];
    bool equal = true;
    equal &= g.base().is_equal(f.base());
    equal &= g.side_b().is_equal(f.side_b());
    equal &= g.side_c().is_equal(f.side_c());
    
    if (equal)
      return true;
  }

  return false;
}

bool ComFunDomain::add_face(const Face& f)
{
  // check if it is in the set
  const bool is_in = seen_before(f);
  if (is_in)
    return false;

  ms_faces.push_back(f);
  return true;
}

bool ComFunDomain::abc_bca_rule(const Face& f)
{
  // Here we also check the Face_2s braid
  if (f.braid_order() < 0)
  {
    std::cout << "Found a non-braiding Face!" << std::endl;
    kill_switch = true;
    return false;
  }

/*  if (f.m_word_bc.get_isom_class() != IsomClass::Elliptic)
  {
    std::cout << "Found non-eliptic vertex" << std::endl;
    kill_switch = true;
    return false;
  }*/

  const Word abc = f.base() * f.side_b() * f.side_c();
  if (m_123.is_equal(abc))
    return true;

  const Word bca = f.side_b() * f.side_c() * f.base();
  if (m_123.is_equal(bca))
    return true;

  return false;
}

bool ComFunDomain::match_2_face(const Face& f)
{
  // there a four transforms of F(x; y, z) to check
  // F(y; x, z)
  const Face temp_face_bac(f.side_b(),
                           f.base(),
                           f.side_c(), true, false);
  if (abc_bca_rule(temp_face_bac)) // check it passese the abc/bca test...
    return add_face(temp_face_bac); // ... and try to add it to the list

  // F(y; z, x)
  const Face temp_face_bca(f.side_b(),
                           f.side_c(),
                           f.base(), true, false);
  if (abc_bca_rule(temp_face_bca)) // check it passese the abc/bca test...
    return add_face(temp_face_bca); // ... and try to add it to the list

  // F(z; x, y)
  const Face temp_face_cab(f.side_c(),
                           f.base(),
                           f.side_b(), true, false);
  if (abc_bca_rule(temp_face_cab)) // check it passese the abc/bca test...
    return add_face(temp_face_cab); // ... and try to add it to the list

  // F(z; y, x)
  const Face temp_face_cba(f.side_c(),
                           f.side_b(),
                           f.base(), true, false);
  if (abc_bca_rule(temp_face_cba)) // check it passese the abc/bca test...
    return add_face(temp_face_cba); // ... and try to add it to the list

  if (kill_switch)
    return false;

  std::cout << "No valid matched 2-Face" << std::endl;
  return false;
}

bool ComFunDomain::match_C_line(const Face& f)
{
  // there are two options for F(a,;b,c)
  // F(a; a^-1ba, a^-1ca)
  const Face temp_face_1(f.base(),
                         conjugate(f.side_b(), f.base().invert()),
                         conjugate(f.side_c(), f.base().invert()),
                         false, true);
  if (abc_bca_rule(temp_face_1)) // check it passese the abc/bca test...
    return add_face(temp_face_1); // ... and try to add it to the list

  // F(a; aba^-1, aca^-1)
  const Face temp_face_2(f.base(),
                         conjugate(f.side_b(), f.base()),
                         conjugate(f.side_c(), f.base()),
                         false, true);
  if (abc_bca_rule(temp_face_2)) // check it passese the abc/bca test...
    return add_face(temp_face_2); // ... and try to add it to the list

  if (kill_switch)
    return false;

  std::cout << "No valid matched C-line" << std::endl;
  return false;
}

bool ComFunDomain::build_f_domain(const unsigned int max,
                                  const bool verbose)
{
  bool changed = true;
  unsigned int i = 0;
  while (changed)
  {
    changed = false;
    //std::vector<Face> prev_Face_2s = ms_faces;
    //std::vector<Face>::iterator cur_Face_2;
    //for (cur_Face_2 = prev_Face_2s.begin(); cur_Face_2 != prev_Face_2s.end(); ++cur_Face_2)
    size_t upto = ms_faces.size();
    for (size_t i = 0; i < upto; ++i)  
    {
      // try to add the matched two-Face to fdom
      if (!ms_faces[i].is_face_matched())
      {
        const bool matched_2 = match_2_face(ms_faces[i]);
        changed |= matched_2;
        if (matched_2)
        {
//          std::cout << ms_faces[i].as_string() << std::endl;
          ms_faces[i].set_face_matched(true);
//          std::cout << ms_faces[i].as_string() << std::endl;
        }
      }

      if (!ms_faces[i].is_side_paired())
      {
        const bool side_paired =  match_C_line(ms_faces[i]);
        changed |= side_paired;
        if (side_paired)
          ms_faces[i].set_side_paired(true);
      }
    }

    if (kill_switch)
      return false;

    if (verbose)
      std::cout << (changed ? "CONTINUE" : "STOP") << " - (" 
                << ms_faces.size() << ")" << std::endl;

    if (i > max)
    {
      std::cout << "WARNING: DID NOT TERMINATE" << std::endl;
      return false;
    }
    ++i;
  }
  return true;
}