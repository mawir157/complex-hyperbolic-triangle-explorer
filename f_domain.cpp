#include "f_domain.h"

Face::Face(const Word& w_1, const Word& w_2, const Word& w_3,
           const bool face, const bool side) :
    m_word_side_a(w_1)
  , m_word_side_b(w_2)
  , m_word_side_c(w_3)
  , m_face_matched(face)
  , m_side_paired(side)
  , m_braid_order(get_braid_relation(w_2, w_3))
{}

bool Face::operator<(const Face& f) const
{
  return (m_word_side_a < f.m_word_side_a) &&
         (m_word_side_b < f.m_word_side_b) &&
         (m_word_side_c < f.m_word_side_c);
}

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
  {}

bool ComFunDomain::seen_before(const Face& f) const
{
  for (unsigned int i = 0; i < ms_faces.size(); ++i)
  {
    Face g = ms_faces[i];
    bool equal = true;
    equal &= g.m_word_side_a.is_equal(f.m_word_side_a);
    equal &= g.m_word_side_b.is_equal(f.m_word_side_b);
    equal &= g.m_word_side_c.is_equal(f.m_word_side_c);
    
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

bool ComFunDomain::abc_bca_rule(const Face& f) const
{
  // Here we also check the faces braid
  if (f.m_braid_order < 0)
  {
    std::cout << "Found an non-braiding face!" << std::endl;
    return false;
  }

  const Word abc = f.m_word_side_a * f.m_word_side_b * f.m_word_side_c;
  if (m_123.is_equal(abc))
    return true;

  const Word bca = f.m_word_side_b * f.m_word_side_c * f.m_word_side_a;
  if (m_123.is_equal(bca))
    return true;

  return false;
}

bool ComFunDomain::match_2_face(const Face& f)
{
  // there a four transforms of F(x; y, z) to check
  // F(y; x, z)
  const Face temp_face_bac(f.m_word_side_b,
                           f.m_word_side_a,
                           f.m_word_side_c, true, false);
  if (abc_bca_rule(temp_face_bac)) // check it passese the abc/bca test...
    return add_face(temp_face_bac); // ... and try to add it to the list

  // F(y; z, x)
  const Face temp_face_bca(f.m_word_side_b,
                           f.m_word_side_c,
                           f.m_word_side_a, true, false);
  if (abc_bca_rule(temp_face_bca)) // check it passese the abc/bca test...
    return add_face(temp_face_bca); // ... and try to add it to the list

  // F(z; x, y)
  const Face temp_face_cab(f.m_word_side_c,
                           f.m_word_side_a,
                           f.m_word_side_b, true, false);
  if (abc_bca_rule(temp_face_cab)) // check it passese the abc/bca test...
    return add_face(temp_face_cab); // ... and try to add it to the list

  // F(z; y, x)
  const Face temp_face_cba(f.m_word_side_c,
                           f.m_word_side_b,
                           f.m_word_side_a, true, false);
  if (abc_bca_rule(temp_face_cba)) // check it passese the abc/bca test...
    return add_face(temp_face_cba); // ... and try to add it to the list

  std::cout << "No valid matched 2-face" << std::endl;
  return false;
}

bool ComFunDomain::match_C_line(const Face& f)
{
  // there are two options for F(a,;b,c)
  // F(a; a^-1ba, a^-1ca)
  const Face temp_face_1(f.m_word_side_a,
                         conjugate(f.m_word_side_b, f.m_word_side_a.invert()),
                         conjugate(f.m_word_side_c, f.m_word_side_a.invert()),
                         false, true);
  if (abc_bca_rule(temp_face_1)) // check it passese the abc/bca test...
    return add_face(temp_face_1); // ... and try to add it to the list

  // F(a; aba^-1, aca^-1)
  const Face temp_face_2(f.m_word_side_a,
                         conjugate(f.m_word_side_b, f.m_word_side_a),
                         conjugate(f.m_word_side_c, f.m_word_side_a),
                         false, true);
  if (abc_bca_rule(temp_face_2)) // check it passese the abc/bca test...
    return add_face(temp_face_2); // ... and try to add it to the list

  std::cout << "No valid matched C-line" << std::endl;
  return false;
}

void ComFunDomain::build_f_domain(const unsigned int max)
{
  bool changed = true;
  unsigned int i = 0;
  while (changed)
  {
    changed = false;
    //std::vector<Face> prev_faces = ms_faces;
    //std::vector<Face>::iterator cur_face;
    //for (cur_face = prev_faces.begin(); cur_face != prev_faces.end(); ++cur_face)
    size_t upto = ms_faces.size();
    for (size_t i = 0; i < upto; ++i)  
    {
      // try to add the matched two-face to fdom
      //if (!cur_face->m_face_matched)
      if (!ms_faces[i].m_face_matched)
      {
        //const bool matched_2 = match_2_face(*cur_face);
        const bool matched_2 = match_2_face(ms_faces[i]);
        changed |= matched_2;
        if (matched_2)
        {
//          std::cout << ms_faces[i].as_string() << std::endl;
          ms_faces[i].m_face_matched = true;
//          std::cout << ms_faces[i].as_string() << std::endl;
        }
      }

      //if (!cur_face->m_side_paired)
      if (!ms_faces[i].m_side_paired)
      {
        //const bool side_paired =  match_C_line(*cur_face);
        const bool side_paired =  match_C_line(ms_faces[i]);
        changed |= side_paired;
        if (side_paired)
          ms_faces[i].m_side_paired = true;
      }
    }

    std::cout << (changed ? "CONTINUE" : "STOP") << " - (" 
              << ms_faces.size() << ")" << std::endl;

    if (i > max)
      break;
    ++i;
  }
  return;
}