#include "f_domain.h"

Face::Face(const Word& w_1, const Word& w_2, const Word& w_3,
           const bool face, const bool side) :
    m_word_side_a(w_1)
  , m_word_side_b(w_2)
  , m_word_side_c(w_3)
  , m_face_matched(face)
  , m_side_paired(side)
  {}

bool Face::operator<(const Face& f) const
{
  return (m_word_side_a < f.m_word_side_a) &&
         (m_word_side_b < f.m_word_side_b) &&
         (m_word_side_c < f.m_word_side_c);
}

ComFunDomain::ComFunDomain(const Word& centre) :
    m_123(centre)
  {}

bool ComFunDomain::seen_before(const Face& f) const
{
  return false;
}

bool ComFunDomain::add_face(const Face& f)
{
  // check if it is in the set
  const bool is_in = seen_before(f);
  if (is_in)
    return false;

  ms_faces.insert(f);
  return true;
}

bool ComFunDomain::abc_bca_rule(const Face& f) const
{
  const Word abc = f.m_word_side_a * f.m_word_side_b * f.m_word_side_c;
  if (m_123.is_equal(&abc))
    return true;

  const Word bca = f.m_word_side_b * f.m_word_side_c * f.m_word_side_a;
  if (m_123.is_equal(&bca))
    return true;

  return false;
}

bool ComFunDomain::match_2_face(const Face& f)
{
  // there a four transforms of F(x; y, z) to check
  // F(y; x, z)
  const Face temp_face_bac(f.m_word_side_b,
                           f.m_word_side_a,
                           f.m_word_side_c, false, true);
  if (abc_bca_rule(temp_face_bac)) // check it passese the abc/bca test...
    return add_face(temp_face_bac); // ... and try to add it to the list

  // F(y; z, x)
  const Face temp_face_bca(f.m_word_side_b,
                           f.m_word_side_c,
                           f.m_word_side_a, false, true);
  if (abc_bca_rule(temp_face_bca)) // check it passese the abc/bca test...
    return add_face(temp_face_bca); // ... and try to add it to the list

  // F(z; x, y)
  const Face temp_face_cab(f.m_word_side_c,
                           f.m_word_side_a,
                           f.m_word_side_b, false, true);
  if (abc_bca_rule(temp_face_cab)) // check it passese the abc/bca test...
    return add_face(temp_face_cab); // ... and try to add it to the list

  // F(z; y, x)
  const Face temp_face_cba(f.m_word_side_c,
                           f.m_word_side_b,
                           f.m_word_side_a, false, true);
  if (abc_bca_rule(temp_face_cba)) // check it passese the abc/bca test...
    return add_face(temp_face_cba); // ... and try to add it to the list

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
}

void ComFunDomain::build_f_domain()
{
  // check if it is in the set
  std::set<Face> prev_faces = ms_faces;

  std::set<Face>::iterator cur_face;
  bool changed = false;
  for (cur_face = prev_faces.begin(); cur_face != prev_faces.end(); ++cur_face)
  {
    // try to add the matched two-face to fdom
    if (!cur_face->m_face_matched)
      {
      const bool matched_2 = match_2_face(*cur_face);
      changed |= matched_2;
      if (matched_2)
        cur_face->m_face_matched = true;
    }

    if (!cur_face->m_side_paired)
    {
      const bool side_paired =  match_C_line(*cur_face);
      changed |= side_paired;
      if (side_paired)
        cur_face->m_side_paired = true;
    }
  }
  return;
}