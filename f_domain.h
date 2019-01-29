#pragma once

#include "matrixfns.h"
#include "word.h"

class Face
{
  public:
    Face(const Word& w_1, const Word& w_2, const Word& w_3,
         const bool face = false, const bool side = false);
    const Word m_word_side_a;
    const Word m_word_side_b;
    const Word m_word_side_c;
    bool m_face_matched;
    bool m_side_paired;

  private:




  public:
    bool operator<(const Face& f) const;
};

class ComFunDomain
{
  public:
    ComFunDomain(const Word& centre);
    void build_f_domain();
    bool add_face(const Face& f);

  private:
    const Word m_123;
    std::set<Face> ms_faces;
    bool seen_before(const Face& f) const;

    bool abc_bca_rule(const Face& f) const;

    bool match_2_face(const Face& f);
    bool match_C_line(const Face& f);
    bool populate_pyramid(const Face& f);
};

