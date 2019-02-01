#pragma once

#include "matrixfns.h"
#include "word.h"
#include "getbraidrelns.h"

class Face
{
  public:
    Face(const Word& w_1, const Word& w_2, const Word& w_3,
         const bool face = false, const bool side = false);
    std::string as_string() const;
    Word base() const { return m_word_side_a; }
    Word side_b() const { return m_word_side_b; }
    Word side_c() const { return m_word_side_c; }
    int braid_order() const { return m_braid_order; }
    bool is_face_matched() const { return m_face_matched; }
    void set_face_matched(const bool b) { m_face_matched = b; }
    bool is_side_paired() const { return m_side_paired; }
    void set_side_paired(const bool b) { m_side_paired = b; }

  private:
    const Word m_word_side_a;
    const Word m_word_side_b;
    const Word m_word_side_c;
    bool m_face_matched;
    bool m_side_paired;
    const int m_braid_order;
    const Word m_word_bc;
};

class Face_3
{
  public:
    std::vector<Face> m_v;
  private:
};

class ComFunDomain
{
  public:
    ComFunDomain(const Word& centre);
    bool build_f_domain(const unsigned int max = 25, const bool verbose = false);
    bool add_face(const Face& f);
    size_t face_count() const { return ms_faces.size(); }

  private:
    const Word m_123;
    bool kill_switch;
    std::vector<Face> ms_faces;
    bool seen_before(const Face& f) const;

    bool abc_bca_rule(const Face& f);

    bool match_2_face(const Face& f);
    bool match_C_line(const Face& f);
    bool populate_pyramid(const Face& f);
};
