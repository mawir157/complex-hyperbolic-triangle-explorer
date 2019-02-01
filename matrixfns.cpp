#include "matrixfns.h"

mat_sig::mat_sig(unsigned int p, unsigned int n, unsigned int g) :
    m_pos(p)
  , m_nul(n)
  , m_neg(g) {}

void mat_sig::print_sig() const
{
  std::cout << "(" << m_pos << "," << m_neg << ")" <<std::endl;
}

void print_e_structure(const CompMat3& matrix, const std::string& Name, const double tol)
{
  arma::cx_vec eigvals = arma::eig_gen( matrix );
  std::cout << "eigenvalues of " << Name << " = (" 
            << std::real(eigvals[0]) << ", " 
            << std::real(eigvals[1]) << ", "
            << std::real(eigvals[2]) << ")\n";
  unsigned int pos_evals = 0;
  unsigned int neg_evals = 0;
  for (size_t i = 0; i < eigvals.size(); ++i)
  {
    if (std::abs(std::real(eigvals[i])) > TOL)
    {
      if (std::real(eigvals[i]) > 0.0)
        pos_evals += 1;
      else
        neg_evals += 1;
    }
  }
  std::cout << Name << " has signature (" << pos_evals << ", " << neg_evals << ")\n";
  std::cout << Name << " has Goldman trace formula " << goldman(matrix) << "\n";
}

mat_sig get_mat_sig(const CompMat3& matrix, const double tol)
{
  arma::cx_vec eigvals = arma::eig_gen( matrix );
  unsigned int pos_evals = 0;
  unsigned int nul_evals = 0;
  unsigned int neg_evals = 0;
  for (size_t i = 0; i < eigvals.size(); ++i)
  {
    if (std::abs(std::real(eigvals[i])) > TOL)
    {
      if (std::real(eigvals[i]) > 0.0)
        ++pos_evals;
      else
        ++neg_evals;
    }
    else
      ++nul_evals;
  }
  mat_sig ret(pos_evals, nul_evals, neg_evals);
  return ret;
}

point get_neg_evec(const CompMat3& M, const CompMat3& H)
{
  arma::cx_vec eigval;
  arma::cx_mat eigvec;

  arma::eig_gen(eigval, eigvec, M);

  CompMat3 D = arma::trans(eigvec) * H * eigvec;

  size_t neg_index = 0;
  bool seen_neg = false;
  for (size_t i = 0; i < 3; ++i)
  {
    double norm = std::real(D(i, i));
    if (norm < 0)
    {
      if (seen_neg) { std::cout << "WARNING: There are two negative eigen-vectors!\n"; }

      neg_index = i;
      seen_neg = true;
    }
  }

  if (!seen_neg) { std::cout << "WARNING: There are no negative eigen-vectors\n"; }

  return eigvec.col(neg_index);
}

double goldman(const CompMat3& matrix)
{
  comp_d trace = arma::trace(matrix);
  double ab_tr_sq = std::abs(trace) * std::abs(trace);
  return ab_tr_sq * ab_tr_sq - 8.0 * std::real(trace * trace * trace) + 
         18.0 * ab_tr_sq - 27.0;
}

IsomClass mat_iso_class(const CompMat3& matrix)
{
  if (is_id(matrix))
    return IsomClass::Identity;

  double g_man = goldman(matrix);

  // if g_man is very small we're either parabolic or a reflection
  if (std::abs(g_man) < TOL)
  {
    arma::cx_vec eigvals = arma::eig_gen( matrix );
    size_t mod_1_count = 0;
    for (size_t i = 0; i < eigvals.size(); ++i)
      if ((std::abs(eigvals[i]) - 1.0) < TOL)
        mod_1_count += 1;

//    std::cout << mod_1_count << "\n";
    if (mod_1_count == 3)
    {
//      std::cout << "E(" 
//                << std::abs(std::real(eigvals[0])) << ", " 
//                << std::abs(std::real(eigvals[1])) << ", "
//                << std::abs(std::real(eigvals[2])) << ")\n";
      return IsomClass::Reflection;
    }
    else
    {
//      std::cout << "P(" 
//                << std::abs(std::real(eigvals[0])) << ", " 
//                << std::abs(std::real(eigvals[1])) << ", "
//                << std::abs(std::real(eigvals[2])) << ")\n";
      return IsomClass::Parabolic;
    }
  }

  if (g_man > 0)
    return IsomClass::Loxodromic;

  return IsomClass::Elliptic;
}

bool is_id(const CompMat3& matrix)
{
  CompMat3 id = arma::eye<CompMat3>(3,3);
  if (arma::approx_equal(matrix, id, "absdiff", TOL))
    return true;

  if (arma::approx_equal(matrix, omega * id, "absdiff", TOL))
    return true;

  if (arma::approx_equal(matrix, omega * omega * id, "absdiff", TOL))
    return true;

  return false;
}

bool is_close(const CompMat3& A, const CompMat3& B)
{
  if (arma::approx_equal(A, B, "absdiff", TOL))
    return true;

  if (arma::approx_equal(A, omega * B, "absdiff", TOL))
    return true;

  if (arma::approx_equal(A, omega * omega * B, "absdiff", TOL))
    return true;

  return false;  
}

int mat_order(const CompMat3& matrix, const int max_order)
{
  CompMat3 temp = matrix;
  for (int i = 1; i < max_order; ++i)
  {
    if (is_id(temp))
      return i;

    temp *= matrix;
  }

  return -2;
}

int mat_order_alt(const CompMat3& matrix, const int max_order)
{
  arma::cx_vec eigval = arma::eig_gen( matrix );
  const double eval_1 = (std::arg(eigval[0]) + 2.0 * PI) / (2.0 * PI);
  const double eval_2 = (std::arg(eigval[1]) + 2.0 * PI) / (2.0 * PI);
  const double eval_3 = (std::arg(eigval[2]) + 2.0 * PI) / (2.0 * PI);

  for (int i = 1; i < max_order; ++i)
  {
    double t1 = frac_part(i * eval_1);
    double t2 = frac_part(i * eval_2);
    double t3 = frac_part(i * eval_3);

    if ((near_val(t1, 0.0)       && near_val(t2, 0.0)       && near_val(t3, 0.0)) ||
        (near_val(t1, 1.0 / 3.0) && near_val(t2, 1.0 / 3.0) && near_val(t3, 1.0 / 3.0)) ||
        (near_val(t1, 2.0 / 3.0) && near_val(t2, 2.0 / 3.0) && near_val(t3, 2.0 / 3.0)))
    {
      return i;
    }
  }
  return -1;
}

// we use this as an alternative to std::abs as it is faster, due to no sqrt
inline double jmt_abs_fast(comp_d z)
{
  return std::abs(std::max(std::real(z), std::imag(z)));
}

inline bool near_val(const double d, const double v, const double tol)
{
  return (std::abs(d - v) < tol);
}

inline double frac_part(const double d)
{
  return std::abs(d - std::round(d));
}

bool eq_trace(const comp_d& tr_1, const comp_d& tr_2)
{

  if (jmt_abs_fast(tr_1 - tr_2) < TOL)
    return true;

  if (jmt_abs_fast(tr_1 - omega * tr_2) < TOL)
    return true;

  if (jmt_abs_fast(tr_1 - omega * omega * tr_2) < TOL)
    return true;

  return false;
}

Generator inverse(const Generator gen)
{
    switch (gen)
    {
      case Generator::R1: return Generator::E1; break;
      case Generator::R2: return Generator::E2; break;
      case Generator::R3: return Generator::E3; break;
      case Generator::E1: return Generator::R1; break;
      case Generator::E2: return Generator::R2; break;
      case Generator::E3: return Generator::R3; break;
      case Generator::ID: return Generator::ID; break;
    }
}

point get_polar(const point& p1, const point& p2, const CompMat3& H)
{
  point polar_vector;
  return polar_vector;
}

point normalize_point(const point& p, const int index)
{
  size_t size_index = 0;
  if ((index <  0) || (index > p.size()))
    size_index = 0;
  else
    size_index = (size_t)index;

  point norm_point(3);
  for (size_t i = 0; i < norm_point.size(); ++i)
    norm_point[i] = p[i] / p[size_index];

  return norm_point;
}
