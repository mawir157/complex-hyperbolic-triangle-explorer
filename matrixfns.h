#pragma once

#include <armadillo>
#include <complex>
#include <vector>
#include <iostream> 
#include <string>
#include <cmath>
#include <chrono>

typedef std::complex<double> comp_d;
typedef arma::Mat<comp_d> CompMat3;
typedef arma::cx_dvec point;

struct mat_sig
{
  mat_sig(unsigned int p, unsigned int n, unsigned int g);

  unsigned int m_pos;
  unsigned int m_nul;
  unsigned int m_neg;
};

static const bool VERBOSE = false;
static const double PI = 3.14159265358979323;
static const double TOL = 1e-6;
static const comp_d omega(-1.0 / 2.0, std::sqrt(3.0) / 2.0);

enum IsomClass { Identity, Elliptic, Reflection, Parabolic, Loxodromic };
enum Generator { ID, R1, E1, R2, E2, R3, E3 };

Generator inverse(const Generator gen);
//////////////////////////////////////////////////////////////////////////////
// 
// MATRIX FUNCTIONS
//
// Is a matrix the identity ?
bool is_id(const CompMat3& matrix);
// Are two matrices close ?
bool is_close(const CompMat3& A, const CompMat3& B);
// Get the isometry class of a matrix
IsomClass mat_iso_class(const CompMat3& matrix);
// Get the order of a matrix 
int mat_order(const CompMat3& matrix, const int max_order=1000);
int mat_order_alt(const CompMat3& matrix, const int max_order=1000);
// Print the eigen-structure of ta matrix (for debugging)
void print_e_structure(const CompMat3& matrix, const std::string& Name, const double tol=TOL);
mat_sig get_mat_sig(const CompMat3& matrix, const double tol=TOL);
// Returns the unique -ve eigen vector of a matrix wrt Hermitian form H
point get_neg_evec(const CompMat3& M, const CompMat3& H);
// The Goldman trace value of a matrix
double goldman(const CompMat3& matrix);


// various numerical nonsense function
inline double jmt_abs_fast(comp_d z);
bool eq_trace(const comp_d& tr_1, const comp_d& tr_2);
inline double jmt_abs_fast(comp_d z);
inline bool near_val(const double d, const double v, const double tol=TOL);
inline double frac_part(const double d);

point get_polar(const point& p1, const point& p2, const CompMat3& H);
point normalize_point(const point& p, const int index = -1);
inline comp_d h_n(const point& p1, const CompMat3& H, const point& p2)
{
  return (arma::trans(p1) * H * p2).eval()(0, 0);
}