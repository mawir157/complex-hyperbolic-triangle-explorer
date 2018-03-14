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

static const bool VERBOSE = false;
static const double PI = 3.14159265358979323;
static const double TOL = 1e-6;
static const comp_d omega(-1.0 / 2.0, std::sqrt(3.0) / 2.0);

bool is_id(const CompMat3& matrix);

enum IsomClass { Identity, Elliptic, Reflection, Parabolic, Loxodromic };
enum Generator { ID, R1, E1, R2, E2, R3, E3 };
Generator inverse(const Generator gen);

void print_e_structure(const CompMat3& matrix, const std::string Name, const double tol=TOL);
arma::cx_dvec get_neg_evec(const CompMat3& M, const CompMat3& H);
double goldman(const CompMat3& matrix);
IsomClass mat_iso_class(const CompMat3& matrix);
bool is_id(const CompMat3& matrix);
bool is_close(const CompMat3& A, const CompMat3& B);
int mat_order(const CompMat3& matrix, const int max_order=1000);
int mat_order_alt(const CompMat3& matrix, const int max_order=1000);
inline double jmt_abs_fast(comp_d z);
bool eq_trace(comp_d tr_1, comp_d tr_2);
inline double jmt_abs_fast(comp_d z);
inline bool near_val(const double d, const double v, const double tol=TOL);
inline double frac_part(const double d);
