//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#include <cmath>
#include <random>
#include <functional>
#include <thread>
#include "amath583.hpp"
#include "omp.h"

double ompTwoNorm(const Vector& x) {
  double norm = 0.0, sum = 0.0;
  size_t i;
  size_t num_threads = omp_get_max_threads();
  size_t schedule_length = x.num_rows() / num_threads;
  #pragma omp parallel for num_threads(num_threads) schedule(dynamic, schedule_length) default(none) reduction(+:sum) private(i) shared(x)
  {
    for( i = 0; i < x.num_rows(); ++i){
      sum += x(i)*x(i);
    }
  }
  norm = sqrt(sum);
  return norm;
}

void ompMatvec(const COOMatrix& A, const Vector& x, Vector& y) { A.ompMatvec(x,y); }
void ompMatvec(const CSRMatrix& A, const Vector& x, Vector& y) { A.ompMatvec(x,y); }

void matvec(const COOMatrix& A, const Vector& x, Vector& y) { A.matvec(x, y); }
void matvec(const CSRMatrix& A, const Vector& x, Vector& y) { A.matvec(x, y); }

void randomize(Vector& x) {
  static std::default_random_engine             generator;
  static std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  static auto                                   dice = std::bind(distribution, generator);

  for (int i = 0; i < x.num_rows(); ++i) {
    x(i) = dice();
  }
}

void zeroize(Vector& x) {
  for (size_t i = 0; i < x.num_rows(); ++i) {
    x(i) = 0;
  }
}

double twoNorm(const Vector& x) {
  double sum = 0.0;
  for (size_t i = 0; i < x.num_rows(); ++i) {
    sum += x(i) * x(i);
  }
  return sqrt(sum);
}

void piscetize(COOMatrix& A, size_t xpoints, size_t ypoints) {
  assert(A.num_rows() == A.num_cols());
  assert(xpoints * ypoints == A.num_rows());

  A.clear();

  for (size_t j = 0; j < xpoints; j++) {
    for (size_t k = 0; k < ypoints; k++) {
      size_t jrow = j * ypoints + k;

      if (j != 0) {
        size_t jcol = (j - 1) * ypoints + k;
        A.push_back(jrow, jcol, -1.0);
      }
      if (k != 0) {
        size_t jcol = j * ypoints + (k - 1);
        A.push_back(jrow, jcol, -1.0);
      }

      A.push_back(jrow, jrow, 4.0);

      if (k != ypoints - 1) {
        size_t jcol = j * ypoints + (k + 1);
        A.push_back(jrow, jcol, -1.0);
      }
      if (j != xpoints - 1) {
        size_t jcol = (j + 1) * ypoints + k;
        A.push_back(jrow, jcol, -1.0);
      }
    }
  }
}

void piscetize(CSRMatrix& A, size_t xpoints, size_t ypoints) {
  assert(A.num_rows() == A.num_cols());
  assert(xpoints * ypoints == A.num_rows());

  A.clear();
  A.openForPushBack();

  for (size_t j = 0; j < xpoints; j++) {
    for (size_t k = 0; k < ypoints; k++) {
      size_t jrow = j * ypoints + k;

      if (j != 0) {
        size_t jcol = (j - 1) * ypoints + k;
        A.push_back(jrow, jcol, -1.0);
      }
      if (k != 0) {
        size_t jcol = j * ypoints + (k - 1);
        A.push_back(jrow, jcol, -1.0);
      }

      A.push_back(jrow, jrow, 4.0);

      if (k != ypoints - 1) {
        size_t jcol = j * ypoints + (k + 1);
        A.push_back(jrow, jcol, -1.0);
      }
      if (j != xpoints - 1) {
        size_t jcol = (j + 1) * ypoints + k;
        A.push_back(jrow, jcol, -1.0);
      }
    }
  }
  A.closeForPushBack();
}

/*  
 * Driver helper function 
 * user_dim : user given matrix dimension
 * isCSR : when set to true will run for CSR matrix, otherwise it will run for COO
*/
void driver_helper(size_t user_dim, bool isCSR)
{

  unsigned long dim = user_dim * user_dim;
  double seq_time, par_time, seq_norm, par_norm;
  size_t num_cores = std::thread::hardware_concurrency();
  size_t thread_count = 0;
  
  #pragma omp parallel reduction(+:thread_count)
  thread_count += 1;

  Vector x(dim), y_seq(dim), y_par(dim);
  randomize(x);

  size_t times_seq = 100000 / user_dim < 100 ? 100 : 10000 / user_dim  ;
  size_t times_omp = times_seq * 4;

  if (isCSR){
    //std::cout << "Running CSR\n";
    CSRMatrix A(dim, dim);
    piscetize(A, user_dim, user_dim);

    /* Sequential */
    seq_time = timed_matvec(A, x, y_seq, times_seq, false);
    
    /* Parallel */
    par_time = timed_matvec(A, x, y_par, times_omp, true);
  }
  else
  {
    //std::cout << "Running COO\n";
    COOMatrix A(dim, dim);
    piscetize(A, user_dim, user_dim);

    /* Sequential */
    seq_time = timed_matvec(A, x, y_seq, times_seq, false);

    /* Parallel */
    par_time = timed_matvec(A, x, y_par, times_omp, true);
  }

  seq_norm = twoNorm(y_seq);
  par_norm = ompTwoNorm(y_par);  

  /* output */
  double speed_up = 0;
  if (par_time > 0) speed_up = seq_time / par_time;
  double norm_diff = std::abs(seq_norm - par_norm);
  
  std::cout << std::setprecision(15) << num_cores << "\t" << thread_count << "\t" << dim << "\t";
  std::cout << seq_time << "\t" << par_time << "\t" << speed_up << "\t" << norm_diff   << std ::endl;
  //std::cout << std::setprecision(30) << seq_norm << "\t" << par_norm << std::endl;
}