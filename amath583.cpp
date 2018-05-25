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
#include "amath583.hpp"
#include "omp.h"

double ompTwoNorm(const Vector& x) {
  double norm = 0.0, sum = 0.0;
  int nthread = 4;
  #pragma omp parallel for num_threads(nthread) reduction(+:norm) shared(x)
  {
    for(size_t i = 0; i < x.num_rows(); ++i){
      sum += x(i)*x(i);
    }
  }
  norm = std::sqrt(sum);
  return norm;
}

//void ompMatvec(const COOMatrix& A, const Vector& x, Vector& y) { A.ompMatvec(x,y) }
//void ompMatvec(const CSRMatrix& A, const Vector& x, Vector& y) { A.ompMatvec(x,y) }

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
