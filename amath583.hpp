//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2018
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine
//

#ifndef __AMATH583_HPP
#define __AMATH583_HPP

#include "Vector.hpp"
#include "COOMatrix.hpp"
#include "CSRMatrix.hpp"

double ompTwoNorm(const Vector& x);
double twoNorm(const Vector& x);
void randomize(Vector& x);
void matvec(const COOMatrix& A, const Vector& x, Vector& y);
void ompMatvec(const COOMatrix& A, const Vector& x, Vector& y);
void matvec(const CSRMatrix& A, const Vector& x, Vector& y);
void ompMatvec(const CSRMatrix& A, const Vector& x, Vector& y);

void piscetize(COOMatrix& A, size_t xpoints, size_t ypoints);
void piscetize(CSRMatrix& A, size_t xpoints, size_t ypoints);

void zeroize(Vector& x);

#endif // __AMATH583_HPP
