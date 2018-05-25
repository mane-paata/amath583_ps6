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
#include "Timer.hpp"

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

/*  Templatized timed sequential matvec (COO and CSR) */
template <typename MatrixType>
double seq_matvec(const MatrixType& A, const Vector& x, Vector& y, size_t times){
	double time_elapsed = 0.0;
	Timer t;
	for(int j = 0; j < times; ++j){
		zeroize(y);
		t.start();
		matvec(A, x, y);
		t.stop();
		time_elapsed += t.elapsed();
	}
	return time_elapsed/times;
}


/*  Templatized timed parallel ompMatvec (COO and CSR) */
template <typename MatrixType>
double omp_matvec(const MatrixType& A, const Vector& x, Vector& y, size_t times){
	double time_elapsed = 0.0;
	Timer t;
	for(int j = 0; j < times; ++j){
		zeroize(y);
		t.start();
		ompMatvec(A, x, y);
		t.stop();
		time_elapsed += t.elapsed();
	}
	return time_elapsed/times;
}

#endif // __AMATH583_HPP
