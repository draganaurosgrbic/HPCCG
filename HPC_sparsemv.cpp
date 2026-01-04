
//@HEADER
// ************************************************************************
// 
//               HPCCG: Simple Conjugate Gradient Benchmark Code
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// BSD 3-Clause License
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
// 
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/////////////////////////////////////////////////////////////////////////

// Routine to compute matrix vector product y = Ax where:
// First call exchange_externals to get off-processor values of x

// A - known matrix 
// x - known vector
// y - On exit contains Ax.

/////////////////////////////////////////////////////////////////////////

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <string>
#include <cmath>
#include "HPC_sparsemv.hpp"

#include <stdint.h>

void csr_spmv(csr_t *m, const double *x, double *y) {
  const double * __restrict__ nz = (const double *)__builtin_assume_aligned(m->nz, 64);
  const int64_t * __restrict__ col_ind = (const int64_t *)__builtin_assume_aligned(m->col_ind, 64);
  const int64_t * __restrict__ row_ptr = (const int64_t *)__builtin_assume_aligned(m->row_ptr, 64);

  const double * __restrict__ x_ptr = (const double *)__builtin_assume_aligned(x, 64);
  double * __restrict__ y_ptr = (double *)__builtin_assume_aligned(y, 64);

  #pragma omp parallel for schedule(static)
  for (uint64_t i = 0; i < m->m; ++i) {
    double sum = 0.0;
    for (int64_t j = row_ptr[i]; j < row_ptr[i+1]; ++j) {
      sum += nz[j] * x_ptr[col_ind[j]];
    }
    y_ptr[i] = sum;
  }
}

void ellpack8_spmv(ellpack8_t *m, const double * __restrict__ x, double * __restrict__ y) {
  uint64_t nrows = m->m;

  const auto * __restrict__ nz = (const ellpack8_row_nz_t *)__builtin_assume_aligned(m->nz, 64);
  const auto * __restrict__ col_ind = (const ellpack8_row_ind_t *)__builtin_assume_aligned(m->col_ind, 64);
  
  const double * __restrict__ x_ptr = (const double *)__builtin_assume_aligned(x, 64);
  double * __restrict__ y_ptr = (double *)__builtin_assume_aligned(y, 64);

  #pragma omp parallel for schedule(static)
  for (uint64_t i = 0; i < nrows; i++) {
    y_ptr[i] = nz[i].val[0] * x_ptr[col_ind[i].col[0]] +
               nz[i].val[1] * x_ptr[col_ind[i].col[1]] +
               nz[i].val[2] * x_ptr[col_ind[i].col[2]] +
               nz[i].val[3] * x_ptr[col_ind[i].col[3]] +
               nz[i].val[4] * x_ptr[col_ind[i].col[4]] +
               nz[i].val[5] * x_ptr[col_ind[i].col[5]] +
               nz[i].val[6] * x_ptr[col_ind[i].col[6]] +
               nz[i].val[7] * x_ptr[col_ind[i].col[7]];
  }
}

void ellpack7_spmv(ellpack7_t *m, const double * __restrict__ x, double * __restrict__ y) {
  uint64_t nrows = m->m;

  const double * __restrict__ nz0 = (const double *)__builtin_assume_aligned(m->nz[0], 64);
  const double * __restrict__ nz1 = (const double *)__builtin_assume_aligned(m->nz[1], 64);
  const double * __restrict__ nz2 = (const double *)__builtin_assume_aligned(m->nz[2], 64);
  const double * __restrict__ nz3 = (const double *)__builtin_assume_aligned(m->nz[3], 64);
  const double * __restrict__ nz4 = (const double *)__builtin_assume_aligned(m->nz[4], 64);
  const double * __restrict__ nz5 = (const double *)__builtin_assume_aligned(m->nz[5], 64);
  const double * __restrict__ nz6 = (const double *)__builtin_assume_aligned(m->nz[6], 64);

  const int64_t * __restrict__ ind0 = (const int64_t *)__builtin_assume_aligned(m->col_ind[0], 64);
  const int64_t * __restrict__ ind1 = (const int64_t *)__builtin_assume_aligned(m->col_ind[1], 64);
  const int64_t * __restrict__ ind2 = (const int64_t *)__builtin_assume_aligned(m->col_ind[2], 64);
  const int64_t * __restrict__ ind3 = (const int64_t *)__builtin_assume_aligned(m->col_ind[3], 64);
  const int64_t * __restrict__ ind4 = (const int64_t *)__builtin_assume_aligned(m->col_ind[4], 64);
  const int64_t * __restrict__ ind5 = (const int64_t *)__builtin_assume_aligned(m->col_ind[5], 64);
  const int64_t * __restrict__ ind6 = (const int64_t *)__builtin_assume_aligned(m->col_ind[6], 64);

  const double * __restrict__ x_ptr = (const double *)__builtin_assume_aligned(x, 64);
  double * __restrict__ y_ptr = (double *)__builtin_assume_aligned(y, 64);

  #pragma omp parallel for schedule(static)
  for (uint64_t i = 0; i < nrows; i++) {
    y_ptr[i] = nz0[i] * x_ptr[ind0[i]] +
               nz1[i] * x_ptr[ind1[i]] +
               nz2[i] * x_ptr[ind2[i]] +
               nz3[i] * x_ptr[ind3[i]] +
               nz4[i] * x_ptr[ind4[i]] +
               nz5[i] * x_ptr[ind5[i]] +
               nz6[i] * x_ptr[ind6[i]];
  }
}


void ellpack7_tiled_spmv(const ellpack7_tiled_t *matrix, const double * __restrict__ x, double * __restrict__ y) {
  uint64_t nrows = matrix->m;

  const double * __restrict__ nz = (const double *)__builtin_assume_aligned(matrix->nz, 64);
  const int64_t * __restrict__ col_ind = (const int64_t *)__builtin_assume_aligned(matrix->col_ind, 64);

  const double * __restrict__ x_ptr = (const double *)__builtin_assume_aligned(x, 64);
  double * __restrict__ y_ptr = (double *)__builtin_assume_aligned(y, 64);

  #pragma omp parallel for schedule(static)
  for (uint64_t j = 0; j < nrows; j += 8) {
    uint64_t i = j * 7;
    
    #pragma omp simd
    for (uint64_t row = 0; row < 8; ++row) {
      y_ptr[j + row] = 
        nz[i + row]      * x_ptr[col_ind[i + row]] + 
        nz[i + row + 8]  * x_ptr[col_ind[i + row + 8]] +  
        nz[i + row + 16] * x_ptr[col_ind[i + row + 16]] +  
        nz[i + row + 24] * x_ptr[col_ind[i + row + 24]] +
        nz[i + row + 32] * x_ptr[col_ind[i + row + 32]] +
        nz[i + row + 40] * x_ptr[col_ind[i + row + 40]] +
        nz[i + row + 48] * x_ptr[col_ind[i + row + 48]];
    }
  }
}


int HPC_sparsemv(HPC_Sparse_Matrix *A, const double * const x, double * const y) {
  if (A->selected_format == "tiled" && A->ell7_tiled != NULL) {
      ellpack7_tiled_spmv(A->ell7_tiled, x, y);
      return 0;
  } 
  else if (A->selected_format == "ell8" && A->ell8 != NULL) {
      ellpack8_spmv(A->ell8, x, y);
      return 0;
  }
  else if (A->selected_format == "ell7" && A->ell7 != NULL) {
      ellpack7_spmv(A->ell7, x, y);
      return 0;
  }
  else if (A->selected_format == "csr" && A->csr != NULL) {
      csr_spmv(A->csr, x, y);
      return 0;
  }

  const int nrow = (const int) A->local_nrow;

#ifdef USING_OMP
#pragma omp parallel for schedule(static)
#endif
  for (int i=0; i< nrow; i++)
    {
      double sum = 0.0;
      const double * const cur_vals = 
     (const double * const) A->ptr_to_vals_in_row[i];

      const int    * const cur_inds = 
     (const int    * const) A->ptr_to_inds_in_row[i];

      const int cur_nnz = (const int) A->nnz_in_row[i];

      for (int j=0; j< cur_nnz; j++)
          sum += cur_vals[j]*x[cur_inds[j]];
      y[i] = sum;
    }
  return(0);
}
