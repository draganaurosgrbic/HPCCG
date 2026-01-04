
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

#ifndef HPC_SPARSE_MATRIX_H
#define HPC_SPARSE_MATRIX_H

// These constants are upper bounds that might need to be changes for 
// pathological matrices, e.g., those with nearly dense rows/columns.

const int max_external = 100000;
const int max_num_messages = 500;
const int max_num_neighbors = max_num_messages;

#include <stdint.h>
#include <string>

typedef struct csr_t {
  uint64_t m;
  uint64_t n;
  uint64_t non_zeros;
  int64_t *row_ptr;
  int64_t *col_ind;
  double *nz;
} csr_t;

typedef struct ellpack8_row_nz_t { double val[8]; } ellpack8_row_nz_t;
typedef struct ellpack8_row_ind_t { int64_t col[8]; } ellpack8_row_ind_t;

typedef struct ellpack8_t {
  uint64_t m;
  uint64_t n;
  ellpack8_row_nz_t *nz;
  ellpack8_row_ind_t *col_ind;
} ellpack8_t;

typedef struct ellpack7_t {
  uint64_t m;
  uint64_t n;
  double *nz[7];
  int64_t *col_ind[7];
} ellpack7_t;

#define ELLPACK7_BLOCK_SIZE 8
typedef struct ellpack7_tiled_t {
  uint64_t m;
  uint64_t n;
  uint64_t non_zeros;
  double *nz;
  int64_t *col_ind;
  int total_elements_count;
} ellpack7_tiled_t;

typedef struct vector_t {
  uint64_t length;
  double *vals;
} vector_t;

struct HPC_Sparse_Matrix_STRUCT {
  char   *title;
  int start_row;
  int stop_row;
  int total_nrow;
  long long total_nnz;
  int local_nrow;
  int local_ncol;  // Must be defined in make_local_matrix
  int local_nnz;
  int  * nnz_in_row;
  double ** ptr_to_vals_in_row;
  int ** ptr_to_inds_in_row;
  double ** ptr_to_diags;

  csr_t            *csr;
  ellpack8_t       *ell8;
  ellpack7_t       *ell7;
  ellpack7_tiled_t *ell7_tiled;
  std::string selected_format;

#ifdef USING_MPI
  int num_external;
  int num_send_neighbors;
  int *external_index;
  int *external_local_index;
  int total_to_be_sent;
  int *elements_to_send;
  int *neighbors;
  int *recv_length;
  int *send_length;
  double *send_buffer;
#endif

  double *list_of_vals;   //needed for cleaning up memory
  int *list_of_inds;      //needed for cleaning up memory

};
typedef struct HPC_Sparse_Matrix_STRUCT HPC_Sparse_Matrix;


void destroyMatrix(HPC_Sparse_Matrix * &A);

#ifdef USING_SHAREDMEM_MPI
#ifndef SHAREDMEM_ALTERNATIVE
void destroySharedMemMatrix(HPC_Sparse_Matrix * &A);
#endif
#endif

#endif

