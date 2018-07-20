cimport numpy as np
import cython
from cython.parallel cimport parallel
from libcpp cimport bool
from libc.stdlib cimport malloc, free
cimport openmp

ctypedef double complex cdouble
ctypedef float complex cfloat

cdef extern from "csr_matvec.h":
    void csr_matvec[I,T1,T2,T3](const bool,const I,const I[],const I[],const T1[],
                              const T2,const T3[],T3[]) nogil

ctypedef fused index:
  np.int32_t
  np.int64_t

ctypedef fused T1:
  float
  double
  float complex
  double complex

ctypedef fused T2:
  float
  double
  float complex
  double complex

ctypedef fused T3:
  float
  double
  float complex
  double complex

@cython.boundscheck(False)   
def _csr_matvec(bool overwrite_y, index[:] Ap, index[:] Aj,
                  T1[:] Ax, T2 alpha, T3[:] Xx, T3[:] Yx):
  cdef index nr = Yx.shape[0]

  if T1 is double and T2 is double and T3 is double:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is double and T2 is double and T3 is cdouble:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is double and T2 is cdouble and T3 is cdouble:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is cdouble and T2 is cdouble and T3 is cdouble:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is cdouble and T2 is double and T3 is cdouble:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is float and T2 is float and T3 is float:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is float and T2 is float and T3 is cfloat:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is float and T2 is cfloat and T3 is cfloat:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is cfloat and T2 is cfloat and T3 is cfloat:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  elif T1 is cfloat and T2 is float and T3 is cfloat:
    with nogil:
      csr_matvec(overwrite_y,nr,&Ap[0],&Aj[0],&Ax[0],alpha,&Xx[0],&Yx[0])
  else:
    raise ValueError("imcompatbile types")

