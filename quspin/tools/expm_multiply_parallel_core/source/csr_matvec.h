#ifndef __CSR_MATVEC_H__
#define __CSR_MATVEC_H__


#if defined(_OPENMP)
#include "csrmv_merge.h"
#include <iostream>
template<typename I, typename T1,typename T2,typename T3>
void inline csr_matvec(const bool overwrite_y,
						const I n,
						const I Ap[],
						const I Aj[],
						const T1 Ax[],
						const T2 a,
						const T3 x[],
							  I rco[],
							  T3 vco[],
							  T3 y[])
{
	csrmv_merge(overwrite_y,n,Ap,Aj,Ax,a,x,rco,vco,y);
}

template<typename I, typename T1,typename T2,typename T3>
void inline csr_matvec(const bool overwrite_y,
						const I n,
						const I Ap[],
						const I Aj[],
						const T1 Ax[],
						const T2 a,
						const T3 x[],
							  T3 y[])
{
	int nthread = omp_get_max_threads();
	I * rco = new I[nthread];
	T3 * vco = new T3[nthread];

	for(int tid=0;tid<nthread;tid++){
		rco[tid] = 0;
		vco[tid] = 0;
	}

	#pragma omp parallel
	{
		csrmv_merge(overwrite_y,n,Ap,Aj,Ax,a,x,rco,vco,y);
	}

	delete [] rco;
	delete [] vco;
	rco = NULL;
	vco = NULL;
}

#else
template<typename I, typename T1,typename T2,typename T3>
void csr_matvec(const bool overwrite_y,
				const I n,
				const I Ap[],
				const I Aj[],
				const T1 Ax[],
				const T2 a,
				const T3 x[],
					  I rco[],
					  T3 vco[],
					  T3 y[])
{
	if(overwrite_y){
		for(I k = 0; k<n; k++){
			T3 sum = 0;
			for(I jj = Ap[k]; jj < Ap[k+1]; jj++){
				sum += Ax[jj] * x[Aj[jj]];
			}
			y[k] = a*sum;
		}
	}else{
		for(I k = 0; k<n; k++){
			T3 sum = 0;
			for(I jj = Ap[k]; jj < Ap[k+1]; jj++){
				sum += Ax[jj] * x[Aj[jj]];
			}
			y[k] += a*sum;
		}
	}

}

template<typename I, typename T1,typename T2,typename T3>
void csr_matvec(const bool overwrite_y,
				const I n,
				const I Ap[],
				const I Aj[],
				const T1 Ax[],
				const T2 a,
				const T3 x[],
					  T3 y[])
{
	if(overwrite_y){
		for(I k = 0; k<n; k++){
			T3 sum = 0;
			for(I jj = Ap[k]; jj < Ap[k+1]; jj++){
				sum += Ax[jj] * x[Aj[jj]];
			}
			y[k] = a*sum;
		}
	}else{
		for(I k = 0; k<n; k++){
			T3 sum = 0;
			for(I jj = Ap[k]; jj < Ap[k+1]; jj++){
				sum += Ax[jj] * x[Aj[jj]];
			}
			y[k] += a*sum;
		}
	}

}
#endif

#endif