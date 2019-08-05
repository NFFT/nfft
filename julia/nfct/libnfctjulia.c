#include "config.h"

#include <stdio.h>

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "nfft3.h"
#include "infft.h"

nfct_plan* jnfct_alloc(void){
	nfct_plan* p = nfft_malloc(sizeof(nfct_plan));
	return p;
}

void jnfct_init(nfct_plan* p, int d, int* N, int M, int* n, int m, unsigned int f1, unsigned int f2){
    nfct_init_guru(p,d,N,M,n,m,f1,f2);
}

double* jnfct_set_x(nfct_plan* p, double* X){
	int M = p->M_total;
	int d = p->d;
	int r,c;
	for (r = 0; r < M; r++)
		for (c = 0; c < d; c++)
			p->x[d*r+c] = X[d*r+c];
	nfct_precompute_one_psi(p);
	return p->x;
}

// setting Fourier coefficients and returning pointer for access by Julia
double* jnfct_set_fhat(nfct_plan* p,double* f_hat){
	int n = p->N_total;
	int k;
	for (k=0;k<n;k++)
		p->f_hat[k] = f_hat[k];
	return p->f_hat;
}

// setting values and returning pointer for access by Julia
double* jnfct_set_f(nfct_plan* p,double* f){
	int M = p->M_total;
	int j;
	for (j=0;j<M;j++)
		p->f[j] = f[j];
	return p->f;
}

// nfct trafo, return pointer to values for access by Julia if pointer isn't set
double* jnfct_trafo(nfct_plan* p){
	nfct_trafo(p);
	return p->f;
}

// nfct adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double* jnfct_adjoint(nfct_plan* p){
	nfct_adjoint(p);
	return p->f_hat;
}

// nfct trafo, return pointer to values for access by Julia if pointer isn't set
double* jnfct_trafo_direct(nfct_plan* p){
	nfct_trafo_direct(p);
	return p->f;
}

// nfct adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double* jnfct_adjoint_direct(nfct_plan* p){
	nfct_adjoint_direct(p);
	return p->f_hat;
}

// nfct plan finalizer
void jnfct_finalize(nfct_plan* p){
	nfct_finalize(p);
	nfft_free(p);
}
