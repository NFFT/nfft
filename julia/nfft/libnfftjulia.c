#include "config.h"

#include <stdio.h>

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "nfft3.h"

nfft_plan* jnfft_alloc(void) {
	nfft_plan* p = nfft_malloc(sizeof(nfft_plan));
	return p;
}

void jnfft_init(nfft_plan* p, int d, int* N, int M, int* n, int m, unsigned int f1, unsigned int f2){
    nfft_init_guru(p,d,N,M,n,m,f1,f2);
}

double* jnfft_set_x(nfft_plan* p, double* X){
	int M = p->M_total;
	int d = p->d;
	int r,c;
	for (r = 0; r < M; r++)
		for (c = 0; c < d; c++)
			p->x[d*r+c] = X[d*r+c];
	nfft_precompute_one_psi(p);
	return p->x;
}

// setting Fourier coefficients and returning pointer for access by Julia
double _Complex* jnfft_set_fhat(nfft_plan* p,double _Complex* f_hat){
	int n = p->N_total;
	int k;
	for (k=0;k<n;k++)
		p->f_hat[k] = f_hat[k];
	return p->f_hat;
}

// setting values and returning pointer for access by Julia
double _Complex* jnfft_set_f(nfft_plan* p,double _Complex* f){
	int M = p->M_total;
	int j;
	for (j=0;j<M;j++)
		p->f[j] = f[j];
	return p->f;
}

// nfft trafo, return pointer to values for access by Julia if pointer isn't set
double _Complex* jnfft_trafo(nfft_plan* p){
	nfft_trafo(p);
	return p->f;
}

// nfft adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double _Complex* jnfft_adjoint(nfft_plan* p){
	nfft_adjoint(p);
	return p->f_hat;
}

// nfft trafo, return pointer to values for access by Julia if pointer isn't set
double _Complex* jnfft_trafo_direct(nfft_plan* p){
	nfft_trafo_direct(p);
	return p->f;
}

// nfft adjoint, return pointer to coefficients for access by Julia if pointer isn't set
double _Complex* jnfft_adjoint_direct(nfft_plan* p){
	nfft_adjoint_direct(p);
	return p->f_hat;
}

// nfft plan finalizer
void jnfft_finalize(nfft_plan* p){
	nfft_finalize(p);
	nfft_free(p);
}
