#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
  #include <complex.h>
#endif

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "fastsum.h"
#include "kernels.h"
#include "infft.h"

fastsum_plan* jfastsum_alloc(){
	fastsum_plan* p = nfft_malloc(sizeof(fastsum_plan));
	return p;
}
// c wird von Julia als Float64-Pointer übergeben

void jfastsum_init( fastsum_plan* p, int d, char* s, double* c, unsigned int f, int n, int ps, float eps_I, float eps_B ){
  C (*kernel)(R, int, const R *);
  
	if ( strcmp(s, "gaussian") == 0 )
		kernel = gaussian;
	else if ( strcmp(s, "multiquadric") == 0 )
		kernel = multiquadric;
	else if ( strcmp(s, "inverse_multiquadric") == 0 )
		kernel = inverse_multiquadric;
	else if ( strcmp(s, "logarithm") == 0 )
		kernel = logarithm;
	else if ( strcmp(s, "thinplate_spline") == 0 )
		kernel = thinplate_spline;
	else if ( strcmp(s, "one_over_square") == 0 )
		kernel = one_over_square;
	else if ( strcmp(s, "one_over_modulus") == 0 )
		kernel = one_over_modulus;
	else if ( strcmp(s, "one_over_x") == 0 )
		kernel = one_over_x;
	else if ( strcmp(s, "inverse_multiquadric3") == 0 )
		kernel = inverse_multiquadric3;
	else if ( strcmp(s, "sinc_kernel") == 0 )
		kernel = sinc_kernel;
	else if ( strcmp(s, "cosc") == 0 )
		kernel = cosc;
	else if ( strcmp(s, "cot") == 0 )
		kernel = kcot;
	else if ( strcmp(s, "one_over_cube") == 0 )
		kernel = one_over_cube;
	else if ( strcmp(s, "log_sin") == 0 )
		kernel = log_sin;
	else if ( strcmp(s, "laplacian_rbf") == 0 )
		kernel = laplacian_rbf;
	else {
		s = "multiquadric";
		kernel = multiquadric;
	}
	
	fastsum_init_guru_kernel( p, d, kernel, c, f | STORE_PERMUTATION_X_ALPHA, n, ps, eps_I, eps_B);
}

double* jfastsum_set_x( fastsum_plan* p, double* x, int N, int nn, int m ){	
	fastsum_init_guru_source_nodes( p, N, nn, m );
	
	int d = p -> d;
	
	if ( p -> permutation_x_alpha == NULL ) {
	
		for ( int k = 0; k < N; k++ )
			for ( int t = 0; t < d; t++)
				p -> x[k*d+t] = x[k*d+t];
			
	} else {
		
		for ( int k = 0; k < N; k++ )
			for ( int t = 0; t < d; t++)
				p -> x[k*d+t] = x[p->permutation_x_alpha[k]*d+t];
		
	}
	
	fastsum_precompute_source_nodes( p );
		
	return p -> x;
}

double* jfastsum_set_y( fastsum_plan* p, double* y, int M, int nn, int m ){
	fastsum_init_guru_target_nodes( p, M, nn, m );
	
	int d = p -> d;
	
	for ( int j = 0; j < M; j++ )
		for ( int t = 0; t < d; t++ )
			p -> y[j*d+t] = y[j*d+t];
		
	fastsum_precompute_target_nodes( p );
		
	return p -> y;
}

double _Complex* jfastsum_set_alpha( fastsum_plan* p, double _Complex* alpha ){
	int N = p -> N_total;
	
	for ( int k = 0; k < N; k++ )
		if ( p -> permutation_x_alpha == NULL )
			p -> alpha[k] = alpha[k];
		else
			p -> alpha[k] = alpha[p->permutation_x_alpha[k]];
			
	return p -> alpha;
}


double _Complex* jfastsum_trafo( fastsum_plan* p ){
	fastsum_trafo( p );
	return p -> f;
}

double _Complex* jfastsum_exact( fastsum_plan* p ){
	fastsum_exact( p );
	return p -> f;
}

void jfastsum_finalize( fastsum_plan* p ){
	fastsum_finalize_source_nodes( p );
	fastsum_finalize_target_nodes( p );
	fastsum_finalize_kernel( p );
	return;
}