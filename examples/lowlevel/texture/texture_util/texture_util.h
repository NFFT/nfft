#include<stdio.h>

#define MAX_LINE 1000

void read_h(int *N1_ptr, double **h_phi_ptr, double **h_theta_ptr, FILE * in,
						FILE * out);

void read_r(int N1, int *N2_ptr, double **r_ptr, FILE * in, FILE * out);

void read_grid(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
							 double **h_theta_ptr, double **r_ptr, FILE * h_in, FILE * r_in,
							 FILE * out);

void read_x(int *N1_ptr, int *N2_ptr, complex ** x_ptr, FILE * in,
						FILE * out);

void read_samples(int *N1_ptr, int *N2_ptr, double **h_phi_ptr,
									double **h_theta_ptr, double **r_ptr, complex ** x_ptr,
									FILE * h_in, FILE * r_in, FILE * x_in, FILE * out);

void read_omega(int *N_ptr, complex ** omega_ptr, FILE * in, FILE * out);

void write_r(int N2, double *r, FILE * out);

void write_h(int N1, double *h_phi, double *h_theta, FILE * out);

void write_omega(int N, complex * omega, FILE * out);

void write_x(int N1, int N2, complex * x, FILE * out);

double equidist(double start, double end, int i, int n, int incl);

const char *grid_descr[2];

void calculate_grid(int h_phi_count, int h_theta_count, int r_phi_count,
										int r_theta_count, double *h_phi, double *h_theta,
										double *r, int type);

void error(const char *msg);

const char *omega_policy_descr[2];

void init_omega(complex * omega, int N, int policy);

void mult_error(int N1, int N2, complex *x, double min_err, double max_err);

void block_mult_error(int N1, int N2, complex *x, double min_err,
											double max_err);
