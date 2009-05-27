#include <assert.h>
#include <nfft3.h>
#include <nfft3util.h>
#include "data_nfft3util.h"

double normalise_phi(double angle)
{
	while (angle >= PI) {
		angle -= 2 * PI;
	}
	while (angle < -PI) {
		angle += 2 * PI;
	}
	assert(angle < PI);
	return angle;
}

double normalise_theta(double angle)
{
	while (angle > PI) {
		angle -= 2 * PI;
	}
	while (angle < 0) {
		angle += 2 * PI;
	}
	assert(angle <= PI);
	return angle;
}

double deg2rad(double angle)
{
	return angle / 360.0 * 2 * PI;
}

double invert_phi(double phi)
{
	return (phi + PI);
}

double invert_theta(double theta)
{
	return (PI - theta);
}

void normalise_h(int N1, double *h_phi, double *h_theta)
{
	int i;

	for (i = 0; i < N1; i++) {
		h_phi[i] = normalise_phi(h_phi[i]);
		h_theta[i] = normalise_theta(h_theta[i]);
	}
}

void normalise_r(int N2, double *r)
{
	int j;

	for (j = 0; j < N2; j++) {
		r[2 * j] = normalise_phi(r[2 * j]);
		r[2 * j + 1] = normalise_theta(r[2 * j + 1]);
	}
}

void expand_r(int N2, double *r)
{
	int j;

	for (j = 0; j < N2 / 2; j++) {
		r[2 * j + N2] = invert_phi(r[2 * j]);
		r[2 * j + 1 + N2] = invert_theta(r[2 * j + 1]);
	}
}

void expand_h(int N1, double *h_phi, double *h_theta)
{
	int i;

	for (i = 0; i < N1 / 2; i++) {
		h_phi[i + N1 / 2] = invert_phi(h_phi[i]);
		h_theta[i + N1 / 2] = invert_theta(h_theta[i]);
	}
}
