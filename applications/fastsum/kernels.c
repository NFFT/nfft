/*! \file kernels.c
 *  \brief File with predefined kernels for the fast summation algorithm.
 */

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <float.h>

/** 
 * \addtogroup applications_fastsum
 * \{
 */

complex double gaussian(double x, int der, const double *param)    /* K(x)=exp(-x^2/c^2) */
{
  double c=param[0];
  double value=0.0;

  switch (der)
  {
    case  0 : value=exp(-x*x/(c*c)); break;
    case  1 : value=-2.0*x/(c*c)*exp(-x*x/(c*c)); break;
    case  2 : value=2.0*exp(-x*x/(c*c))*(-c*c+2.0*x*x)/(c*c*c*c); break;
    case  3 : value=-4.0*x*exp(-x*x/(c*c))*(-3.0*c*c+2.0*x*x)/(c*c*c*c*c*c); break;
    case  4 : value=4.0*exp(-x*x/(c*c))*(3.0*c*c*c*c-12.0*c*c*x*x+4.0*x*x*x*x)/(c*c*c*c*c*c*c*c); break;
    case  5 : value=-8.0*x*exp(-x*x/(c*c))*(15.0*c*c*c*c-20.0*c*c*x*x+4.0*x*x*x*x)/pow(c,10.0); break;
    case  6 : value=8.0*exp(-x*x/(c*c))*(-15.0*c*c*c*c*c*c+90.0*x*x*c*c*c*c-60.0*x*x*x*x*c*c+8.0*x*x*x*x*x*x)/pow(c,12.0); break;
    case  7 : value=-16.0*x*exp(-x*x/(c*c))*(-105.0*c*c*c*c*c*c+210.0*x*x*c*c*c*c-84.0*x*x*x*x*c*c+8.0*x*x*x*x*x*x)/pow(c,14.0); break;
    case  8 : value=16.0*exp(-x*x/(c*c))*(105.0*c*c*c*c*c*c*c*c-840.0*x*x*c*c*c*c*c*c+840.0*x*x*x*x*c*c*c*c-224.0*x*x*x*x*x*x*c*c+16.0*x*x*x*x*x*x*x*x)/pow(c,16.0); break;
    case  9 : value=-32.0*x*exp(-x*x/(c*c))*(945.0*c*c*c*c*c*c*c*c-2520.0*x*x*c*c*c*c*c*c+1512.0*x*x*x*x*c*c*c*c-288.0*x*x*x*x*x*x*c*c+16.0*x*x*x*x*x*x*x*x)/pow(c,18.0); break;
    case 10 : value=32.0*exp(-x*x/(c*c))*(-945.0*pow(c,10.0)+9450.0*x*x*c*c*c*c*c*c*c*c-12600.0*x*x*x*x*c*c*c*c*c*c+5040.0*x*x*x*x*x*x*c*c*c*c-720.0*x*x*x*x*x*x*x*x*c*c+32.0*pow(x,10.0))/pow(c,20.0); break;
    case 11 : value=-64.0*x*exp(-x*x/(c*c))*(-10395.0*pow(c,10.0)+34650.0*x*x*c*c*c*c*c*c*c*c-27720.0*x*x*x*x*c*c*c*c*c*c+7920.0*x*x*x*x*x*x*c*c*c*c-880.0*x*x*x*x*x*x*x*x*c*c+32.0*pow(x,10.0))/pow(c,22.0); break;
    case 12 : value=64.0*exp(-x*x/(c*c))*(10395.0*pow(c,12.0)-124740.0*x*x*pow(c,10.0)+207900.0*x*x*x*x*c*c*c*c*c*c*c*c-110880.0*x*x*x*x*x*x*c*c*c*c*c*c+23760.0*x*x*x*x*x*x*x*x*c*c*c*c-2112.0*pow(x,10.0)*c*c+64.0*pow(x,12.0))/pow(c,24.0); break;
    default : value=0.0;
  }

  return value;
}

complex double multiquadric(double x, int der, const double *param)    /* K(x)=sqrt(x^2+c^2) */
{
  double c=param[0];
  double value=0.0;

  switch (der)
  {
    case  0 : value=sqrt(x*x+c*c); break;
    case  1 : value=1.0/(sqrt(x*x+c*c))*x; break;
    case  2 : value=c*c/sqrt(pow(x*x+c*c,3.0)); break;
    case  3 : value=-3.0*x*c*c/sqrt(pow(x*x+c*c,5.0)); break;
    case  4 : value=3.0*c*c*(4.0*x*x-c*c)/sqrt(pow(x*x+c*c,7.0)); break;
    case  5 : value=-15.0*x*c*c*(4.0*x*x-3.0*c*c)/sqrt(pow(x*x+c*c,9.0)); break;
    case  6 : value=45.0*c*c*(8.0*x*x*x*x-12.0*x*x*c*c+c*c*c*c)/sqrt(pow(x*x+c*c,11.0)); break;
    case  7 : value=-315.0*x*c*c*(8.0*x*x*x*x-20.0*x*x*c*c+5.0*c*c*c*c)/sqrt(pow(x*x+c*c,13.0)); break;
    case  8 : value=315.0*c*c*(64.0*x*x*x*x*x*x-240.0*x*x*x*x*c*c+120.0*x*x*c*c*c*c-5.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,15.0)); break;
    case  9 : value=-2835.0*x*c*c*(64.0*x*x*x*x*x*x-336.0*x*x*x*x*c*c+280.0*x*x*c*c*c*c-35.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,17.0)); break;
    case 10 : value=14175.0*c*c*(128.0*x*x*x*x*x*x*x*x-896.0*x*x*x*x*x*x*c*c+1120.0*x*x*x*x*c*c*c*c-280.0*x*x*c*c*c*c*c*c+7.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,19.0)); break;
    case 11 : value=-155925.0*x*c*c*(128.0*x*x*x*x*x*x*x*x-1152.0*x*x*x*x*x*x*c*c+2016.0*x*x*x*x*c*c*c*c-840.0*x*x*c*c*c*c*c*c+63.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,21.0)); break;
    case 12 : value=467775.0*c*c*(1260.0*x*x*c*c*c*c*c*c*c*c-21.0*pow(c,10.0)+512.0*pow(x,10.0)-5760.0*x*x*x*x*x*x*x*x*c*c+13440.0*x*x*x*x*x*x*c*c*c*c-8400.0*x*x*x*x*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,23.0)); break;
    default : value=0.0;
  }

  return value;
}

complex double inverse_multiquadric(double x, int der, const double *param)    /* K(x)=1/sqrt(x^2+c^2) */
{
  double c=param[0];
  double value=0.0;

  switch (der)
  {
    case  0 : value=1.0/sqrt(x*x+c*c); break;
    case  1 : value=-1.0/(sqrt(pow(x*x+c*c,3.0)))*x; break;
    case  2 : value=(2.0*x*x-c*c)/sqrt(pow(x*x+c*c,5.0)); break;
    case  3 : value=-3.0*x*(2.0*x*x-3.0*c*c)/sqrt(pow(x*x+c*c,7.0)); break;
    case  4 : value=3.0*(8.0*x*x*x*x-24.0*x*x*c*c+3.0*c*c*c*c)/sqrt(pow(x*x+c*c,9.0)); break;
    case  5 : value=-15.0*x*(8.0*x*x*x*x-40.0*x*x*c*c+15.0*c*c*c*c)/sqrt(pow(x*x+c*c,11.0)); break;
    case  6 : value=45.0*(16.0*x*x*x*x*x*x-120.0*x*x*x*x*c*c+90.0*x*x*c*c*c*c-5.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,13.0)); break;
    case  7 : value=-315.0*x*(16.0*x*x*x*x*x*x-168.0*x*x*x*x*c*c+210.0*x*x*c*c*c*c-35.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,15.0)); break;
    case  8 : value=315.0*(128.0*x*x*x*x*x*x*x*x-1792.0*x*x*x*x*x*x*c*c+3360.0*x*x*x*x*c*c*c*c-1120.0*x*x*c*c*c*c*c*c+35.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,17.0)); break;
    case  9 : value=-2835.0*x*(128.0*x*x*x*x*x*x*x*x-2304.0*x*x*x*x*x*x*c*c+6048.0*x*x*x*x*c*c*c*c-3360.0*x*x*c*c*c*c*c*c+315.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,19.0)); break;
    case 10 : value=14175.0*(256.0*pow(x,10.0)-5760.0*x*x*x*x*x*x*x*x*c*c+20160.0*x*x*x*x*x*x*c*c*c*c-16800.0*x*x*x*x*c*c*c*c*c*c+3150.0*x*x*c*c*c*c*c*c*c*c-63.0*pow(c,10.0))/sqrt(pow(x*x+c*c,21.0)); break;
    case 11 : value=-155925.0*x*(256.0*pow(x,10.0)-7040.0*x*x*x*x*x*x*x*x*c*c+31680.0*x*x*x*x*x*x*c*c*c*c-36960.0*x*x*x*x*c*c*c*c*c*c+11550.0*x*x*c*c*c*c*c*c*c*c-693.0*pow(c,10.0))/sqrt(pow(x*x+c*c,23.0)); break;
    case 12 : value=467775.0*(231.0*pow(c,12.0)+190080.0*x*x*x*x*x*x*x*x*c*c*c*c-16632.0*x*x*pow(c,10.0)-295680.0*x*x*x*x*x*x*c*c*c*c*c*c+138600.0*x*x*x*x*c*c*c*c*c*c*c*c+1024.0*pow(x,12.0)-33792.0*pow(x,10.0)*c*c)/sqrt(pow(x*x+c*c,25.0)); break;
    default : value=0.0;
  }

  return value;
}

complex double logarithm(double x, int der, const double *param)    /* K(x)=log |x| */
{
  double value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=log(fabs(x)); break;
    case  1 : value=(x<0 ? -1 : 1)/fabs(x); break;
    case  2 : value=-1/(x*x); break;
    case  3 : value=2.0*(x<0 ? -1 : 1)/pow(fabs(x),3.0); break;
    case  4 : value=-6.0/(x*x*x*x); break;
    case  5 : value=24.0*(x<0 ? -1 : 1)/pow(fabs(x),5.0); break;
    case  6 : value=-120.0/(x*x*x*x*x*x); break;
    case  7 : value=720.0*(x<0 ? -1 : 1)/pow(fabs(x),7.0); break;
    case  8 : value=-5040.0/(x*x*x*x*x*x*x*x); break;
    case  9 : value=40320.0*(x<0 ? -1 : 1)/pow(fabs(x),9.0); break;
    case 10 : value=-362880.0/pow(x,10.0); break;
    case 11 : value=3628800.0*(x<0 ? -1 : 1)/pow(fabs(x),11.0); break;
    case 12 : value=-39916800.0/pow(x,12.0); break;
    case 13 : value=479001600.0/pow(x,13.0); break;
    case 14 : value=-6227020800.0/pow(x,14.0); break;
    case 15 : value=87178291200.0/pow(x,15.0); break;
    case 16 : value=-1307674368000.0/pow(x,16.0); break;
    case 17 : value=20922789888000.0/pow(x,17.0); break;
    default : value=0.0;
  }

  return value;
}

complex double thinplate_spline(double x, int der, const double *param)    /* K(x) = x^2 log |x| */
{
  double value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=x*x*log(fabs(x)); break;
    case  1 : value=2.0*x*log(fabs(x))+x; break;
    case  2 : value=2.0*log(fabs(x))+3.0; break;
    case  3 : value=2.0/x; break;
    case  4 : value=-2.0/(x*x); break;
    case  5 : value=4.0/(x*x*x); break;
    case  6 : value=-12.0/(x*x*x*x); break;
    case  7 : value=48.0/(x*x*x*x*x); break;
    case  8 : value=-240.0/(x*x*x*x*x*x); break;
    case  9 : value=1440.0/(x*x*x*x*x*x*x); break;
    case 10 : value=-10080.0/(x*x*x*x*x*x*x*x); break;
    case 11 : value=80640.0/(x*x*x*x*x*x*x*x*x); break;
    case 12 : value=-725760.0/pow(x,10.0); break;
    default : value=0.0;
  }

  return value;
}

complex double one_over_square(double x, int der, const double *param)    /* K(x) = 1/x^2 */
{
  double value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=1.0/(x*x); break;
    case  1 : value=-2.0/(x*x*x); break;
    case  2 : value=6.0/(x*x*x*x); break;
    case  3 : value=-24.0/(x*x*x*x*x); break;
    case  4 : value=120.0/(x*x*x*x*x*x); break;
    case  5 : value=-720.0/(x*x*x*x*x*x*x); break;
    case  6 : value=5040.0/(x*x*x*x*x*x*x*x); break;
    case  7 : value=-40320.0/(x*x*x*x*x*x*x*x*x); break;
    case  8 : value=362880.0/pow(x,10.0); break;
    case  9 : value=-3628800.0/pow(x,11.0); break;
    case 10 : value=39916800.0/pow(x,12.0); break;
    case 11 : value=-479001600.0/pow(x,13.0); break;
    case 12 : value=6227020800.0/pow(x,14.0); break;
    default : value=0.0;
  }

  return value;
}

complex double one_over_modulus(double x, int der, const double *param)    /* K(x) = 1/|x| */
{
  double value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=1.0/fabs(x); break;
    case  1 : value=-1/x/fabs(x); break;
    case  2 : value=2.0/pow(fabs(x),3.0); break;
    case  3 : value=-6.0/(x*x*x)/fabs(x); break;
    case  4 : value=24.0/pow(fabs(x),5.0); break;
    case  5 : value=-120.0/(x*x*x*x*x)/fabs(x); break;
    case  6 : value=720.0/pow(fabs(x),7.0); break;
    case  7 : value=-5040.0/(x*x*x*x*x*x*x)/fabs(x); break;
    case  8 : value=40320.0/pow(fabs(x),9.0); break;
    case  9 : value=-362880.0/(x*x*x*x*x*x*x*x*x)/fabs(x); break;
    case 10 : value=3628800.0/pow(fabs(x),11.0); break;
    case 11 : value=-39916800.0/pow(x,11.0)/fabs(x); break;
    case 12 : value=479001600.0/pow(fabs(x),13.0); break;
    default : value=0.0;
  }

  return value;
}

complex double one_over_x(double x, int der, const double *param)    /* K(x) = 1/x */
{
  double value=0.0;

  (void)param;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=1.0/x; break;
    case  1 : value=-1.0/(x*x); break;
    case  2 : value=2.0/(x*x*x); break;
    case  3 : value=-6.0/(x*x*x*x); break;
    case  4 : value=24.0/(x*x*x*x*x); break;
    case  5 : value=-120.0/(x*x*x*x*x*x); break;
    case  6 : value=720.0/(x*x*x*x*x*x*x); break;
    case  7 : value=-5040.0/(x*x*x*x*x*x*x*x); break;
    case  8 : value=40320.0/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-362880.0/pow(x,10.0); break;
    case 10 : value=3628800.0/pow(x,11.0); break;
    case 11 : value=-39916800.0/pow(x,12.0); break;
    case 12 : value=479001600.0/pow(x,13.0); break;
    default : value=0.0;
  }

  return value;
}

complex double inverse_multiquadric3(double x, int der, const double *param)    /* K(x) = 1/sqrt(x^2+c^2)^3 */
{
  double c=param[0];
  double value=0.0;

  switch (der)
  {
    case  0 : value=1.0/(sqrt(pow(x*x+c*c,3.0))); break;
    case  1 : value=-3.0/sqrt(pow(x*x+c*c,5.0))*x; break;
    case  2 : value=3.0*(4.0*x*x-c*c)/sqrt(pow(x*x+c*c,7.0)); break;
    case  3 : value=-15.0*x*(4.0*x*x-3.0*c*c)/sqrt(pow(x*x+c*c,9.0)); break;
    case  4 : value=45.0*(8.0*x*x*x*x-12.0*x*x*c*c+c*c*c*c)/sqrt(pow(x*x+c*c,11.0)); break;
    case  5 : value=-315.0*x*(8.0*x*x*x*x-20.0*x*x*c*c+5.0*c*c*c*c)/sqrt(pow(x*x+c*c,13.0)); break;
    case  6 : value=315.0*(64.0*x*x*x*x*x*x-240.0*x*x*x*x*c*c+120.0*x*x*c*c*c*c-5.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,15.0)); break;
    case  7 : value=-2835.0*x*(64.0*x*x*x*x*x*x-336.0*x*x*x*x*c*c+280.0*x*x*c*c*c*c-35.0*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,17.0)); break;
    case  8 : value=14175.0*(128.0*x*x*x*x*x*x*x*x-896.0*x*x*x*x*x*x*c*c+1120.0*x*x*x*x*c*c*c*c-280.0*x*x*c*c*c*c*c*c+7.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,19.0)); break;
    case  9 : value=-155925.0*x*(128.0*x*x*x*x*x*x*x*x-1152.0*x*x*x*x*x*x*c*c+2016.0*x*x*x*x*c*c*c*c-840.0*x*x*c*c*c*c*c*c+63.0*c*c*c*c*c*c*c*c)/sqrt(pow(x*x+c*c,21.0)); break;
    case 10 : value=467775.0*(512.0*pow(x,10.0)-5760.0*x*x*x*x*x*x*x*x*c*c+13440.0*x*x*x*x*x*x*c*c*c*c-8400.0*x*x*x*x*c*c*c*c*c*c+1260.0*x*x*c*c*c*c*c*c*c*c-21.0*pow(c,10.0))/sqrt(pow(x*x+c*c,23.0)); break;
    case 11 : value=-6081075.0*x*(512.0*pow(x,10.0)-7040.0*x*x*x*x*x*x*x*x*c*c+21120.0*x*x*x*x*x*x*c*c*c*c-18480.0*x*x*x*x*c*c*c*c*c*c+4620.0*x*x*c*c*c*c*c*c*c*c-231.0*pow(c,10.0))/sqrt(pow(x*x+c*c,25.0)); break;
    case 12 : value=42567525.0*(1024.0*pow(x,12.0)+27720.0*x*x*x*x*c*c*c*c*c*c*c*c+33.0*pow(c,12.0)-2772.0*x*x*pow(c,10.0)-73920.0*x*x*x*x*x*x*c*c*c*c*c*c+63360.0*x*x*x*x*x*x*x*x*c*c*c*c-16896.0*pow(x,10.0)*c*c)/sqrt(pow(x*x+c*c,27.0)); break;
    default : value=0.0;
  }

  return value;
}

complex double sinc_kernel(double x, int der, const double *param)    /* K(x) = sin(cx)/x */
{
  double c=param[0];
  double value=0.0;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=sin(c*x)/x; break;
    case  1 : value=(cos(c*x)*c*x-sin(c*x))/(x*x); break;
    case  2 : value=-(sin(c*x)*c*c*x*x+2.0*cos(c*x)*c*x-2.0*sin(c*x))/(x*x*x); break;
    case  3 : value=-(cos(c*x)*c*c*c*x*x*x-3.0*sin(c*x)*c*c*x*x-6.0*cos(c*x)*c*x+6.0*sin(c*x))/(x*x*x*x); break;
    case  4 : value=(sin(c*x)*c*c*c*c*x*x*x*x+4.0*cos(c*x)*c*c*c*x*x*x-12.0*sin(c*x)*c*c*x*x-24.0*cos(c*x)*c*x+24.0*sin(c*x))/(x*x*x*x*x); break;
    case  5 : value=(cos(c*x)*c*c*c*c*c*x*x*x*x*x-5.0*sin(c*x)*c*c*c*c*x*x*x*x-20.0*cos(c*x)*c*c*c*x*x*x+60.0*sin(c*x)*c*c*x*x+120.0*cos(c*x)*c*x-120.0*sin(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+6.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x-30.0*sin(c*x)*c*c*c*c*x*x*x*x-120.0*cos(c*x)*c*c*c*x*x*x+360.0*sin(c*x)*c*c*x*x+720.0*cos(c*x)*c*x-720.0*sin(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=-(cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-7.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-42.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+210.0*sin(c*x)*c*c*c*c*x*x*x*x+840.0*cos(c*x)*c*c*c*x*x*x-2520.0*sin(c*x)*c*c*x*x-5040.0*cos(c*x)*c*x+5040.0*sin(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+8.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-56.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-336.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+1680.0*sin(c*x)*c*c*c*c*x*x*x*x+6720.0*cos(c*x)*c*c*c*x*x*x-20160.0*sin(c*x)*c*c*x*x-40320.0*cos(c*x)*c*x+40320.0*sin(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=(cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-9.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-72.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+504.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+3024.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x-15120.0*sin(c*x)*c*c*c*c*x*x*x*x-60480.0*cos(c*x)*c*c*c*x*x*x+181440.0*sin(c*x)*c*c*x*x+362880.0*cos(c*x)*c*x-362880.0*sin(c*x))/pow(x,10.0); break;
    case 10 : value=-(sin(c*x)*pow(c,10.0)*pow(x,10.0)+10.0*cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-90.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-720.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+5040.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+30240.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x-151200.0*sin(c*x)*c*c*c*c*x*x*x*x-604800.0*cos(c*x)*c*c*c*x*x*x+1814400.0*sin(c*x)*c*c*x*x+3628800.0*cos(c*x)*c*x-3628800.0*sin(c*x))/pow(x,11.0); break;
    case 11 : value=-(cos(c*x)*pow(c,11.0)*pow(x,11.0)-11.0*sin(c*x)*pow(c,10.0)*pow(x,10.0)-110.0*cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+990.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+7920.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-55440.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-332640.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+1663200.0*sin(c*x)*c*c*c*c*x*x*x*x+6652800.0*cos(c*x)*c*c*c*x*x*x-19958400.0*sin(c*x)*c*c*x*x-39916800.0*cos(c*x)*c*x+39916800.0*sin(c*x))/pow(x,12.0); break;
    case 12 : value=(sin(c*x)*pow(c,12.0)*pow(x,12.0)+12.0*cos(c*x)*pow(c,11.0)*pow(x,11.0)-132.0*sin(c*x)*pow(c,10.0)*pow(x,10.0)-1320.0*cos(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+11880.0*sin(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+95040.0*cos(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-665280.0*sin(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-3991680.0*cos(c*x)*c*c*c*c*c*x*x*x*x*x+19958400.0*sin(c*x)*c*c*c*c*x*x*x*x+79833600.0*cos(c*x)*c*c*c*x*x*x-239500800.0*sin(c*x)*c*c*x*x-479001600.0*cos(c*x)*c*x+479001600.0*sin(c*x))/pow(x,13.0); break;
    default : value=0.0;
  }

  return value;
}

complex double cosc(double x, int der, const double *param)    /* K(x) = cos(cx)/x */
{
  double c=param[0];
  double value=0.0;
  double sign;

  if (x<0) sign=-1.0; else sign=1.0;
  x=fabs(x);

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value=cos(c*x)/x; break;
    case  1 : value=-(sin(c*x)*c*x+cos(c*x))/(x*x); break;
    case  2 : value=(-cos(c*x)*c*c*x*x+2.0*sin(c*x)*c*x+2.0*cos(c*x))/(x*x*x); break;
    case  3 : value=(sin(c*x)*c*c*c*x*x*x+3.0*cos(c*x)*c*c*x*x-6.0*sin(c*x)*c*x-6.0*cos(c*x))/(x*x*x*x); break;
    case  4 : value=(cos(c*x)*c*c*c*c*x*x*x*x-4.0*sin(c*x)*c*c*c*x*x*x-12.0*cos(c*x)*c*c*x*x+24.0*sin(c*x)*c*x+24.0*cos(c*x))/(x*x*x*x*x); break;
    case  5 : value=-(sin(c*x)*c*c*c*c*c*x*x*x*x*x+5.0*cos(c*x)*c*c*c*c*x*x*x*x-20.0*sin(c*x)*c*c*c*x*x*x-60.0*cos(c*x)*c*c*x*x+120.0*sin(c*x)*c*x+120.0*cos(c*x))/(x*x*x*x*x*x); break;
    case  6 : value=-(cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-6.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-30.0*cos(c*x)*c*c*c*c*x*x*x*x+120.0*sin(c*x)*c*c*c*x*x*x+360.0*cos(c*x)*c*c*x*x-720.0*sin(c*x)*c*x-720.0*cos(c*x))/(x*x*x*x*x*x*x); break;
    case  7 : value=(sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+7.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-42.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-210.0*cos(c*x)*c*c*c*c*x*x*x*x+840.0*sin(c*x)*c*c*c*x*x*x+2520.0*cos(c*x)*c*c*x*x-5040.0*sin(c*x)*c*x-5040.0*cos(c*x))/(x*x*x*x*x*x*x*x); break;
    case  8 : value=(cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-8.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-56.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+336.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x+1680.0*cos(c*x)*c*c*c*c*x*x*x*x-6720.0*sin(c*x)*c*c*c*x*x*x-20160.0*cos(c*x)*c*c*x*x+40320.0*sin(c*x)*c*x+40320.0*cos(c*x))/(x*x*x*x*x*x*x*x*x); break;
    case  9 : value=-(sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+9.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-72.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-504.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+3024.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x+15120.0*cos(c*x)*c*c*c*c*x*x*x*x-60480.0*sin(c*x)*c*c*c*x*x*x-181440.0*cos(c*x)*c*c*x*x+362880.0*sin(c*x)*c*x+362880.0*cos(c*x))/pow(x,10.0); break;
    case 10 : value=-(cos(c*x)*pow(c,10.0)*pow(x,10.0)-10.0*sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-90.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+720.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+5040.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-30240.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-151200.0*cos(c*x)*c*c*c*c*x*x*x*x+604800.0*sin(c*x)*c*c*c*x*x*x+1814400.0*cos(c*x)*c*c*x*x-3628800.0*sin(c*x)*c*x-3628800.0*cos(c*x))/pow(x,11.0); break;
    case 11 : value=(sin(c*x)*pow(c,11.0)*pow(x,11.0)+11.0*cos(c*x)*pow(c,10.0)*pow(x,10.0)-110.0*sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x-990.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x+7920.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x+55440.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x-332640.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x-1663200.0*cos(c*x)*c*c*c*c*x*x*x*x+6652800.0*sin(c*x)*c*c*c*x*x*x+19958400.0*cos(c*x)*c*c*x*x-39916800.0*sin(c*x)*c*x-39916800.0*cos(c*x))/pow(x,12.0); break;
    case 12 : value=(cos(c*x)*pow(c,12.0)*pow(x,12.0)-12.0*sin(c*x)*pow(c,11.0)*pow(x,11.0)-132.0*cos(c*x)*pow(c,10.0)*pow(x,10.0)+1320.0*sin(c*x)*c*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x*x+11880.0*cos(c*x)*c*c*c*c*c*c*c*c*x*x*x*x*x*x*x*x-95040.0*sin(c*x)*c*c*c*c*c*c*c*x*x*x*x*x*x*x-665280.0*cos(c*x)*c*c*c*c*c*c*x*x*x*x*x*x+3991680.0*sin(c*x)*c*c*c*c*c*x*x*x*x*x+19958400.0*cos(c*x)*c*c*c*c*x*x*x*x-79833600.0*sin(c*x)*c*c*c*x*x*x-239500800.0*cos(c*x)*c*c*x*x+479001600.0*sin(c*x)*c*x+479001600.0*cos(c*x))/pow(x,13.0); break;
    default : value=0.0;
  }
  value*=pow(sign,der);

  return value;
}

complex double cot(double x, int der, const double *param)   /* K(x) = cot(cx) */
{
  double c=param[0];
  double value=0.0;

  if (fabs(x)<DBL_EPSILON) value=0.0;
  else switch (der)
  {
    case  0 : value = 1.0/tan(c * x); break;
    case  1 : value = -(1.0 + pow(1.0/tan(c * x), 2.0)) * c; break;
    case  2 : value = 2.0 * 1.0/tan(c * x) * (1.0 + pow(1.0/tan(c * x), 2.0)) * c * c; break;
    case  3 : value = -2.0 * (1.0 + pow(1.0/tan(c * x), 2.0)) * pow(c, 3.0) * (1.0 + 3.0 * pow(1.0/tan(c * x), 2.0)); break;
    case  4 : value = 8.0 * (1.0 + pow(1.0/tan(c * x), 2.0)) * pow(c, 4.0) * 1.0/tan(c * x) * (2.0 + 3.0 * pow(1.0/tan(c * x), 2.0)); break;
    case  5 : value = -0.8e1 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.5e1) * (0.15e2 * pow(1.0/tan(c * x), 0.2e1) + 0.15e2 * pow(1.0/tan(c * x), 0.4e1) + 0.2e1); break;
    case  6 : value = 0.16e2 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.6e1) * 1.0/tan(c * x) * (0.60e2 * pow(1.0/tan(c * x), 0.2e1) + 0.45e2 * pow(1.0/tan(c * x), 0.4e1) + 0.17e2); break;
    case  7 : value = -0.16e2 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.7e1) * (0.525e3 * pow(1.0/tan(c * x), 0.4e1) + 0.315e3 * pow(1.0/tan(c * x), 0.6e1) + 0.231e3 * pow(1.0/tan(c * x), 0.2e1) + 0.17e2); break;
    case  8 : value = 0.128e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.8e1) * 1.0/tan(c * x) * (0.630e3 * pow(1.0/tan(c * x), 0.4e1) + 0.315e3 * pow(1.0/tan(c * x), 0.6e1) + 0.378e3 * pow(1.0/tan(c * x), 0.2e1) + 0.62e2); break;
    case  9 : value = -0.128e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.9e1) * (0.6615e4 * pow(1.0/tan(c * x), 0.6e1) + 0.2835e4 * pow(1.0/tan(c * x), 0.8e1) + 0.5040e4 * pow(1.0/tan(c * x), 0.4e1) + 0.1320e4 * pow(1.0/tan(c * x), 0.2e1) + 0.62e2); break;
    case 10 : value = 0.256e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.10e2) * 1.0/tan(c * x) * (0.37800e5 * pow(1.0/tan(c * x), 0.6e1) + 0.14175e5 * pow(1.0/tan(c * x), 0.8e1) + 0.34965e5 * pow(1.0/tan(c * x), 0.4e1) + 0.12720e5 * pow(1.0/tan(c * x), 0.2e1) + 0.1382e4); break;
    case 11 : value = -0.256e3 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.11e2) * (0.467775e6 * pow(1.0/tan(c * x), 0.8e1) + 0.155925e6 * pow(1.0/tan(c * x), 0.10e2) + 0.509355e6 * pow(1.0/tan(c * x), 0.6e1) + 0.238425e6 * pow(1.0/tan(c * x), 0.4e1) + 0.42306e5 * pow(1.0/tan(c * x), 0.2e1) + 0.1382e4); break;
    case 12 : value = 0.1024e4 * (0.1e1 + pow(1.0/tan(c * x), 0.2e1)) * pow(c, 0.12e2) * 1.0/tan(c * x) * (0.1559250e7 * pow(1.0/tan(c * x), 0.8e1) + 0.467775e6 * pow(1.0/tan(c * x), 0.10e2) + 0.1954260e7 * pow(1.0/tan(c * x), 0.6e1) + 0.1121670e7 * pow(1.0/tan(c * x), 0.4e1) + 0.280731e6 * pow(1.0/tan(c * x), 0.2e1) + 0.21844e5); break;
    default : value=0.0;
  }

  return value;
}

/* \} */

/* kernels.c */
