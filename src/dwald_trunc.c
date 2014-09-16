#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double dwald_trunc(double t, double lambda, double alpha, double v, double d)
{
  double d;
  
  d = alpha * sqrt( lambda / (2 * pi * pow(t, 3) * (lambda * t * v + 1)) ) *
      1 / pnorm(d / sqrt(v), 0, 1, 1, 0) *
      exp( - (lambda * pow(d * t - alpha, 2)) / 
      (2 * t * (lambda * t * v + 1)) ) *
      pnorm( (lambda * alpha * v + d, 0, 1, 1, 0) / 
      (sqrt(lambda * t * pow(v, 2) + v) ) );

  return d;
}