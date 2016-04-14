// Driftrate as Gamma - Wald Distribution

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <gsl/gsl_specfunc.h>

double laguerrel_d(double n, double a, double x)
{
  double Lf;
  double c1;
  double c2;
  double c3;

    if (n+a+1 <= 0)
    {
        if (n+a+1 == 0)
        {
            c1 = digamma(1);
        } else {
            if (is_int(n+a+1))
            {
                c1 = neg_int_gamma(n+a+1);
            } else
            {
                c1 = neg_gamma(n+a+1);
            }
        }
    }else {
        c1 = r_gamma(n+a+1);
    }
    
    if (n+1 <= 0)
    {
        if (n+1 == 0)
        {
            c2 = digamma(1);
        } else
        {
            if (is_int(n+1))
            {
                c2 = neg_int_gamma(n+1);
            } else
            {
                c2 = neg_gamma(n+1);
            }
        }
    }else
    {
        c2 = r_gamma(n+1);
    }
    
    if (a+1 <= 0) {
        if (a+1 == 0) {
            c3 = digamma(1);
        } else {
            if (is_int(a+1)) {
                c3 = neg_int_gamma(a+1);
            } else {
                c3 = neg_gamma(a+1);
            }
        }
    }else {
        c3 = r_gamma(a+1);
    }

  Lf = c1/(c2*c3) * gsl_sf_hyperg_1F1(-n,a+1,x);

  return Lf;
}

SEXP  laguerrel(SEXP n, SEXP a, SEXP x) {
  double d;
  SEXP value;
    
  d = laguerrel_d(REAL(n)[0], REAL(a)[0], REAL(x)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}    
