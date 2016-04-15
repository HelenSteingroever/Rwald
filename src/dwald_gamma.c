// Driftrate as Gamma - Wald Distribution

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <gsl/gsl_specfunc.h>

int factorial(int n)
{
    if (n<=1)
        return(1);
    else
        n=n*factorial(n-1);
    return(n);
}

double r_gamma(double x)
{
  return gammafn(x); //exp(gamma(x));
}

double phi(double x)
{
  double c = 0;
    for (int i=1; i<=10; i++) {
    c = c + 1./i;
  }
  return(c);
}

double neg_gamma(double x)
{
    return(-M_PI / (-x*r_gamma(-x)*sin(M_PI*(-x))));
}

double neg_int_gamma(double x)
{
    return(pow(-1.,-x) / factorial(-x)) * (phi(-x) + digamma(1));
}

int is_int(double f)
{
	return floorf(f)==f;
}

double LaguerreL(double n, double a, double x) 
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

double erf(double x)
{
  double e;
  
  e = 2 * pnorm(x * sqrt(2.), 0,1, 1,0) - 1;
  
  return e;
}

double dwald_gamma_d_log(double t, double alpha, double tau, double kappa)
{
    double d;
    
    if (alpha == 1 && tau == 1 && kappa == 1 )
    {
        d = (- 1. / (2.*t)) - log(2.) -2*log(t);
    }
    else if (tau == 1)
    {
        d = log(alpha) - ((2.*alpha*kappa-1.)/(2.*t*pow(kappa,2))) + 
        log(1. + erf( (alpha*kappa-1.) / (kappa*sqrt(2.*t)) )) -
        log(2.) - 2*log(t) - log(kappa);
    }
    else
    {
        double L1; 
        double L2; 
        double L3; 
        double L4; 
        double C1; 
        double C2; 
                                                   
        L1 = LaguerreL(-(.5)*tau+.5, .5, (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
        
        L2 = LaguerreL(-(.5)*tau+.5, 3./2., (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
        
        L3 = LaguerreL(-(.5)*tau, .5, (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
        
        L4 = LaguerreL(-(.5)*tau, 3./2., (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
        
        C1 = sin((.5)*M_PI*tau)*r_gamma(-(.5)*tau+3./2.);
        
        C2 = sqrt(2.)*pow(alpha,3)*pow(kappa,3)*sqrt(t);
        
        d = log(-(1.0/16.0)*pow(2.,((.5)*tau+.5))*alpha*exp(-(.5)*pow(alpha,2)/t)*pow(kappa,(-tau-3.))*
        pow(t,(-(.5)*tau-7.0/2.0))*M_PI*
        
        (
         
         C1*
         L1*
         C2 +
         
         C1*
         L1*
         sqrt(2.)*alpha*pow(kappa,3)*tau*pow(t,(3./2.)) -
         
         C1*
         L2*
         C2 -
         
         C1*
         L1*
         sqrt(2.)*alpha*pow(kappa,3)*pow(t,(3./2.)) -
         
         3.*C1*
         L1*
         sqrt(2.)*pow(alpha,2)*pow(kappa,2)*sqrt(t) -
         
         C1*
         L1*
         sqrt(2.)*pow(kappa,2)*tau*pow(t,(3./2.)) +
         
         3.*C1*
         L2*
         sqrt(2.)*pow(alpha,2)*pow(kappa,2)*sqrt(t)-
         
         2.*
         L3*
         cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*pow(alpha,2)*pow(kappa,3)*t +
         
         2.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
         L4*
         pow(alpha,2)*pow(kappa,3)*t +
         
         C1*
         L1*
         sqrt(2.)*pow(kappa,2)*pow(t,(3./2.)) -
         
         2.*
         L3*
         cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*pow(kappa,3)*pow(t,2) +
         
         3.*C1*
         L1*
         sqrt(2.)*alpha*kappa*sqrt(t) -
         
         3.*C1*
         L2*
         sqrt(2.)*alpha*kappa*sqrt(t) +
         
         4.*
         L3*
         cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*alpha*pow(kappa,2)*t -
         
         4.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
         L4*
         alpha*pow(kappa,2)*t -
         
         C1*
         L1*
         sqrt(2.)*sqrt(t) +
         
         C1*
         L2*
         sqrt(2.)*sqrt(t) -
         
         2.*
         L3*
         cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*kappa*t +
         
         2.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
         L4*
         kappa*t
         
         ))-
        
        log(r_gamma(tau)) -log(C1) - log(cos((.5)*M_PI*tau)) - log(r_gamma(-(.5)*tau+2.));
    }
    return(d);
}

double dwald_gamma_d(double t, double alpha, double tau, double kappa, int give_log)
{
    double d;

    if(give_log) 
      d = dwald_gamma_d_log(t, alpha, tau, kappa);
    else
    {
      if (alpha == 1 && tau == 1 && kappa == 1 )
      {
          d = exp(- 1. / (2.*t)) / (2. * pow(t,2));
      }
      else if (tau == 1)
      {
          d = alpha*exp(-(2.*alpha*kappa-1.)/(2.*t*pow(kappa,2)))* (1. + erf( (alpha*kappa-1.) / (kappa*sqrt(2.*t)) )) / (2.*pow(t,2)*kappa);
      }
      else
      {
          double L1; 
          double L2; 
          double L3; 
          double L4;
          double C1; 
          double C2; 
                                                     
          L1 = LaguerreL(-(.5)*tau+.5, .5, (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
          
          L2 = LaguerreL(-(.5)*tau+.5, 3./2., (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
          
          L3 = LaguerreL(-(.5)*tau, .5, (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
          
          L4 = LaguerreL(-(.5)*tau, 3./2., (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
          
          C1 = sin((.5)*M_PI*tau)*r_gamma(-(.5)*tau+3./2.);
          
          C2 = sqrt(2.)*pow(alpha,3)*pow(kappa,3)*sqrt(t);
          
          d = -(1.0/16.0)*pow(2.,((.5)*tau+.5))*alpha*exp(-(.5)*pow(alpha,2)/t)*pow(kappa,(-tau-3.))*
          pow(t,(-(.5)*tau-7.0/2.0))*M_PI*
          
          (
           
           C1*
           L1*
           C2 +
           
           C1*
           L1*
           sqrt(2.)*alpha*pow(kappa,3)*tau*pow(t,(3./2.)) -
           
           C1*
           L2*
           C2 -
           
           C1*
           L1*
           sqrt(2.)*alpha*pow(kappa,3)*pow(t,(3./2.)) -
           
           3.*C1*
           L1*
           sqrt(2.)*pow(alpha,2)*pow(kappa,2)*sqrt(t) -
           
           C1*
           L1*
           sqrt(2.)*pow(kappa,2)*tau*pow(t,(3./2.)) +
           
           3.*C1*
           L2*
           sqrt(2.)*pow(alpha,2)*pow(kappa,2)*sqrt(t)-
           
           2.*
           L3*
           cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*pow(alpha,2)*pow(kappa,3)*t +
           
           2.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
           L4*
           pow(alpha,2)*pow(kappa,3)*t +
           
           C1*
           L1*
           sqrt(2.)*pow(kappa,2)*pow(t,(3./2.)) -
           
           2.*
           L3*
           cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*pow(kappa,3)*pow(t,2) +
           
           3.*C1*
           L1*
           sqrt(2.)*alpha*kappa*sqrt(t) -
           
           3.*C1*
           L2*
           sqrt(2.)*alpha*kappa*sqrt(t) +
           
           4.*
           L3*
           cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*alpha*pow(kappa,2)*t -
           
           4.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
           L4*
           alpha*pow(kappa,2)*t -
           
           C1*
           L1*
           sqrt(2.)*sqrt(t) +
           
           C1*
           L2*
           sqrt(2.)*sqrt(t) -
           
           2.*
           L3*
           cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*kappa*t +
           
           2.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
           L4*
           kappa*t
           
           )/
          
          (r_gamma(tau)*C1*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.));
      }
    }
    return(d);
}
     
SEXP  dwald_gamma(SEXP t, SEXP alpha, SEXP tau, SEXP kappa, SEXP give_log) {
  double d;
  SEXP value;
    
  d = dwald_gamma_d(REAL(t)[0], REAL(alpha)[0], REAL(tau)[0], REAL(kappa)[0], LOGICAL(give_log)[0]);
    
  PROTECT(value = allocVector(REALSXP, 1));
  REAL(value)[0] = d;
  UNPROTECT(1);
  return value;
}    
