#include <R.h> 
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>


double F77_SUB(cpnorm)(q, mu, sigma)
double *q, *mu, *sigma;
{
double p;
int one, zero;
one=1;
zero=0;
p = pnorm(*q, *mu, *sigma, one, zero);
return p;
}

double F77_SUB(cqnorm)(p, mu, sigma)
double *p, *mu, *sigma;
{
double q;
int one, zero;
one=1;
zero=0;
q = qnorm(*p, *mu, *sigma, one, zero);
return q;
}

/*
void F77_SUB(setseed)(iseed1,iseed2)
unsigned int *iseed1, *iseed2;
{
set_seed(*iseed1, *iseed2);
}
*/
void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(crunif)(void) { return unif_rand(); }

double F77_SUB(crnorm)(void) { return norm_rand(); } /* added subroutine */

/* Add the Rprintf interface to Fortran */
void F77_SUB(ifprintf)(i, flag, cr)
int *i, *flag, *cr;
{ 
  if ((*i)>0) {
    if ((*flag)==0) {
	  if ((*cr)==0) {
        Rprintf("Working on iteration %d...", *i);
	  } else {
        Rprintf("Working on iteration %d...\n", *i);
      }	  
	} else {
	  if ((*cr)==0) {
        Rprintf("  working on subject %d...", *i);
	  } else {
        Rprintf("  working on subject %d...\n", *i);
      }
	}
  } else {
    Rprintf("\n");
  }
}
