#include <math.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "util/numerical_random.h"

#define MAX_FACT_FOR_DOUBLE 170 //max allowed for fact of double
double G_factorial[MAX_FACT_FOR_DOUBLE+1];


float gammln(float xx)
//Returns the value ln[Γ(xx)] for xx > 0.
{
    //Internal arithmetic will be done in double precision, a nicety that you can omit if ﬁve-ﬁgure
    //accuracy is good enough.
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
         24.01409824083091,-1.231739572450155,
         0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

float poissonHighPrecision(float xm)
// Returns as a ﬂoating-point number an integer value that is a random deviate drawn from a
// Poisson distribution of mean xm, using ran1(idum) as a source of uniform random deviates.
{
    //cout << xm << endl;
    float gammln(float xx);
    float ran2();
    static float sq,alxm,g,oldm=(-1.0);               //oldm is a ﬂag for whether xm has changed
    float em,t,y;                                     //     since last call.
    if (xm < 12.0) {                         //Use direct method.
         if (xm != oldm) {
              oldm=xm;
              g=exp(-xm);                    //If xm is new, compute the exponential.
         }
         em = -1;
         t=1.0;
         do {                               //Instead of adding exponential deviates it is equiv-
              ++em;                          //     alent to multiply uniform deviates. We never
              t *= ran2();               //     actually have to take the log, merely com-
         } while (t > g);                   //      pare to the pre-computed exponential.
    } else {                               //  Use rejection method.
         if (xm != oldm) {                  // If xm has changed since the last call, then pre-
              oldm=xm;                        //    compute some functions that occur below.
              sq=sqrt(2.0*xm);
              alxm=log(xm);
              g=xm*alxm-gammln(xm+1.0);
              //The function gammln is the natural log of the gamma function, as given in §6.1.
         }
         do {
              do {                          // y is a deviate from a Lorentzian comparison func-
                   y=tan(PI*ran2());    //      tion.
               em=sq*y+xm;                //em is y, shifted and scaled.
          } while (em < 0.0);             //Reject if in regime of zero probability.
          em=floor(em);                   //The trick for integer-valued distributions.
          t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
          //The ratio of the desired distribution to the comparison function; we accept or
          //reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
          //that t never exceeds 1.
      } while (ran2() > t);
  }
  return em;
}

//**********************************************************
// LOW PRECISION POISSON
//**********************************************************
//**********************************************************

void inline allocFactorial(){
	G_factorial[0] = 1;
	for(int i=1;i<=(MAX_FACT_FOR_DOUBLE);i++)
		G_factorial[i] = G_factorial[i-1] * i;

// 	for(int i=0;i<=(MAX_FACT_FOR_DOUBLE);i++)
// 		cout << i << " " << G_factorial[i] << endl;
}


int poissonLowPrecision(double mut)
//does not support lamda>1930
{
  int m;
  double soma, x;

  soma = 0;
  m = 0;
  x = ran2();
  while( soma<x )
    {
      if( m<=60 ){
	soma += pow(mut,(double)m)*exp(-mut)/G_factorial[m];

	}
	if( m>60 ){
		soma += exp(-mut+(double)m)*pow((mut/(double)m),(double)m)*pow((2*(double)m*PI),-0.5);
	}
       if( m>2000 )
 	soma = 1.0;
      m++;
    }
    m--;
  return(m);
}

