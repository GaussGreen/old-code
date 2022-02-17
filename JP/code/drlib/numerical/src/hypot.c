#include "imsl_inc.h"
/* Part of MS CRT - get error when linking if this is here due to duplicate
   symbols */
#if !defined( _MSC_VER )

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
#define ABS2(_A) ((_A) < 0.0 ? -(_A) : (_A))
#define MAX2(_A, _B) ((_A) < (_B) ? (_B) : (_A))
#define MIN2(_A), _B) ((_A) > (_B) ? (_B) : (_A))
*/

double hypot(double a, double b)
{
	double p,r,s,t,u;
	if (fabs(a) > fabs(b)) 
	   p = fabs(a);
	else
	   p = fabs(b);
	if ( p == 0.0) return(p);
/*	r = (MIN2(ABS2(a),ABS2(b))/p); */
	if (fabs(a) < fabs(b))
	   r = fabs(a)/p;
	else
	   r = fabs(b)/p;
/*	printf("r=%g\n",r); */
	r *= r;
	while (1) {
		t = 4.0 + r;
/*		printf(" ... t=%g\n",t); */
		if (t == 4.0) return(p);
		s = r/t;
		u = 1.0 + 2.0*s;
		p *= u;
		r = (s/u)*(s/u)*r;
	}
}
#endif
