/* ==========================================================================
   FILE_NAME:	num_f_pythag.c

   PURPOSE:     returns a^2 + b^2 preventing numerical overflows   
   ========================================================================== */

#include        "utallhdr.h"
#include		<math.h"

#define NRANSI

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software koV219. */
