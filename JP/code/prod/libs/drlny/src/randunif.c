/************************************************************************
 * Module:	DRL
 * Submodule:	RAND
 * Function:	Random Number Generation and Simulation
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <float.h>
#include <stdlib.h>


#if defined(CLIB)
#include "goodrnd.h"
#else
#endif

#include "drlrand.h"


static	double	drlRan2(long *idum);

/*f---------------------------------------------------------------------
 * Simulation of uniform distribution on an interval [0,1].
 * "idum" is the seed of the random number generator and should not
 * be change between successive calls (can be NULL, if which case a default
 * of 0 is used).
 */

DLL_EXPORT(int)
DrlRandUnif(long *idum, double *rnd)
{
#if defined (CLIB_RAN)
static	long	idum1 = -1L;
	idum = (idum != NULL ? idum : &idum1);
	return GtoRandomUniform((idum != NULL ? idum : &idum1), rnd);
#else
	*rnd = drlRan2(idum);
	return(SUCCESS);
#endif
}



/*f---------------------------------------------------------------------
 * Simulation of uniform distribution on an interval [x_low,x_high].
 * "idum" is the seed of the random number generator and should not
 * be change between successive calls (can be NULL, if which case a default
 * of 0 is used).
 */

DLL_EXPORT(double)
DrlDoubleRand(long *idum, double xLow, double xHigh)
{
static	long	idum1;
static	double	rnd;
	DrlRandUnif((idum != NULL ? idum : &idum1), &rnd);
	return ((xLow) + ((xHigh)-(xLow))*((double) rnd));
}



/*f---------------------------------------------------------------------
 * Returns gaussian deviates. The seed "idum" should not
 * be changed between consecutive function calls.
 */

DLL_EXPORT(double)
DrlGaussSimul(long *idum)
{
#ifdef	CLIB
static	double	rnd;
static	int iset=0;
static	double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			DrlRandUnif(idum, &rnd);
			v1=2.0*rnd-1.0;
			DrlRandUnif(idum, &rnd);
			v2=2.0*rnd-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
#endif
}


/*f---------------------------------------------------------------------
 * Simulation of uniform distribution on an interval [x_low,x_high].
 * "idum" is the seed of the random number generator and should not
 * be change between successive calls (can be NULL).
 */

DLL_EXPORT(int)
DrlIntRand(long *idum, int iLow, int iHigh)
{
static	long	idum1;
static	double	rnd;
	DrlRandUnif((idum != NULL ? idum : &idum1), &rnd);
	return (int) floor (((iLow) + ((iHigh)-(iLow))*((double) rnd)) + 0.5);
}



/*----------------------------------------------------------------------
 * Private copy of NRC routine.
 */



#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static	double
drlRan2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



