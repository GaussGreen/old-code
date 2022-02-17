/* file: SciMClib.h */

#ifndef SCIMCLIB_H
#define SCIMCLIB_H

#if !defined(max)
#define max( a, b ) ((a) > (b) ? (a) : (b))
#endif
#if !defined(min)
#define min( a, b ) ((a) < (b) ? (a) : (b))
#endif

#include <sys/timeb.h>

// the random seed will depend on a seconds clock rather
// than a millisec clock if the following line is commented out 
#define MILLISECONDCLOCK	  

/* Uniform Random Deviates */
static void StdSetURDSeed(
	int newSeed,
	int *iURDseed
	);
static void SciSetURDSeed(
	int newSeed,
	int *iURDseed
	);
static double ranStd(
	int *iURDseed
	);
static double ran1(
	int *iURDseed
	);
static double ran2(
	int *iURDseed
	);
static double ran3(
	int *iURDseed
	);


/* Normal Deviates */
static double randNormBM(
	double (*URDfunc)(int *iURDseed),
	int *iURDseed
	);
static double randNormICNCore(
	double URD
	);
static double randNormICN(
	double (*URDfunc)(int *iURDseed),
	int *iURDseed
	);

static double randPoisson(
	double (*URDfunc)(int *iURDseed),
	int *iURDseed,
	double xm
	);
static double gammln(
	double xx
	);
static int bincoeff(
	int n,
	int k
	);

/**********************************************************************/

static void StdSetURDSeed(
	int newSeed,
	int *iURDseed
	)
{

#ifdef MILLISECONDCLOCK
	struct timeb tstruct;
	if (newSeed < 0) {
		ftime(&tstruct);
		*iURDseed = tstruct.millitm;
	}
#else
    time_t timer;
	if (newSeed < 0) {
		timer = time(0);
		*iURDseed = (int)timer;
	}
#endif
	else {
		if (newSeed == 0) {
			*iURDseed = 1;
		}
		else
		*iURDseed = abs(newSeed);
	}
	srand((unsigned)*iURDseed);
	return;
}

/**********************************************************************/

static void SciSetURDSeed(
	int newSeed,
	int *iURDseed
	)
{
#ifdef MILLISECONDCLOCK
	struct timeb tstruct;
	if (newSeed < 0) {
		ftime(&tstruct);
		*iURDseed = -1*tstruct.millitm;
	}
#else
    time_t timer;
	if (newSeed < 0) {
		timer = time(0);
		*iURDseed = -1*((int)timer);
	}
#endif
	else {
		if (newSeed ==0){
			*iURDseed = -1;
		}
		else
		*iURDseed = -1*abs(newSeed);
	}
	return;
}

/**********************************************************************/


static double ranStd(
	int *iURDseed
	)
{
	/* system standard random number generator */
	double rconst = 1./(double)RAND_MAX;

	return rconst*rand();
}

/**********************************************************************/


static double ran1(
	int *iURDseed
	)
{
	/* This routine, ran1, is based on the routine ran1 from the book 
	"Numerical recipes in C" (Cambridge University Press), Copyright (C)
	1986-1992 by Numerical Recipes Software. Used by permission. Use 
	of this routine other than as an integral part of SciMC requires 
	an additional license from Numerical Recipes Software. Further
	distribution in any form is prohibited */
	/* period = O(10^9) */

	#define IA 16807
	#define IM 2147483647
	#define AM (1.0/IM)
	#define IQ 127773
	#define IR 2836
	#define IR2 3791
	#define mNTAB 32
	#define mNDIV (1+(IM-1)/mNTAB)
	#define EPS 1.2e-7
	#define RNMX (1.0-EPS)

	int j;
	long k;
	static long iy = 0;
	static long iv[mNTAB];
	double temp;

	if(*iURDseed <= 0 || !iy){				/*Initializer*/
		if(-*iURDseed < 1) *iURDseed=1;		/*prevent idum=0*/
		else *iURDseed = -*iURDseed;
		for (j=mNTAB+7; j>=0; j--) {		/*8 warm ups then shuffle*/
			k = *iURDseed/IQ;
			*iURDseed=IA*(*iURDseed-k*IQ)-k*IR;
			if (*iURDseed<0) *iURDseed += IM;
			if (j<mNTAB) iv[j] = *iURDseed;
		}
		iy=iv[0];
	}
	k = *iURDseed/IQ;
	*iURDseed = IA*(*iURDseed-k*IQ)-k*IR;
	if (*iURDseed < 0) *iURDseed += IM;
	j = iy/mNDIV;
	iy = iv[j];
	iv[j] = *iURDseed;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

/**********************************************************************/


static double ran2(
	int *iURDseed
	)
{
	/* This routine, ran2, is based on the routine ran2 from the book 
	"Numerical recipes in C" (Cambridge University Press), Copyright (C)
	1986-1992 by Numerical Recipes Software. Used by permission. Use 
	of this routine other than as an integral part of SciMC requires 
	an additional license from Numerical Recipes Software. Further
	distribution in any form is prohibited */
	/* period = O(10^18); based on L'Ecuyer with shuffle */

	#define IM1 2147483563
	#define IM2 2147483399
	#define AM1 (1.0/IM1)
	#define IMM1 (IM1 -1)
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

	int j;
	long k;
	static long idum2=123456789;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if(*iURDseed <= 0){					/*Initializer*/
		if(-*iURDseed < 1) *iURDseed=1; /*prevent idum=0*/
		else *iURDseed = -*iURDseed;
		idum2 = *iURDseed;
		for (j=NTAB+7;j>=0;j--) {		/*8 warm ups then shuffle*/
			k = *iURDseed/IQ1;
			*iURDseed=IA1*(*iURDseed-k*IQ1)-k*IR1;
			if (*iURDseed < 0) *iURDseed += IM1;
			if (j < NTAB) iv[j] = *iURDseed;
		}
		iy=iv[0];
	}
	k = *iURDseed/IQ1;
	*iURDseed = IA1*(*iURDseed-k*IQ1)-k*IR1;
	if (*iURDseed < 0) *iURDseed += IM1;
	k = idum2/IQ2;
	idum2 = IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j = iy/NDIV;
	iy = iv[j]-idum2;
	iv[j] = *iURDseed;
	if (iy < 1) iy += IMM1;
	if ((temp=AM1*iy) > RNMX) return RNMX;
	else return temp;
}

/**********************************************************************/


static double ran3(
	int *iURDseed
	)
{
	/* This routine, ran3, is based on the routine ran3 from the book 
	"Numerical recipes in C" (Cambridge University Press), Copyright (C)
	1986-1992 by Numerical Recipes Software. Used by permission. Use 
	of this routine other than as an integral part of SciMC requires 
	an additional license from Numerical Recipes Software. Further
	distribution in any form is prohibited */
	/* For sequence period see Knuth vol II SS 3.2 */

	#define MBIG 1000000000
	#define MSEED 161803398
	#define MZ 0
	#define FAC (1.0/MBIG)

	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*iURDseed  <= 0 || iff == 0) {    /*Initialize*/
		iff=1;
		mj=labs(MSEED-labs(*iURDseed ));
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			if (mk<MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++) {
			for(i=1;i<=55;i++) {
				ma[i] -=ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		}
		inext=0;
		inextp=31;
		*iURDseed=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

/**********************************************************************/

static double randNormBM(
	double (*URDfunc)(int *iURDseed),
	int *iURDseed
	)
{
	/* normal distrubuted random deviates via Box-Mueller
       algorithm with reuse of extras. URDfunc is a pointer
       to a function returning random numbers uniformly
       distributed on the interval (0,1) */

	double r1,u1,u2,rad,q;
	static double r2;
	static int gotone=0;

	if (gotone==0) {
		rad = 2.;
		while (rad>1.) {
			u1 = 2.*URDfunc(iURDseed) - 1.;
			u2 = 2.*URDfunc(iURDseed) - 1.;
			rad = u1*u1+u2*u2;
		}
		q = sqrt(-2.*log(rad)/rad);
		r1 = u1*q;
		r2 = u2*q;
		gotone = 1;
		return r1;
	}
	else {
		gotone = 0;
		return r2;
	}
}


/**********************************************************************/

static double randNormICNCore(
 double URD
 )
{
    /* normal distrubuted random deviates via an
       Inverse Cumulative Normal Expansion.  The
       expansion coefficients are from the CNDEV
       algorithm by Moro and Boyle */

    double ax,y,r,x;
    static double a[4] = { 2.50662823884,   -18.61500062529,
                          41.39119773534,   -25.44106049637};
    static double b[4] = {-8.47351093090,    23.08336743743,
                         -21.06224101826,     3.13082909833};
    static double c[9] = {0.3374754822726147, 0.9761690190917186,
                          0.1607979714918209, 0.0276438810333863,
                          0.0038405729373609, 0.0003951896511919,
                          0.0000321767881768, 0.0000002888167364,
                          0.0000003960315187};

    x = URD - 0.5;
    ax = min(fabs(x),0.4999999999);
    if (ax < 0.38) {
        y = ax*ax;
        r = ax*(a[0]+y*(a[1]+y*(a[2]+y*a[3])))
                /(1.+y*(b[0]+y*(b[1]+y*(b[2]+y*b[3]))));
        }
    else {
        y = log(-log(0.5-ax));
        r = c[0]+y*(c[1]+y*(c[2]+y*(c[3]+y*(c[4]+
                 y*(c[5]+y*(c[6]+y*(c[7]+y*c[8])))))));
        }
    if (x<0) {r = -r;}
    return r;
}

/**********************************************************************/


static double randNormICN(
	double (*URDfunc)(int *iURDseed),
	int *iURDseed
	)
{
 /* normal distrubuted random deviates via an
    Inverse Cumulative Normal Expansion */

    double URD;

    URD = URDfunc(iURDseed);
    return randNormICNCore(URD);
}

/**********************************************************************/


static double randPoisson(
	double (*URDfunc)(int *iURDseed),
	int *iURDseed,
	double xm
	)
{

	/* Poisson distributed random deviates via Num Rec p.295 */
	/* Uses log of the gamma function */

	#define PI 3.1415926543
	static double sq, alxm, g, oldm=(-1.);
	double em,t,y;

	if (xm<12.0) {
		if (xm != oldm) {
			oldm = xm;
			g = exp(-xm);
		}
		em = -1.;
		t = 1.;
		do {++em; t *= URDfunc(iURDseed);} while (t>g);
		} 
	else {
		if (xm != oldm) {
			oldm = xm;
			sq = sqrt(2.*xm);
			alxm = log(xm);
			g = xm*alxm - gammln(xm+1.);
		}
	do {
		do {y = tan(PI*URDfunc(iURDseed)); em = sq*y+xm;} while (em<0.);
			em = floor(em);
			t = 0.9*(1. + y*y)*exp(em*alxm - gammln(em + 1.) - g);
		} while (URDfunc(iURDseed)>t);
	}
	return em;
}

/**********************************************************************/


static double gammln(
	double xx
	)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.120865097386617e-2,-0.5395239384953e-5};
	int j;

	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


static int bincoeff(
	int n,
	int k
	)
{
	/* binomial coeff C(n+k,k) */

	int i;
	double b=1.;

	for (i=1; i<=k; i++) {
		b = b*(n+i)/i;
	}
	return (int)b;
}

#endif
