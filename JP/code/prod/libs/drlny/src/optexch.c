/****************************************************************
 * Module:	DRL
 * Submodule:	OPTIO
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>	
#include <float.h>	
#include <ctype.h>	

#include "drlmath.h"		/* DrlCumNorm */

#include "drloptio.h"		/* prototype consistency */



/*f--------------------------------------------------------------
 * Options formulae : Margrabe spread option.
 *
 * <br><br>
 * Margrabe analytics:
 * computes and returns price and sensitivities of a lognormal
 * option
 * E[(+-(S_2 - S_1), 0)]
 * with $t_2&gt;= t_1$.
 * The argument <i> what</i> defines the type of sensitivity
 * to be computed: 
 * <i> "P"</i> for premium,
 * <i> "D1"</i> for delta w.r.t asset 1
 * <i> "D2"</i> for delta w.r.t asset 2
 */

DLL_EXPORT(int)
DrlMargrabe(
	double t1,	/* (I) time to expiration # 1*/
	double t2,	/* (I) time to expiration (t1 < t2) # 2 */
	double s1,	/* (I) underlying # 1 */
	double s2,	/* (I) underlying # 2 */
	double vol1,	/* (I) base volatility # 1*/
	double vol2,	/* (I) base volatility # 2*/
	double volf,	/* (I) fwd volatility # 2*/
	double rho,	/* (I) correlation */
	char *callPut,	/* (I) "C" for call, "P" for put */
	char *what,	/* (I) see below */
	double *retVal)	/* (O) output */
{
static	char	routine[] = "DrlMargrabe";
	int	status = FAILURE;
	double	volc, vol12, m, d1, d2, texp, vega, v;

#undef	N
#define	N(x)		DrlCumNorm(x)
#undef	N1
#define	N1(x)		(exp(-0.5*(x)*(x))*0.3989422804)
#undef	SQR
#define	SQR(x)		((x)*(x))

	ASSERT_OR_DONE(t2 >= 0e0);
	ASSERT_OR_DONE(t2 >= t1);
	ASSERT_OR_DONE(s1 > 0e0);
	ASSERT_OR_DONE(s2 > 0e0);
	ASSERT_OR_DONE(vol1 > 0e0);
	ASSERT_OR_DONE(vol2 > 0e0);
	ASSERT_OR_DONE(volf > 0e0);
	ASSERT_OR_DONE((rho >= -1e0) && (rho <= 1e0));



	if (t1 > 0e0) {
	    /**
	     **  t1 > 0
	     **/

	    vol12 = sqrt(SQR(vol2) - SQR(volf) * (t2 - t1)/t2);
	    volc = sqrt((SQR(vol1)*t1 -
		2e0 * rho * vol1 * sqrt(t1*(SQR(vol2)*t2 - SQR(volf)*(t2-t1)))
		+ SQR(vol2)*t2) / t2);
	    texp = t2;

	    d1 = (log(s2/s1) + 5e-1 * SQR(volc) * texp) / (volc * sqrt(texp));
	    d2 = d1 - volc * sqrt(texp);

	    switch (toupper(callPut[0])) {
	    case 'C':
		/*
		 * Call
		 */
		switch (toupper(what[0])) {
		case 'P':
		    v = s2 * N(d1) - s1 * N(d2);
		    break;
		case 'D':
		    switch (toupper(what[1])) {
			case '1':
			    v = N(d1);
			    break;
			case '2':
			    v = - N(d2);
			    break;
			default:
			    GtoErrMsg("%s: bad arg `%s'\n", routine, what);
			    goto done;
		    }
		    break;

		case 'G':
		    switch (toupper(what[1])) {
			case '1':
			    v = N1(d1) / (s2 * volc * sqrt(texp));
			    break;
			case '2':
			    v = N1(d2) / (s1 * volc * sqrt(texp));
			    break;
			case '3':
			    v = - (s2/s1) * N1(d1) / (s2 * volc * sqrt(texp));
			    break;
			default:
			    GtoErrMsg("%s: bad arg `%s'\n", routine, what);
			    goto done;
		    }
		    break;

		case 'V':
		    vega =  s2 * sqrt(texp) * N1(d1);
		    switch (toupper(what[1])) {
			case '1':
			    v = (vol1 * t1 + rho * sqrt(t1*t2) * vol12)
				/ volc;
			    break;
			case '2':
			    v = (vol2 * t2 +
				rho * sqrt(t1*t2) * vol12 * vol1 / vol2)
				/ volc;
			    break;
			case 'R':
			    v = (sqrt(t1*t2) * vol1 * vol12 / vol2)
				/ volc;
			    break;
			default:
			    GtoErrMsg("%s: bad arg `%s'\n", routine, what);
			    goto done;
		    }
		    break;

		default:
		    GtoErrMsg("%s: bad arg `%s'\n", routine, what);
		    goto done;
		}
		break;

	    default:
		GtoErrMsg("%s: bad option type\n", routine);
		goto done;
	    }
	} else {
	    /**
	     **  t1 <= 0
	     **/
	    texp = t2;
	    volc = vol2;
	    d1 = (log(s2/s1) + 5e-1 * SQR(volc) * texp) / (volc * sqrt(texp));
	    d2 = d1 - volc * sqrt(texp);

	    switch (toupper(what[0])) {
	    case 'C':
	    case 'c':
		/*
		 * Call
		 */

		switch (what[0]) {
		case 'p': /* price */
		case 'P':
			v = s2 * N(d1) - s1 * N(d2);
			break;
		case 'd': /* delta */
		case 'D':
			switch (what[1]) {
			case '\0':
				v = N(d1);
				break;
			case 'v':
			case 'V':
				v = - d2 * N1(d1) / volc;
				break;
			case 't':
			case 'T':
				v = - d2 * N1(d1) / (2.0 * texp);
				break;
			default:
				goto done;
				break;
			}
			break;
		case 'g': /* gamma */
		case 'G':
			v = N1(d1) / (s2 * volc * sqrt(texp));
			break;
		case 'v': /* vega */
		case 'V':
			switch (what[1]) {
			case '\0':
				v = s2 * sqrt(texp) * N1(d1);
				break;
			case 'v': /* wisoo */
			case 'V':
				v = s2 * sqrt(texp) * d1 * d2 * N1(d1) / volc;
				break;
			default:
				goto done;
				break;
			}
			break;
		case 'T':
		case 't':
			v = - s2 * N1(d1) * volc / (2.0 * sqrt(texp));
			break;
		default:
			goto done;
			break;
		}
		break;
#ifdef	_SKIP
	    case 'P':
	    case 'p':
		/*
		 * Put
		 */

		switch (what[0]) {
		case 'p': /* price */
		case 'P':
			v = k * N(-d2) - p * N(-d1);
			break;
		case 'd': /* delta */
		case 'D':
			switch (what[1]) {
			case '\0':
				v = N(d1)-1;
				break;
			case 'v':
			case 'V':
				v = - d2 * N1(d1) / vol;
				break;
			case 't':
			case 'T':
				v = - d2 * N1(d1) / (2.0 * texp);
				break;
			}
			break;
		case 'g': /* gamma */
		case 'G':
			v = N1(d1) / (p * vol * sqrt(texp));
			break;
		case 'v': /* vega */
		case 'V':
			switch (what[1]) {
			case '\0':
				v = p * sqrt(texp) * N1(d1);
				break;
			case 'v': /* wisoo */
			case 'V':
				v = p * sqrt(texp) * d1 * d2 * N1(d1) / vol;
				break;
			default:
				goto done;
				break;
			}
			break;
		case 'T':
		case 't':
			v = - p * N1(d1) * vol / (2.0 * sqrt(texp));
			break;
		case 'Z':
		case 'z':
			v = - d2 * N1(d1) / (2.0 * texp);
			break;
		default:
			goto done;
			break;
		}
		break;
#endif
	    default:
		goto done;
	    }
	}

	*retVal = v;
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
#undef	N
#undef	N1
#undef	SQR
}



#undef	N
#define	N(x)		DrlCumNorm(x)
static	const	double	twopi = 6.28318530;
static	double	m1, m2, s1, s2, rho, stk;
static	double	iab(double a, double b, double stk);
static	double	G(double y);
static	double	trapzd(double (*func)(double), double a, double b, int n);
static	double	qtrap(double (*func)(double), double a, double b);
static	double	drlSpreadOptionR(double texp, double p1, double p2, double vol1,
			double vol2, double cor, double strike, char *callPut);

/*f--------------------------------------------------------------
 * Options formulae : spread option (with non zero strike).
 *
 * <br><br>
 * Option on the forward spread of two assets.
 */

DLL_EXPORT(int)
DrlSpreadOption(
	double texp,	/* (I) time to exp. */
	double p1,	/* (I) asset 1 */
	double p2,	/* (I) asset 2 */
	double vol1,	/* (I) vol 1 */
	double vol2,	/* (I) vol 2 */
	double cor,	/* (I) correlation */
	double strike,	/* (I) strike */
	char *callPut,	/* (I) see below */
	char *what,	/* (I) see below */
	double *retVal)	/* (O) output */
{
static	char	routine[] = "DrlSpreadOption";
	int	status = FAILURE;

	double	xu, xd, xm, dx, vu, vd, vm, v;

	switch (what[0]) {
	case 'p': /* price */
	case 'P':
		v = drlSpreadOptionR(texp, p1, p2, vol1, vol2, cor, strike, callPut);
		break;
	case 'd': /* delta */
	case 'D':
		switch (what[1]) {
		case '1':
			dx = p1 * 1e-3;
			xd = p1 - dx;
			if (xd <= 1e-4) xd = 1e-4;
			xm = xd + dx;
			xu = xm + dx;
			vd = drlSpreadOptionR(texp, xd, p2, vol1, vol2, cor, strike, callPut);
			vu = drlSpreadOptionR(texp, xu, p2, vol1, vol2, cor, strike, callPut);
			v = (vu - vd) / (xu - xd);
			break;
		case '2':
			dx = p2 * 1e-3;
			xd = p2 - dx;
			if (xd <= 1e-4) xd = 1e-4;
			xm = xd + dx;
			xu = xm + dx;
			vd = drlSpreadOptionR(texp, p1, xd, vol1, vol2, cor, strike, callPut);
			vu = drlSpreadOptionR(texp, p1, xu, vol1, vol2, cor, strike, callPut);
			v = (vu - vd) / (xu - xd);
			break;
		default:
			goto done;
		}
		break;
	case 'g': /* gamma */
	case 'G':
		switch (what[1]) {
		case '1':
			dx = p1 * 1e-3;
			xd = p1 - dx;
			if (xd <= 1e-4) xd = 1e-4;
			xm = xd + dx;
			xu = xm + dx;
			vd = drlSpreadOptionR(texp, xd, p2, vol1, vol2, cor, strike, callPut);
			vm = drlSpreadOptionR(texp, xm, p2, vol1, vol2, cor, strike, callPut);
			vu = drlSpreadOptionR(texp, xu, p2, vol1, vol2, cor, strike, callPut);
			v = (vu - 2e0 * vm + vd) / (dx * dx);
			break;
		case '2':
			dx = p2 * 1e-3;
			xd = p2 - dx;
			if (xd <= 1e-4) xd = 1e-4;
			xm = xd + dx;
			xu = xm + dx;
			vd = drlSpreadOptionR(texp, p1, xd, vol1, vol2, cor, strike, callPut);
			vm = drlSpreadOptionR(texp, p1, xm, vol1, vol2, cor, strike, callPut);
			vu = drlSpreadOptionR(texp, p1, xu, vol1, vol2, cor, strike, callPut);
			v = (vu - 2e0 * vm + vd) / (dx * dx);
			break;
		case '3':
			v  = drlSpreadOptionR(texp, p1,           p2,           vol1, vol2, cor, strike, callPut);
			vu = drlSpreadOptionR(texp, p1*(1.+1e-3), p2,           vol1, vol2, cor, strike, callPut);
			vd = drlSpreadOptionR(texp, p1,           p2*(1.+1e-3), vol1, vol2, cor, strike, callPut);
			vm = drlSpreadOptionR(texp, p1*(1.+1e-3), p2*(1.+1e-3), vol1, vol2, cor, strike, callPut);
			v = (vm + v - vu - vd) / (p1 * p2 * 1e-6);
			break;
		default:
			goto done;
		}
		break;
	case 'v': /* vega */
	case 'V':
		switch (what[1]) {
		case '1':
			dx = 1e-4;
			xd = vol1;
			xu = vol1 + dx;
			vd = drlSpreadOptionR(texp, p1, p2, xd, vol2, cor, strike, callPut);
			vu = drlSpreadOptionR(texp, p1, p2, xu, vol2, cor, strike, callPut);
			v = (vu - vd) / dx;
			break;
		case '2':
			dx = 1e-4;
			xd = vol2;
			xu = vol2 + dx;
			vd = drlSpreadOptionR(texp, p1, p2, vol1, xd, cor, strike, callPut);
			vu = drlSpreadOptionR(texp, p1, p2, vol1, xu, cor, strike, callPut);
			v = (vu - vd) / dx;
			break;
		case '3':	/* rho */
			dx = 1e-2;
			xm = rho;
			if (xm >=  .98) xm =  .98;
			if (xm <= -.98) xm = -.98;
			xd = xm - dx;
			xu = xm + dx;
			vd = drlSpreadOptionR(texp, p1, p2, vol1, vol2, xd, strike, callPut);
			vu = drlSpreadOptionR(texp, p1, p2, vol1, vol2, xu, strike, callPut);
			v = (vu - vd) / (xu - xd);
			break;
		default:
			goto done;
		}
		break;
	case 't': /* theta */
	case 'T':
		vd = drlSpreadOptionR(texp          , p1, p2, vol1, vol2, cor, strike, callPut);
		vu = drlSpreadOptionR(texp+1e-4, p1, p2, vol1, vol2, cor, strike, callPut);
		v = (vu - vd) / 1e-4;
		break;
	default:
		goto done;
	}

	*retVal = v;
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}

/*---------------------------------------------------------------
 * Quasi-Exact Spread Option Formula
 */

static	double
drlSpreadOptionR(
	double texp,
	double p1,
	double p2,
	double vol1,
	double vol2,
	double cor,
	double strike,
	char *callPut)
{
static	char	routine[] = "drlSpreadOptionR";
	double	u, ymax;
	int	m;


	ymax = 5.0;

	/* test for maturity = 0 */
	if (texp <  1e-8) {
		switch (callPut[0]) {
		case 'C':
		case 'c':
			return(MAX(p1-p2-strike, 0.));
		case 'P':
		case 'p':
			return(MAX(-(p1-p2-strike), 0.));
		default:
		GtoErrMsg("%s: bad option type.\n", routine);
			return DBL_MAX;
		}
	}


	switch (callPut[0]) {
	case 'C':
	case 'c':
		/* call */
		m1 = log(p1) - vol1*vol1*0.5*texp;
		m2 = log(p2) - vol2*vol2*0.5*texp;
		s1 = vol1 * sqrt(texp);
		s2 = vol2 * sqrt(texp);
		rho = cor;
		stk = strike;
		break;
	case 'P':
	case 'p':
		/* put */
		m1 = log(p2) - vol2*vol2*0.5*texp;
		m2 = log(p1) - vol1*vol1*0.5*texp;
		s1 = vol2 * sqrt(texp);
		s2 = vol1 * sqrt(texp);
		rho = cor;
		stk = -strike;
		break;
	default:
		GtoErrMsg("%s: bad option type.\n", routine);
		return DBL_MAX;


	}

	/*
	 * Computes E[(exp(m1+s1*X) - exp(m2+s2*Y) - stk)+]
	 * where (X,Y) are std normal with correlation rho.
	 */
	/*for (m=1; m<=8; m++) {
		u = trapzd(G, -ymax, ymax, m);
		u /= sqrt(twopi);
		printf("m=%2d   u=%10.6f\n", m, u);
	}*/
	u = qtrap(G, -ymax, ymax);
	u /= sqrt(twopi);



	return(u);
}


/*---------------------------------------------------------------
 *
 */

static double
iab(double a, double b, double E)
{
	double	iabval;

	double	d1, d2;

	if (E >= 0) {
		d1 = ((b - log(E)) / a) + a;
		d2 = (b - log(E)) / a;
		iabval =  exp(a*a*0.5+b)*N(d1) - E*N(d2);
	} else {
		iabval =  exp(a*a*0.5 + b) - E;
	}

	return iabval;
}


static double
G(double y)
{
	return iab(s1*sqrt(1.-rho*rho), m1+rho*s1*y, exp(m2+s2*y)+stk)
			* exp(-y*y*0.5);
}


/*---------------------------------------------------------------
 *
 */


#define FUNC(x) ((*func)(x))
static double
trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
		/*printf("trapz: n=%2d  s=%10.6f\n", n, s);*/
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		/*printf("trapz: n=%2d  s=%10.6f\n", n, s);*/
		return s;
	}
}
#undef FUNC

#define	EPS	1.0e-7
#define	JMAX	9

static double
qtrap(double (*func)(double), double a, double b)
{
	int	j;
	double	s, olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s = trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	return s;
}
#undef EPS
#undef JMAX
#undef	N
