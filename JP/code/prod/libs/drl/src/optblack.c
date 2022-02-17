/****************************************************************
 * Module:	DRL
 * Submodule:	OPTIO
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */
#include <math.h>	
#include <string.h>
#include <float.h>	

#include "drlmath.h"		/* DrlCumNorm */

#include "drloptio.h"		/* prototype consistency */


#undef	ISFLAG
#define	ISFLAG(b)		(!strncmp(what,b,strlen(b)))


/*f--------------------------------------------------------------
 * Options formulae : Black analytics.
 *
 *
 * <br><br>
 * Computes and returns price and sensitivities of a lognormal
 * option  E[max(+-(S - K), 0)].
 * The argument <i> what</i> defines the type of sensitivity
 * to be computed: 
 * <i> "P"</i> for premium,
 * <i> "D"</i> for delta
 * <i> "DV"</i> for d-del-v
 * <i> "DT"</i> for d-del-t,
 * <i> "G"</i> for gamma 
 * <i> "V"</i> for vega 
 * <i> "VV"</i> for wisoo
 * <i> "T"</i> for theta
 * <br> Return SUCCESS/FAILURE.
 */

int
DrlBlack(
	double texp,	/* (I) time to expiration */
	double p,	/* (I) exp. value of underlying */
	double vol,	/* (I) volatility */
	double k,	/* (I) strike */
	char *callPut,	/* (I) "C" for call, "P" for put */
	char *what,	/* (I) see below */
	double *retVal)	/* (O) return value */
{
static	char	routine[] = "DrlBlack";
	int	status = FAILURE;

	double	v, d1, d2;

#undef	N
#undef	N1
#define	N(x)		DrlCumNorm(x)
#define	N1(x)		(exp(-0.5*(x)*(x))*0.3989422804)

	if (texp < 0e0) {
		DrlErrMsg("%s: texp (%lf) < 0.\n", routine, texp);
		goto done;
	}
	if (p <= 0e0) {
		DrlErrMsg("%s: underlying (%lf) <= 0.\n", routine, p);
		goto done;
	}
	if (vol < 0e0) {
		DrlErrMsg("%s: volatility (%lf) < 0.\n", routine, vol);
		goto done;
	}
	if (k < 0e0) {
		DrlErrMsg("%s: strike (%lf) < 0.\n", routine, k);
		goto done;
	}
	if (texp == 0e0) {
	}

	d1 = (log(p/k) + 0.5 * vol * vol *texp) / (vol * sqrt(texp));
	d2 = d1 - vol * sqrt(texp);

	switch (callPut[0]) {
	case 'C':
	case 'c':
		/*
		 * Call
		 */

		switch (what[0]) {
		case 'p': /* price */
		case 'P':
			v = p * N(d1) - k * N(d2);
			break;
		case 'd': /* delta */
		case 'D':
			switch (what[1]) {
			case '\0':
				v = N(d1);
				break;
			case 'v':
			case 'V':
				v = - d2 * N1(d1) / vol;
				break;
			case 't':
			case 'T':
				v = - d2 * N1(d1) / (2.0 * texp);
				break;
			default:
				goto done;
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
			}
			break;
		case 'T':
		case 't':
			v = - p * N1(d1) * vol / (2.0 * sqrt(texp));
			break;
		default:
			goto done;
		}
		break;

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
		}
		break;
	default:
		goto done;
	}

	*retVal = v;
	status = SUCCESS;
done:
	if (status != SUCCESS)
		DrlErrMsg("%s: failed.\n", routine);
	return(status);
#undef	N
#undef	N1
}
/*f--------------------------------------------------------------
 * Options formulae : implied Black volatility.
 *
 * <br><br>
 * See <i> DrlBlack</i>.
 * <br> Return SUCCESS/FAILURE.
 */

int
DrlBlackImplVol(
	double texp,	/* (I) time to expiration (y) */
	double fwdp,	/* (I) future price */
	double prem,	/* (I) premium */
	double stk,	/* (I) strike */
	char *cp,	/* (I) "C" for call, "P" for put */
	char *what,	/* (I) see <i> Black</i> */
	double *retVal)	/* (O) return value */
{
#define ITMAX	100
#define EPS	DBL_EPSILON
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define	BLPRICE(z)	(DrlBlack(texp, fwdp, (z), stk, cp, what,\
			 &blackValue), blackValue - prem)

static	char	routine[] = "BlackIV";
	int	status = FAILURE;

	double	x1 = 1e-6,
		x2 = 1.0,
		tol = DBL_EPSILON,
		blackValue;

	/* patch routine here */
	int iter;
	double	a=x1,b=x2,c=x2,d,e,min1,min2;
	double	fa = BLPRICE(a),
		fb = BLPRICE(b),
		fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		goto done;
	}

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) {
			*retVal = b;
			status = SUCCESS;
			goto done;
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=BLPRICE(b);
	}
done:
	if (status != SUCCESS)
		DrlErrMsg("%s: failed with inputs "
			"te=%lf fw=%lf pr=%lf stk=%lf cp=%c what=%c.\n",
			routine, texp, fwdp, prem, stk, cp[0], what[0]);
	return (status);

#undef	ITMAX
#undef	EPS
#undef	SIGN
#undef	BLPRICE
}




/*f--------------------------------------------------------------
 * Options formulae : binary option analytics.
 *
 * <br><br>
 * Computes and returns price and sensitivities of a lognormal
 * binary option {S > K}>
 * The argument <i> what</i> defines the type of sensitivity
 * to be computed: 
 * <i> "P"</i> for premium,
 * <i> "D"</i> for delta,
 * <i> "G"</i> for gamma,
 * <i> "V"</i> for vega,
 * <i> "T"</i> for theta.
 * <br> Return SUCCESS/FAILURE.
 */

int
DrlBinary(
	double t,	/* (I) time to expiration */
	double s,	/* (I) expected value of underlying */
	double v,	/* (I) volatility */
	double k,	/* (I) strike */
	char *what,	/* (I) see below */
	double *retVal)	/* (O) return value */
{
static	char	routine[] = "DrlBinary";
	int	status = FAILURE;
	double	d1, d2;

#undef	N0
#undef	N1
#undef	N2
#define	N0(x)	DrlCumNorm(x)
#define	N1(x)	(exp(-(x)*(x)*0.5e0) * 0.398942280401e0)
#define	N2(x)	(-(x) * exp(-(x)*(x)*0.5e0) * 0.398942280401e0)

	/*
	 */
	if (t <= 1e-6) {
	    if (ISFLAG("P")) {
		*retVal = (s >= k ? 1. : 0.);
	    } else {
		*retVal = 0;
	    }
	} else {
	    /* t>1e-6 */

	    d1 = (log(s/k) + v*v*t*0.5) / (v * sqrt(t));
	    d2 = (log(s/k) - v*v*t*0.5) / (v * sqrt(t));

	    if (ISFLAG("P")) {
		*retVal = N0(d2);

	    } else if (ISFLAG("D")) {
		*retVal = N1(d2) / (s*v*sqrt(t));

	    } else if (ISFLAG("G")) {
		*retVal = - (N1(d2) / (s*v*sqrt(t))) * (1. + d2 / (v*sqrt(t))) / s;
	    } else if (ISFLAG("V")) {
		*retVal = - d1 * N1(d2) / v;
	    } else if (ISFLAG("T")) {
		*retVal = - d1 * N1(d2) * 0.5 / t;
	    } else {
		DrlErrMsg("%s: bad option type.\n", routine);
		goto done;
	    }
	}

	status = SUCCESS;
done:
	return (status);
#undef	N0
#undef	N1
#undef	N2
}


