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

#include "drlmath.h"		/* DrlCumNorm */

#include "drloptio.h"		/* prototype consistency */


#define	N(x)		DrlCumNorm(x)
#define	N1(x)		(exp(-0.5*(x)*(x))*0.3989422804)

/*f--------------------------------------------------------------
 * Options formulae : option with normal distribution.
 *
 * <br><br>
 * The argument <i> what</i> defines the type of sensitivity
 * to be computed: 
 * <i> "P"</i> for premium,
 * <i> "D"</i> for delta (${\partial v / \partial S}$),
 * <i> "G"</i> for gamma (${\partial^2 v / \partial S^2}$),
 * <i> "V"</i> for vega (${\partial v / \partial \sigma}$),
 * <i> "T"</i> for theta  ($-{\partial v / \partial t}$).
 * <br> Returns SUCCESS/FAILURE.
 */


DLL_EXPORT(int)
DrlNormOption(
	double texp,	/* (I) time to expiration */
	double p,	/* (I) expected underlying */
	double vol,	/* (I) normal volatility */
	double k,	/* (I) strike */
	char *cp,	/* (I) "C" for call, "P" for put */
	char *what,	/* (I) see below */
	double *retVal)	/* (O) output */
{
static	char	routine[] = "DrlNormOption";
	int	status  = FAILURE;
	double	v, d1;
	double	scaledVol;

	if (texp < 0e0) {
		GtoErrMsg("%s: time to expiration (%lf) < 0e0.\n",
			routine, texp);
		goto done;
	}

	if (vol < 0e0) {
		GtoErrMsg("%s: normal volatility (%lf) < 0e0.\n",
			routine, vol);
		goto done;
	}

	/* Floor scaled vol at some very small value just so we
	 * can use the same pricing/sensitivity equations without
	 * having to handle zero scaled vol as a special case. */
	scaledVol = MAX(vol * sqrt(texp), DBL_EPSILON);

	d1 = (p-k)/scaledVol;

	switch (cp[0]) {
	case 'C':
	case 'c':
		/*
		 * Call
		 */

		switch (what[0]) {
		case 'p': /* price */
		case 'P':
			v = (p-k)*N(d1) + scaledVol*N1(d1);
			break;
		case 'd': /* delta */
		case 'D':
			v = N(d1);
			break;
		case 'g': /* gamma */
		case 'G':
			v = N1(d1)/scaledVol;
			break;
		case 'v': /* vega */
		case 'V':
			v = sqrt(texp) * N1(d1);
			break;
		case 'T':
		case 't':
			v = vol * N1(d1) / MAX(2.0 * sqrt(texp), DBL_EPSILON);
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
			v = -(p-k)*N(-d1) + scaledVol*N1(-d1);
			break;
		case 'd': /* delta */
		case 'D':
			v = N(d1)-1;
			break;
		case 'g': /* gamma */
		case 'G':
			v = N1(d1)/scaledVol;
			break;
		case 'v': /* vega */
		case 'V':
			v = sqrt(texp) * N1(d1);
			break;
		case 'T':
		case 't':
			v = vol * N1(d1) / MAX(2.0 * sqrt(texp), DBL_EPSILON);
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
		GtoErrMsg("%s: failed.\n", routine);
	return(SUCCESS);
}


/*f--------------------------------------------------------------
 * Options formulae : normal distribution implied volatility.
 *
 * <br><br>
 * Computes implied volatity for normal options.
 * See <i> DrlNormOption</i>.
 * <br> Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DrlNormOptionImplVol(
	double texp,	/* (I) time to expiration */
	double fp,	/* (I) expected underlying */
	double prem,	/* (I) premium */
	double k,	/* (I) strike */
	char *cp,	/* (I) see DrlNormOption */
	char *what,	/* (I) see DrlNormOption */
	double *retVal)	/* (O) output (implied volatility) */
{
static	char	routine[] = "DrlNormOptionImplVol";
#define ITMAX	100
#define EPS	DBL_EPSILON
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define	BLPRICE(z)	(DrlNormOption(texp, fp, (z), k, cp, what, &value),\
				value - prem)

	double	x1 = 1e-6,
		x2 = 1.0,
		value,
		tol = DBL_EPSILON;

	/* patch routine here */
	int	iter;
	double	a=x1,b=x2,c=x2,d,e,min1,min2;
	double	fa = BLPRICE(a),
		fb = BLPRICE(b),
		fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		GtoErrMsg("%s: cannot solve.\n", routine);
		return(FAILURE);
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
			return(SUCCESS);
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
	GtoErrMsg("%s: too many iterations.\n", routine);
	return(FAILURE);

}
