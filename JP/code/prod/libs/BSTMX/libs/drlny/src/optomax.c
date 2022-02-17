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
#include <ctype.h>
#include <string.h>

#include "nvarnorm.h"		/* C Analytics */
#include "drlmath.h"		/* DrlCumNorm */

#include "drloptio.h"		/* prototype consistency */


static	double	D2(double z, double k, double s, double t)
{
	return ((log(z/k) - .5*s*s*t) / (s*sqrt(t)));
}

/* rounds to zero */
#define	RZERO(x)	{x = (fabs(x) <= 1e-8 ? 1e-8 : x);}


/*f--------------------------------------------------------------
 * Options formulae : max/min of two assets.
 *
 * <br><br>
 * Computes the option on the maximum or minimum of two
 * assets.
 * The argument <i> cp</i> defines the type of option
 * to be computed: 
 * <i> "CA"</i> for call on max
 * ($E\left[ \max(\max(S_1,S_2) - K), 0) \right]$),
 * <i> "CB"</i> for call on min
 * ($E\left[ \max(\min(S_1,S_2) - K), 0) \right]$),
 * <i> "PA"</i> for put on max
 * ($E\left[ \max(K - \max(S_1,S_2)), 0) \right]$),
 * <i> "PB"</i> for put on min
 * ($E\left[ \max(K - \min(S_1,S_2)), 0) \right]$).
 * <br> Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DrlOpMax2Sec(
	double t,	/* (I) time to expiration */
	double s1,	/* (I) exp asset 1 */
	double s2,	/* (I) exp asset 2 */
	double sigma1,	/* (I) vol of asset 1 */
	double sigma2,	/* (I) vol of asset 2 */
	double rho,	/* (I) correlation bet $S_1$ and $S_2$ */
	double k,	/* (I) strike */
	char *cp,	/* (I) see below */
	char *what,	/* (I) "P" for premium (only available) */
	double *retVal)	/* (O) output */
{
static	char	routine[] = "DrlOpMax2Sec";
	int	status = FAILURE;

	char	cType, mType;
	double	v1, v2, v12, r12, r112, r212;

#define	N(x)		DrlCumNorm(x)
#define	N1(x)		(exp(-0.5*(x)*(x))*0.3989422804)
#define	N2(a,b,cor)	DrlCumBiNorm(a,b,cor)	

	/* */

	if ((cp[0] == '\0') || (cp[0] == '\0')) {
		GtoErrMsg("%s: bad option type.\n", routine);
		goto done;
	}

	cType = toupper(cp[0]);
	mType = toupper(cp[1]);

	/*
	 *
	 */

	RZERO(t);
	RZERO(s1);
	RZERO(s2);
	RZERO(sigma1);
	RZERO(sigma2);
	RZERO(k);
	rho = (fabs(1.-rho) <= 1e-8 ?  1.-1e-8 : rho);
	rho = (fabs(1.+rho) <= 1e-8 ? -1.+1e-8 : rho);

	/*
	 *
	 */
	v1 = sigma1;
	v2 = sigma2;
	r12 = rho;
	v12	= sqrt(v1*v1 + v2*v2 - 2.*r12*v1*v2);
	r112	= (v1 - r12*v2) / v12;
	r212	= (v2 - r12*v1) / v12;


	if ((cType == 'C') && (mType == 'A')) {
	/*
	 * Call on maximum
	 */
	*retVal  = - k * (1. - N2(-D2(s1,k,v1,t), -D2(s2,k,v2,t), r12));
	*retVal += + s1 * N2(-D2(s2/s1, 1., v12, t), -D2(k/s1, 1., v1, t), r112);
	*retVal += + s2 * N2(-D2(s1/s2, 1., v12, t), -D2(k/s2, 1., v2, t), r212);

	} else if ((cType == 'C') && (mType == 'B')) {
	/*
	 * Call on minimum
	 */
	*retVal  = - k * N2(D2(s1,k,v1,t), D2(s2,k,v2,t), r12);
	*retVal += s1 * N2(D2(s2/s1, 1., v12, t), -D2(k/s1, 1., v1, t), -r112);
	*retVal += s2 * N2(D2(s1/s2, 1., v12, t), -D2(k/s2, 1., v2, t), -r212);

	} else if ((cType == 'P') && (mType == 'B')) {
	/*
	 * Put on minimum
	 */
	*retVal  = k * (1. - N2(D2(s1/k,1.,v1,t), D2(s2/k,1.,v2,t), r12));
	*retVal += - s1 * N2(D2(s2/s1, 1., v12, t), D2(k/s1, 1., v1, t), r112);
	*retVal += - s2 * N2(D2(s1/s2, 1., v12, t), D2(k/s2, 1., v2, t), r212);

	} else if ((cType == 'P') && (mType == 'A')) {
	/*
	 * Put on maximum
	 */
	*retVal  = k * N2(-D2(s1,k,v1,t), -D2(s2,k,v2,t), r12);
	*retVal += - s1 * N2(-D2(s2/s1, 1., v12, t), D2(k/s1, 1., v1, t), -r112);
	*retVal += - s2 * N2(-D2(s1/s2, 1., v12, t), D2(k/s2, 1., v2, t), -r212);

	} else {
		GtoErrMsg("%s: bad option type.\n", routine);
		goto done;
	}

	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);

#undef	N
#undef	N1
#undef	N2
}


/*---------------------------------------------------------------
 * Static vriables for CumTriNorm
 */

#define	DIM_MAX		10

static	int	dim,
		tmpI[DIM_MAX];

static	double	t,
		v[DIM_MAX+1],
		s[DIM_MAX+1],
		r[DIM_MAX+1][DIM_MAX+1],
		vcomp[DIM_MAX+1][DIM_MAX+1],
		rcomp[DIM_MAX+1][DIM_MAX+1][DIM_MAX+1];




#define	N3	DrlCumTriNorm
static	double	eps = 1e-7;


/*---------------------------------------------------------------
 * 3D cumulative normal.
 *
 * <br><br>
 * Convenience function for $N_3$ (3-dimensional
 * cumulative normal distribution) calling
 * the <i> C Analytics Library</i>.
 */


DLL_EXPORT(double)
DrlCumTriNorm(
	double x1, double x2, double x3,
	double r12, double r13, double r23)
{
	long	dimension = 3,
		ifail, errCode;
static	long	inf[3];
	double	a[3], b[3], sig[3],
		prob, bound;

	inf[0] = 1L;
	inf[1] = 1L;
	inf[2] = 1L;

	a[0] = x1;
	a[1] = x2;
	a[2] = x3;

	b[0] = 1e2;	/* should not be actually used */
	b[1] = 1e2;
	b[2] = 1e2;

	sig[0] = r12;
	sig[1] = r13;
	sig[2] = r23;

	bound = 1e-6;
	
	errCode = GtoMultiNormalCum(dimension, inf, a, b, sig,
		eps, &prob, &bound, &ifail);

	/*printf("CumTriNorm: x=(%lf, %lf, %lf) r=(%lf, %lf, %lf) N3=%lf\n",
		x1, x2, x3, r12, r13, r23, prob);*/



	/*if (ifail != 0 || errCode == FAILURE) {*/
	if (ifail != 0) {
	    GtoErrMsg("CumTriNorm: x=(%lf, %lf, %lf) "
		"r=(%lf, %lf, %lf) N3=%lf\n",
		x1, x2, x3, r12, r13, r23, prob);
	     GtoErrMsg("            errCode = %d    IFAULT = %d\n",
		errCode, ifail);

		return(1e10);
	} else {
		return(prob);
	}
}



/*f--------------------------------------------------------------
 * Options formulae : max of 3 assets.
 *
 * <br><br>
 * Option on the maximum of 3 assets:
 * call on maximum
 * max(max(S_1,S_2,S_3) - K, 0), 
 * or put on maximum
 * max(K - (S_1,S_2,S_3), 0).
 * Argument list:
 * <br>
 * <br>[t] time to expiration in years,
 * <br>[sVec] array of length 3 containing the expected values
 * of the assets,
 * <br>[sigmaVec] 
 * array of length 3 containing the volatilities
 * <br>[rhoVec]
 * array of length 3 containing the correlations: 
 * rhoVec[0] = rho_{12},
 * rhoVec[1] = rho_{13},
 * rhoVec[2] = rho_{23},
 * <br>[strike] strike of the option,
 * <br>[cp] <i> "CA"</i> for call/max, <i> "PA"</i> for put/max,
 * <br>[what] <i> "P"</i> for premium.
 * <br>
 */

DLL_EXPORT(int)
DrlOpMax3Sec(
	double tExp,		/* (I) time to exp */
	double *sVec,		/* (I) exp. assets */
	double *sigmaVec,	/* (I) vector of vol. */
	double *rhoVec,		/* (I) vector of corr. */
	double strike,		/* (I) strike */
	char *cp,		/* (I) see below */
	char *what,		/* (I) see below */
	double *retVal)		/* (O) output */
{
static	char	routine[] = "DrlOpMax3Sec";
	int	status = FAILURE;

	int	i, j, k;
	char	cType, mType;
const	int	DIM = 3;


#define	BTERM(i,j,k,l,c1,c2,c3)\
			(s[i] * N3(\
			- (c1)*D2(s[j], s[i], vcomp[i][j], tExp),\
			- (c2)*D2(s[k], s[i], vcomp[i][k], tExp),\
			- (c3)*D2(s[l], s[i], vcomp[i][l], tExp),\
			(c1)*(c2)*rcomp[i][j][k],\
			(c1)*(c3)*rcomp[i][j][l],\
			(c2)*(c3)*rcomp[i][k][l]))
#define	CNZERO(x)		{x = (fabs(x) <  1e-8 ? 1e-8 : x);}
#define	CNNEG(x)		{x = ((x) < 1e-8 ? 1e-8 : (x));}
#define	CNLE1(x)		{x = (fabs(1.-x) <= 1e-8 ?  1.-1e-8 : x);\
				 x = (fabs(1.+x) <= 1e-8 ? -1.+1e-8 : x);}



	/*
	 *
	 */
	if ((cp[0] == '\0') || (cp[1] == '\0')) {
		GtoErrMsg("%s: bad option type.\n", routine);
		goto done;
	}

	cType = toupper(cp[0]);
	mType = toupper(cp[1]);

	/*
	 * Check the values
	 */

	CNZERO(tExp);
	CNZERO(sVec[0]);
	CNZERO(sVec[1]);
	CNZERO(sVec[2]);
	CNZERO(sigmaVec[0]);
	CNZERO(sigmaVec[1]);
	CNZERO(sigmaVec[2]);
	CNZERO(strike);
	CNLE1(rhoVec[0]);
	CNLE1(rhoVec[1]);
	CNLE1(rhoVec[2]);

	/*
	 * Calculate here composite vol and correlations
	 */

	s[0] = strike;
	s[1] = sVec[0];
	s[2] = sVec[1];
	s[3] = sVec[2];

	/* volatility */
	v[0] = 0;
	v[1] = sigmaVec[0];
	v[2] = sigmaVec[1];
	v[3] = sigmaVec[2];

	/* corrrelation */
	r[1][2] = r[2][1] = rhoVec[0];
	r[1][3] = r[3][1] = rhoVec[1];
	r[2][3] = r[3][2] = rhoVec[2];
	
	for(i=0; i<=DIM;i++) r[0][i] = r[i][0] = 0e0;
	for(i=1; i<=DIM;i++) r[i][i] = 1e0;

	/* composite volatility */
	for(i=0; i<=DIM;i++)
	for(j=0; j<=DIM;j++)
		vcomp[i][j] = sqrt(v[i]*v[i]+v[j]*v[j]-2e0*v[i]*v[j]*r[i][j]);

	/* composite correlation */
	for(j=0; j<=DIM;j++)
	for(k=j+1; k<=DIM;k++) {
		i=0;
		rcomp[i][j][k] = rcomp[i][k][j] = r[j][k];

		for(i=1; i<=DIM;i++) {
			rcomp[i][j][k] = rcomp[i][k][j] = 
				(v[j]*v[k]*r[j][k] - v[i]*v[j]*r[i][j] 
				- v[i]*v[k]*r[i][k] + v[i]*v[i])/
				(vcomp[i][j] * vcomp[i][k]);

		}
	}



	/*
	 *
	 */
	if (strcmp(what, "P") || strcmp(what, "P")) {
		GtoErrMsg("%s: sensitivities N/A\n", routine);
	}

	if ((cType == 'C') && (mType == 'A')) {
		/*
		 * Call on maximum
		 */

		*retVal  =  BTERM(0,1,2,3,  1e0, 1e0, 1e0) - s[0];
		*retVal +=  BTERM(1,0,2,3,  1e0, 1e0, 1e0);
		*retVal +=  BTERM(2,0,1,3,  1e0, 1e0, 1e0);
		*retVal +=  BTERM(3,0,1,2,  1e0, 1e0, 1e0);


	} else if ((cType == 'P') && (mType == 'A')) {
	/*
	 * Put on maximum
	 */
		*retVal  =  BTERM(0,1,2,3,  1e0, 1e0, 1e0);
		*retVal += -BTERM(1,0,2,3, -1e0, 1e0, 1e0);
		*retVal += -BTERM(2,0,1,3, -1e0, 1e0, 1e0);
		*retVal += -BTERM(3,0,1,2, -1e0, 1e0, 1e0);

	} else if ((cType == 'C') && (mType == 'B')) {
	/*
	 * Call on minimum
	 */
		*retVal  = -BTERM(0,1,2,3, -1e0, -1e0, -1e0);
		*retVal +=  BTERM(1,0,2,3,  1e0, -1e0, -1e0);
		*retVal +=  BTERM(2,0,1,3,  1e0, -1e0, -1e0);
		*retVal +=  BTERM(3,0,1,2,  1e0, -1e0, -1e0);

	} else if ((cType == 'P') && (mType == 'B')) {
	/*
	 * Put on minimum
	 */
		*retVal  =  BTERM(0,1,2,3, -1e0, -1e0, -1e0) - s[0];
		*retVal +=  BTERM(1,0,2,3, -1e0, -1e0, -1e0);
		*retVal +=  BTERM(2,0,1,3, -1e0, -1e0, -1e0);
		*retVal +=  BTERM(3,0,1,2, -1e0, -1e0, -1e0);

	} else {
		GtoErrMsg("%s: bad option type\n", routine);
		goto done;
	}

	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);

#undef	BTERM
#undef	CNZERO
#undef	CNNEG
#undef	CNLE1
}


/*---------------------------------------------------------------
 * 
 */

static	double
BTerm(int idx)
{
	int	j, k;
static	double	aVec[DIM_MAX],
		rVec[(DIM_MAX)*(DIM_MAX)];
	double	retVal;

	for(j=0; j<=dim;j++) tmpI[j] = j;
	tmpI[0] = idx;
	tmpI[idx] = 0;


	for(j=1; j<=dim;j++) {
		aVec[j-1] = -D2(s[j], s[idx], vcomp[idx][tmpI[j]], t);

	}

	for(j=2; j<=dim;j++)
	for(k=1; k<=j-1;k++) {
		rVec[((j-2)*(j-1))/2+k-1] = rcomp[idx][tmpI[j]][tmpI[k]];
	}

	retVal = DrlCumMultiNorm(dim, aVec, rVec);

#if defined(DEBUG)
	for(j=0; j<=dim-1;j++)
		printf("BTerm %d  [%d] a=%lf  r=%lf\n",
			idx, j, aVec[j], rVec[j]);
	printf("BTerm %d  %lf\n", idx, retVal);
#endif

	return(retVal);
}



/*---------------------------------------------------------------
 *  DO NOT USE
 */

DLL_EXPORT(double)
DrlOpMaxNSec(
	double texp,		/* (I) time to expiration */
	int dimension,
	double *f,		/* (I) f[0]=strike, f[1]..f[n]=underlying */
	double *vol,		/* (I) v[1]..v[n]=vols */
	double *rhovec,		/* (I) rho[i][j] = vrho[(i-1)*i/2+j-1]
				    (lower triang by row starting 0) */
	char *cp,		/* (I) "C"=call, "P"=put, "A"=max, "B"=min */
	char *what,		/* (I) option type */
	double *retVal)		/* (O) output */
{
static	char	routine[] = "DrlOpMaxNSec";
	int	status = FAILURE;

	int	i, j, k;
	char	cType, mType;


	/*
	 *
	 */
	if ((cp[0] == '\0') || (cp[0] == '\0')) {
		GtoErrMsg("%s: bad options type.\n", routine);
		goto done;
	}

	cType = toupper(cp[0]);
	mType = toupper(cp[1]);

	/*
	 * Check valied values
	 */
	if ((dimension < 2) || (dimension > DIM_MAX-1)) {
		GtoErrMsg("%s: bad # of dimensions.\n", routine);
		goto done;
	}

	/*
	 *
	 */
	t = texp;
	dim = dimension;

	for(i=0; i<=dim;i++) {
		v[i] = vol[i];
		s[i] = f[i];
	}
	v[0] = 0e0;	/* vol of stk ! */


	for(i=2; i<=dim;i++)
	for(j=1; j<=i-1;j++) {
		r[i][j] = r[j][i] = rhovec[((i-2)*(i-1))/2 +j -1];
	}
	for(i=0; i<=dim;i++) r[0][i] = r[i][0] = 0e0;
	for(i=1; i<=dim;i++) r[i][i] = 1e0;


	for(i=0; i<=dim;i++)
	for(j=0; j<=dim;j++)
		vcomp[i][j] = sqrt(v[i]*v[i]+v[j]*v[j]-2e0*v[i]*v[j]*r[i][j]);

	for(j=0; j<=dim;j++)
	for(k=j+1; k<=dim;k++) {
		i=0;
		rcomp[i][j][k] = rcomp[i][k][j] = r[j][k];

		for(i=1; i<=dim;i++) 
			rcomp[i][j][k] = rcomp[i][k][j] = 
				(v[j]*v[k]*r[j][k] - v[i]*v[j]*r[i][j] 
				- v[i]*v[k]*r[i][k] + v[i]*v[i])/
				(vcomp[i][j] * vcomp[i][k]);

		}

#define	BTERM(i,j,k,l)	(s[i] * N3(\
				-D2(s[j], s[i], vcomp[i][j], t),\
				-D2(s[k], s[i], vcomp[i][k], t),\
				-D2(s[l], s[i], vcomp[i][l], t),\
				rcomp[i][j][k],\
				rcomp[i][j][l],\
				rcomp[i][k][l]))


	/*
	 *
	 */

	if ((cType == 'C') && (mType == 'A')) {
	/*
	 * Call on maximum
	 */
		*retVal  = - s[0] * (1.0 - BTerm(0));
		for(i=1; i<=dim;i++) {
			*retVal += s[i] * BTerm(i);
		}


	} else if ((cType == 'C') && (mType == 'B')) {
	/*
	 * Call on minimum
	 */
		goto done;

	} else if ((cType == 'P') && (mType == 'B')) {
	/*
	 * Put on minimum
	 */
		goto done;

	} else if ((cType == 'P') && (mType == 'A')) {
	/*
	 * Put on maximum
	 */
		goto done;
	} else {
		goto done;
	}

	status = SUCCESS;

done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}

