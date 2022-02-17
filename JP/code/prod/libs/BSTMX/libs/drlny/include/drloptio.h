/****************************************************************
 * Module:	DRL
 * Submodule:	OPTION - Option Formulas
 * File:	drloptio.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef _drloptio_H
#define _drloptio_H

#include "drlstd.h"

#ifdef	_SKIP
#include <string.h>
/*---------------------------------------------------------------
 * Useful macros
 */

#define	If(a,b,c)		((a) ? (b) : (c))
#ifndef	MAX
#define	MAX(a,b)		((a)>=(b) ? (a) : (b))
#endif
#ifndef	MIN
#define	MIN(a,b)		((a)>=(b) ? (b) : (a))
#endif
#define	SQR(a)			((a)*(a))
#define	ISFLAG(b)		(!strncmp(what,b,strlen(b)))

#define	DRL_ERROR_VALUE	(1e15)
#define	INFTY			(1e15)

#define	CNZERO(x)		(x = (fabs(x) <  1e-8 ? 1e-8 : x))
#define	CNNEG(x)		(x = ((x) < 1e-8 ? 1e-8 : (x)))
#define	CNLE1(x)		(x = (fabs(1.-x) <= 1e-8 ?  1.-1e-8 : x),\
				 x = (fabs(1.+x) <= 1e-8 ? -1.+1e-8 : x))

#endif



/*---------------------------------------------------------------
 * Functions
 */

/*
 * Gaussian Integrals
 */
#ifdef	_SKIP
extern	DLL_EXPORT(double)	DrlCumNorm(double x);
extern	DLL_EXPORT(double)	DrlDenNorm(double x);
extern	DLL_EXPORT(double)	DrlCumBiNorm(double a, double b, double p);
extern	DLL_EXPORT(double)	DrlCumTriNorm(double x1, double x2, double x3,
			double r12, double r13, double r23);
extern	DLL_EXPORT(double)	DrlCumMultiNorm(int dim, double *a, double *r);
#endif


/*
 * Options on one underlying
 */

extern	DLL_EXPORT(int)	DrlBinary(double t, double s, double v,
				double k, char *what,
				double *retVal);

extern	DLL_EXPORT(int)	DrlNormOption(double texp, double p,
				double vol, double k,
				char *cp, char *what,
				double *retVal);
extern	DLL_EXPORT(int)	DrlNormOptionImplVol(double texp,
				double fp, double prem,
				double k, char *cp, char *what,
				double *retVal);


extern	DLL_EXPORT(int)	DrlBlack(double texp, double p, double vol, double k,
				char *cp, char *what, double *retVal);
extern	DLL_EXPORT(int)	DrlBlackImplVol(double texp, double fwdp, double prem,
				double stk, char *cp, char *what,
				double *retVal);


/*
 * Options on two underlyings
 */

extern	DLL_EXPORT(int)	DrlOpMax2Sec(double t, double s1, double s2,
				double sigma1, double sigma2, double rho,
				double k, char *cp, char *what,
				double *retVal);

extern	DLL_EXPORT(int)	DrlSpreadOption(double texp, double p1, double p2,
				double vol1, double vol2, double cor,
				double strike, char *callPut, char *what,
				double *retVal);

extern	DLL_EXPORT(int)	DrlMargrabe(double t1, double t2,
				double s1, double s2,
				double vol1, double vol2, double volf,
				double rho,
				char *callPut, char *what, double *retVal);


extern	DLL_EXPORT(int)	DrlOptionDoubleKO(double expiration,
				double forward2, double strike, double vol2,
				double spot1, double forward1, double vol1,
				double upperBarrier, double lowerBarrier,
				double correlation, char   optionType,
				double *premium);

/*
 * Option on several assets
 */
extern	DLL_EXPORT(double) DrlCumTriNorm(
				double x1, double x2, double x3,
				double r12, double r13, double r23);

extern	DLL_EXPORT(int)	DrlOpMax3Sec(double t, double *sVec,
				double *sigmaVec, double *rhoVec,
				double strike,
				char *cp, char *what,
				double *retVal);

/*
 * Special use
 */


extern DLL_EXPORT(int) ConvexityC_BS2Q (
    double    S,         /* (I) Annualized volatility  */
    double    QLeft,     /* (I) Q left                 */
    double    QRight,    /* (I) Q right                */
    double    FwdSh,     /* (I) Fwd shift              */
    double    *CC);      /* (O) Constant               */

extern	DLL_EXPORT(int)	DrlOptBS2Q(
	double Y,                /* (I) Fwd yield              */
	double K,                /* (I) Strike                 */
	double T,                /* (I) Option expiration      */
	double S,                /* (I) Annualized volatility  */
	char   CoP,              /* (I) Call or put            */
	double QLeft,            /* (I) Q left                 */
	double QRight,           /* (I) Q right                */
	double FwdSh,            /* (I) Fwd shift              */
	double *value);		/* (O) option value */

extern	DLL_EXPORT(int)	DrlOptBS2QImplVol(
	double Y,                /* (I) Fwd yield              */
	double K,                /* (I) Strike                 */
	double T,                /* (I) Option expiration      */
	double P,                /* (I) Price of option        */
	char   CoP,              /* (I) Call or put            */
	double QLeft,            /* (I) Q left                 */
	double QRight,           /* (I) Q right                */
	double FwdSh,            /* (I) Fwd shift              */
	double VolGuess,         /* (I) Initial vol guess      */
	double *implVol);	/* (O) implied volatility */

extern	DLL_EXPORT(int)	DrlATMTo2QVol(
	double texp,	/* (I) time to expiration */
	double fwd,	/* (I) exp. value of underlying */
	double vol,	/* (I) volatility */
	KVolType vType, /* (I) LOGVOL, NORMVOL */
	double strike,	/* (I) strike */
	double QLeft,	/* (I) Q left */
	double QRight,	/* (I) Q right */
	double FwdSh,	/* (I) Fwd shift */
	char *callPut,	/* (I) "C" for call, "P" for put */
	double *retVal);/* (O) return value */



#endif
