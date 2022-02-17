/****************************************************************
 * Module:	DRL
 * Submodule:	OPTION - Option Formulas
 * File:	drloptil.h
 * Function:	Wrappers for options functions.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef _drloptil_H
#define	_drloptil_H
#include "drlstd.h"


extern	DLL_EXPORT(int)	DrlBinaryL(double *t, double *s, double *v,
	double *k, char *what, double *premium);

extern	DLL_EXPORT(int)	DrlBlackL(double *t, double *p, double *vol,
	double *k, char *callPut, char *what,
	double *premium);

extern	DLL_EXPORT(int)	DrlNormDenL( double *inRate, double *outRate);
extern	DLL_EXPORT(int)	DrlCumNormL( double *inRate, double *outRate);
extern	DLL_EXPORT(int)	DrlCumBiNormL(double *a, double *b, double *r,
	double *retVal);
extern	DLL_EXPORT(int)	DrlCumTriNormL(double *s, double *rho, double *retVal);


extern	DLL_EXPORT(int)	DrlNormOptionL(double *t, double *p, double *vol,
	double *k, char *cp, char *what, double *premium);
extern	DLL_EXPORT(int)	DrlNormOptionImplVolL(double *t, double *p, double *prem,
	double *k, char *cp, char *what, double *vol);


extern	DLL_EXPORT(int)	DrlOpMax2SecL(double *t, double *s1, double *s2,
	double *vol1, double *vol2, double *rho,
	double *k, char *cp, char *what,
	double *retVal);

extern	DLL_EXPORT(int)	DrlOpMax3SecL(double *t, double *s, double *vol,
	double *rho, double *k, char *cp,
	char *what, double *retVal);
 
extern	DLL_EXPORT(int)	DrlSpreadOptionL(double *t, double *p1, double *p2,
	double *vol1, double *vol2, double *corr,
	double *k, char *callPut, char *what,
	double *premium);

extern	DLL_EXPORT(int)	DrlMargrabeL(
	double *t1L,	/*  1 F (I) time to expiration # 1*/
	double *t2L,	/*  2 F (I) time to expiration (t1 < t2) # 2 */
	double *s1L,	/*  3 F (I) underlying # 1 */
	double *s2L,	/*  4 F (I) underlying # 2 */
	double *vol1L,	/*  5 F (I) base volatility # 1*/
	double *vol2L,	/*  6 F (I) base volatility # 2*/
	double *volfL,	/*  7 F (I) fwd volatility # 2*/
	double *rhoL,	/*  8 F (I) correlation */
	char *callPutL,	/*  9 C (I) \"C\" for call, \"P\" for put */
	char *whatL,	/* 10 C (I) see below */
	double *retValL);/*     (O) */

extern	DLL_EXPORT(int)	DrlDoubleKOOptionL(
	double	*expiration,	/*  1 'F' (I) */
	double *forward2,	/*  2 'F' (I) payoff parameters*/
	double *strike,		/*  3 'F' (I) */
	double *vol2,		/*  4 'F' (I) */
	double *spot1,		/*  5 'F' (I) barrier variable parameters*/
	double *forward1,	/*  6 'F' (I) */
	double *vol1,		/*  7 'F' (I) */
	double *upperBarrier,	/*  8 'F' (I) */
	double *lowerBarrier,	/*  9 'F' (I) */
	double *correlation, 	/* 10 'F' (I) corr var1 var2*/
	char   *optionType,	/* 11 'C' (I) C for call and P for put*/
	double *premium);	/*    'F' (O) output*/



#endif
