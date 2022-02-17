/****************************************************************
 * Module:	DRL
 * Submodule:	OPTIO
 * File:	
 * Function:	Wrapper
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <float.h>

#include "drloptio.h"
#include "drlmath.h"

#include "drloptil.h"		/* prototype consistency */

/*-----------------------------------------------------------------------
 * LIL wrapper for DrlBinary.
 */

DLL_EXPORT(int)
DrlBinaryL(
	double *t,
	double *s,
	double *v,
	double *k,
	char *what,
	double *premium)
{

	if (	((int) t[0] != 1) ||
		((int) s[0] != 1) ||
		((int) v[0] != 1) ||
		((int) k[0] != 1) ||
		((int) what[0]	!= 1))
	{
		return(FAILURE);
	}

	return DrlBinary(t[1], s[1],
		v[1],
		k[1],
		&what[1],
		&premium[1]);
}


/*---------------------------------------------------------------
 * LIL wrapper for DrlBlack.
 */

DLL_EXPORT(int)
DrlBlackL(
	double	*t,
	double	*p,
	double	*vol,
	double	*k,
	char	*callPut,
	char	*what,
	double	*premium)
{
	if (	((int) t[0]	!= 1) ||	
		((int) p[0]	!= 1) ||
		((int) vol[0]	!= 1) ||
		((int) k[0]	!= 1) ||
		((int) callPut[0] != 1) ||
		((int) what[0]	!= 1))
	{
		return(FAILURE);
	}

	return DrlBlack(t[1], p[1], vol[1], k[1],
			callPut+1, what+1,
			&premium[1]);
}


/*---------------------------------------------------------------
 * LIL wrapper.
 */

DLL_EXPORT(int)
DrlCumBiNormL(
	double *a,
	double *b,
	double *r,
	double *retVal)
{
	if (	((int) a[0]	!= 1) ||	
		((int) b[0]	!= 1) ||
		((int) r[0]	!= 1))
	{
		return(FAILURE);
	}
	retVal[1] = DrlCumBiNorm(a[1], b[1], r[1]);
	return(SUCCESS) ;
}


/*----------------------------------------------------------------------
 *
 */

DLL_EXPORT(int)
DrlCumNormL(
             double    *inRate,     /* (I) input rate */
             double    *outRate)
{
	static char routine[] = "DrlCumNormL";

	if (inRate[0] != 1) {
		GtoErrMsg("%s: Number of cells for single element "
			"argument is not one\n", routine);
		return(FAILURE);
	}


	outRate[1] = DrlCumNorm(inRate[1]);

	return (SUCCESS);
}

/*----------------------------------------------------------------------
 *
 */

DLL_EXPORT(int)
DrlNormDenL(
             double    *inRate,     /* (I) input rate */
             double    *outRate)
{
	static char routine[] = "DrNormDenL";

	if (inRate[0] != 1) {
		GtoErrMsg("%s: Number of cells for single element "
			"argument is not one\n", routine);
		return(FAILURE);
	}


	outRate[1] = DrlDenNorm(inRate[1]);

	return (SUCCESS);
}


/*f---------------------------------------------------------------------
 * Wrapper for normal option formula.
 */

DLL_EXPORT(int)
DrlNormOptionL(
	double *t,		/*  1, F (I) time to exp */
	double *p,		/*  2, F (I) asset */
	double *vol,		/*  3, F (I) volatility (normal) */
	double *k,		/*  4, F (I) strike */
	char *cp,		/*  5, C (I) option type (C,P) */
	char *what,		/*  6, C (I) value type (P,D,G,V,T) */
	double *premium)	/*       (O) return value */
{

	if (	((int) t[0]	!= 1) ||	
		((int) p[0]	!= 1) ||
		((int) vol[0]	!= 1) ||
		((int) k[0]	!= 1) ||
		((int) cp[0]	!= 1) ||
		((int) what[0]	!= 1))
	{
		return(FAILURE);
	}

	return DrlNormOption(t[1], p[1], vol[1], k[1],
		&cp[1],
		&what[1],
		&premium[1]);
}


/*f---------------------------------------------------------------------
 * Wrapper for normal option implied volatility.
 */

DLL_EXPORT(int)
DrlNormOptionImplVolL(
	double *t,		/*  1, F (I) time to expiration */
	double *p,		/*  2, F (I) asset */
	double *prem,		/*  3, F (I) option premium */
	double *k,		/*  4, F (I) strike */
	char *cp,		/*  5, C (I) option type (C,P) */
	char *what,		/*  6, C (I) value type (P,D,G,etc) */
	double *vol)		/*       (O) implied vol */
{

	if (	((int) t[0]	!= 1) ||	
		((int) p[0]	!= 1) ||
		((int) prem[0]	!= 1) ||
		((int) k[0]	!= 1) ||
		((int) cp[0]	!= 1) ||
		((int) what[0]	!= 1) ||
		((int) vol[0]	!= 1))
	{
		return(FAILURE);
	}

	return DrlNormOptionImplVol(
		t[1], p[1], prem[1], k[1],
		cp+1, what+1,
		&vol[1]);
}


/*f---------------------------------------------------------------------
 * Wrapper for option on the max/min of 2 assets using
 * analytical formulas.
 */

DLL_EXPORT(int)
DrlOpMax2SecL(
	double *t,		/*  1, F (I) time to exp  */
	double *s1,		/*  2, F (I) asset 1 */
	double *s2,		/*  3, F (I) asset 2 */
	double *vol1,		/*  4, F (I) vol asset 1 */
	double *vol2,		/*  5, F (I) vol asset 2 */
	double *rho,		/*  6, F (I) correlation  */
	double *k,		/*  7, F (I) strike */
	char *cp,		/*  8, C (I) option type (CA,CB,PA,PB) */
	char *what,		/*  9, C (I) 'P'remuim */
	double *retVal)		/*       (O) return value */
{

	if (	((int) t[0]	!= 1) ||	
		((int) s1[0]	!= 1) ||
		((int) s2[0]	!= 1) ||
		((int) vol1[0]	!= 1) ||
		((int) vol2[0]	!= 1) ||
		((int) rho[0]	!= 1) ||
		((int) k[0]	!= 1) ||
		((int) cp[0]	!= 1) ||
		((int) what[0]	!= 1))
	{
		return(FAILURE);
	}


	return DrlOpMax2Sec(t[1], s1[1], s2[1],
			vol1[1], vol2[1], rho[1], k[1],
			cp+1, what+1,
			&retVal[1]);

}


/*f---------------------------------------------------------------------
 * Wrapper for option on the max/min of 3 assets.
 */

DLL_EXPORT(int)
DrlOpMax3SecL(
	double *t,		/*  1, F (I)  time to exp */
	double *s,		/*  2, F (I)  array of assets [3] */
	double *vol,		/*  3, F (I)  array of volatilities [3] */
	double *rho,		/*  4, F (I)  array of correlations [3] */
	double *k,		/*  5, F (I)  strike */
	char *cp,		/*  6, C (I)  option type (CA,CB,PA,PB) */
	char *what,		/*  7, C (I)  'P'remium */
	double *retVal)		/*       (O)  return value */
{

	if (	((int) t[0]	!= 1) ||	
		((int) s[0]	!= 3) ||
		((int) vol[0]	!= 3) ||
		((int) rho[0]	!= 3) ||
		((int) k[0]	!= 1) ||
		((int) cp[0]	!= 1) ||
		((int) what[0]	!= 1))
	{
		return(FAILURE);
	}

	return DrlOpMax3Sec(t[1], s+1, vol+1, rho+1,
			k[1],
			cp+1, what+1,
			&retVal[1]);
}



/*----------------------------------------------------------------------
 *
 */

DLL_EXPORT(int)
DrlCumTriNormL(
	double *s,		/*  1, F (I)  */
	double *rho,		/*  2, F (I)  */
	double *retVal)		/*       (O)  */
{

	if (	((int) s[0]	!= 3) ||
		((int) rho[0]	!= 3) )
	{
		return(FAILURE);
	}

	retVal[1] = DrlCumTriNorm(
		s[1], s[2], s[3],
		rho[1], rho[2], rho[3]);

	return(SUCCESS) ;
}


/*f----------------------------------------------------------------------
 * Wrapper for Margrabe formula.
 */


DLL_EXPORT(int)
DrlMargrabeL(
	double *t1L,	/*  1 F (I) time to expiration # 1*/
	double *t2L,	/*  2 F (I) time to expiration (t1 < t2) # 2 */
	double *s1L,	/*  3 F (I) underlying # 1 */
	double *s2L,	/*  4 F (I) underlying # 2 */
	double *vol1L,	/*  5 F (I) base volatility # 1*/
	double *vol2L,	/*  6 F (I) base volatility # 2*/
	double *volfL,	/*  7 F (I) fwd volatility # 2*/
	double *rhoL,	/*  8 F (I) correlation */
	char *callPutL,	/*  9 C (I) "C" for call, "P" for put */
	char *whatL,	/* 10 C (I) see below */
	double *retValL)/*    F (O) */
{
static	char	routine[] = "MargrabeL";
	int	status = FAILURE;

	WRAP_CHECK_SCALAR(t1L);
	WRAP_CHECK_SCALAR(t2L);
	WRAP_CHECK_SCALAR(s1L);
	WRAP_CHECK_SCALAR(s2L);
	WRAP_CHECK_SCALAR(vol1L);
	WRAP_CHECK_SCALAR(vol2L);
	WRAP_CHECK_SCALAR(volfL);
	WRAP_CHECK_SCALAR(rhoL);
	WRAP_CHECK_SCALAR(callPutL);
	WRAP_CHECK_SCALAR(whatL);
	WRAP_CHECK_SCALAR(retValL);

	return DrlMargrabe(
		t1L[1],
		t2L[1],
		s1L[1],
		s2L[1],
		vol1L[1],
		vol2L[1],
		volfL[1],
		rhoL[1],
		&callPutL[WRAP_STR_IDX(1)],
		&whatL[WRAP_STR_IDX(1)],
		&retValL[1]);
done:
	return(status);
}




/*f----------------------------------------------------------------------
 * Wrapper for spread option formula (using numerical)
 * integration.
 */

DLL_EXPORT(int)
DrlSpreadOptionL(
	double *t,		/*  1, F (I) time to exp */
	double *p1,		/*  2, F (I) asset 1 */
	double *p2,		/*  3, F (I) asset 2 */
	double *vol1,		/*  4, F (I) vol asset 1 */
	double *vol2,		/*  5, F (I) vol asset 2 */
	double *cor,		/*  6, F (I) correlation */
	double *k,		/*  7, F (I) strike */
	char *callPut,		/*  8, C (I) C,P */
	char *what,		/*  9, C (I) P,D1,D2,G1,G2,G3,V1, etc.. */
	double *premium)	/*       (O) return value */
{
	if (	((int) t[0]	!= 1) ||	
		((int) p1[0]	!= 1) ||
		((int) p2[0]	!= 1) ||
		((int) vol1[0]	!= 1) ||
		((int) vol2[0]	!= 1) ||
		((int) cor[0]	!= 1) ||
		((int) k[0]	!= 1) ||
		((int) callPut[0]	!= 1) ||
		((int) what[0]	!= 1))
	{
		return(FAILURE);
	}

	return DrlSpreadOption(t[1], p1[1], p2[1],
		vol1[1], vol2[1], cor[1], k[1], callPut+1, what+1,
		&premium[1]);

	return(SUCCESS);
}


/*-----------------------------------------------------------------------
 * Author:     Bing-Le Wu
 */

DLL_EXPORT(int)
DrlDoubleKOOptionL(
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
	double *premium)	/*    'F' (O) output*/
{
static	char	routine[] = "DrlDoubleKOOptionL";

	if (    ((int) expiration[0] !=1) ||
		((int) forward2[0] != 1) ||
		((int) strike[0] != 1) ||
		((int) vol2[0] != 1) ||
		((int) spot1[0] != 1)||
		((int) forward1[0] != 1) ||
		((int) vol1[0] != 1) ||
		((int) upperBarrier[0] != 1) ||
		((int) lowerBarrier[0] != 1) ||
		((int) correlation[0] != 1) ||
		((int) optionType[0] != 1))
	{
		return(FAILURE);
	}

	if (expiration[1] <= 0e0)
	{
		GtoErrMsg("%s: expiration is less than zero.\n", routine);
		return(FAILURE);
	}

	if (vol2[1] <= 0)
	{
		GtoErrMsg("%s: vol for pay-off variable is not positive.\n",
			routine);
		return(FAILURE);
	}
	
	
	if(vol1[1]<=0) 
	{        
		GtoErrMsg("%s: vol for knock-out variable is not positive.\n",
			routine);
		return(FAILURE);
	}
	if(upperBarrier[1]<=lowerBarrier[1])
	{
		GtoErrMsg("%s: upper barrier is not higher than "
			" lower barrier.\n", routine);
		return(FAILURE);
	}
	if(upperBarrier[1]<=spot1[1])
	{        
		GtoErrMsg("%s: upper barrier is not higher than spot.\n",
			routine);
		return(FAILURE);
	}
	if(spot1[1]<=lowerBarrier[1])
	{ 
		GtoErrMsg("%s: lower barrier is not lower than spot.\n",
			routine);
		return(FAILURE); 
	}


	if (DrlOptionDoubleKO(
		expiration[1], forward2[1], strike[1], vol2[1], 
		spot1[1], forward1[1], vol1[1], upperBarrier[1], 
		lowerBarrier[1],  correlation[1],
		optionType[1], &premium[1])
		!= SUCCESS)
			return(FAILURE);

	return(SUCCESS) ;

}


