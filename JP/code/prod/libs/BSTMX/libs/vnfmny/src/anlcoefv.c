/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>

#include "convert.h"

#include "drlsort.h"

#define	_vnfm_SOURCE
#include "vnfmanly.h"

/* Now are static routine */

	int	VnfmA(VnfmData *that, double S, double T, double *A);
	int	VnfmQ(VnfmData *that, double S, double T, double *y, double *Q);
	int	VnfmB(VnfmData *that, double S,
			double tMat, int freq, double *y, double *B);

	double	H(double S1, double S2, double T, double beta1, double beta2);



#define	__DEBUG__
#undef	__DEBUG__

#define xDrKey

/*c-@CDOC(idxn="Model Specification",typn=comm)-----------------------
 * <h3>Model specification</h3>
 * The model is 
 * an $N$-factor model with exponential factors:
 * <br>
 * {dR(t,T) / R(t,T)} &=&
 *    \sum_{i=0}^N \alpha_i\sigma_i(t) e^{-\beta_i(T-t)} dW_i(t)\\
 * dW_i(t) dW_j(t) &=& \rho_{ij}(t) dt
 * <br>
 * where
 * <br>
 * <br> $(\beta_i)_{i=1..N}$ are the mean-reversion coefficients,
 * <br> $(\alpha_i)_{i=1..N}$ are the weight coefficients,
 * <br> $(\sigma_i(t))_{i=1..N}$ are the spot volatilities
 * of the factors which are assumed to be piecewise constant,
 * <br> $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$ are the spot correlations
 * of the factors which are assumed to be piecewise constant.
 * <br>
 * %The time line  is defined by
 * %$0=t_0 < t_1 < ... t_{M-1}$.
 * %<center>
 * %<!-- (400,50)(-10,-25)
 * %\put(  0, 0){\vector( 1, 0){300}}
 * %\put( 00,-2){\line( 0, 1){4}}
 * %\put( 00,-10){\makebox(0,0){$t_0=0$}}
 * %\put( 40,-2){\line( 0, 1){4}}
 * %\put( 40,-10){\makebox(0,0){$t_1$}}
 * %\put( 80,-2){\line( 0, 1){4}}
 * %\put( 80,-10){\makebox(0,0){$t_2$}}
 * %\put(200,-2){\line( 0, 1){4}}
 * %\put(200,-10){\makebox(0,0){$t_{i}$}}
 * %\put(240,-2){\line( 0, 1){4}}
 * %\put(240,-10){\makebox(0,0){$t_{i+1}$}}
 * %-->
 * %</center>
 * %If $x$ is any piecewise constant time-dependent quantity,
 * %(for example $\sigma(t)$, $\rho(t)$),
 * %then $x_i$ designates the (constant) value of $x(t)$ of the
 * %interval $t_i &lt;= t < t_{i+1}$.
 * %(or $t_i &lt;= t < +\infty$ if $i=N-1$).
 * % 
 * %The matrix of correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$ is
 * %stored as a vector <i> rho[k]</i> with the matrix index
 * %$(i+1,j+1)$ with $1&lt;= i<j &lt;= N$
 * %correponding to the vector index
 * %$$k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$$
 * %with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * %The matrix of factor covariances $J_{ij}(S)$ defined by
 * %<br>
 * %J_{ij}(S) &=&integral_0^S
 * % 		\alpha_i\sigma_i(t) \alpha_j\sigma_j(t) \rho_{ij}(t)
 * %		e^{-(\beta_i+\beta_j) (S-t)} dt.
 * %<br>
 * %for $1&lt;= i&lt;= j &lt;= N$
 * %is stored as a vector <i> J[k]</i> with the matrix index
 * %$(i+1,j+1)$ with $1&lt;= i&lt;= j &lt;= N$
 * %correponding to the vector index
 * %$$k = {N(N+1)/ 2} - {(N-i)(N-i+1)/ 2} + (j-i)$$
 * %with  $0&lt;= k &lt;= {N(N+1)/ 2}-1$.
 */

/*--------------------------------------------------------------
 * Macros used for some simple integrals:
 * <br>
 *	K(a,b,\alpha) &=& integral_a^b e^{-\alpha*t} dt \\
 *	L(a,\alpha)   &=& integral_0^a e^{-\alpha*t} dt \\
 *	H(a,b,T, \beta_i, beta_j) &=& integral_a^b (1-e^{-\beta_i*(T-t)})(1-e^{-\beta_j*(T-t)}) dt \\
 * <br>
 */

#ifndef	_vnfm_SOURCE	/* now in header file */
#error	should be defined for source code
#endif



/*f-@CDOC(catn="Low Level Analytical Routines")--------------------
 * Given the model timeline $t_0<...<t_{n-1}$ and $x$, computes
 * $$j \;=\; \sup\{i \;\mbox{such that}\; t_i <= x\}.$$
 * If $x < t_0$, returns $t=0$ with error code 1 (0 otherwise).
 */

int
VnfmFloorIdx(VnfmData *that, double x, int *j)
{
	register int jl, ju, jm;

	/* searches for the closest lower index in the timeline
	 * for a given maturity time (ACT/365F).
	 * to avoid rounding errors, 1 hour (1/365/60 = 45.66210046e-6)
	 * is added
	 */
	x+= 45.66210046e-6;


	jl = (-1);
	ju = that->fNDates;
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= that->fTime[jm])
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
	if (*j == -1) {
		(*j)++;
		return(1);
	}
	return(0);
}


int
VnfmDateFloorIdx(VnfmData *that, TDate x, int *j)
{
	register int jl, ju, jm;


	jl = (-1);
	ju = that->fNDates;
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= that->fDate[jm])
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
	if (*j == -1) {
		(*j)++;
		return(1);
	}
	return(0);
}




/*f-------------------------------------------------------------
 * Update all cached coefficients.
 *                                                             
 * <br><br>
 * Computes all internal coefficients and forward rates used
 * in the closed form analytics. This function assumes that the
 * data structure "that" contains all the parameters.
 * It should be called once <i> before</i> any of the volatility or 
 * correlation computation routines are called.
 * It has to be called again anytime the parameters
 * are changed in the data structure.
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
VnfmComputeCoeff(
	VnfmData *that)		/* (I/O) model parameters */
{
	int	errCode = 1;
	int	idx;

extern	int	DrSecurity(unsigned int c);

#ifdef DrKey
	/*DR_SECURITY_CHECK();*/
	if (DrSecurity((unsigned int) 2) != SUCCESS)
		goto done;
#endif

	for (idx=0; idx<=that->fNf-1; idx++) {
		if (IS_ALMOST_ZERO(that->fBeta[idx])) that->fBeta[idx] = 1e-8;
	}


	if (VnfmComputeFwdRates(that) != 0) goto done;
	if (VnfmComputeJ(that) != 0) goto done;

#ifdef	__DEBUG__
	VnfmPrintCoeff(that, vnfmFpLog);
#endif
	/* make it through OK */
	errCode = 0;
done:
	if (errCode != 0) {
		GtoErrMsg("VnfmComputeCoeff: failed\n");
	}
	return(errCode);
}


/*f-------------------------------------------------------------
 * Update cached coefficients for J integrals.
 *                                                             
 * <br><br>
 * Computes and stores internally the volatility integrals
 * <blockquote>
 * J_{ij}(S) =integral_0^S
 * 		alpha_i sigma_i(t) alpha_j sigma_j(t) rho_{ij}(t)
 *		e^{-(beta_i+beta_j) (S-t)} dt.
 * </blockquote>
 * Should not be called directly in applications
 * (use <i> VnfmComputeCoeff</i> instead.
 * Returns 0 if OK.
 */

int
VnfmComputeJ(
	VnfmData *that)		/* (I) model parameters */
{
	int	i, j,
		nDim = NDIM,
		jIdx, rIdx,
		n;
	double	lambda, dt;

	/*
	 * for 0 <= i < j <= n-1,
	 * corr[i][j] = fRho[ n(n-1)/2 - (n-i)*(n-i-1)/2 + (j-i-1) ]
	 * for 0 <= i <= j <= n-1,
	 * J[i][j] = J[ n(n+1)/2 - (n-i)*(n-i+1)/2 + (j-i) ]
	 */

	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);
	    lambda = BETA[i] + BETA[j];

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmComputeJ(%3d): i=%d j=%d  jIdx=%2d  rIdx=%2d\n",
		__LINE__, i, j, jIdx, rIdx);
#endif
	    n = 0;
	    that->fJ[jIdx][n] = 0e0;

	    for (n=1; n<=NDATES-1; n++) {
		dt = TT[n] - TT[n-1];
		that->fJ[jIdx][n] =
		    ALPHA[i] * ALPHA[j] *
			SIGMA[i][n-1] * SIGMA[j][n-1] *
			(i == j ? 1e0 : RHO[rIdx][n-1]) *
			L(dt, lambda)
		    + exp(-dt*lambda) * that->fJ[jIdx][n-1];
	    }
	}

	return(0);
}

/*f-------------------------------------------------------------
 * Update cached coefficients for J integrals (partial recalculation).
 *                                                             
 * <br><br>
 * Same as <i> VnfmComputeJ</i>, but updates the $J$ intergals
 * only between (volatility) timeline indices
 * <i> idxStart</i> and <i> idxEnd</i>.
 * Returns 0 if OK.
 */

int
VnfmComputeJPartial(
	VnfmData *that,		/* (I) model parameters */
	int idxStart,		/* (I) first index */
	int idxEnd)		/* (I) last index */
{
register int	i, j,
		nDim = NDIM,
		jIdx, rIdx,
		n;
	double	lambda, dt;

	/* */
	idxStart = MAX(idxStart, 1);
	idxEnd   = MIN(idxEnd, NDATES-1);

	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);
	    lambda = BETA[i] + BETA[j];

	    n = 0;
	    that->fJ[jIdx][n] = 0e0;

	    for (n=idxStart; n<=idxEnd; n++) {
		dt = TT[n] - TT[n-1];
		that->fJ[jIdx][n] =
		    ALPHA[i] * ALPHA[j] *
			SIGMA[i][n-1] * SIGMA[j][n-1] *
			(i == j ? 1e0 : RHO[rIdx][n-1]) *
			L(dt, lambda)
		    + exp(-dt*lambda) * that->fJ[jIdx][n-1];
	    }
	}

	return(0);
}


/*f-------------------------------------------------------------
 * Compute J integrals.
 *                                                             
 * <br><br>
 * Computes and returns the volatility integrals
 * <br>
 * J_{ij}(S_1,S_2,T) &=&integral_{S_1}^{S_2}
 *		\alpha_i \alpha_j
 * 		\sigma_i(t) \sigma_j(t) \rho_{ij}(t)
 *		e^{-(\beta_i+\beta_j) (T-t)} dt.
 * <br>
 * where $S_1$ is the option start, $S_2$ the option expiration,
 * $T$ the rate reset and $0&lt;= S_1 &lt;= S_2 &lt;= T$.
 * The routine interpolates the integrals stored internally,
 * and <i> should always be called after</i> <i> VnfmComputeCoeff</i>
 * has been called.
 * Returns 0 iff OK.
 */

int
VnfmJ(
	VnfmData *that,		/* (I) model parameters */
	double S1,		/* (I) option start */
	double S2,		/* (I) option expiration */
	double T,		/* (I) rate reset */
	double *J0)		/* (O) J[i][j] stored as vector (see above) */
{

	int	i, j,
		nDim = NDIM,
		jIdx, rIdx,
		m, n;
	double	lambda;
	/* find closest but lower */
	DrlDoubleArrayFloorIdx(TT, NDATES, S1, &m);
	DrlDoubleArrayFloorIdx(TT, NDATES, S2, &n);

	/*
	 * for 0 <= i < j <= n-1,
	 * corr[i][j] = fRho[ n(n-1)/2 - (n-i)*(n-i-1)/2 + (j-i-1) ]
	 * for 0 <= i <= j <= n-1,
	 * J[i][j] = J[ n(n+1)/2 - (n-i)*(n-i+1)/2 + (j-i) ]
	 */

	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);
	    lambda = BETA[i] + BETA[j];

	    J0[jIdx] = 
		  exp(-(T - TT[n])*lambda) * that->fJ[jIdx][n] 
	    	+ ALPHA[i] * ALPHA[j] *
		    SIGMA[i][n] * SIGMA[j][n] *
		    (i == j ? 1e0 : RHO[rIdx][n]) *
		    exp(-(T-S2)*lambda) * L(S2-TT[n], lambda)
		- exp(-(T-TT[m])*lambda) * that->fJ[jIdx][m]
	    	- ALPHA[i] * ALPHA[j] *
		    SIGMA[i][m] * SIGMA[j][m] *
		    (i == j ? 1e0 : RHO[rIdx][m]) *
		    exp(-(T-S1)*lambda)* L(S1-TT[m], lambda);
	}

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmJ(%3d): S1=%lf S2=%lf T=%lf m=%2d n=%2d\n",
		__LINE__, S1, S2, T, m, n);
	DrlFPrintf(vnfmFpLog, "VnfmJ     :");
	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) 
		DrlFPrintf(vnfmFpLog, " J[%2d]=%lf ", JIDX(i,j), J0[JIDX(i,j)]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif

	return(0);
}


/*f-------------------------------------------------------------
 * Computes and returns volatility integrals.
 *                                                             
 * <br><br>
 * Calculate
 * <blockquote>
 * I_{ij}(S_1,S_2,T) = (1 / \beta_i \beta_j) integral_{S_1}^{S_2}
 *		alpha_i alpha_j
 * 		sigma_i(t) sigma_j(t) rho_{ij}(t)
 *		(1-e^{-beta_i (T-t)})(1-e^{-beta_j (T-t)}) dt.
 * </blockquote>
 * which is used in computing the accumulated vol between $S_1$ and $S_$ for
 * a ON rate averaging over $S_1$ and $T$ resetting at $T$. 
 * $0&lt;= S_1 &lt;= S_2 &lt;= T$.
 * The routine interpolates the integrals stored internally,
 * and <i> should always be called after</i> <i> VnfmComputeCoeff</i>
 * has been called.
 * Returns 0 iff OK.
 */

int
VnfmI(
	VnfmData *that,		/* (I) model parameters */
	double S1,		/* (I) option start */
	double S2,		/* (I) option expiration */
	double T,		/* (I) rate reset */
	double *I0)		/* (O) I[i][j] stored as vector (see above) */
{

	int	i, j, l,
		nDim = NDIM,
		jIdx, rIdx,
		m, n;
	double	beta_i, beta_j;

	double	h;

	/* find closest but lower */
	DrlDoubleArrayFloorIdx(TT, NDATES, S1, &m);
	DrlDoubleArrayFloorIdx(TT, NDATES, S2, &n);

	/*
	 * for 0 <= i < j <= n-1,
	 * corr[i][j] = fRho[ n(n-1)/2 - (n-i)*(n-i-1)/2 + (j-i-1) ]
	 * for 0 <= i <= j <= n-1,
	 * I[i][j] = I[ n(n+1)/2 - (n-i)*(n-i+1)/2 + (j-i) ]
	 */

	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);

	    I0[jIdx] = 0e0;

	    if (n == m)
	    {
		/* integrate from S1 to S2 using constant
		 * vol at index n 
		 */
	     	I0[jIdx] = 
	    	 	ALPHA[i] * ALPHA[j] *
		    	SIGMA[i][n] * SIGMA[j][n] *
		    	(i == j ? 1e0 : RHO[rIdx][n]) *
			H(S1, S2, T, BETA[i], BETA[j]);

	    }
	    else
	    {
		/* Integral from S1 to TT[m+1] */
		if (!IS_ALMOST_ZERO(TT[m+1] - S1))
	     	    I0[jIdx] = 
	    	 	ALPHA[i] * ALPHA[j] *
		    	SIGMA[i][m] * SIGMA[j][m] *
		    	(i == j ? 1e0 : RHO[rIdx][m]) *
			H(S1, TT[m+1], T, BETA[i], BETA[j]);

		/* Integral from TT[n] to S2 */
		if (!IS_ALMOST_ZERO(TT[n] - S2))
	     	    I0[jIdx] += 
	    	 	ALPHA[i] * ALPHA[j] *
		    	SIGMA[i][n] * SIGMA[j][n] *
		    	(i == j ? 1e0 : RHO[rIdx][n]) *
			H(TT[n], S2, T, BETA[i], BETA[j]);


		/* from TT[m+1] to TT[n] */
		for (l=m+1; l<n; l++)
	     	    I0[jIdx] += 
	    	 	ALPHA[i] * ALPHA[j] *
		    	SIGMA[i][l] * SIGMA[j][l] *
		    	(i == j ? 1e0 : RHO[rIdx][l]) *
			H(TT[l], TT[l+1], T, BETA[i], BETA[j]);

	    }
	}

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmJ(%3d): S1=%lf S2=%lf T=%lf m=%2d n=%2d\n",
		__LINE__, S1, S2, T, m, n);
	DrlFPrintf(vnfmFpLog, "VnfmI     :");
	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) 
		DrlFPrintf(vnfmFpLog, " I[%2d]=%lf ", JIDX(i,j), I0[JIDX(i,j)]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif

	return(0);
}



/*f-------------------------------------------------------------
 * Calculate the B coefficients.
 *                                                             
 * <br><br>
 * Computes and returns the B coefficients for coupon yield volatility.
 * "S" is the time to reset of the option,
 * "tMat" the forward maturity of the bond,
 * "freq" its frequency (1,2,4,12).\\
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>).
 * Returns 0 iff OK.
 * <b>WARNING: the B in the code is the B of the note
 * multiplied by the $1/(1-Z)$ factor. </b>
 */

int
VnfmB(
	VnfmData *that,	/* (I) model parameters */
	double S,	/* (I) expiration */
	double tMat,	/* (I) fwd bond maturity (ROUNDED TO NEAREST!) */
	int freq,	/* (I) bond frequency */
	double *yield,	/* (O) par yield */
	double *B)	/* (O) B[i](S,tMat) */
{
	/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*
	 *	New Version:					*
	 *	We do simple stubs to make it consistent	*
	 *	with simple rate				*
	 *xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

	int	nDim = NDIM,
		m, n, idxC, nC, idxF;
	double	T,
		A[VNFM_NDIMMAX], An[VNFM_NDIMMAX],/* A coefficients */
		z, zn,				/* forward zero */
		b[VNFM_NDIMMAX],		/* temporary B1, B2 storage */
		dur,				/* duration */
		dcf,				/* accrual */
		y,				/* pay yield */
		dt;
#define	ONE_DAY	2.73972602e-3

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmB: ");
#endif

	dt      = 1e0/freq;
	nC      = (int) floor((tMat-ONE_DAY)/dt) + 1;	/* # coupons */
	dur     = 0e0;
	for (idxF=0; idxF<=nDim-1; idxF++) b[idxF] = 0e0;

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "S=%lf tMat=%lf freq=%d : nC=%d\n",
		S, tMat, freq, nC);
#endif


	for (idxC=1; idxC<=nC; idxC++) {
		/* Coupon payment date */
		T = tMat - (nC - idxC) * dt + S;

		/* compute A coefficients */
		VnfmA(that, S, T, A);

		/* */
		DrlDoubleArrayFloorIdx(ZTT, NZDATES, S, &n);
		DrlDoubleArrayFloorIdx(ZTT, NZDATES, T, &m);

		/* compute forward zero */
		z = (that->fZero[m] / that->fZero[n])
			*exp( - RATE[m] * (T - ZTT[m])
			      + RATE[n] * (S - ZTT[n]) );

		/* last coupon: add principal */
		if (idxC == nC) {
			for (idxF=0; idxF<=nDim-1; idxF++) An[idxF] = A[idxF];
			zn  = z;
		}

		/* 1s coupon: case of a front stub */
		if ((idxC == 1) && (T - dt < S)) {
			dcf = (T - S);
		} else {
			dcf = dt;
		}


#ifdef	__DEBUG__
		DrlFPrintf(vnfmFpLog, "VnfmB(%3d): cpn %2d/%2d "
			"S=%7.4f T=%7.4f (%7.4f)  n=%2d m=%2d dcf=%lf "
			"z=%lf zrcc=%lf zr1c=%lf\n",
			__LINE__, idxC, nC, S, T, T-S,
			n, m, dcf, z, -log(z) / (T-S),
			pow(z, -1./(T-S))-1.);
		DrlFPrintf(vnfmFpLog, "          :");
		for (idxF=0; idxF<=nDim-1; idxF++)
			DrlFPrintf(vnfmFpLog, " A[%2d]=%lf ", idxF, A[idxF]);
		DrlFPrintf(vnfmFpLog, "\n");
#endif

		/* compute duration and beta-duration */
		dur += z * dcf;
		for (idxF=0; idxF<=nDim-1; idxF++) 
			b[idxF]  += A[idxF] * z * dcf;
	}

	/* Compute par yield: simple stub */
	y = (1e0 - zn) / dur;

	*yield = y;

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmB(%3d):", __LINE__);
	for (idxF=0; idxF<=nDim-1; idxF++)
		DrlFPrintf(vnfmFpLog, " B[%2d]=%lf ", idxF, b[idxF]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif

	for (idxF=0; idxF<=nDim-1; idxF++) {
		b[idxF] = (y * b[idxF] + An[idxF] * zn) / dur;
		B[idxF] = b[idxF];
	}


#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog,
		"VnfmB(%3d): S=%lf T=%lf freq=%3d y=%lf dur=%lf acc=%lf\n",
		__LINE__, S, tMat, freq, y, dur, (S-T));
#endif



	return(0);
#undef	__DEBUG__
}

/*f-------------------------------------------------------------
 * Calculate the A coefficients.
 *                                                             
 * <br><br>
 * Compute and returns coefficient A (used for computation
 * of the MM rates volatility) defined by
 * <blockquote>
 * A_j(S,T) = integral_S^T r^{1-2q}(t) e^{-beta_j*(t-S)} dt.
 * </blockquote>
 * where $q$ is the backboneq with cevPower = $1-2q$
 * (0=lognormal, 0.5=normal).
 *
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>).
 * Returns 0 iff OK.
 */

int
VnfmA(
	VnfmData *that,		/* (I) model parameters */
	double S,		/* (I) forward rate start */
	double T,		/* (I) forward rate end */
	double *A)		/* (O) A_i(S,T) */
{
	int	i, j, m, n,
		nDim = NDIM;

	DrlDoubleArrayFloorIdx(ZTT, NZDATES, S, &n);
	DrlDoubleArrayFloorIdx(ZTT, NZDATES, T, &m);

	for (j=0; j<=nDim-1; j++) {

	    A[j] = 0e0;

	    for (i=n+1; i<m; i++) {
		A[j] += FN(RATE[i],that->fBackBoneQ) *
				K(ZTT[i]-S, ZTT[i+1]-S, BETA[j]);
	    }

	    if (n+1 <= m) {
		A[j] += FN(RATE[n],that->fBackBoneQ) *
				L((ZTT[n+1]-S), BETA[j]);
				/* L(..) same as K(0e0, ZTT[n+1]-S, BETA[j]) */
		A[j] += FN(RATE[m],that->fBackBoneQ) *
				K(ZTT[m]-S, T-S, BETA[j]);
	    }
	    else {
		A[j] += FN(RATE[n],that->fBackBoneQ) *
				L((T-S), BETA[j]);
				/* L(..) same as K(0e0, T-S, BETA[j]) */
	    }

	}

#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmA(%3d): S=%8.4f T=%lf n=%2d m=%2d\n",
		__LINE__, S, T, m, n);
	DrlFPrintf(vnfmFpLog, "VnfmA     :");
	for (j=0; j<=nDim-1; j++) DrlFPrintf(vnfmFpLog, "%lf  ", A[j]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif

#undef	FN

	return(0);
}


/*f-------------------------------------------------------------
 * Calculate the Q coefficients.
 *                                                             
 * <br><br>
 * Compute and returns Q coefficients for a simple rate volatility
 * defined by
 * <blockquote>
 * Q_j(S,T) = {1 / log Z(S,T)} A_j(S,T).
 * </blockquote>
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>).\\
 * <b> WARNING:</b> does the adjustment from cc to  simple rate volatility.
 * Returns 0 iff OK.
 */

int
VnfmQ(
	VnfmData *that,		/* (I) model parameters */
	double S,		/* (I) forward rate start */
	double T,		/* (I) forward rate end */
	double *yield,		/* (O) simple yield */
	double *Q)		/* (O) Q_(S,T) */
{
	double	logz,		/* log of fwd zero */
		adj;		/* vol adjustment */
	int	j, m, n,
		nDim = NDIM;

	/*
	 * Compute A coefficients (same as VnfmA,
	 * but included here to save a function call)
	 */
	VnfmA(that, S, T, Q);

	/*
	 * Compute the forward zero Z(S,T). DO NOT USE linear
	 * interp on the rate, but constant forward (to be consistent
	 * with computation of A
	 */
	DrlDoubleArrayFloorIdx(ZTT, NZDATES, T, &m);
	DrlDoubleArrayFloorIdx(ZTT, NZDATES, S, &n);

	logz = log(that->fZero[m] / that->fZero[n])
		- RATE[m] * (T - ZTT[m])
		+ RATE[n] * (S - ZTT[n]);


	/*
	 * Adjustment for continuous compounding rate vol
	 * to simple rate vol:
	 *	vol_sc = vol_cc * (1+(T-S)*Ysc) * log(1+(T-S)*Ysc)
	 *				/ ((T-S)*Ysc)
	 *	       = vol_cc * (T-S)* Ycc * exp((T-S)*Ycc)
	 *				/ (exp((T-S)*Ycc - 1)
	 *
	 * adj = -logz / (1. - exp(logz));
	 */

	/* 
	 * 1+RT = 1/z
	 * dR = 1/(TZ) * d(log(z))
	 * Q = A / (TZ)
	 */
	*yield = (1. / exp(logz) - 1.) / (T-S);

	for (j=0; j<=nDim-1; j++) {
		Q[j] = Q[j] / (T-S) / exp(logz);
	}


#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "VnfmQ(%3d): z=%lf  zrsc=%lf  adj=%lf\n",
		__LINE__, exp(logz), pow(exp(logz),-1/(T-S))-1., adj);
	DrlFPrintf(vnfmFpLog, "VnfmQ     :");
	for (j=0; j<=nDim-1; j++) DrlFPrintf(vnfmFpLog, "%lf  ", Q[j]);
	DrlFPrintf(vnfmFpLog, "\n");
#endif
	return(0);
}


/*f-------------------------------------------------------------
 * Calculate the QB coefficients.
 *                                                             
 * <br><br>
 * Computes and returns the QB coefficients for yield volatility.
 * "tExp" is the time to reset of the option,
 * "rateReset" the reset of the rate,
 * "rateMat" the forward maturity of the rate,
 * "rateFreq" its frequency (1,2,4,12).\\
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>)
 * and calls <i> VnfmQ</i> or <i> VnfmB</i>.
 * Returns 0 iff OK.
 */

int
VnfmQBCoeff(
	VnfmData *that,		/* (I) model parameters */
	double S,		/* (I) option expiration */
	double rateReset,	/* (I) rate reset */
	double rateMat,		/* (I) rate forward maturity */
	int rateFreq,		/* (I) rate frequency (0,1,2,4,12) */
	double *yield,		/* (O) forward rate */
	double *QB)		/* (O) coefficients QB */
{
	int	idxF,
		nDim = NDIM;


#define	ONE_WEEK	0.019230769e0

	/* compute the Q coefficients */
	if (rateMat > ONE_WEEK) {
	    if (rateFreq > 0) {
		VnfmB(that, rateReset, rateMat, rateFreq, yield, QB);
	    } else {
		VnfmQ(that, rateReset, rateReset+rateMat, yield, QB);
	    }
	} else {
	    VnfmQ(that, rateReset, rateReset+rateMat, yield, QB);
	}

	if (!IS_ALMOST_ZERO(rateReset-S)) {
	    for (idxF=0; idxF<=nDim-1;idxF++)
		QB[idxF] *= exp(-(rateReset-S)*BETA[idxF]);
	}

	return(SUCCESS);
}


/*f-------------------------------------------------------------
 * Update cached rates.
 *                                                             
 * <br><br>
 * Compute and store forward rates and zero coupons prices.
 * Should not be called directly in applications
 * (use <i> VnfmComputeCoeff</i> instead.
 * Returns 0 iff OK.
 */

int
VnfmComputeFwdRates(
	VnfmData *that)		/* (I/O) model parameters */
{
static	const char	routine[] = "VnfmComputeFwdRates";
	int		status = FAILURE;
	int		i, nDaysToSpot;
	double		zStart, zEnd;

	GtoErrMsgOn();
	/*
	 *
	 */
	nDaysToSpot = that->fZcCurve->fBaseDate - REFDATE;

	/* compute forward rates */
	that->fZero[0] = 1e0;
	for (i=1; i<=NZDATES-1; i++) {
		/*
		 * Compute forward rate: fwdRate[i] is the 
		 * continuous compounding fwd rate on period [i, i+1].
		 * Defined by
		 *	exp(-fwdRate[i]*dt) = Z[i+1] / Z[i]
		 * where
		 *	dt = t(i+1)-t(i)
		 */

/*
		if (GtoDiscountDate((that->fZcCurve)->fArray[i-2].fDate + nDaysToSpot,
				that->fZcCurve,	
        GTO_LINEAR_INTERP, &zStart) != SUCCESS)
			goto done;
*/

		zStart = that->fZero[i-1];
#ifndef	VNFM_V5X
		if (GtoDiscountDate((that->fZcCurve)->fArray[i-1].fDate
				+ nDaysToSpot,
				that->fZcCurve,
				GTO_LINEAR_INTERP, &zEnd) != SUCCESS)
			goto done;
#else
		if (GtoDiscountDate(that->fDate[i] + nDaysToSpot,
				that->fZcCurve,
				GTO_LINEAR_INTERP,
				&zEnd) != SUCCESS)
			goto done;
#endif	/*VNFM_V5X*/


		that->fZero[i] = zEnd;
		that->fRate[i-1] = - log(zEnd / zStart) /
						(ZTT[i] - ZTT[i-1]);
	}

	that->fRate[NZDATES-1] = that->fRate[NZDATES-2];


#ifdef	__DEBUG__
	DrlFPrintf(vnfmFpLog, "%s(%d):\n", routine, __LINE__);

	DrlFPrintf(vnfmFpLog, "\n");
	DrlFPrintf(vnfmFpLog, "\t%4d t=%9.6f  %s  fRate=%7.4f z=%12.10f\n",
			    0,
			    ZTT[0],
			    GtoFormatDate(REFDATE),
			    that->fRate[0]*100e0,
			    that->fZero[0]);
  
	for (i=1; i<=NZDATES-1; i++) {
		DrlFPrintf(vnfmFpLog, "\t%4d t=%9.6f  %s  fRate=%7.4f z=%12.10f\n",
			i,
			ZTT[i],
			GtoFormatDate((that->fZcCurve)->fArray[i-1].fDate),
			that->fRate[i]*100e0,
			that->fZero[i]
		);
	}
#endif

	status = SUCCESS;
done:
	if (status!= SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f-------------------------------------------------------------
 * Calculate the G coefficients.
 *                                                             
 * <br><br>
 * Compute and returns $G$ coefficients for a simple rate volatility
 * defined by
 * <blockquote>
 * G(S,T) = {log Z(S,T)^{q-1} / (T-S)^{q}}.
 * </blockquote>
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>).\\
 * Returns 0 iff OK.
 */

int
VnfmG(
	VnfmData *that,		/* (I) model parameters */
	double S,		/* (I) forward rate start */
	double T,		/* (I) forward rate end */
	double *G)		/* (O) Q_(S,T) */
{
	double	logz;		/* log of fwd zero */
	int	m, n;


	/*
	 * Compute the forward zero Z(S,T). DO NOT USE linear
	 * interp on the rate, but constant forward (to be consistent
	 * with computation of A
	 */
	DrlDoubleArrayFloorIdx(ZTT, NZDATES, T, &m);
	DrlDoubleArrayFloorIdx(ZTT, NZDATES, S, &n);

	logz = log(that->fZero[m] / that->fZero[n])
		- RATE[m] * (T - ZTT[m])
		+ RATE[n] * (S - ZTT[n]);

	*G = pow(logz/(T-S), that->fBackBoneQ)/logz;

	return(0);
}


#define MIN_BETA 	1e-4
double	H(
	double	S1,
	double	S2,
	double	T,
	double	beta1,
	double	beta2)
{
	double	h;

	/* 
	 * Set miminum for beta to aviod a mysterious 
	 * numerical problem that result in h < 0 in 
	 * some cases.
	 */
	beta1 = MAX(beta1, MIN_BETA);
	beta2 = MAX(beta2, MIN_BETA);

        h = ( S2 - S1
       	    - K(T-S2, T-S1, beta1)
       	    - K(T-S2, T-S1, beta2)
       	    + K(T-S2, T-S1, beta1+beta2))
	    / (beta1 * beta2);

	return h;
}
