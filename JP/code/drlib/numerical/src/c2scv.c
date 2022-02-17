#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  C2SCV/DC2SCV (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 2, 1989

    Purpose:    Compute a smooth cubic spline from noisy data using
                cross-validation to estimate the smoothing parameter.

    Usage:      CALL C2SCV (NDATA, XDATA, FDATA, IEQUAL, BREAK, CSCOEF,
                            WEIGHT, WK, YWK, IWK)

    Arguments:
       NDATA  - Number of data points.  (Input)
                NDATA must be at least 3.
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
                XDATA must be distinct.
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       IEQUAL - A flag alerting the subroutine that the data is
                equally spaced.  (Input)
       BREAK  - Array of length NDATA containing the breakpoints
                for the piecewise cubic representation.  (Output)
       CSCOEF - Matrix of size 4 by NDATA containing the local
                coefficients of the cubic pieces.  (Output)
       WEIGHT - Array of lenght NDATA containing the relative
                weights.  (Input)
       WK     - Work array of length 7*(NDATA+2).
       YWK    - Work array of length NDATA.
       IWK    - Work array of length NDATA.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* IEQUAL is not used , but leave the calling sequence intact. */
#ifdef ANSI
void imsl_c2scv(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mint *iequal, Mfloat break_[], Mfloat *cscoef,
	   Mfloat weight[], Mfloat *wk, Mfloat ywk[], Mint iwk[])
#else
void imsl_c2scv(ndata, xdata, fdata, iequal, break_, cscoef,
	   weight, wk, ywk, iwk)
	Mint            *ndata;
	Mfloat           xdata[], fdata[];
	Mint            *iequal;
	Mfloat           break_[], *cscoef, weight[], *wk, ywk[];
	Mint             iwk[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
   
#define WK(I_,J_)	(wk+(I_)*(*ndata + 1-(0)+1)+(J_))
	Mint             i, j;
        Mint             num_weights_zero = 0;
	Mfloat           amach4, avar, avdf, avh, delta, err, gf1, gf2,
	                gf3, gf4, p, q, r1, r2, r3, r4, stat[6], var;
	static Mfloat    ratio = 2.0;
	static Mfloat    tau = 1.618033989;



	imsl_e1psh("IMSL_C2SCV");

	var = -F_ONE;

        /* CHECK WEIGHTS */
        for (i = 1; i <= *ndata; i++) {
                if (weight[i-1] == F_ZERO) num_weights_zero++;
                if (weight[i - 1] < F_ZERO) {
                        imsl_e1sti(1, i-1);
                        imsl_e1str(1, weight[i - 1]);

/*                        imsl_ermes(4, 7, "All elements of the argument WEIGHT must be greater than or equal to zero, but WEIGHT(%(i1)) = %(r1).");
*/
                        imsl_e1stl(1, "X");
                        imsl_ermes(IMSL_FATAL, IMSL_NEGATIVE_WEIGHTS);
                        goto L_9000;
                }
        }
        if (num_weights_zero == *ndata){

/*                        imsl_ermes(5, 10, "At least one element of WEIGHT must be greater than zero, but WEIGHT(I) = 0.0 for I =0,...,NDATA-1. ");
*/
                        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_NO_POS_ELMNT);
                              
                        goto L_9000;
                }

	/* Check for sorted array */
	for (i = 2; i <= *ndata; i++) {
		if (xdata[i - 2] >= xdata[i - 1]) {
			/* Check that XDATA values are distinct */
			if (xdata[i - 2] == xdata[i - 1]) {
				j = i - 1;
				imsl_e1sti(1, j-1);
				imsl_e1sti(2, i-1);
				imsl_e1str(1, xdata[i - 1]);

/*				imsl_ermes(4, 2, "Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_DUPLICATE_XDATA_VALUES);
				goto L_9000;
			} else {
				goto L_20;
			}
		}
	}
	/*
	 * Data is already sorted.  Move to BREAK, DFWK  and CSCOEF.
	 */
	scopy(*ndata, xdata, 1, break_, 1);
	scopy(*ndata, fdata, 1, cscoef, 4);
	goto L_60;
	/* Set initial permutation */
L_20:
	for (i = 1; i <= *ndata; i++) {
		iwk[i - 1] = i;
	}
	/* Find sorting permutation */
	imsl_svrgp(*ndata, xdata, break_, iwk);
	/* Apply permutation */
	for (i = 1; i <= *ndata; i++) {
		*CSCOEF(i - 1, 0) = fdata[iwk[i - 1] - 1];
	}
	/* Check the XDATA values are distinct */
	for (i = 2; i <= *ndata; i++) {
		if (break_[i - 2] == break_[i - 1]) {
			imsl_e1sti(1, iwk[i - 2]-1);
			imsl_e1sti(2, iwk[i - 1]-1);
			imsl_e1str(1, break_[i - 1]);

/*			imsl_ermes(4, 2, "Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
*/
                        imsl_ermes(IMSL_FATAL, IMSL_DUPLICATE_XDATA_VALUES);
			goto L_9000;
		}
	}

L_60:
	amach4 = imsl_amach(4);

	imsl_c3scv(break_, &avh, iwk, weight, &avdf, ndata, ywk, cscoef, wk,
		   WK(3, 0));

	if (imsl_n1rty(0) != 0)
		goto L_9000;
	avar = var;
	/*
	 * Find local minimum of GCV or the expected mean square error
	 */
	r1 = F_ONE;
	r2 = ratio * r1;
	/* CALL SPFIT1 */
	imsl_c4scv(break_, &avh, weight, ndata, &r2, &p, &q, &gf2, &avar,
		   stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));
L_70:
	imsl_c4scv(break_, &avh, weight, ndata, &r1, &p, &q, &gf1, &avar,
		   stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));
	if (gf1 < gf2) {
		/* Exit if P zero */
		if (p <= F_ZERO)
			goto L_120;
		r2 = r1;
		gf2 = gf1;
		r1 /= ratio;
		goto L_70;
	}
	r3 = ratio * r2;
L_80:
	imsl_c4scv(break_, &avh, weight, ndata, &r3, &p, &q, &gf3, &avar,
		   stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));
	if (gf3 < gf2) {
		/* Exit if Q zero */
		if (q <= F_ZERO)
			goto L_120;
		r2 = r3;
		gf2 = gf3;
		r3 *= ratio;
		goto L_80;
	}
	r2 = r3;
	gf2 = gf3;
	delta = (r2 - r1) / tau;
	r4 = r1 + delta;
	r3 = r2 - delta;
	imsl_c4scv(break_, &avh, weight, ndata, &r3, &p, &q, &gf3, &avar,
		   stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));
	imsl_c4scv(break_, &avh, weight, ndata, &r4, &p, &q, &gf4, &avar,
		   stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));

	/*
	 * Golden section search for local minimum
	 */
L_90:
	if (gf3 < gf4) {
		r2 = r4;
		gf2 = gf4;
		r4 = r3;
		gf4 = gf3;
		delta /= tau;
		r3 = r2 - delta;
		imsl_c4scv(break_, &avh, weight, ndata, &r3, &p, &q, &gf3, &avar,
		       stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));
		goto L_100;
	}
	r1 = r3;
	gf1 = gf3;
	r3 = r4;
	gf3 = gf4;
	delta /= tau;
	r4 = r1 + delta;
	imsl_c4scv(break_, &avh, weight, ndata, &r4, &p, &q, &gf4, &avar,
		   stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));
L_100:
	err = (r2 - r1) / (r1 + r2);
	if (err * err + F_ONE > F_ONE && err > amach4)
		goto L_90;
	r1 = (r1 + r2) * F_HALF;
	/*
	 * Calculate spline coefficients
	 */

	imsl_c4scv(break_, &avh, weight, ndata, &r1, &p, &q, &gf1, &avar,
		   stat, ywk, cscoef, wk, WK(3, 0), WK(5, 0), WK(6, 0));
	/* CALL SPCOF1 */
L_120:
	imsl_c5scv(break_, &avh, weight, ndata, &p, &q, ywk, cscoef, WK(5, 0),
		   WK(6, 0));
	/*
	 * Optionally calculate standard error estimates
	 */
	if (var < F_ZERO) {
		avar = stat[5];
		var = avar / (avdf * avdf);
	}
	scopy(6, stat, 1, WK(0, 0), 1);
	*WK(0, 5) = stat[5] / (avdf * avdf);
	*WK(0, 6) = avdf * avdf;

	*CSCOEF(*ndata - 1, 0) = F_ZERO;
	*CSCOEF(*ndata - 1, 1) = F_ZERO;
	*CSCOEF(*ndata - 1, 2) = F_ZERO;
	*CSCOEF(*ndata - 1, 3) = F_ZERO;
	scopy(*ndata, ywk, 1, cscoef, 4);
	sscal(*ndata, F_TWO, CSCOEF(0, 2), 4);
	sscal(*ndata, F_SIX, CSCOEF(0, 3), 4);
	*CSCOEF(*ndata - 1, 1) = F_ZERO;
	*CSCOEF(*ndata - 1, 2) = F_ZERO;
	*CSCOEF(*ndata - 1, 3) = F_ZERO;
L_9000:
	imsl_e1pop("IMSL_C2SCV");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C3SCV/DC3SCV (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 27, 1989

    Purpose:    Compute a smooth cubic spline from noisy data using
                cross-validation to estimate the smoothing parameter.

    Usage:      CALL C3SCV (XDATA, AVH, IWK, DFCOPY, AVDY, NDATA,
                            YDATA, CSCOEF, R, T)

    Arguments:  (See CSSCV)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* IWK is not used, but leave the calling sequence intact. */
#ifdef ANSI
void imsl_c3scv(Mfloat xdata[], Mfloat *avh, Mint iwk[], Mfloat dfcopy[], Mfloat *avdy, 
                Mint *ndata, Mfloat ydata[],Mfloat *cscoef, Mfloat *r, Mfloat *t)
#else
void imsl_c3scv(xdata, avh, iwk, dfcopy, avdy, ndata, ydata,
	   cscoef, r, t)
	Mfloat           xdata[], *avh;
	Mint             iwk[];
	Mfloat           dfcopy[], *avdy;
	Mint            *ndata;
	Mfloat           ydata[], *cscoef, *r, *t;
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
   
#define R(I_,J_)	(r+(I_)*(*ndata + 1-(0)+1)+(J_))
#define T(I_,J_)	(t+(I_)*(*ndata + 1-(0)+1)+(J_))
	Mint             i;
	Mfloat           e, f, g, h;


	imsl_e1psh("IMSL_C3SCV ");
	/* Check NDATA */
	if (*ndata <= 2) {
		imsl_e1sti(1, *ndata);

/*		imsl_ermes(5, 1, "The number of data points must be at least 3 while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_3_PTS);
		goto L_9000;
	}
	/* Get average XDATA spacing in AVH */
	g = F_ZERO;
	for (i = 1; i <= (*ndata - 1); i++) {
		h = xdata[i] - xdata[i - 1];
		g += h;
	}
	*avh = g / (*ndata - 1);
	/* Scale relative weights */
	g = imsl_sdot(*ndata, dfcopy, 1, dfcopy, 1);
	*avdy = sqrt(g / *ndata);
	sscal(*ndata, F_ONE / *avdy, dfcopy, 1);
	/* Initialize H,F--- */
	h = (xdata[1] - xdata[0]) / *avh;
	/* F = (FDATA(2)-FDATA(1))/H */
	f = (*CSCOEF(1, 0) - *CSCOEF(0, 0)) / h;
	/* Calculate A,T,R--- */
	for (i = 2; i <= (*ndata - 1); i++) {
		g = h;
		h = (xdata[i] - xdata[i - 1]) / *avh;
		e = f;
		/* F = (FDATA(I+1)-FDATA(I))/H */
		f = (*CSCOEF(i, 0) - *CSCOEF(i - 1, 0)) / h;
		ydata[i - 1] = f - e;
		*T(0, i) = F_TWO * (g + h) / F_THREE;
		*T(1, i) = h / F_THREE;
		*R(2, i) = dfcopy[i - 2] / g;
		*R(0, i) = dfcopy[i] / h;
		*R(1, i) = -dfcopy[i - 1] / g - dfcopy[i - 1] / h;
	}
	/* Calculate C = R'*R--- */
	*R(1, *ndata) = F_ZERO;
	*R(2, *ndata) = F_ZERO;
	*R(2, *ndata + 1) = F_ZERO;
	for (i = 2; i <= (*ndata - 1); i++) {
		*CSCOEF(i - 1, 1) = *R(0, i) ** R(0, i) + *R(1, i) ** R(1, i) + *R(2, i) *
			*R(2, i);
		*CSCOEF(i - 1, 2) = *R(0, i) ** R(1, i + 1) + *R(1, i) ** R(2, i + 1);
		*CSCOEF(i - 1, 3) = *R(0, i) ** R(2, i + 2);
	}

L_9000:
	imsl_e1pop("IMSL_C3SCV ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C4SCV/DC4SCV (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 27, 1989

    Purpose:    Compute a smooth cubic spline from noisy data using
                cross-validation to estimate the smoothing parameter.

    Usage:      CALL C4SCV (XDATA, AVH, DFCOPY, NDATA, RHO, P, Q, FUN,
                            VAR, STAT, YDATA, CSCOEF, R, T, U, V)

    Arguments:  (See CSSCV)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* VAR is not used, but leave the calling sequence intact. */
#ifdef ANSI
void imsl_c4scv(Mfloat xdata[], Mfloat *avh, Mfloat dfcopy[], Mint *ndata,
                Mfloat *rho, Mfloat *p, Mfloat *q, Mfloat *fun,
	        Mfloat *var, Mfloat stat[], Mfloat ydata[], Mfloat *cscoef,
                Mfloat *r, Mfloat *t, Mfloat u[], Mfloat v[])
#else
void imsl_c4scv(xdata, avh, dfcopy, ndata, rho, p, q, fun,
	   var, stat, ydata, cscoef, r, t, u, v)
	Mfloat           xdata[], *avh, dfcopy[];
	Mint            *ndata;
	Mfloat          *rho, *p, *q, *fun, *var, stat[], ydata[], *cscoef,
	               *r, *t, u[], v[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
   
#define R(I_,J_)	(r+(I_)*(*ndata + 1-(0)+1)+(J_))
#define T(I_,J_)	(t+(I_)*(*ndata + 1-(0)+1)+(J_))
	Mint             i;
	Mfloat           e, f, g, h, rho1;

	/*
	 * Use P and Q instead of RHO to prevent overflow or underflow
	 */
	rho1 = F_ONE + *rho;
	*p = *rho / rho1;
	*q = F_ONE / rho1;
	if (rho1 == F_ONE)
		*p = F_ZERO;
	if (rho1 == *rho)
		*q = F_ZERO;
	/*
	 * Rational cholesky decomposition of P*C + Q*T
	 */
	f = F_ZERO;
	g = F_ZERO;
	h = F_ZERO;
	*R(0, 0) = F_ZERO;
	*R(0, 1) = F_ZERO;
	for (i = 2; i <= (*ndata - 1); i++) {
		*R(2, i - 2) = g ** R(0, i - 2);
		*R(1, i - 1) = f ** R(0, i - 1);
		*R(0, i) = F_ONE / (*p ** CSCOEF(i - 1, 1) + *q ** T(0, i) - f ** R(1, i - 1) -
				  g ** R(2, i - 2));
		f = *p ** CSCOEF(i - 1, 2) + *q ** T(1, i) - h ** R(1, i - 1);
		g = h;
		h = *p ** CSCOEF(i - 1, 3);
	}
	/* Solve for U */
	u[0] = F_ZERO;
	u[1] = F_ZERO;
	for (i = 2; i <= (*ndata - 1); i++) {
		u[i] = ydata[i - 1] - *R(1, i - 1) * u[i - 1] - *R(2, i - 2) * u[i - 2];
	}
	u[*ndata] = F_ZERO;
	u[*ndata + 1] = F_ZERO;
	for (i = *ndata - 1; i >= 2; i--) {
		u[i] = *R(0, i) * u[i] - *R(1, i) * u[i + 1] - *R(2, i) * u[i + 2];
	}
	/* Calculate residual vector V */
	e = F_ZERO;
	h = F_ZERO;
	for (i = 1; i <= (*ndata - 1); i++) {
		g = h;
		h = (u[i + 1] - u[i]) / ((xdata[i] - xdata[i - 1]) / *avh);
		v[i] = dfcopy[i - 1] * (h - g);
		e += v[i] * v[i];
	}
	v[*ndata] = dfcopy[*ndata - 1] * (-h);
	e += v[*ndata] * v[*ndata];
	/*
	 * Calculate upper three bands of inverse matrix
	 */
	*R(0, *ndata) = F_ZERO;
	*R(1, *ndata) = F_ZERO;
	*R(0, *ndata + 1) = F_ZERO;
	for (i = *ndata - 1; i >= 2; i--) {
		g = *R(1, i);
		h = *R(2, i);
		*R(1, i) = -g ** R(0, i + 1) - h ** R(1, i + 1);
		*R(2, i) = -g ** R(1, i + 1) - h ** R(0, i + 2);
		*R(0, i) += -g ** R(1, i) - h ** R(2, i);
	}
	/* Calculate trace */
	f = F_ZERO;
	g = F_ZERO;
	h = F_ZERO;
	for (i = 2; i <= (*ndata - 1); i++) {
		f += *R(0, i) ** CSCOEF(i - 1, 1);
		g += *R(1, i) ** CSCOEF(i - 1, 2);
		h += *R(2, i) ** CSCOEF(i - 1, 3);
	}
	f += F_TWO * (g + h);
	/* Calculate statistics */
	stat[0] = *p;
	stat[1] = f ** p;
	stat[2] = *ndata * e / (f * f);
	stat[3] = e ** p ** p / *ndata;
	stat[5] = e ** p / f;
	stat[4] = stat[5] - stat[3];
	*fun = stat[2];
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C5SCV/DC5SCV (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 27, 1989

    Purpose:    Compute a smooth cubic spline from noisy data using
                cross-validation to estimate the smoothing parameter.

    Usage:      CALL C5SCV (XDATA, AVH, DFCOPY, NDATA, P, Q, YDATA,
                            CSCOEF, U, V)

    Arguments:  (See CSSCV)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c5scv(Mfloat xdata[], Mfloat *avh, Mfloat dfcopy[], Mint *ndata, Mfloat *p,
                Mfloat *q, Mfloat ydata[], Mfloat *cscoef, Mfloat u[], Mfloat v[])
#else
void imsl_c5scv(xdata, avh, dfcopy, ndata, p, q, ydata, cscoef,
	   u, v)
	Mfloat           xdata[], *avh, dfcopy[];
	Mint            *ndata;
	Mfloat          *p, *q, ydata[], *cscoef, u[], v[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
   
	Mint             i;
	Mfloat           h, qh;

	/* Calculate YDATA */
	qh = *q / (*avh ** avh);
	/*
	 * Maybe done as follows: SSCAL(NDATA,P,DFCOPY,1)
	 * SHPROD(NDATA,DFCOPY,1,V,1,DFCOPY,1) SCOPY(NDATA,CSCOEF,4,YDATA,1)
	 * SAXPY(NDATA,-1.,DFCOPY,1,YDATA,1) SSCAL(NDATA,QH,U,1)
	 */
	for (i = 1; i <= *ndata; i++) {
		ydata[i - 1] = *CSCOEF(i - 1, 0) - *p * dfcopy[i - 1] * v[i];
		u[i] *= qh;
	}
	/* Calculate C */
	for (i = 1; i <= (*ndata - 1); i++) {
		h = xdata[i] - xdata[i - 1];
		*CSCOEF(i - 1, 3) = (u[i + 1] - u[i]) / (F_THREE * h);
		*CSCOEF(i - 1, 1) = (ydata[i] - ydata[i - 1]) / h - (h ** CSCOEF(i - 1, 3) +
								  u[i]) * h;
		*CSCOEF(i - 1, 2) = u[i];
	}
	return;
}				/* end of function */
