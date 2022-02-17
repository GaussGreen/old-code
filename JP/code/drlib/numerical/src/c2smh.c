#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  C2SMH/DC2SMH (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 15, 1986

    Purpose:    Compute a smooth cubic spline from noisy data.

    Usage:      CALL C2SMH (NDATA, XDATA, FDATA, WEIGHT, SMPAR, BREAK,
                            CSCOEF, WK, IWK)

    Arguments:  (See CSSMH)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c2smh(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mfloat *smpar, Mfloat break_[],
	   Mfloat *cscoef, Mfloat wk[], Mint iwk[])
#else
void imsl_c2smh(ndata, xdata, fdata, weight, smpar, break_,
	   cscoef, wk, iwk)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], weight[], *smpar, break_[], *cscoef,
	                wk[];
	Mint             iwk[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
	Mint             i, ifdata, ir, ir1, ir2, it, iu, iv, iweigh;
        Mint             num_weights_zero = 0;



	imsl_e1psh("IMSL_C2SMH ");
	/* Check NDATA */
	if (*ndata < 2) {
		imsl_e1sti(1, *ndata);

/*		imsl_ermes(5, 4, "The number of data points must be 2 or more while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_2_PTS);
	}
	/* Check SMPAR */
	if (*smpar < F_ZERO) {
		imsl_e1str(1, *smpar);

/*		imsl_ermes(5, 2, "The smoothing parameter must be non-negative while SMPAR = %(r1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_SMPAR_VALUE);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Partition workspace */
	ifdata = 1;
	iweigh = ifdata + *ndata;
	ir = iweigh + *ndata;
	ir1 = ir + *ndata + 1;
	ir2 = ir1 + *ndata;
	it = ir2 + *ndata + 2;
	iu = it + *ndata;
	iv = iu + *ndata + 2;
	/*
	 * Sort XDATA into BREAK, FDATA into CSCOEF(1,I) and re-arrange
	 * WEIGHT into WK(IWEIGH+I).
	 */
	for (i = 1; i <= *ndata; i++) {
		iwk[i - 1] = i;
	}

	imsl_c1sor(*ndata, xdata, fdata, break_, cscoef, 4, iwk);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	for (i = 1; i <= *ndata; i++) {
                if (weight[i-1] == F_ZERO) num_weights_zero++;
		if (weight[iwk[i - 1] - 1] < F_ZERO) {
			imsl_e1sti(1, iwk[i - 1]-1);
			imsl_e1str(1, weight[iwk[i - 1] - 1]);
/*                        imsl_ermes(4, 3, "All elements of the argument WEIGHT must be greater than or equal to zero, but WEIGHT(%(i1)) = %(r1).");
*/
                        imsl_e1stl(1, "X");
                        imsl_ermes(IMSL_FATAL, IMSL_NEGATIVE_WEIGHTS);
			goto L_9000;
		} else {
			wk[iweigh + i - 2] = weight[iwk[i - 1] - 1];
		}
	}
        if (num_weights_zero == *ndata){

/*                        imsl_ermes(5, 10, "At least one element of WEIGHT must be greater than zero, but WEIGHT(I) = 0.0 for I =0,...,NDATA-1. ");
*/
                          imsl_ermes(IMSL_TERMINAL,
			  IMSL_SPLINE_NO_POS_ELMNT);
                              
                        goto L_9000;
                }
	/*
	 * We need to keep original FDATA so use workspace.
	 */
	scopy(*ndata, cscoef, 4, &wk[ifdata - 1], 1);

	imsl_c3smh(ndata, &wk[ifdata - 1], &wk[iweigh - 1], smpar, break_,
	       cscoef, &wk[ir - 1], &wk[ir1 - 1], &wk[ir2 - 1], &wk[it - 1],
		   &wk[iu - 1], &wk[iv - 1], iwk);

L_9000:
	;
	imsl_e1pop("IMSL_C2SMH ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C3SMH/DC3SMH (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 3, 1986

    Purpose:    Compute a smooth cubic spline from noisy data.

    Usage:      CALL C3SMH (NDATA, FDATA, WEIGHT, SMPAR, BREAK, CSCOEF,
                            R, R1, R2, T, U, V, IWK)

    Arguments:
       NDATA  - Number of data points.  (Input)
                It must be at least 2.
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       WEIGHT - Array of length NDATA containing the relative weights
                of the data points.  (Input)
                Recommended values are the standard deviations of the
                ordinates FDATA.
       SMPAR  - A nonnegitive number which controls the smoothing.
                (Input)
                The spline function S returned is such that
                the sum from I=1 to NDATA of
                ((S(XDATA(I)-FDATA(I))/FWEIGH(I))**2
                is less than or equal to SMPAR.  It is recommended that
                SMPAR lie in the confidence interval of this sum, i.e.,
                NDATA-SQRT(2*NDATA) .LE. SMPAR .LE. NDATA+SQRT(2*NDATA).
       BREAK  - Array of length NDATA containing the breakpoints
                of the piecewise cubic representation.  (Output)
       CSCOEF - Matrix of size 4 by NDATA containing the local
                coefficients of the cubic pieces.  (Output)
       R      - Work vector of length NDATA.
       R1     - Work vector of length NDATA.
       R2     - Work vector of length NDATA+2.
       T      - Work vector of length NDATA.
       U      - Work vector of length NDATA+2.
       V      - Work vector of length NDATA.
       IWK    - Work vector of length NDATA.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c3smh(Mint *ndata, Mfloat fdata[], Mfloat weight[], Mfloat *smpar, Mfloat break_[], Mfloat *cscoef,
	   Mfloat r[], Mfloat r1[], Mfloat r2[], Mfloat t[], Mfloat u[], Mfloat v[], Mint iwk[])
#else
void imsl_c3smh(ndata, fdata, weight, smpar, break_, cscoef,
	   r, r1, r2, t, u, v, iwk)
	Mint            *ndata;
	Mfloat           fdata[], weight[], *smpar, break_[], *cscoef, r[],
	                r1[], r2[], t[], u[], v[];
	Mint             iwk[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
	Mint             _l0, _l1, _l2, i, icount, itmax;
	Mfloat           _f0, _f1, a[2], b0, b1, e, eps, f,
	                f2, g, h, p, sse, wk[14];


	imsl_e1psh("IMSL_C3SMH ");
	/*
	 * If SMAPR is zero, call CSDEC to get natural spline fit.
	 */
	if (*smpar == F_ZERO) {
                _l0 = 2;
                _l1 = 2;
                _f0 = F_ZERO;
                _f1 = F_ZERO;
		imsl_c2dec(ndata, break_, fdata,&_l0, &_f0, &_l1,
			  &_f1, break_, cscoef, iwk);
		goto L_9000;
	}
	icount = 0;
	/*
	 * See if SMPAR implies linear fit. The weights for FNLSQ are not the
	 * same as the weights for CSSMH, they must be inverted and squared.
	 */
	for (i = 1; i <= *ndata; i++) {
		t[i - 1] = imsl_fi_power(F_ONE / weight[i - 1], 2);
	}
	/* Compute linear least squares */
        _l0 = 1;
        _l1 = 1;
        _l2 = 1;
	imsl_f2lsq(imsl_c4smh, &_l0, &_l1, ndata, break_, fdata, &_l2,
		   t, a, &sse, wk);
	if (*smpar >= sse) {
		b0 = a[0];
		b1 = a[1];
		/* V = new F values for CSINT */
		for (i = 1; i <= *ndata; i++) {
			v[i - 1] = b0 + break_[i - 1] * b1;
		}
		imsl_c2int(ndata, break_, v, break_, cscoef, iwk);
		goto L_9000;
	}
	/* Set up working areas */
	itmax = 50;
	eps = imsl_amach(4);
	r[0] = F_ZERO;
	r[1] = F_ZERO;
	r1[*ndata - 1] = F_ZERO;
	r2[*ndata] = F_ZERO;
	r2[*ndata + 1] = F_ZERO;
	u[0] = F_ZERO;
	u[1] = F_ZERO;
	u[*ndata] = F_ZERO;
	u[*ndata + 1] = F_ZERO;
	p = F_ZERO;
	h = break_[1] - break_[0];
	f2 = -*smpar;
	f = (fdata[1] - fdata[0]) / h;
	for (i = 2; i <= (*ndata - 1); i++) {
		g = h;
		h = break_[i] - break_[i - 1];
		e = f;
		f = (fdata[i] - fdata[i - 1]) / h;
		*CSCOEF(i - 1, 0) = f - e;
		t[i - 1] = (g + h) * (F_TWO / F_THREE);
		r2[i] = weight[i - 2] / g;
		r[i] = weight[i] / h;
		r1[i - 1] = -weight[i - 1] / g - weight[i - 1] / h;
	}
	for (i = 2; i <= (*ndata - 1); i++) {
		*CSCOEF(i - 1, 1) = imsl_fi_power(r[i], 2) + imsl_fi_power(r1[i - 1], 2) + imsl_fi_power(r2[i], 2);
		*CSCOEF(i - 1, 2) = r[i] * r1[i] + r1[i - 1] * r2[i + 1];
		*CSCOEF(i - 1, 3) = r[i] * r2[i + 2];
	}
	/* Next iteration */
L_50:
	;
	for (i = 2; i <= (*ndata - 1); i++) {
		r1[i - 2] = f * r[i - 1];
		r2[i - 2] = g * r[i - 2];
		r[i] = F_ONE / (p ** CSCOEF(i - 1, 1) + t[i - 1] - f * r1[i - 2] -
			      g * r2[i - 2]);
		u[i] = *CSCOEF(i - 1, 0) - r1[i - 2] * u[i - 1] - r2[i - 2] * u[i - 2];
		f = p ** CSCOEF(i - 1, 2) + (break_[i] - break_[i - 1]) / F_THREE -
			h * r1[i - 2];
		g = h;
		h = *CSCOEF(i - 1, 3) * p;
	}
	for (i = *ndata - 1; i >= 2; i--) {
		u[i] = r[i] * u[i] - r1[i - 1] * u[i + 1] - r2[i] * u[i + 2];
	}
	e = F_ZERO;
	h = F_ZERO;
	/* Compute U and accumulate E */
	for (i = 1; i <= (*ndata - 1); i++) {
		g = h;
		h = (u[i + 1] - u[i]) / (break_[i] - break_[i - 1]);
		v[i - 1] = (h - g) * imsl_fi_power(weight[i - 1], 2);
		e += v[i - 1] * (h - g);
	}
	g = -h * imsl_fi_power(weight[*ndata - 1], 2);
	v[*ndata - 1] = g;
	e += -g * h;
	g = f2;
	f2 = p * e * p;
	icount += 1;
	if (icount >= itmax) {
                imsl_e1sti(1,itmax);

/*		imsl_ermes(3, 1, "The maximum number of iterations has been reached.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_MAX_ITERATIONS_REACHED);
		goto L_100;
	}
if (*smpar == 40.0){
    }
	if (fabs(f2 - *smpar) <= *smpar * 100.0 * eps)
		goto L_100;
	f = F_ZERO;
	h = (v[1] - v[0]) / (break_[1] - break_[0]);
	for (i = 2; i <= (*ndata - 1); i++) {
		g = h;
		h = (v[i] - v[i - 1]) / (break_[i] - break_[i - 1]);
		g = h - g - r1[i - 2] * r[i - 1] - r2[i - 2] * r[i - 2];
		f += g * r[i] * g;
		r[i] = g;
	}
	h = e - p * f;
	if (h != F_ZERO) {
		/*
		 * Update the LaGrange multiplier P for the next iteration
		 */
		p += (*smpar - f2) / ((sqrt(*smpar / e) + p) * h);
		goto L_50;
	}
	/*
	 * If E .LE. S, compute the coefficients and return.
	 */
L_100:
	;
	for (i = 1; i <= *ndata; i++) {
		*CSCOEF(i - 1, 0) = fdata[i - 1] - p * v[i - 1];
		*CSCOEF(i - 1, 2) = u[i];
	}

	for (i = 1; i <= (*ndata - 1); i++) {
		h = break_[i] - break_[i - 1];
		*CSCOEF(i - 1, 3) = (*CSCOEF(i, 2) - *CSCOEF(i - 1, 2)) / (F_THREE *
									 h);
		*CSCOEF(i - 1, 1) = (*CSCOEF(i, 0) - *CSCOEF(i - 1, 0)) / h - (h *
				 *CSCOEF(i - 1, 3) + *CSCOEF(i - 1, 2)) * h;
	}
	/* Define the rest of CSCOEF */
	*CSCOEF(*ndata - 1, 1) = F_ZERO;
	*CSCOEF(*ndata - 1, 2) = F_ZERO;
	*CSCOEF(*ndata - 1, 3) = F_ZERO;
	/* Scale for correct normalization */
	sscal(*ndata - 1, F_TWO, CSCOEF(0, 2), 4);
	sscal(*ndata - 1, F_SIX, CSCOEF(0, 3), 4);

L_9000:
	imsl_e1pop("IMSL_C3SMH ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C4SMH/DC4SMH (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 3, 1986

    Purpose:    Compute a smooth cubic spline from noisy data.

    Usage:      C4SMH(K, X)

    Arguments:
       K      - Basis function number.  (Input)
       X      - Data point where function is to be evaluated.  (Input)
       C4SMH  - Function value.  (Output)
                C4SMH is equal to X.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* K is not used, but leave the calling sequence intact. */
#ifdef ANSI
Mfloat imsl_c4smh(Mint k, Mfloat x)
#else
Mfloat imsl_c4smh(k, x)
	Mint            k;
	Mfloat          x;
#endif
{
	Mfloat           c4smh_v;


	c4smh_v = x;
	return (c4smh_v);
}				/* end of function */
