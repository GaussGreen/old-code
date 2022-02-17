#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mint PROTO(l_ismax,(Mint, Mfloat[], Mint));
/*  -----------------------------------------------------------------------
    IMSL Name:  B2VLS/DB2VLS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 13, 1986

    Purpose:    Compute the variable knot B-spline least squares to
                given data.

    Usage:      CALL B2VLS (NDATA, XDATA, FDATA, WEIGHT, KORDER, NCOEF,
                            XGUESS, XKNOT, BSCOEF, SSQ, IWK, WK)

    Arguments:
       NDATA  - Number of data points.  (Input)
                NDATA must be at least 2.
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       WEIGHT - Array of length NDATA containing the weights.  (Input)
       KORDER - Order of the spline.  (Input)
                KORDER must be less than or equal to NDATA.
       NCOEF  - Number of B-spline coefficients.  (Input)
                NCOEF must be less than or equal to NDATA.
       XGUESS - Array of length NCOEF+KORDER containing the initial
                guess of knots.  (Input)
                XGUESS must be nonincreasing.
       XKNOT  - Array of length NCOEF+KORDER containing the
                (nondecreasing) knot sequence.  (Output)
       BSCOEF - Array of length NCOEF containing the B-spline
                representation.  (Output)
       SSQ    - The square root of the sum of the squares of the error.
       IWK    - Work array of length NDATA.
       WK     - Work array of length NCOEF*(6+2*KORDER)+
                KORDER*(7-KORDER)+3*NDATA+3.
    Remark:
       Informational errors
       Type Code
         3   12 The knots found to be optimal are stacked more than
                KORDER.  This indicates fewer knots will produce the
                same error sum of squares.  The knots have been
                separated slightly.
         4   9  The multiplicity of the knots in XGUESS cannot exceed
                the order of the spline.
         4   10 XGUESS must be nondecreasing.

    GAMS:       K1a1a1; K1b

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b2vls(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mint *korder, Mint *ncoef,
	   Mfloat xguess[], Mfloat xknot[], Mfloat bscoef[], Mfloat *ssq, Mint iwk[], Mfloat wk[])
#else
void imsl_b2vls(ndata, xdata, fdata, weight, korder, ncoef,
	   xguess, xknot, bscoef, ssq, iwk, wk)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], weight[];
	Mint            *korder, *ncoef;
	Mfloat           xguess[], xknot[], bscoef[], *ssq;
	Mint             iwk[];
	Mfloat           wk[];
#endif
{
	Mint             i, ibrk, idd, ifdata, ippcf, irwk, istep, ittt,
	                iweigh, ixdata, mult;


	imsl_e1psh("IMSL_B2VLS ");
	/* Check KORDER */
	if (*korder < 1) {
		imsl_e1sti(1, *korder);

/*		imsl_ermes(5, 16, "The order of the spline must be at least 1 while KORDER = %(i1) is given.");
*/
		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
		goto L_9000;
	}
	/* Check NCOEF */
	if (*ncoef <= *korder) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *korder);

/*		imsl_ermes(5, 17, "The number of coefficients must be larger than the order of the spline while NCOEF = %(i1) and KORDER = %(i2) are given.");
*/
		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS);
		goto L_9000;
	}
	/* Check NDATA */
	if (*ncoef > *ndata) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *ndata);

/*		imsl_ermes(5, 18, "The number of coefficients must be at less than or equal to the number of data points while NCOEF = %(i1) and NDATA = %(i2) are given.");
*/
		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS_2); 
		goto L_9000;
	}
	/* Check initial knot sequence XGUESS */
	mult = 1;
	for (i = 2; i <= (*ncoef + *korder); i++) {
		if (xguess[i - 1] == xguess[i - 2]) {
			mult += 1;
			if (mult > *korder) {
				imsl_e1sti(1, i-1 - mult + 1);
				imsl_e1sti(2, i-1);
				imsl_e1str(1, xguess[i - 1]);
				imsl_e1sti(3, *korder);

/*				imsl_ermes(4, 9, "The knots XGUESS(%(i1)) through XGUESS(%(i2)) are all equal to %(r1).  The multiplicity of the knots must not exceed KORDER = %(i3).");
*/
				imsl_ermes(IMSL_FATAL,
				IMSL_XGUESS_MULTIPLICITY);
				goto L_9000;
			}
		} else if (xguess[i - 1] < xguess[i - 2]) {
			imsl_e1sti(1, i - 2);
			imsl_e1sti(2, i - 1);
			imsl_e1str(1, xguess[i - 2]);
			imsl_e1str(2, xguess[i - 1]);

/*			imsl_ermes(4, 10, "The knot XGUESS(%(i1)) = %(r1) and XGUESS(%(i2)) = %(r2).  The knots must be nondecreasing.");
*/
			imsl_ermes(IMSL_FATAL,IMSL_XGUESS_NOT_INCREASING);
			goto L_9000;
		} else {
			mult = 1;
		}
	}
	/* Check arguments */
	imsl_b3lsq(ndata, korder, xguess, weight, ncoef);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Partition workspace */
	idd = 1;
	ibrk = idd + 3 * (*ncoef + *korder);
	ippcf = ibrk + *ncoef - *korder + 2;
	/*
	 * STEP also holds BRNEW which, if K=1 is of length NCOEF+1.
	 */
	istep = ippcf + *korder * (*ncoef - *korder + 1);
	ittt = istep + *ncoef + 1;
	ixdata = ittt + *ncoef + *korder;
	ifdata = ixdata + *ndata;
	iweigh = ifdata + *ndata;
	irwk = iweigh + *ndata;
	/*
	 * Sort XDATA into XSORT and FDATA into FSORT. Set initial
	 * permutation
	 */
	for (i = 1; i <= *ndata; i++) {
		iwk[i - 1] = i;
	}
	imsl_svrgp(*ndata, xdata, &wk[ixdata - 1], iwk);
	/*
	 * NOTE: A check for unique XDATA was removed 11/89. Rearrange FDATA
	 * and WEIGHT.
	 */
	wk[ifdata - 1] = fdata[iwk[0] - 1];
	wk[iweigh - 1] = weight[iwk[0] - 1];
	for (i = 2; i <= *ndata; i++) {
		wk[ifdata + i - 2] = fdata[iwk[i - 1] - 1];
		wk[iweigh + i - 2] = weight[iwk[i - 1] - 1];
	}

	imsl_b3vls(ndata, &wk[ixdata - 1], &wk[ifdata - 1], &wk[iweigh - 1],
		   korder, ncoef, xguess, xknot, bscoef, ssq, &wk[idd - 1], &wk[ibrk - 1],
	      &wk[ippcf - 1], &wk[istep - 1], &wk[ittt - 1], &wk[irwk - 1]);
L_9000:
	imsl_e1pop("IMSL_B2VLS ");
	return;
}				/* end of function */

/* Structured by FOR_STRUCT, v0.2, on 08/28/90 at 17:13:45
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  B3VLS/DB3VLS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 13, 1986

    Purpose:    Compute the variable knot B-spline least squares to
                given data.

    Usage:      CALL B3VLS (NDATA, XDATA, FDATA, WEIGHT, KORDER, NCOEF,
                            XGUESS, XKNOT, BSCOEF, SSQ, DD, BREAK,
                            PPCOEF, STEP, TTT, WK)

    Arguments:
       NDATA  - Number of data points (must be at least 2).  (Input)
       XDATA  - Array of length NDATA containing the data point
                abscissae.  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       WEIGHT - Array of length NDATA containing the weights.  (Input)
       KORDER - Order of the spline.  (Input)
                KORDER must be less than or equal to NDATA.
       NCOEF  - Number of spline coefficients.  (Input)
                It cannot be greater than NDATA.
       XGUESS - Vector of length NCOEF+KORDER containing the intial
                guess of knots.  (Input)
       XKNOT  - Array of length NCOEF+KORDER containing the
                (nondecreasing) knot sequence.  (Output)
       BSCOEF - Array of length NCOEF containing the B-spline
                representation.  (Output)
       SSQ    - The square root of the sum of the squares of the error.
                (Output)
       DD     - Work array of size 3*(NCOEF+KORDER).
       BREAK  - Work array of length NCOEF-KORDER+2 to hold a breakpoint
                vector.
       PPCOEF - Work array of size KORDER*(NCOEF-KORDER+1) to hold
                piecewise polynomial coefficients.
       STEP   - Work arrary of length NCOEF to hols step size.
       TTT    - Work array of length NCOEF+KORDER that holds a
                copy of XKNOT.
       WK     - Work vector of length KORDER*(NCOEF+3).

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b3vls(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mint *korder, Mint *ncoef,
           Mfloat xguess[], Mfloat xknot[], Mfloat bscoef[], Mfloat *ssq, Mfloat *dd, Mfloat break_[],
           Mfloat *ppcoef, Mfloat step[], Mfloat ttt[], Mfloat wk[])
#else

void imsl_b3vls(ndata, xdata, fdata, weight, korder, ncoef,
	   xguess, xknot, bscoef, ssq, dd, break_, ppcoef, step, ttt, wk)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], weight[];
	Mint            *korder, *ncoef;
	Mfloat           xguess[], xknot[], bscoef[], *ssq, *dd, break_[],
	               *ppcoef, step[], ttt[], wk[];
#endif
{
#define DD(I_,J_)	(dd+(I_)*(3)+(J_))
#define PPCOEF(I_,J_)	(ppcoef+(I_)*(*korder)+(J_))
	Mint             _d_l, _d_m, _do0, _do1, i, icount, ido,
	                iii, iiii, j, kk, knot, m1, m2, m3,
	                maxfn, mult, nf, nloop, numgt, numlt;
	Mfloat           a, b, curerr, eps, eps2, scale1,
	                scale2, ss, ssnew, ssorig, ssqold, steps, tol,
	                x, xacc, xguss, znum;


	imsl_e1psh("IMSL_B3VLS");
	/*
	 * SCALE DATA SO THE PROBLEM IS SOLVED ON THE CLOSED INTERVAL (0,1).
	 */
	scale1 = xdata[*ndata - 1] - xdata[0];
	scale2 = xdata[0];
	/* SCALE XDATA */
	for (i = 1; i <= *ndata; i++) {
		xdata[i - 1] = (xdata[i - 1] - scale2) / scale1;
	}
	/* SCALE XGUESS */
	for (i = 1; i <= (*ncoef + *korder); i++) {
		xguess[i - 1] = (xguess[i - 1] - scale2) / scale1;
	}

	ssqold = imsl_amach(2);
	znum = (Mfloat) (*ncoef - *korder);
	if (znum == F_ZERO)
		znum = F_ONE;
	eps = imsl_amach(4);
	imsl_scopy(*ncoef + *korder, xguess, 1, xknot, 1);
	/*
	 * TRY TO GET A BETTER INITIAL KNOT SEQUENCE
	 */
	imsl_b4lsq(ndata, xdata, fdata, weight, korder, xknot, ncoef, bscoef,
		   &wk[0], &wk[*korder ** ncoef]);
	ssorig = imsl_b6vls(ndata, xdata, fdata, weight, korder, xknot, ncoef,
			    bscoef, ppcoef, break_, wk);
	nloop = 15;
	imsl_b4vls(&nloop, ndata, xdata, fdata, weight, korder, xknot, ncoef,
		   bscoef, ppcoef, break_, wk, step, dd);
	ssnew = imsl_b6vls(ndata, xdata, fdata, weight, korder, xknot, ncoef,
			   bscoef, ppcoef, break_, wk);
	/*
	 * IF NEW KNOT SEQUNCE DOES NOT GIVE A SMALLER SUM OF SQUARES THEN
	 * USE THE ORIGINAL KNOT SEQUENCE.
	 */
	if (ssorig < ssnew) {
		imsl_scopy(*ncoef + *korder, xguess, 1, xknot, 1);
	}
	/*
	 * NOW MAKE SURE THAT THE SMALLEST INTERIOR KNOT IS GREATER THAN THE
	 * MINIMUM XDATA VALUE, AND THAT THE LARGEST INTERIOR KNOT IS LESS
	 * THAN THE LARGEST XDATA VALUE.
	 */
	eps2 = sqrt(imsl_amach(4));
	/*
	 * FIND THE NUMBER OF KNOTS OUTSIDE THE THE INTERVAL CONTAINING THE
	 * XDATA VALUES.
	 */
	numlt = 0;
	numgt = 0;
	for (i = *korder + 1; i <= *ncoef; i++) {
		if (xknot[i - 1] <= (xdata[0] + *ncoef * eps2))
			numlt += 1;
		if (xknot[i - 1] >= (xdata[*ndata - 1] - *ncoef * eps2))
			numgt += 1;
	}
	/*
	 * MOVE THE KNOTS THAT WERE OUTSIDE OF THE XDATA VALUES SO THAT WILL
	 * BE INSIDE THE XDATA VALUES.
	 */
	for (i = 1; i <= numlt; i++) {
		xknot[*korder + i - 1] = xdata[0] + i * (eps2);
	}
	for (i = 1; i <= numgt; i++) {
		xknot[*ncoef - i] = xdata[*ndata - 1] - i * (eps2);
	}

	tol = 100.0 * eps * (xknot[*ncoef] - xknot[*korder - 1]) / (F_TWO * znum);
	/* LOOP 70 IS THE NUMBER OF ITERATIONS */
	for (iii = 1; iii <= 2; iii++) {
		for (knot = *korder + 1; knot <= *ncoef; knot++) {
			step[knot - 1] = xknot[knot] - xknot[knot - 2];
		}
		/*
		 * LOOP 60 IS THE NUMBER OF PASSES THROUGH THE KNOTS FOR EACH
		 * ITERATION.
		 */
		for (iiii = 1; iiii <= 2; iiii++) {
			if (((iiii + iii + 1) / 2) * 2 == (iiii + iii + 1)) {
				m1 = *korder + 1;
				m2 = *ncoef;
				m3 = 1;
			} else {
				m1 = *ncoef;
				m2 = *korder + 1;
				m3 = -1;
			}
			/*
			 * STORE COPY OF XKNOT IN TOP ROW OF DD TO BE USED IN
			 * THE ACCELERATION ROUTINE.
			 */
			imsl_scopy(*ncoef + *korder, xknot, 1, DD(0, 0), 3);
			for (icount = 1; icount <= 2; icount++) {
				for (knot = m1, _do0 = DOCNT(m1, m2, _do1 = m3); _do0 > 0; knot += _do1, _do0--) {
					steps = step[knot - 1] / (F_TWO + imsl_fi_power(F_TWO, iii));
					if (steps > tol) {
						ido = 0;
						a = xknot[knot - 2];
						if (knot == *korder + 1)
							a = xdata[0] + eps2;
						b = xknot[knot];
						if (knot == *ncoef)
							b = xdata[*ndata - 1] - eps2;
						xguss = xknot[knot - 1];
						xacc = imsl_f_max(100 * eps, steps) / (imsl_fi_power(F_TWO, iii +
									2));
						maxfn = 10 + (iii - 1) * 4;
						/*
						 * THIS ROUTINE USES REVERSE
						 * COMMUN.
						 */
				L_70:
						imsl_b7vls(&ido, ssq, &a, &b, &xguss, &steps,
							 &xacc, &maxfn, &x);
						if (ido == 1) {
							xknot[knot - 1] = x;
							imsl_b4lsq(ndata, xdata, fdata, weight, korder,
								   xknot, ncoef, bscoef, &wk[0], &wk[*korder ** ncoef]);
							*ssq = imsl_b6vls(ndata, xdata, fdata, weight,
									  korder, xknot, ncoef, bscoef, ppcoef,
								break_, wk);
							goto L_70;
						}
						xknot[knot - 1] = x;
					}
				}
				for (knot = *korder + 1; knot <= *ncoef; knot++) {
					step[knot - 1] = fabs(xknot[knot - 1] - *DD(knot - 1, 0));
				}
				j = l_ismax(*ncoef - *korder, &step[*korder],
					       1) + *korder;
				ss = step[j - 1] / F_FOUR;
				if (ss <= tol)
					goto L_120;
				imsl_scopy(*ncoef + *korder, xknot, 1, DD(0, icount), 3);
			}
			curerr = fabs(ssqold - *ssq);
			if (curerr <= 100.0 * eps ** ssq)
				goto L_120;
			ssqold = *ssq;
			/*
			 * WE NOW TRY TO ACCELLARATE THE KNOT VECTOR USING
			 * POLNOMIAL INTERPOLATION.
			 */
			imsl_b8vls(dd, korder, ncoef, ttt);
			imsl_b4lsq(ndata, xdata, fdata, weight, korder, ttt, ncoef,
				   bscoef, &wk[0], &wk[*korder ** ncoef]);
			*ssq = imsl_b6vls(ndata, xdata, fdata, weight, korder, ttt,
					  ncoef, bscoef, ppcoef, break_, wk);
			/*
			 * NOW THE OLD AND NEW VALUES FOR SSQ ARE COMPARED.
			 * IF THE OLD VALUE IS SMALLER THEN NO ACCELLARATION
			 * OCURRS.
			 */
			if (*ssq < ssqold) {
				imsl_scopy(*ncoef + *korder, ttt, 1, xknot, 1);
			} else {
				*ssq = ssqold;
			}
		}
L_120:
		;
	}
	/* CHECK FINAL KNOT SEQUENCE XKNOT */
	imsl_scopy(*ncoef + *korder, xknot, 1, ttt, 1);
	nf = 0;
L_130:
	mult = 1;
	for (i = 2; i <= (*ncoef + *korder); i++) {
		if (ttt[i - 1] == ttt[i - 2]) {
			mult += 1;
			if (mult > *korder) {

/*				imsl_ermes(3, 12, "THE KNOTS FOUND TO BE OPTIMAL ARE STACKED MORE THAN KORDER.  THIS INDICATES FEWER KNOTS WILL PRODUCE THE SAME ERROR SUM OF SQUARES.  THE KNOTS HAVE BEEN SEPARATED SLIGHTLY.");
*/
                                imsl_ermes(IMSL_WARNING,
				IMSL_OPT_KNOTS_STACKED_1);
				for (kk = *korder + 1; kk <= (*ncoef + 1); kk++) {
					ttt[kk - 1] += (Mfloat) (kk) * eps;
				}
				eps *= F_TEN;
				nf += 1;
				if (nf < 4)
					goto L_130;

/*				imsl_ermes(4, 8, "THE KNOTS FOUND TO BE OPTIMAL ARE STACKED MORE THAN KORDER.  THIS INDICATES FEWER KNOTS WILL PRODUCE THE SAME ERROR SUM OF SQUARES.");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_OPT_KNOTS_STACKED_2);
				/*
				 * KNOTS WEREN'T SEPARATED, USE ORIGINAL
				 * KNOTS AND SSQ.
				 */
				*ssq = sqrt(*ssq);
				goto L_9000;
			}
		} else {
			mult = 1;
		}
	}
	/* STACK KNOTS AT ENDPOINTS */
	if (nf != 0)
		imsl_scopy(*ncoef + *korder, ttt, 1, xknot, 1);
	if (*korder - 1 > 0)
		sset(*korder - 1, xknot[*korder - 1], &xknot[0], 1);
	if (*korder - 1 > 0)
		sset(*korder - 1, xknot[*ncoef], &xknot[*ncoef + 1], 1);
	imsl_b4lsq(ndata, xdata, fdata, weight, korder, xknot, ncoef, bscoef,
		   &wk[0], &wk[*korder ** ncoef]);
	*ssq = imsl_b6vls(ndata, xdata, fdata, weight, korder, xknot, ncoef,
			  bscoef, ppcoef, break_, wk);
	*ssq = sqrt(*ssq);
	/*
	 * RESCALE DATA TO THE ORIGINAL INTERVAL.
	 */
L_9000:
	for (i = 1; i <= (*ncoef + *korder); i++) {
		xguess[i - 1] = xguess[i - 1] * scale1 + scale2;
		xknot[i - 1] = xknot[i - 1] * scale1 + scale2;
	}
	imsl_e1pop("IMSL_B3VLS");
	return;
}				/* end of function */

/*----------------------------------------------------------------------- */

/*  IMSL Name:  B4VLS/DB4VLS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 13, 1986

    Purpose:    Compute the variable knot B-spline least squares to
                given data.

    Usage:      CALL B4VLS (NLOOP, NDATA, XDATA, FDATA, WEIGHT, KORDER,
                            XKNOT, NCOEF, BSCOEF, PPCOEF, BREAK, WK,
                            BRNEW, COEFG)

    Arguments:
       NLOOP  - Number of of iterations.  (Input)
       NDATA  - Number of data points (must be at least 2).  (Input)
       XDATA  - Array of length NDATA containing the data point
                abscissae.  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       WEIGHT - Array of length NDATA containing the weights.  (Input)
       KORDER - Order of the spline, KORDER .LE. NDATA.  (Input)
       XKNOT  - Array of length NCOEF+KORDER containing the
                (nondecreasing) knot sequence.  (Output)
       NCOEF  - Number of spline coefficients.  (Input)
                It cannot be greater than NDATA.
       BSCOEF - Array of length NCOEF containing the spline
                representation.  (Output)
       PPCOEF - Array of length KORDER*(NCOEF+KORDER-1) containing the
                piecewise polynomial coefficients.  (Output)
       BREAK  - Array of length NCOEF-KORDER+2 containing the
                breakpoints of the PP representation.  (Output)
       WK     - Work arrary of length KORDER*(NCOEF+3).
       BRNEW  - Array of length NCOEF-KORDER+2 containing the new
                breakpoints of the PP representation.  (Output)
       COEFG  - Work array of length 2*(NCOEF-KORDER+1).

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b4vls(Mint *nloop, Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mint *korder,
	   Mfloat xknot[], Mint *ncoef, Mfloat bscoef[], Mfloat ppcoef[], Mfloat break_[], Mfloat wk[], 
           Mfloat brnew[], Mfloat coefg[])
#else
void imsl_b4vls(nloop, ndata, xdata, fdata, weight, korder,
	   xknot, ncoef, bscoef, ppcoef, break_, wk, brnew, coefg)
	Mint            *nloop, *ndata;
	Mfloat           xdata[], fdata[], weight[];
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], ppcoef[], break_[], wk[], brnew[], coefg[];
#endif
{
	Mint             _l0, i, nintv;


	for (i = 1; i <= *nloop; i++) {
		imsl_b3cpp(korder, xknot, ncoef, bscoef, &nintv, break_, ppcoef,
		  &wk[0], &wk[*korder], &wk[*korder * 2], &wk[*korder * 3]);
		break_[nintv] = xknot[*ncoef];
                _l0 = *ncoef-*korder+1;
		imsl_b5vls(break_, ppcoef, &nintv, korder, brnew, &_l0, coefg);
		imsl_scopy(*ncoef - *korder + 1, brnew, 1, &xknot[*korder - 1],
			   1);
		imsl_b4lsq(ndata, xdata, fdata, weight, korder, xknot, ncoef,
			   bscoef, &wk[0], &wk[*korder ** ncoef]);
	}
	return;
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B5VLS/DB5VLS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 13, 1986

    Purpose:    Compute the variable knot B-spline least squares to
                given data.

    Usage:      CALL B5VLS (BREAK, COEF, L, K, BRKNEW, LNEW, COEFG)

    Arguments:
       BREAK  - Work array of length L to hold a breakpoint vector
       COEF   - Work array of size K*L to hold piecewise polynomial
                coefficients.
       L      - Length of BREAK.  (Input)
       K      - Order of the PP.  (Input)
       BRKNEW - Vector of length LNEW containing the new breakpoint
                vector.  (Output)
       LNEW   - Length of BRKNEW.  (Input)
       COEFG  - Coefficient array of size 2*L which holds the piecewise
                linear coefficients of the function which BRKNEW wil
                be equidistributed.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b5vls(Mfloat break_[], Mfloat *coef, Mint *l, Mint *k, Mfloat brknew[], Mint *lnew, Mfloat *coefg)
#else
void imsl_b5vls(break_, coef, l, k, brknew, lnew, coefg)
	Mfloat           break_[], *coef;
	Mint            *l, *k;
	Mfloat           brknew[];
	Mint            *lnew;
	Mfloat          *coefg;
#endif
{
#define COEF(I_,J_)	(coef+(I_)*(*k)+(J_))
#define COEFG(I_,J_)	(coefg+(I_)*(2)+(J_))
	Mint             i, j;
	Mfloat           dif, difprv, oneovk, step, stepi;


	brknew[0] = break_[0];
	brknew[*lnew] = break_[*l];
	/* If G is constant, BRKNEW is uniform. */
	if (*l > 1) {
		/*
		 * Construct the continuous piecewise linear function G.
		 */
		oneovk = F_ONE / (Mfloat) (*k);
		*COEFG(0, 0) = F_ZERO;
		difprv = fabs(*COEF(1, *k - 1) - *COEF(0, *k - 1)) / (break_[2] -
								 break_[0]);
		for (i = 2; i <= *l; i++) {
			dif = fabs(*COEF(i - 1, *k - 1) - *COEF(i - 2, *k - 1)) /
				(break_[i] - break_[i - 2]);
			*COEFG(i - 2, 1) = pow(dif + difprv, oneovk);
			*COEFG(i - 1, 0) = *COEFG(i - 2, 0) + *COEFG(i - 2, 1) * (break_[i - 1] -
							     break_[i - 2]);
			difprv = dif;
		}
		*COEFG(*l - 1, 1) = pow(F_TWO * difprv, oneovk);
		/* STEP = G(B)/LNEW */
		step = (*COEFG(*l - 1, 0) + *COEFG(*l - 1, 1) * (break_[*l] -
					 break_[*l - 1])) / (Mfloat) (*lnew);
		/* If G is constant, BRKNEW is uniform. */
		if (step > F_ZERO) {
			/*
			 * FOR I=2,...,LNEW, construct BRKNEW(I) = A +
			 * G**(-1)(STEPI), with STEPI = (I-1)*STEP . This
			 * requires inversion of the P.linear function G .
			 * For this, J is found so that G(BREAK(J)) .LE.
			 * STEPI .LE. G(BREAK(J+1)) and then BRKNEW(I) =
			 * BREAK(J) + (STEPI-G(BREAK(J)))/DG(BREAK(J)) . The
			 * midpoint is chosen if DG(BREAK(J)) = 0 .
			 */
			j = 1;
			for (i = 2; i <= *lnew; i++) {
				stepi = (Mfloat) (i - 1) * step;
		L_20:
				if (j != *l) {
					if (stepi > *COEFG(j, 0)) {
						j += 1;
						goto L_20;
					}
				}
				if (*COEFG(j - 1, 1) != F_ZERO) {
					brknew[i - 1] = break_[j - 1] + (stepi - *COEFG(j - 1, 0)) /
						*COEFG(j - 1, 1);
				} else {
					brknew[i - 1] = (break_[j - 1] + break_[j]) / F_TWO;
				}
			}
			goto L_9000;
			/* If G is constant, BRKNEW is uniform. */
		}
	}
	step = (break_[*l] - break_[0]) / (Mfloat) (*lnew);
	for (i = 2; i <= *lnew; i++) {
		brknew[i - 1] = break_[0] + (Mfloat) (i - 1) * step;
	}
L_9000:
	return;
}				/* end of function */

/*----------------------------------------------------------------------- */

/*  IMSL Name:  B6VLS/DB6VLS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 13, 1986

    Purpose:    Compute the variable knot B-spline least squares to
                given data.

    Usage:      B6VLS(NDATA, XDATA, FDATA, WEIGHT, KORDER, XKNOT, NCOEF,
                      BSCOEF, PPCOEF, BREAK, WK)

    Arguments:
       NDATA  - Number of data points (must be at least 2).  (Input)
       XDATA  - Array of length NDATA containing the data point
                abscissae.  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       WEIGHT - Array of length NDATA containing the weights.  (Input)
       KORDER - Order of the spline, KORDER .LE. NDATA.  (Input)
       XKNOT  - Array of length NCOEF+KORDER containing the
                (nondecreasing) knot sequence.  (Output)
       NCOEF  - Number of spline coefficients.  (Input)
                It cannot be greater than NDATA.
       BSCOEF - Array of length NCOEF containing the spline
                representation.  (Output)
       PPCOEF - Array of length KORDER*(NCOEF+KORDER-1) containing the
                piecewise polynomial coefficients.  (Output)
       BREAK  - Array of length NCOEF-KORDER+2 containing the
                breakpoints of the PP representation.  (Output)
       WK     - Work array of length KORDER*(NCOEF+3).
       B6VLS  - The sum of the squares of the error.  (Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_b6vls(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mint *korder, Mfloat xknot[],
	   Mint *ncoef, Mfloat bscoef[], Mfloat ppcoef[], Mfloat break_[], Mfloat wk[])
#else
Mfloat imsl_b6vls(ndata, xdata, fdata, weight, korder, xknot,
	   ncoef, bscoef, ppcoef, break_, wk)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], weight[];
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], ppcoef[], break_[], wk[];
#endif
{
	Mint             _l0, j, nintv;
	Mfloat           b6vls_v, s;
	double          a, b, c, d;

        _l0 = 0;
	a = F_ZERO;
	if (*ndata / *korder <= *ncoef - *korder) {
		for (j = 1; j <= *ndata; j++) {
			b = weight[j - 1];
			c = fdata[j - 1];
			d = imsl_b3der(&_l0, &xdata[j - 1], korder, xknot, ncoef,
			    bscoef, &wk[0], &wk[*korder], &wk[*korder * 2]);
			a += b * imsl_fi_power(c - d, 2);
		}
	} else {

		imsl_b3cpp(korder, xknot, ncoef, bscoef, &nintv, break_, ppcoef,
		  &wk[0], &wk[*korder], &wk[*korder * 2], &wk[*korder * 3]);
		for (j = 1; j <= *ndata; j++) {
			b = weight[j - 1];
			c = fdata[j - 1];
			d = imsl_ppder(0, xdata[j - 1], *korder, nintv, break_, ppcoef);
			a += b * imsl_fi_power(c - d, 2);
		}
	}
	s = (Mfloat) (a);
	b6vls_v = s;
	return (b6vls_v);
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B7VLS/DB7VLS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 13, 1986

    Purpose:    Find the minimum of a smooth univariate function within
                an interval.  (Reverse communication)

    Usage:      CALL B7VLS (IDO, F, A, B, XGUESS, STEP, XACC, MAXFN, X)

    Arguments:
       IDO    - Parameter for reverse communication.  (Input/output)
                User has to set it to 0 before calling this subroutine,
                and evaluate the function when it is 1.
       F      - Function value evaluated at X.  (Input)
       A      - Lower bound of the interval in which the minimum has
                to be found.  (Input)
       B      - Upper bound of the interval in which the minimum has
                to be found.  (Input)
                B has to be greater than A.
       XGUESS - An initial guess of the minimum.  (Input)
       STEP   - An order of magnitude estimate of the required
                change in X.  (Input)
       XACC   - The required absolute accuracy in the final
                value of X.  (Input)
                On a normal return there are points on either side of X
                within a distance XACC at which F is no less than F at X.
       MAXFN  - A limit, which must be set to a positive
                integer, on the number of function evaluations.  (Input)
       X      - The scalar variable of the calculation.  (Output)
                It gives the least calculated value of F.

    Remark:
       Informational errors
       Type Code
         3   15 The final value of X is at a bound.  The minimum
                is probably beyond the bound.
         3   14 Computer rounding errors prevent further
                refinement of X.
         4   11 The number of function evaluations has exceeded MAXFN.

    Keywords:   Unconstrained univariate minimization; Safeguarded
                quadratic interpolation

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b7vls(Mint *ido, Mfloat *f, Mfloat *a, Mfloat *b, Mfloat *xguess, Mfloat *step, Mfloat *xacc, Mint *maxfn, Mfloat *x)
#else
void imsl_b7vls(ido, f, a, b, xguess, step, xacc, maxfn, x)
	Mint            *ido;
	Mfloat          *f, *a, *b, *xguess, *step, *xacc;
	Mint            *maxfn;
	Mfloat          *x;
#endif
{
	static Mint      info, is, jump, nf;
	static Mfloat    acc, bl, bound, bu, da, db, dc, fa, fb, fc, fd,
	                fx, h, rho, sigma, st, temp, tol, x0, xa, xb, xc,
	                xd, xx;


	imsl_e1psh("IMSL_B7VLS ");

	jump = *ido + 1;
	if (jump == 2)
		goto L_100;

	*ido = 1;
	*x = *xguess;
	nf = 0;
	goto L_90;
L_20:
	x0 = *x;
	xb = *x;
	fb = fx;
	info = 0;
	/*
	 * ALTER ANY UNSUITABLE INITIAL PARAMETERS, INCLUDING INCREASING THE
	 * STEP IF NECESSARY IN ORDER THAT THE COMPUTER PRECISION GIVES A
	 * CHANGE IN X.
	 */
	if (*x < *a || *x > *b) {
		info = 6;
		goto L_310;
	}
	acc = imsl_f_max(F_ZERO, *xacc);
	bl = *x - *a;
	bu = *b - *x;
	bound = imsl_f_max(bl, bu);
	if (bound <= acc) {
		info = 2;
		goto L_310;
	}
	st = *step;
L_30:
	if (st < F_ZERO) {
		if (fabs(st) > bl)
			st = -bl;
	} else {
		if (fabs(st) > bu)
			st = bu;
	}
	if (st == F_ZERO) {
		st = 0.01 * bound;
		if (bl > bu)
			st = -st;
	}
	if (fabs(st) < acc)
		st = sign(acc, st);
	goto L_50;
L_40:
	st *= F_FIVE;
L_50:
	xa = xb + F_HALF * fabs(st);
	xc = xb + fabs(st);
	if (xa <= xb)
		goto L_40;
	if (xc <= xa)
		goto L_40;

	if (*x + st < *a || *x + st > *b) {
		if (info == 2) {
			goto L_310;
		} else {
			info = 2;
			st = -st;
			goto L_30;
		}
	}
	/*
	 * CALCULATE THE NEXT TRIAL VALUE OF X FROM XB AND ST.
	 */
L_60:
	is = 1;
L_70:
	*x = xb + st;
	h = xb + 1.5 * st;
	if (h <= *a) {
		*x = *a;
		is = 2;
	} else if (h >= *b) {
		*x = *b;
		is = 2;
	}
	/* CALCULATE THE NEXT VALUE OF F. */
L_80:
	if (nf >= *maxfn) {
		info = 4;
		goto L_310;
	}
L_90:
	nf += 1;
	goto L_330;
L_100:
	fx = *f;
	if (nf <= 1)
		goto L_20;
	/*
	 * REVERSE ST IF THE INITIAL STEP SEEMS TO BE UPHILL.
	 */
	if (nf >= 3)
		goto L_110;
	if (fx < fb)
		goto L_120;
	xc = *x;
	fc = fx;
	if (xb == *a || xb == *b)
		goto L_130;
	st = -st;
	goto L_60;
	/*
	 * ENSURE THAT FB IS THE LEAST CALCULATED VALUE OF F.
	 */
L_110:
	xd = xc;
	fd = fc;
	if (fx >= fb)
		goto L_160;
L_120:
	xc = xb;
	fc = fb;
	xb = *x;
	fb = fx;
	/*
	 * CALCULATE AN EXTRA FUNCTION VALUE IF X IS AT A BOUND.
	 */
	if (is >= 4)
		goto L_180;
	if (is <= 1)
		goto L_150;
	if (is == 3)
		goto L_140;
L_130:
	is = 3;
	h = imsl_f_max(0.9e0 * acc, 0.01e0 * fabs(xb - xc));
	*x = xb + sign(h, xc - xb);
	if (fabs(*x - xc) < fabs(*x - xb))
		*x = F_HALF * (xb + xc);
	temp = (xb - *x) / (xb - xc);
	if (temp <= F_ZERO) {
		info = 3;
		goto L_310;
	}
	goto L_80;
	/*
	 * THIS STAGE IS REACHED WHEN A BRACKET IS FOUND NEAR A BOUND.
	 */
L_140:
	xa = xd;
	fa = fd;
	is = 4;
	goto L_210;
	/*
	 * CALCULATE THE NEXT STEP IN THE SEARCH FOR A BRACKET, TRYING TO
	 * OVERSHOOT THE PREDICTED MINIMUM BY THE FACTOR THREE.
	 */
L_150:
	if (nf <= 2)
		goto L_70;
	rho = (xb - xc) / (xb - xd);
	sigma = (fb - fc) / (fb - fd);
	h = F_NINE;
	if (sigma < rho)
		h = 1.5e0 * (rho - sigma / rho) / (sigma - rho);
	h = imsl_f_min(h, F_NINE);
	h = imsl_f_max(h, F_TWO);
	st *= h;
	goto L_70;
	/*
	 * RETURN IF THE MINIMUM SEEMS TO BE BEYOND A BOUND.
	 */
L_160:
	if (is >= 4)
		goto L_170;
	if (is == 3) {
		info = 3;
		goto L_310;
	}
	/*
	 * A BRACKET HAS BEEN FOUND SO BRANCH TO PREDICT THE MINIMUM.
	 */
	xa = *x;
	fa = fx;
	is = 4;
	goto L_210;
	/*
	 * ENSURE THAT XA, XB, XC AND XD ARE ORDERED MONOTONICALLY.
	 */
L_170:
	xc = *x;
	fc = fx;
L_180:
	temp = (xb - xc) / (xa - xd);
	if (temp > F_ZERO)
		goto L_200;
L_190:
	h = xa;
	xa = xd;
	xd = h;
	h = fa;
	fa = fd;
	fd = h;
	/*
	 * IF THERE ARE THREE CONSECUTIVE EQUAL VALUES OF F, ENSURE THAT FB
	 * IS THE MIDDLE ONE.
	 */
L_200:
	if (fa == fb)
		goto L_210;
	if (fd != fb)
		goto L_210;
	if (fc != fb)
		goto L_210;
	xc = xb;
	xb = *x;
	fc = fb;
	fb = fx;
	goto L_190;
	/*
	 * USE THE MINIMA OF TWO QUADRATICS TO CALCULATE A TOLERANCE ON THE
	 * NEXT
	 */
L_210:
	da = (fb - fa) / (xb - xa);
	db = (fc - fb) / (xc - xb);
	if (is >= 5)
		goto L_220;
	is = 5;
	tol = 0.01e0 * fabs(xa - xc);
	goto L_230;
L_220:
	dc = (fd - fc) / (xd - xc);
	tol = F_ZERO;
	if (db == F_ZERO)
		goto L_230;
	tol = fabs(xa - xc);
	temp = (dc - db) / (xa - xc);
	if (temp >= F_ZERO)
		goto L_230;
	h = F_HALF * fabs(db * ((xd - xb) / (dc - db) + (xa - xc) / (db - da)));
	if (h < tol)
		tol = h;
L_230:
	tol = imsl_f_max(tol, 0.9e0 * acc);
	/*
	 * SET X TO THE VALUE THAT MINIMIZES THE INTERPOLATING QUADRATIC.
	 */
	*x = xb;
	if (da != db)
		*x = F_HALF * (da * (xb + xc) - db * (xa + xb)) / (da - db);
	/* ENSURE THAT ABS(XA-XB).GE.ABS(XB-XC). */
	temp = (xa - xb) / (xb - xc);
	if (temp >= F_ONE)
		goto L_240;
	h = xa;
	xa = xc;
	xc = h;
	h = fa;
	fa = fc;
	fc = h;
	/* TEST FOR CONVERGENCE. */
L_240:
	if (fabs(xa - xb) <= acc)
		goto L_320;
	/*
	 * IF ABS(XA-XB).LE.2.9*ACC, CHOOSE THE NEXT X TO AVOID AN
	 * UNNECESSARY FUNCTION EVALUATION.
	 */
	temp = (xa - xb) / (xb - xc);
	if (temp > F_TEN)
		goto L_260;
	if (fabs(xa - xb) > 2.9 * acc)
		goto L_250;
	*x = F_HALF * (xa + xb);
	if (fabs(*x - xb) <= acc)
		goto L_270;
	*x = 0.67 * xb + 0.33 * xa;
	goto L_270;
	/*
	 * IF (XA-XB)/(XB-XC).LE.10, ENSURE THAT THE DISTANCE FROM X TO XB IS
	 * NORMALLY AT LEAST TOL.
	 */
L_250:
	if (fabs(*x - xb) >= tol)
		goto L_270;
	*x = xb + sign(tol, *x - xb);
	if (fabs(*x - xc) < fabs(*x - xb))
		*x = xb + sign(tol, xa - xb);
	if (fabs(*x - xa) < fabs(*x - xb))
		*x = F_HALF * (xa + xb);
	goto L_270;
	/*
	 * WHEN (XA-XB)/(XB-XC).GT.10, ENSURE THAT X IS IN THE LONGER
	 * INTERVAL, AND TRY TO OVERSHOOT THE MINIMUM.
	 */
L_260:
	h = (*x - xb) * sign(F_THREE, xa - xb);
	h = imsl_f_vmax(3,h, tol, fabs(xb - xc));
	h = imsl_f_min(h, 0.1 * fabs(xa - xb));
	*x = xb + sign(h, xa - xb);
	/*
	 * TEST WHETHER ROUNDING ERRORS MAKE THE NEW X THE SAME AS A PREVIOUS
	 * X.
	 */
L_270:
	if (*x == xb)
		goto L_280;
	if ((xa - *x) / (xa - xb) <= F_ZERO)
		goto L_280;
	if ((*x - xc) / (xa - xb) > F_ZERO)
		goto L_80;
	/*
	 * SEEK A NEW VALUE OF X BETWEEN XA AND XB, BUT ROUNDING ERRORS MAY
	 * NOT ALLOW ONE.
	 */
L_280:
	*x = xa;
L_290:
	xx = F_HALF * (*x + xb);
	temp = (xx - xb) / (xa - xb);
	if (temp <= F_ZERO)
		goto L_300;
	temp = (*x - xx) / (xa - xb);
	if (temp <= F_ZERO)
		goto L_300;
	*x = xx;
	goto L_290;
L_300:
	if (*x != xa)
		goto L_80;
	info = 5;
	/* RETURN FROM THE SUBROUTINE. */
L_310:
	;
	info = 0;
	if (info == 2) {

/*		imsl_ermes(3, 13, "The interval between bounds is too small.  No change in X can be made.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_INTERVAL_TOO_SMALL);
	} else if (info == 6) {

/*		imsl_ermes(5, 19, "The initial guess is out of bound.  Rerun with a new guess.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INIT_GUESS_OUT_BOUNDS);
	} else if (info == 5) {

/*		imsl_ermes(3, 14, "Computer rounding errors prevent further refinement of X.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_ROUNDING_ERRORS_IN_X);
	} else if (info == 3) {

/*		imsl_ermes(3, 15, "The final value for X is at a bound.  The minimum is probably beyond the bound.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_FINAL_VALUE_AT_BOUND);
	} else if (info == 4) {
		imsl_e1sti(1, *maxfn);

/*		imsl_ermes(3, 11, "The maximum number of function evaluations, MAXFN = %(i1), has been exceeded.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_MAX_FCN_EVAL_EXCEEDED);
	}
L_320:
	*x = xb;
	*f = fb;
	*ido = 2;

L_330:
	imsl_e1pop("IMSL_B7VLS ");
	return;
}				/* end of function */

/*  -----------------------------------------------------------------------
    IMSL Name:  B8VLS/DB8VLS (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 13, 1986

    Purpose:    Compute the variable knot B-spline least squares to
                given data.

    Usage:      CALL B8VLS (TT, K, NCOEF, TTT)

    Arguments:
       TT     - Array of length 3*(NCOEF+K) containing the last three
                knot vectors.  (Input)
       K      - Order of the spline.  (Input)
       NCOEF  - Number of spline coefficients.  (Input)
       TTT    - Array of length NCOEF+KORDER containing the new
                (nondecreasing) knot sequence.  (Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b8vls(Mfloat *tt, Mint *k, Mint *ncoef, Mfloat ttt[])
#else
void imsl_b8vls(tt, k, ncoef, ttt)
	Mfloat          *tt;
	Mint            *k, *ncoef;
	Mfloat           ttt[];
#endif
{
#define TT(I_,J_)	(tt+(I_)*(3)+(J_))
	Mint             i;
	Mfloat           right, sleft;
	double          a, b, c, d, dif1, dif2, t1, t2;


	t1 = F_THREE;
	t2 = 12.0e0;
	sleft = *TT(*k - 1, 0);
	right = *TT(*ncoef, 0);
	/*
	 * We now apply quadratic interpolation and evaluate the quadratic at
	 * zero to get new knot value.  If new value is not in the valid
	 * interval, then the most current value is used.  The new knots are
	 * in TTT.
	 */
	for (i = *k + 1; i <= *ncoef; i++) {
		a = *TT(i - 1, 2);
		b = *TT(i - 1, 1);
		c = *TT(i - 1, 0);
		dif1 = (b - a) / (t1 - F_ONE);
		dif2 = (c - b) / (t2 - t1);
		d = a - dif1 + (t1 / (t2 - t1)) * (dif2 - dif1);
		if (sleft < d && d < right) {
			ttt[i - *k - 1] = d;
		} else {
			ttt[i - *k - 1] = a;
		}
	}
	/* We now sort the vector TTT */
	imsl_svrgn(*ncoef - *k, ttt, ttt);
	imsl_scopy(*ncoef - *k, ttt, -1, &ttt[*k], -1);
	sset(*k, sleft, ttt, 1);
	sset(*k, right, &ttt[*ncoef], 1);
	return;
}				/* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  ISMAX (Single precision version)
 
    Computer:   FORC/SINGLE
 
    Revised:    August 9, 1986
 
    Purpose:    Find the smallest index of the component of a
                single-precision vector having maximum value.
 
    Usage:      ISMAX(N, SX, INCX)
 
    Arguments:
       N      - Length of vector X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than zero.
       ISMAX  - The smallest index I such that X(I)  is the maximum of
                X(J) for J=1 to N.  (Output)
                X(I) refers to a specific element of SX. See INCX
                argument description.
 
    Keyword:    Level 1 BLAS
 
    GAMS:       D1a2
 
    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support
 
    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
 
    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.
 
  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mint l_ismax(Mint n, Mfloat sx[], Mint incx)
#else
static Mint l_ismax(n, sx, incx)
        Mint           n;
        Mfloat           sx[];
        Mint           incx;
#endif
{
        Mint            i, ismax_v, ix;
        Mfloat           smax;


        ismax_v = 0;
        if (n >= 1) {
                ismax_v = 1;
                if (n != 1) {
                        if (incx != 1) {
                                /* CODE FOR INCREMENT NOT EQUAL TO 1 */
                                ix = 1;
                                smax = sx[0];
                                ix += incx;
                                for (i = 2; i <= n; i++) {
                                        if (sx[ix - 1] > smax) {
                                                ismax_v = i;
                                                smax = sx[ix - 1];
                                        }
                                        ix += incx;
                                }
                        } else {
                                /* CODE FOR INCREMENT EQUAL TO 1 */
                                smax = sx[0];
                                for (i = 2; i <= n; i++) {
                                        if (sx[i - 1] > smax) {
                                                ismax_v = i;
                                                smax = sx[i - 1];
                                        }
                                }
                        }
                }
        }
        return (ismax_v);
}                               /* end of function */
