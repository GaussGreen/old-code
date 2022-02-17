#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B2DER/DB2DER (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Evaluate the I-th derivative of a spline given its
                B-spline representation.

    Usage:      B2DER(IDERIV, X, KORDER, XKNOT, NCOEF, BSCOEF, AJ, DL,
                      DR)

    Arguments:
       IDERIV - Order of the derivative to be evaluated.  (Input)
                In particular, IDERIV = 0 returns the value of the
                spline.
       X      - Point at which the spline is to be evaluated.  (Input)
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NCOEF+KORDER containing the knot
                sequence.  (Input)
                XKNOT must be nondecreasing.
       NCOEF  - Number of B-spline coefficients.  (Input)
       BSCOEF - Array of length NCOEF containing the B-spline
                coefficients.  (Input)
       AJ     - Array of length KORDER.
       DL     - Array of length KORDER.
       DR     - Array of length KORDER.
       B2DER  - Value of the IDERIV-th derivative of the spline at X.
                (Output)

    Remark:
       Informational errors
       Type Code
         4   4  Multiplicity of the knots cannot exceed the order
                of the spline.
         4   5  The knots must be nondecreasing.

    Keyword:    Differentiate

    GAMS:       E3; K6

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_b2der(Mint *ideriv, Mfloat *x, Mint *korder, Mfloat xknot[],
                 Mint *ncoef, Mfloat bscoef[],
	         Mfloat aj[], Mfloat dl[], Mfloat dr[])
#else
Mfloat imsl_b2der(ideriv, x, korder, xknot, ncoef, bscoef,
	   aj, dl, dr)
	Mint            *ideriv;
	Mfloat          *x;
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], aj[], dl[], dr[];
#endif
{
	Mfloat           b2der_v;


	imsl_e1psh("IMSL_B2DER");
	b2der_v = F_ZERO;
	/* CHECK KORDER */
	if (*korder < 1) {
		imsl_e1sti(1, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
		goto L_9000;
	}
	/* CHECK IDERIV */
	if (*ideriv < 0) {
		imsl_e1sti(1, *ideriv);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_DERIV);
		goto L_9000;
	}
	/* CHECK NCOEF */
	if (*ncoef < *korder) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS);
		goto L_9000;
	}
	/* CHECK REST OF THE ARGUMENTS */
	imsl_b3int(korder, xknot, ncoef);
	/* CHECK FOR ERRORS */
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	/*
	 * COMPUTE VALUE WITHOUT FURTHER CHECKING
	 */
	b2der_v = imsl_b3der(ideriv, x, korder, xknot, ncoef, bscoef, aj,
			     dl, dr);

L_9000:
	;
	imsl_e1pop("IMSL_B2DER");
	return (b2der_v);
}				/* end of function */


/*  -----------------------------------------------------------------------
    IMSL Name:  B3DER/DB3DER (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Calculate the value of the I-th derivative of a spline
                given its B-spline representation.

    Usage:      B3DER(IDERIV, X, KORDER, XKNOT, NCOEF, BSCOEF, AJ,
                      DL, DR)

    Arguments:  (See BSDER)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_b3der(Mint *ideriv, Mfloat *x, Mint *korder, Mfloat xknot[],
                 Mint *ncoef, Mfloat bscoef[],
           Mfloat aj[], Mfloat dl[], Mfloat dr[])
#else
Mfloat imsl_b3der(ideriv, x, korder, xknot, ncoef, bscoef,
	   aj, dl, dr)
	Mint            *ideriv;
	Mfloat          *x;
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], aj[], dl[], dr[];
#endif
{
	Mint             i, ilo, j, jcmax, jcmin, jj, mflag, nknot;
	Mfloat           b3der_v, value;


	value = F_ZERO;
	if (*ideriv >= *korder)
		goto L_9000;
	/*
	 * Find I such that 1 .LE. I .LT. NKNOT and XKNOT(I) .LT. XKNOT(I+1)
	 * and XKNOT(I) .LE. X .LT. XKNOT(I+1) . If no such I can be found, X
	 * lies outside the support of the spline F and VALUE = 0. (The
	 * asymmetry in this choice of I makes F right continuous)
	 */
	nknot = *ncoef + *korder;
	imsl_b4der(xknot, &nknot, x, &i, &mflag);
	if (mflag != 0)
		goto L_9000;
	if (*korder <= 1) {
		value = bscoef[i - 1];
		goto L_9000;
	}
	/*
	 * Store the KORDER B-spline coefficients relevant for the XKNOT
	 * interval (XKNOT(I),XKNOT(I+1)) in AJ(1),...,AJ(KORDER) and compute
	 * DL, DR. Set any of the AJ not obtainable. From input to zero. Set
	 * any XKNOTS not obtainable equal to XKNOT(1) to XKNOT(NKNOT)
	 * appropriately.
	 */
	jcmin = 1;
	if (i < *korder) {
		jcmin = *korder - i + 1;
		for (j = 1; j <= i; j++) {
			dl[j - 1] = *x - xknot[i - j];
		}
		sset(*korder - i, F_ZERO, aj, 1);
		sset(*korder - i, dl[i - 1], &dl[i - 1], 1);
	} else {
		for (j = 1; j <= (*korder - 1); j++) {
			dl[j - 1] = *x - xknot[i - j];
		}
	}

	jcmax = *korder;
	if (*ncoef < i) {
		jcmax = nknot - i;
		for (j = 1; j <= jcmax; j++) {
			dr[j - 1] = xknot[i + j - 1] - *x;
		}
		sset(*korder - jcmax, F_ZERO, &aj[jcmax], 1);
		sset(*korder - jcmax, dr[jcmax - 1], &dr[jcmax - 1], 1);
	} else {
		for (j = 1; j <= (*korder - 1); j++) {
			dr[j - 1] = xknot[i + j - 1] - *x;
		}
	}
	scopy(jcmax - jcmin + 1, &bscoef[i - *korder + jcmin - 1], 1,
		   &aj[jcmin - 1], 1);
	/*
	 * Difference the coefficients IDERIV times.
	 */
	for (j = 1; j <= *ideriv; j++) {
		ilo = *korder - j;
		for (jj = 1; jj <= (*korder - j); jj++) {
			aj[jj - 1] = ((aj[jj] - aj[jj - 1]) / (dl[ilo - 1] + dr[jj - 1])) *
				(float) (*korder - j);
			ilo -= 1;
		}
	}
	/*
	 * Compute value at X in (XKNOT(I),XKNOT(I+1)) of IDERIV-th
	 * derivative, given its relevant B-spline coeffs in
	 * AJ(1),...,AJ(KORDER-IDERIV).
	 */
	for (j = *ideriv + 1; j <= (*korder - 1); j++) {
		ilo = *korder - j;
		for (jj = 1; jj <= (*korder - j); jj++) {
			aj[jj - 1] = (aj[jj] * dl[ilo - 1] + aj[jj - 1] * dr[jj - 1]) /
				(dl[ilo - 1] + dr[jj - 1]);
			ilo -= 1;
		}
	}
	value = aj[0];

L_9000:
	;
	b3der_v = value;
	return (b3der_v);
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B4DER/DB4DER (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Compute LEFT = MAX (I, 1 .LE. NKNOT. .AND. XKNOT(I) .LE.
                X

    Usage:      CALL B4DER (XKNOT, NKNOT, X, LEFT, MFLAG)

    Arguments:
       XKNOT  - The knot sequence.  (Input)
       NKNOT  - Number of knots.  (Input)
       X      - The point whose location in XKNOT is to be found.
                (Input)
       LEFT   - Integer whose value is given below.  (Output)
       MFLAG  - Flag defined below.  (Output)
                 LEFT   MFLAG
                   1     -1      IF                    X .LT.  XKNOT(1)
                   I      0      IF   XKNOT(I)    .LE. X .LE. XKNOT(I+1)
                NKNOT     1      IF  XKNOT(NKNOT) .LT. X
                In particular, MFLAG = 0 is the usual case.  MFLAG .NE. 0
                indicates that X lies outside the half open interval
                XKNOT(1) .LE. X .LT. XKNOT(NKNOT). The asymmetric
                treatment of the interval is due to the decision to make
                all PP functions continuous from the right.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
void imsl_b4der(Mfloat xknot[], Mint *nknot, Mfloat *x, Mint *left, Mint *mflag)
#else
void imsl_b4der(xknot, nknot, x, left, mflag)
	Mfloat           xknot[];
	Mint            *nknot;
	Mfloat          *x;
	Mint            *left, *mflag;
#endif
{
	Mint             ihi, istep, middle;
	static Mint      ilo = 1;



	ihi = ilo + 1;
	if (ihi >= *nknot) {
		if (*x >= xknot[*nknot - 1]) {
			*left = *nknot;
			if (*x == xknot[*nknot - 1]) {
		L_10:
				*left -= 1;
				if (*x == xknot[*left - 1])
					goto L_10;
				*mflag = 0;
			} else {
				*mflag = 1;
			}
			goto L_9000;
		} else if (*nknot <= 1) {
			*mflag = -1;
			*left = 1;
			goto L_9000;
		}
		ilo = *nknot - 1;
		ihi = *nknot;
	}
	if (*x < xknot[ihi - 1]) {
		if (*x >= xknot[ilo - 1]) {
			*mflag = 0;
			*left = ilo;
			goto L_9000;
		}
		/*
		 * Now X .LT. XKNOT(ILO) . Decrease ILO to capture X .
		 */
		istep = 1;
L_20:
		;
		ihi = ilo;
		ilo = ihi - istep;
		if (ilo > 1) {
			if (*x >= xknot[ilo - 1])
				goto L_50;
			istep *= 2;
			goto L_20;
		}
		ilo = 1;
		if (*x < xknot[0]) {
			*mflag = -1;
			*left = 1;
			goto L_9000;
		}
		goto L_50;
	}
	/*
	 * Now X .GE. XKNOT(IHI) . Increase IHI to capture X .
	 */
	istep = 1;
L_30:
	;
	ilo = ihi;
	ihi = ilo + istep;
	if (ihi < *nknot) {
		if (*x < xknot[ihi - 1])
			goto L_50;
		istep *= 2;
		goto L_30;
	}
	if (*x >= xknot[*nknot - 1]) {
		*left = *nknot;
		if (*x == xknot[*nknot - 1]) {
	L_40:
			*left -= 1;
			if (*x == xknot[*left - 1])
				goto L_40;
			*mflag = 0;
		} else {
			*mflag = 1;
		}
		goto L_9000;
	}
	ihi = *nknot;
	/*
	 * Now XKNOT(ILO) .LE. X .LE. XKNOT(IHI) . Narrow the inteval.
	 */
L_50:
	;
	middle = (ilo + ihi) / 2;
	if (middle == ilo) {
		*mflag = 0;
		*left = ilo;
		goto L_9000;
	}
	/*
	 * It is assumed that MIDDLE = ILO in case IHI = ILO+1 .
	 */
	if (*x < xknot[middle - 1]) {
		ihi = middle;
	} else {
		ilo = middle;
	}
	goto L_50;
L_9000:
	;
	return;
}				/* end of function */
