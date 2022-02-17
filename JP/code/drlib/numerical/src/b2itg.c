#include "imsl_inc.h"

#if defined( WIN32 )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*  -----------------------------------------------------------------------
    IMSL Name:  B2ITG/DB2ITG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Evaluate the integral of a spline, given its B-spline
                representation.

    Usage:      B2ITG(A, B, KORDER, XKNOT, NCOEF, BSCOEF, TCOEF,
                      AJ, DL, DR)

    Arguments:
       A      - Lower limit of integration.  (Input)
       B      - Upper limit of integration.  (Input)
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length KORDER+NCOEF containing the knot
                sequence.  (Input)
                XKNOT must be nondecreasing.
       NCOEF  - Number of B-spline coefficients.  (Input)
       BSCOEF - Array of length NCOEF containing the B-spline
                coefficients.  (Input)
       TCOEF  - Work array of length KORDER+1.
       AJ     - Work array of length KORDER+1.
       DL     - Work array of length KORDER+1.
       DR     - Work array of length KORDER+1.
       B2ITG  - Value of the integral of the spline from A to B.
                (Output)

    Remark:
       Informational errors
       Type Code
         3   7  The upper and lower endpoints of integration are equal.
         3   8  The lower limit of integration is less than
                XKNOT(KORDER).
         3   9  The upper limit of integration is greater than
                XKNOT(NCOEF+1).
         4   4  Multiplicity of the knots cannot exceed the order
                of the spline.
         4   5  The knots must be nondecreasing.

    Keyword:    Quadrature

    GAMS:       H2

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_b2itg(Mfloat *a, Mfloat *b, Mint *korder, Mfloat xknot[], Mint *ncoef,
	    Mfloat bscoef[], Mfloat tcoef[], Mfloat aj[], Mfloat dl[], Mfloat dr[])
#else
Mfloat imsl_b2itg(a, b, korder, xknot, ncoef, bscoef, tcoef,
	   aj, dl, dr)
	Mfloat          *a, *b;
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], tcoef[], aj[], dl[], dr[];
#endif
{
	Mfloat           b2itg_v;


	imsl_e1psh("IMSL_B2ITG ");
	/* Check KORDER */
	b2itg_v = F_ZERO;

	if (*korder < 1) {
		imsl_e1sti(1, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
	}
	/* Check NCOEF */
	if (*ncoef < *korder) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS);
	}
	/* Check XKNOT */
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	imsl_b3int(korder, xknot, ncoef);
	if (imsl_n1rty(0) == 0) {
		b2itg_v = imsl_b3itg(a, b, korder, xknot, ncoef, bscoef, tcoef,
				     aj, dl, dr);
	}
L_9000:
	imsl_e1pop("IMSL_B2ITG ");
	return (b2itg_v);
}				/* end of function */

/*  -----------------------------------------------------------------------
    IMSL Name:  B3ITG/DB3ITG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Evaluate the integral of a spline, given its B-spline
                representation.

    Usage:      B3ITG(A, B, KORDER, XKNOT, NCOEF, BSCOEF, TCOEF,
                      AJ, DL, DR)

    Arguments:  (See BSITG)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_b3itg(Mfloat *a, Mfloat *b, Mint *korder, Mfloat xknot[], Mint *ncoef,
            Mfloat bscoef[], Mfloat tcoef[], Mfloat aj[], Mfloat dl[], Mfloat dr[])
#else
Mfloat imsl_b3itg(a, b, korder, xknot, ncoef, bscoef, tcoef,
	   aj, dl, dr)
	Mfloat          *a, *b;
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], tcoef[], aj[], dl[], dr[];
#endif
{
	Mint             _l0, ii, isign, kount, lefta, leftb,
	                mflag;
	Mfloat          b3itg_v, ra, rb, sum, vala,
	                valb;


	imsl_e1psh("IMSL_B3ITG ");

	b3itg_v = F_ZERO;

	/*
	 * Assign RA to MIN(A,B) Assign RB to MAX(A,B)
	 */
	if (*a < *b) {
		ra = *a;
		rb = *b;
		isign = 1;
	} else if (*a > *b) {
		ra = *b;
		rb = *a;
		isign = -1;
	} else {

		imsl_ermes(IMSL_WARNING, IMSL_SPLINE_EQUAL_LIMITS);
		goto L_9000;
	}

	lefta = 0;
	/*
	 * If RA is to the left of XKNOT(KORDER) then reset RA to
	 * XKNOT(KORDER) and there is no need to find LEFTA since it is
	 * KORDER.
	 */
	if (ra == xknot[*korder - 1])
		lefta = *korder;
	if (ra < xknot[*korder - 1]) {
		lefta = *korder;
		ra = xknot[*korder - 1];
		if (isign == 1) {

			imsl_ermes(IMSL_WARNING,
			IMSL_LIMITS_LOWER_TOO_SMALL);
		} else {

			imsl_ermes(IMSL_WARNING,
			IMSL_LIMITS_UPPER_TOO_SMALL);
		}
	}
	/*
	 * If RB is to the right of XKNOT(NCOEF+1) then reset RB to
	 * XKNOT(NCOEF+1) and there is no need to find LEFTB since it is
	 * NCOEF.
	 */
	leftb = 0;
	if (rb == xknot[*ncoef])
		leftb = *ncoef;
	if (rb > xknot[*ncoef]) {
		leftb = *ncoef;
		rb = xknot[*ncoef];
		if (isign == 1) {

			imsl_ermes(IMSL_WARNING, IMSL_LIMITS_UPPER_TOO_BIG);
		} else {

			imsl_ermes(IMSL_WARNING, IMSL_LIMITS_LOWER_TOO_BIG);
		}
	}
	/*
	 * If RA=RB or both RA and RB are outside the range of the B-sline
	 * then the integral is zero.
	 */
	if ((ra <= xknot[*korder - 1] && rb <= xknot[*korder - 1]) ||
	    (ra >= xknot[*ncoef] && rb >= xknot[*ncoef]))
		goto L_9000;
	if (ra == rb)
		goto L_9000;
	/*
	 * Now when we find LEFTA it should be in the closed interval
	 * (KORDER,NCOEF)
	 */
	if (lefta == 0) {
                _l0 =  *ncoef + *korder;
		imsl_b5itg(xknot, &_l0, &ra, &lefta, &mflag);
	}
	kount = 0;
	tcoef[kount] = F_ZERO;
	for (ii = 1; ii <= (lefta - *korder); ii++) {
		tcoef[kount] += bscoef[ii - 1] * (xknot[ii + *korder - 1] -
						  xknot[ii - 1]);
	}
	for (ii = lefta - *korder + 1; ii <= lefta; ii++) {
		kount += 1;
		tcoef[kount] = tcoef[kount - 1] + bscoef[ii - 1] * (xknot[ii + *korder - 1] -
							     xknot[ii - 1]);
	}
	imsl_svcal((*korder + 1), ( F_ONE / (Mfloat) (*korder)),
		   tcoef,  1, aj, 1);

        _l0 = *korder + 1;
	vala = imsl_b4itg(&ra, &_l0, xknot, ncoef, aj, dl,
			  dr, &lefta);
	/*
	 * Now when we find LEFTB it should be in the closed interval
	 * (KORDER,NCOEF)
	 */
	if (leftb == 0) {
                _l0 =  *ncoef + *korder;
		imsl_b5itg(xknot, &_l0, &rb, &leftb, &mflag);
	}
	if (leftb - *korder > lefta) {
		sum = tcoef[kount];
		for (ii = lefta + 1; ii <= (leftb - *korder); ii++) {
			sum += bscoef[ii - 1] * (xknot[ii + *korder - 1] - xknot[ii - 1]);
		}
		kount = 0;
		tcoef[kount] = sum;
	} else if (lefta != leftb) {
		kount = 0;
		tcoef[kount] = tcoef[leftb - lefta];
	}
	if (lefta != leftb) {
		for (ii = leftb - *korder + 1; ii <= leftb; ii++) {
			kount += 1;
			tcoef[kount] = tcoef[kount - 1] + bscoef[ii - 1] * (xknot[ii + *korder - 1] -
							     xknot[ii - 1]);
		}
	}
	imsl_svcal((*korder + 1), ( F_ONE / (Mfloat) (*korder)),
		   tcoef, 1, aj, 1);
        _l0 = *korder + 1;
	valb = imsl_b4itg(&rb, &_l0, xknot, ncoef, aj, dl,
			  dr, &leftb);
	b3itg_v = (Mfloat) (isign) * (valb - vala);
L_9000:
	imsl_e1pop("IMSL_B3ITG ");
	return (b3itg_v);
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B4ITG/DB4ITG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Evaluate the integral of a spline, given its B-spline
                representation.

    Usage:      B4ITG(X, KORDER, XKNOT, NCOEF, AJ, DL, DR, I)

    Arguments:
       X      - Point at which spline is to be evaluated.  (Input)
       KORDER - Order of the B-spline.  (Input)
       XKNOT  - Array of length KORDER+NCOEF containing the knot
                sequence.  (Input)
                It must be nondecreasing.
       NCOEF  - Number of B-spline coefficients.  (Input)
       AJ     - Array of length KORDER-1 containing the relevent
                coefficients for the evaluation.  (Input/Output)
                On output AJ is destroyed.
       DL     - Work array of length KORDER-1 to hold differences.
       DR     - Work array of length KORDER-1 to hold differences.
       I      - Pointer to correct knot interval.  (Input)
       B4ITG  - Value of B-sline at X.  (Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* BCOEF is not used here, but leave the calling sequence intact.*/
#ifdef ANSI
Mfloat imsl_b4itg(Mfloat *x, Mint *korder, Mfloat xknot[], Mint *ncoef, 
                  Mfloat aj[], Mfloat dl[], Mfloat dr[], Mint *i)
#else
Mfloat imsl_b4itg(x, korder, xknot, ncoef, aj, dl, dr, i)
	Mfloat          *x;
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           aj[], dl[], dr[];
	Mint            *i;
#endif
{
	Mint             ilo, j, jj;
	Mfloat           b4itg_v, value;

	value = F_ZERO;
	/*
	 * I is such that 1 .LE. I .LT. NKNOT and XKNOT(I) .LT. XKNOT(I+1)
	 * and XKNOT(I) .LE. X .LT. XKNOT(I+1) .
	 * 
	 * The KORDER B-SPLINE coefficients relevant for the XKNOT interval
	 * (XKNOT(I),XKNOT(I+1)) are in AJ(1),...,AJ(KORDER) then compute DL,
	 * DR. Set any of  AJ not obtainable from input to zero. Set any
	 * XKNOTS not obtainable equal to XKNOT(1) or XKNOT(NKNOT)
	 * appropriately.
	 */
	for (j = 1; j <= (*korder - 1); j++) {
		dl[j - 1] = *x - xknot[*i - j];
	}

	for (j = 1; j <= (*korder - 1); j++) {
		dr[j - 1] = xknot[*i + j - 1] - *x;
	}
	/*
	 * Compute value at X in (XKNOT(I),XKNOT(I+1)) given its relevent
	 * B-spline COEFS in AJ(1),...,AJ(KORDER).
	 */
	for (j = 1; j <= (*korder - 1); j++) {
		ilo = *korder - j;
		for (jj = 1; jj <= (*korder - j); jj++) {
			aj[jj - 1] = (aj[jj] * dl[ilo - 1] + aj[jj - 1] * dr[jj - 1]) /
				(dl[ilo - 1] + dr[jj - 1]);
			ilo -= 1;
		}
	}
	value = aj[0];

	b4itg_v = value;
	return (b4itg_v);
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B5ITG/DB5ITG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Compute LEFT = MAX (I, 1 .LE. NKNOT. .AND. XKNOT(I) .LE.
                X

    Usage:      CALL B5ITG (XKNOT, NKNOT, X, LEFT, MFLAG)

    Arguments:
       XKNOT  - The knot sequence.  (Input)
       NKNOT  - Number of knots.  (Input)
       X      - The point whose location in XKNOT is to be found.
                (Input)
       LEFT   - Integer whose value is given below.  (Output)
       MFLAG  - Flag defined below.  (Output)
                 LEFT   MFLAG
                   1     -1      IF                    X .LT.  XKNOT(1)
                   I      0      IF   XKNOT(I)    .LE. X .LT. XKNOT(I+1)
                NKNOT     1      IF  XKNOT(NKNOT) .LE. X
                In particular, MFLAG = 0 is the usual case.  MFLAG .NE. 0
                indicates that X lies outside the half open interval
                XKNOT(1) .LE. X .LT. XKNOT(NKNOT). The asymmetric
                treatment of the interval is due to the decision to make
                all PP functions continuous from the right.

    Remark:
       This routine is based in INTERV in deBoor, p92-93.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b5itg(Mfloat xknot[], Mint *nknot, Mfloat *x, Mint *left, Mint *mflag)
#else
void imsl_b5itg(xknot, nknot, x, left, mflag)
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
			*mflag = 1;
			*left = *nknot;
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
L_10:
		;
		ihi = ilo;
		ilo = ihi - istep;
		if (ilo > 1) {
			if (*x >= xknot[ilo - 1])
				goto L_30;
			istep *= 2;
			goto L_10;
		}
		ilo = 1;
		if (*x < xknot[0]) {
			*mflag = -1;
			*left = 1;
			goto L_9000;
		}
		goto L_30;
	}
	/*
	 * Now X .GE. XKNOT(IHI) . Increase IHI to capture X .
	 */
	istep = 1;
L_20:
	;
	ilo = ihi;
	ihi = ilo + istep;
	if (ihi < *nknot) {
		if (*x < xknot[ihi - 1])
			goto L_30;
		istep *= 2;
		goto L_20;
	}
	if (*x >= xknot[*nknot - 1]) {
		*mflag = 1;
		*left = *nknot;
		goto L_9000;
	}
	ihi = *nknot;
	/*
	 * Now XKNOT(ILO) .LE. X .LT. XKNOT(IHI) . Narrow the interval
	 */
L_30:
	;
	middle = (ilo + ihi) / 2;
	if (middle == ilo) {
		*mflag = 0;
		*left = ilo;
		goto L_9000;
	}
	/*
	 * It is assumed that middle = ILO in case IHI = ILO+1 .
	 */
	if (*x < xknot[middle - 1]) {
		ihi = middle;
	} else {
		ilo = middle;
	}
	goto L_30;
L_9000:
	;
	return;
}				/* end of function */
