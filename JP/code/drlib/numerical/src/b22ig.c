#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B22IG/DB22IG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Evaluate the integral of a tensor product spline on a
                rectangular domain, given its tensor product B-spline
                representation.

    Usage:      B22IG(A, B, C, D, KXORD, KYORD, XKNOT, YKNOT, NXCOEF,
                      NYCOEF, BSCOEF, WK)

    Arguments:
       A      - Lower limit of the X-variable.  (Input)
       B      - Upper limit of the X-variable.  (Input)
       C      - Lower limit of the Y-variable.  (Input)
       D      - Upper limit of the Y-variable.  (Input)
       KXORD  - Order of the spline in the X-direction.  (Input)
       KYORD  - Order of the spline in the Y-direction.  (Input)
       XKNOT  - Array of length NXCOEF+KXORD containing the knot
                sequence in the X-direction.  (Input)
                XKNOT must be nondecreasing.
       YKNOT  - Array of length NYCOEF+KYORD containing the knot
                sequence in the Y-direction.  (Input)
                YKNOT must be nondecreasing.
       NXCOEF - Number of B-spline coefficients in the X-direction.
                (Input)
       NYCOEF - Number of B-spline coefficients in the Y-direction.
                (Input)
       BSCOEF - Array of length NXCOEF*NYCOEF containing the
                tensor-product B-spline coefficients.  (Input)
                BSCOEF is treated internally as a matrix of size
                NXCOEF by NYCOEF.
       B22IG  - Integral of the spline over the rectangle (A,B) by (C,D).
                (Output)

    Remark:
       Informational errors
       Type Code
         3   1  The lower limit of the X integration is less than
                XKNOT(KXORD).
         3   2  The upper limit of the X integration is greater than
                XKNOT(NXCOEF+1).
         3   3  The lower limit of the Y integration is less than
                YKNOT(KYORD).
         3   4  The upper limit of the Y integration is greater than
                YKNOT(NYCOEF+1).
         4   13 Multiplicity of the knots cannot exceed the order
                of the spline.
         4   14 The knots must be nondecreasing.

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
Mfloat imsl_b22ig(Mfloat *a, Mfloat *b, Mfloat *c, Mfloat *d, Mint *kxord, Mint *kyord,
                  Mfloat xknot[], Mfloat yknot[], Mint *nxcoef,
	          Mint *nycoef, Mfloat *bscoef, Mfloat wk[])
#else
Mfloat imsl_b22ig(a, b, c, d, kxord, kyord, xknot, yknot, nxcoef,
	   nycoef, bscoef, wk)
	Mfloat          *a, *b, *c, *d;
	Mint            *kxord, *kyord;
	Mfloat           xknot[], yknot[];
	Mint            *nxcoef, *nycoef;
	Mfloat          *bscoef, wk[];
#endif
{
#define BSCOEF(I_,J_)	(bscoef+(I_)*(*nxcoef)+(J_))
	Mchar            i1[27], i2[30], i3[32], l1[22], r1[23];
	Mint             iaj, idl, idr, isign1, isign2, itcoef, itmpcf,
	                kord;
	Mfloat           b22ig_v, ra, rb, rc, rd, val;


	imsl_e1psh("IMSL_B22IG ");
	/* Check KXORD */
	val = F_ZERO;
	if (*kxord < 1) {
		imsl_e1sti(1, *kxord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_X);
	}
	/* Check KYORD */
	if (*kyord < 1) {
		imsl_e1sti(1, *kyord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_Y);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check NXCOEF */
	if (*nxcoef < *kxord) {
		imsl_e1sti(1, *nxcoef);
		imsl_e1sti(2, *kxord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_X);
	}
	/* Check NYCOEF */
	if (*nycoef < *kyord) {
		imsl_e1sti(1, *nycoef);
		imsl_e1sti(2, *kyord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_Y);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check XKNOT */
	imsl_b32in("X", kxord, xknot, nxcoef);
	/* Check YKNOT */
	imsl_b32in("Y", kyord, yknot, nycoef);
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	/*
	 * If A=B or C=D then integral is zero.
	 */
	if (*a == *b || *c == *d) {
		val = F_ZERO;
	} else {
		/*
		 * Assign RA to MIN(A,B) Assign RB to MAX(A,B)
		 */
		isign1 = 1;
		isign2 = 1;
		if (*a < *b) {
			ra = *a;
			rb = *b;
			isign1 = 1;
		} else if (*a > *b) {
			ra = *b;
			rb = *a;
			isign1 = -1;
		}
		strcpy(l1, "The left endpoint of ");
		strcpy(r1, "The right endpoint of ");
		strcpy(i1, " integration is less than ");
		strcpy(i2, " integration is greater than ");
		strcpy(i3, "  Integration occurs only from ");

		/*
		 * If RA is to the left of XKNOT(KXORD) then reset RA to
		 * XKNOT(KXORD).
		 */
		if (ra < xknot[*kxord - 1]) {
			ra = xknot[*kxord - 1];
			if (isign1 == 1) {

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_LEFT_ENDPT);
			} else {

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_RIGHT_ENDPT);
			}
		}
		/*
		 * If RB is to the right of XKNOT(NXCOEF+1) then reset RB to
		 * XKNOT(NXCOEF+1).
		 */
		if (rb > xknot[*nxcoef]) {
			rb = xknot[*nxcoef];
			if (isign1 == 1) {

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_RIGHT_ENDPT_1);
			} else {

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_LEFT_ENDPT_1);
			}
		}
		/*
		 * If RA=RB or both RA and RB are outside the range of the
		 * B-sline then the integral is zero.
		 */
		if ((ra <= xknot[*kxord - 1] && rb <= xknot[*kxord - 1]) ||
		    (ra >= xknot[*nxcoef] && rb >= xknot[*nxcoef]))
			goto L_9000;
		if (ra == rb)
			goto L_9000;

		if (*c < *d) {
			rc = *c;
			rd = *d;
			isign2 = 1;
		} else if (*c > *d) {
			rc = *d;
			rd = *c;
			isign2 = -1;
		}
		/*
		 * If RC is to the left of YKNOT(KYORD) then reset RC to
		 * YKNOT(KYORD).
		 */
		if (rc < yknot[*kyord - 1]) {
			rc = yknot[*kyord - 1];
			if (isign2 == 1) {
		strcpy(l1, "The left endpoint of ");
		strcpy(r1, "The right endpoint of ");
		strcpy(i1, " integration is less than ");
		strcpy(i2, " integration is greater than ");
		strcpy(i3, "  Integration occurs only from ");

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_LEFT_ENDPT_2);
			} else {

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_RIGHT_ENDPT_2);
			}
		}
		/*
		 * If RD is to the right of YKNOT(NYCOEF+1) then reset RD to
		 * YKNOT(NYCOEF+1).
		 */
		if (rd > yknot[*nycoef]) {
			rd = yknot[*nycoef];
			if (isign2 == 1) {

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_RIGHT_ENDPT_3);
			} else {

				imsl_ermes(IMSL_WARNING,
				IMSL_SPLINE_LEFT_ENDPT_3);
			}
		}
		/*
		 * If RC=RD or both RC and RD are outside the range of the
		 * B-sline then the integral is zero.
		 */
		if ((rc <= yknot[*kyord - 1] && rd <= yknot[*kyord - 1]) ||
		    (rc >= yknot[*nycoef] && rd >= yknot[*nycoef]))
			goto L_9000;
		if (rc == rd)
			goto L_9000;
		/* Partition workspace */
		kord = imsl_i_max(*kxord, *kyord);
		itmpcf = 1;
		itcoef = itmpcf + *nycoef;
		iaj = itcoef + kord + 1;
		idl = iaj + kord + 1;
		idr = idl + kord + 1;
		val = imsl_b32ig(&ra, &rb, &rc, &rd, kxord, kyord, xknot, yknot,
		   nxcoef, nycoef, bscoef, &wk[itmpcf - 1], &wk[itcoef - 1],
				 &wk[iaj - 1], &wk[idl - 1], &wk[idr - 1]);
		val *= (Mfloat) (isign1 * isign2);
	}
L_9000:
	imsl_e1pop("IMSL_B22IG ");
	b22ig_v = val;
	return (b22ig_v);
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B32IG/DB32IG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Evaluate the integral of a tensor product spline on a
                rectangular domain, given its tensor product B-spline
                representation.

    Usage:      B32IG(A, B, C, D, KXORD, KYORD, XKNOT, YKNOT, NXCOEF,
                      NYCOEF, BSCOEF, TMPCF, TCOEF, AJ, DL, DR)

    Arguments:
       A      - Lower limit of the X variable.  (Input)
       B      - Upper limit of the X variable.  (Input)
       C      - Lower limit of the Y variable.  (Input)
       D      - Upper limit of the Y variable.  (Input)
       KXORD  - Order of the spline in the X direction.  (Input)
       KYORD  - Order of the spline in the Y direction.  (Input)
       XKNOT  - Array of length NXCOEF+KXORD containing the knot
                sequence in the X direction.  (Input)
                It must be nondecreasing.
       YKNOT  - Array of length NYCOEF+KYORD containing the knot
                sequence in the Y direction.  (Input)
                It must be nondecreasing.
       NXCOEF - Number of spline coefficients in the X direction.
                (Input)
       NYCOEF - Number of spline coefficients in the Y direction.
                (Input)
       BSCOEF - Array of length NXCOEF*NYCOEF containing the B-spline
                representation.  (Input)
                BSCOEF is treated internally as a mtrix of size
                NXCOEF by NYCOEF.
       TMPCF  - Work array of length NYCOEF to hold new B-spline
                coefficients in the Y-direction.
       TCOEF  - Work vector of length MAX(KXORD,KYORD)+1.
       AJ     - Work vector of length MAX(KXORD,KYORD)+1.
       DL     - Work vector of length MAX(KXORD,KYORD)+1.
       DR     - Work vector of length MAX(KXORD,KYORD)+1.
       B32IG  - Integral of the spline over the rectangle (A,B)x(C,D).
                (Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_b32ig(Mfloat *a, Mfloat *b, Mfloat *c, Mfloat *d, Mint *kxord, Mint *kyord,
                  Mfloat xknot[], Mfloat yknot[], Mint *nxcoef,
	          Mint *nycoef, Mfloat *bscoef, Mfloat tmpcf[], Mfloat tcoef[], 
                  Mfloat aj[], Mfloat dl[], Mfloat dr[]) 
#else
Mfloat imsl_b32ig(a, b, c, d, kxord, kyord, xknot, yknot, nxcoef,
	   nycoef, bscoef, tmpcf, tcoef, aj, dl, dr)
	Mfloat          *a, *b, *c, *d;
	Mint            *kxord, *kyord;
	Mfloat           xknot[], yknot[];
	Mint            *nxcoef, *nycoef;
	Mfloat          *bscoef, tmpcf[], tcoef[], aj[], dl[], dr[];
#endif
{
#define BSCOEF(I_,J_)	(bscoef+(I_)*(*nxcoef)+(J_))
	Mint             i;
	Mfloat           b32ig_v, val;


	val = F_ZERO;
	/*
	 * If A=B or C=D then integral is zero.
	 */
	if (*a == *b || *c == *d) {
		val = F_ZERO;
	} else {
		/*
		 * Generate temporary Y coefficients Later change this loop
		 * to generate only those coefficients needed.
		 */
		for (i = 1; i <= *nycoef; i++) {
			tmpcf[i - 1] = imsl_b3itg(a, b, kxord, xknot, nxcoef, BSCOEF(i - 1, 0),
						  tcoef, aj, dl, dr);
		}
		if (imsl_n1rty(1) == 3) {

/*imsl_ermes(3, 1, "At least one of the limits of 
integration of the X variable is not in the closed interval (XKNOT(KXORD),
XKNOT(NXCOEF+1)).  The integration region of the X variable is set to the 
intersection of the above interval and the closed interval (A,B)."); */

                        imsl_ermes(IMSL_WARNING, IMSL_SPLINE_LIMITS_X);
/*			imsl_ermes(3, 1, "ERROR TOO LONG.  The integration region of the X variable is set to the
intersection.");*/

		}
		/* Compute answer */
		val = imsl_b3itg(c, d, kyord, yknot, nycoef, tmpcf, tcoef, aj,
				 dl, dr);
		if (imsl_n1rty(1) == 3) {

/*imsl_ermes(3, 2, "At least one of the limits of integration of the Y 
variable is not in the closed interval (YKNOT(KYORD),YKNOT(NYCOEF+1)).  
The integration region of the Y variable is set to the intersection of 
the above interval and the closed interval (C,D).");*/

                        imsl_ermes(IMSL_WARNING, IMSL_SPLINE_LIMITS_Y);
/*			imsl_ermes(3, 1, "ERROR TOO LONG.  The integration region of the Y variable is set to the
intersection.");*/
		}
	}
	b32ig_v = val;
	return (b32ig_v);
}				/* end of function */
