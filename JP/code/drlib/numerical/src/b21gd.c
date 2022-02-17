#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 11/04/91 at 12:24:43
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  B21GD/DB21GD (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 13, 1989

    Purpose:    Evaluate the derivative of a spline on a grid, given its
                B-spline representation.

    Usage:      CALL B21GD (IDERIV, N, XVEC, KORDER, XKNOT, NCOEF,
                            BSCOEF, VALUE, PPCOEF, BREAK, LEFT, H, P,
                            WK)

    Arguments:
       IDERIV - Order of the derivative to be evaluated.  (Input)
                In particular, IDERIV = 0 returns the value of the
                spline.
       N      - Length of vector XVEC.  (Input)
       XVEC   - Array of length N containing the points at which the
                spline is to be evaluated.  (Input)
                XVEC should be strictly increasing.
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NCOEF+KORDER containing the knot
                sequence.  (Input)
                XKNOT must be nondecreasing.
       NCOEF  - Number of B-spline coefficients.  (Input)
       BSCOEF - Array of length NCOEF containing the B-spline
                coefficients.  (Input)
       VALUE  - Array of length N containing the values of the IDERIV-th
                derivative of the spline at the points in XVEC.  (Output)
       PPCOEF - Real array of length KORDER*(NCOEF-KORDER+1).
       BREAK  - Real array of length NCOEF-KORDER+2.
       LEFT   - Integer array of length N.
       H      - Real array of length N.
       P      - Real array of length N.
       WK     - Real array of length (KORDER+3)*KORDER

    Remark:
       Informational error
       Type Code
         4   5  The points in XVEC must be strictly increasing.

    GAMS:       E3; K6

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b21gd (Mint *ideriv, Mint *n, Mfloat xvec[], Mint *korder,
		Mfloat xknot[], Mint *ncoef, Mfloat bscoef[], Mfloat value[], 
		Mfloat *ppcoef, Mfloat break_[], Mint left[], 
		Mfloat h[], Mfloat p[], Mfloat wk[])
#else
void imsl_b21gd (ideriv, n, xvec, korder, xknot, ncoef, bscoef,
                value, ppcoef, break_, left, h, p, wk)
    Mint        *ideriv, *n;
    Mfloat       xvec[];
    Mint        *korder;
    Mfloat       xknot[];
    Mint        *ncoef;
    Mfloat       bscoef[], value[], *ppcoef, break_[];
    Mint         left[];
    Mfloat       h[], p[], wk[];
#endif
{
#define PPCOEF(I_,J_)	(ppcoef+(I_)*(*korder)+(J_))
    Mint         i, incoef, j, nppcf;
    Mfloat       fmm;


    imsl_e1psh ("B21GD");
    /* CHECK N */
    if (*n < 1) {
	imsl_e1sti (1, *n);
        imsl_ermes (IMSL_TERMINAL, IMSL_XVEC_LENGTH);
    }
    /* CHECK IDERIV */
    if (*ideriv < 0) {
	imsl_e1sti (1, *ideriv);

        imsl_ermes(IMSL_TERMINAL, IMSL_IDERIV_NOT_POSITIVE);
    }
    /* CHECK KORDER */
    if (*korder < 1) {
	imsl_e1sti (1, *korder);
        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
    }
    /* Check that XVEC is increasing */
    for (i = 2; i <= *n; i++) {
	if (xvec[i - 2] >= xvec[i - 1]) {
	    imsl_e1sti (1, i - 2);
	    imsl_e1sti (2, i-1);
	    imsl_e1str (1, xvec[i - 2]);
	    imsl_e1str (2, xvec[i - 1]);
            imsl_ermes(IMSL_FATAL, IMSL_XVEC_NOT_INCREASING);
	}
    }
    /* CHECK NCOEF */
    if (*ncoef < *korder) {
	imsl_e1sti (1, *ncoef);
	imsl_e1sti (2, *korder);
        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS);
    }

    /* CHECK FOR ERRORS */
    if (imsl_n1rcd (0) != 0)
	goto L_9000;
    /*
     * COMPUTE VALUE WITHOUT FURTHER CHECKING
     */
    incoef = *ncoef;
    /*
     * CONVERT FROM B-SPLINE REPRESENTATION TO PP-REPRESENTATION
     */
    imsl_b2cpp (korder, xknot, ncoef, bscoef, &nppcf, break_, ppcoef, wk);
    *ncoef = incoef;

    for (i = 1; i <= *n; i++) {
	imsl_p3der (*korder, nppcf, break_, xvec[i - 1], &left[i - 1]);
    }
    for (i = 1; i <= *n; i++) {
	h[i - 1] = xvec[i - 1] - break_[left[i - 1] - 1];
	value[i - 1] = 0.0;
    }

    /*
     * Evaluate jderiv-th derivative of I-th polynomial piece at X
     */
    fmm = *korder - *ideriv;
    for (j = *korder; j >= (*ideriv + 1); j--) {
	for (i = 1; i <= *n; i++) {
	    p[i - 1] = *PPCOEF (left[i - 1] - 1, j - 1);
	    value[i - 1] = (value[i - 1] / fmm) * h[i - 1] + p[i - 1];
	}
	fmm -= 1.0;
    }
L_9000:
    ;
    imsl_e1pop ("B21GD");
    return;
}				/* end of function */
#undef PPCOEF
