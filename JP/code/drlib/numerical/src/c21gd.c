#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 11/04/91 at 10:31:33
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C21GD/DC21GD (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 25, 1989

    Purpose:    Evaluate the derivative of a cubic spline on a grid.

    Usage:      CALL C21GD (IDERIV, N, XVEC, NINTV, BREAK, CSCOEF,
                            VALUE, LEFT, H, P)

    Arguments:  See CS1GD.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	KORDER	4

#ifdef ANSI
void imsl_c21gd (Mint *ideriv, Mint *n, Mfloat xvec[], Mint *nintv,
		Mfloat break_[], Mfloat **cscoef, Mfloat value[],
		Mint left[], Mfloat h[], Mfloat p[])
#else
void imsl_c21gd (ideriv, n, xvec, nintv, break_, cscoef, value,
                left, h, p)
    Mint        *ideriv, *n;
    Mfloat       xvec[];
    Mint        *nintv;
    Mfloat       break_[], cscoef[][4], value[];
    Mint         left[];
    Mfloat       h[], p[];
#endif
{
    Mint         i, j;
    Mfloat       fmm;


    imsl_e1psh ("C21GD ");
    /* IN CASE OF ERRORS */
    for (i = 1; i <= *n; i++) {
	value[i - 1] = 0.0;
    }
    /* Check argument NINTV */
    if (*nintv < 1) {
	imsl_e1sti (1, *nintv);
	imsl_ermes(IMSL_TERMINAL, IMSL_NINTV_NOT_POSITIVE);
    }
    /* Check argument IDERIV */
    if (*ideriv < 0) {
	imsl_e1sti (1, *ideriv);
	imsl_ermes(IMSL_TERMINAL, IMSL_IDERIV_NOT_POSITIVE);
    }
    /* Check argument N */
    if (*n < 1) {
	imsl_e1sti (1, *n);
	imsl_ermes(IMSL_TERMINAL, IMSL_XVEC_LENGTH);
    }
    /* Check for errors */
    if (imsl_n1rty (0) != 0)
	goto L_9000;

    /* Check that XVEC is increasing */
    for (i = 2; i <= *n; i++) {
	if (xvec[i - 2] >= xvec[i - 1]) {
	    imsl_e1sti (1, i - 2);
	    imsl_e1sti (2, i-1);
	    imsl_e1str (1, xvec[i - 2]);
	    imsl_e1str (2, xvec[i - 1]);
	    imsl_ermes(IMSL_FATAL, IMSL_XVEC_NOT_INCREASING);
	    goto L_9000;
	}
    }

    for (i = 1; i <= *n; i++) {
	imsl_p3der (KORDER, *nintv, break_, xvec[i - 1], &left[i - 1]);
    }
    for (i = 1; i <= *n; i++) {
	h[i - 1] = xvec[i - 1] - break_[left[i - 1] - 1];
	value[i - 1] = 0.0;
    }

    /*
     * Evaluate jderiv-th derivative of I-th polynomial piece at X
     */
    fmm = KORDER - *ideriv;
    for (j = KORDER; j >= (*ideriv + 1); j--) {
      for (i = 1; i <= *n; i++) {
	p[i - 1] = *((Mfloat *)cscoef + (j-1) + 4*(left[i - 1] - 1));
	/* 
	 * The following line cause segmentation violations on ANSI platforms.
	 * it has been replaced by the line above.
	 */
/*	p[i - 1] = cscoef[left[i - 1] - 1][j - 1];*/
	value[i - 1] = (value[i - 1] / fmm) * h[i - 1] + p[i - 1];
      }
      fmm -= 1.0;
    }
    
L_9000:
    ;
    imsl_e1pop ("C21GD ");
    return;
}				/* end of function */

#undef KORDER
