#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/26/90 at 14:44:09
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CSFRG/DCSFRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 23, 1985

    Purpose:    Extend a real symmetric matrix defined in its upper
                triangle to its lower triangle.

    Usage:      CALL CSFRG (N, A, LDA)

    Arguments:
       N      - Order of the matrix A.  (Input)
       A      - N by N symmetric matrix of order N to be filled out.
                (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)

    Keywords:   Basic matrix operation; Matrix conversion

    GAMS:       D1b9

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_csfrg(Mint *n, Mfloat *a, Mint *lda)
#else
void imsl_csfrg(n, a, lda)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             i;


#ifdef DOUBLE
	imsl_e1psh("imsl_dcsfrg");
#else
	imsl_e1psh("imsl_csfrg");
#endif

	if (*n <= 0) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "N = %(i1).  The order of A, N, must be greater than 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_A_AND_N_GT_ZERO);
	}
	if (*n > *lda) {
		imsl_e1sti(1, *n);
		imsl_e1sti(2, *lda);

/*		imsl_ermes(5, 2, "N = %(i1) and LDA = %(i2).  The order of A, N, must be less than or equal to the leading dimension of A, LDA.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_LE_LDA);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	/*
	 * Copy upper triangular values to lower triangular values
	 */
	for (i = 1; i <= (*n - 1); i++) {
		scopy(*n - i, A(i, i - 1), *lda, A(i - 1, i), 1);
	}
	/* Exit section */
L_9000:
#ifdef DOUBLE
	imsl_e1pop("imsl_dcsfrg");
#else
	imsl_e1pop("imsl_csfrg");
#endif
	return;
}				/* end of function */
