#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:20:34
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CRGRG/DCRGRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 5, 1985

    Purpose:    Copy a real general matrix.

    Usage:      CALL CRGRG (N, A, LDA, B, LDB)

    Arguments:
       N      - Order of the matrices.  (Input)
       A      - Matrix of order N.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Matrix of order N containing a copy of A.  (Output)
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)

    Keyword:    Basic matrix operation

    GAMS:       D1b8

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_crgrg(Mint n, Mfloat *a, Mint lda, Mfloat *b, Mint ldb)
#else
void imsl_crgrg(n, a, lda, b, ldb)
	Mint             n;
	Mfloat          *a;
	Mint             lda;
	Mfloat          *b;
	Mint             ldb;
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
#define B(I_,J_)	(b+(I_)*(ldb)+(J_))
	Mint             j;


	imsl_e1psh("CRGRG");
	/* Check N */
	if (n < 1) {
		imsl_e1sti(1, n);
/*	        (5, 1, "The argument N = %(i1).  It must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LARGER_N_VALUE_NEEDED);
		goto L_9000;
	}
	/* Check LDA */
	if (lda < n) {
		imsl_e1sti(1, lda);
		imsl_e1sti(2, n);
/*		(5, 2, "The argument LDA = %(i1).  It must be at least as large as N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LARGER_LDA_VALUE_NEEDED);
		goto L_9000;
	}
	/* Check LDB */
	if (ldb < n) {
		imsl_e1sti(1, ldb);
		imsl_e1sti(2, n);
/*		(5, 3, "The argument LDB = %(i1).  It must be at least as large as N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LARGER_LDB_VALUE_NEEDED);
		goto L_9000;
	}
	/* Copy */
	if (lda == n && ldb == n) {
		imsl_scopy(n * n, a, 1, b, 1);
	} else if (lda >= ldb) {
		for (j = 1; j <= n; j++) {
			imsl_scopy(n, A(j - 1, 0), 1, B(j - 1, 0), 1);
		}
	} else {
		for (j = n; j >= 1; j--) {
			imsl_scopy(n, A(j - 1, 0), -1, B(j - 1, 0), -1);
		}
	}

L_9000:
	;
	imsl_e1pop("CRGRG");
	return;
}				/* end of function */
