#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 12/12/89 at 14:01:33
    Options SET: fmt=t
  -----------------------------------------------------------------------
    IMSL Name:  STRSV  (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 10, 1986

    Purpose:    Solve a triangular system, x = T**(-1)*x, where T is a
                a triangular matrix, all single precision.

    Usage:      CALL STRSV (UPLO, TRANS, DIAG, N, A, LDA, X, INCX)

    Arguments:
       UPLO   - Character specifing the storage structure.
                (Input)
                   UPLO              Structure
                'U' or 'u'      Matrix is upper triangular
                'L' or 'l'      Matrix is lower triangular
       TRANS  - Character specifing if the transpose solution is to be
                computed.  (Input)
                   TRANS              Meaning
                'N' or 'n'      Compute x = A**(-1)*x
                'T' or 't'      Compute x = A'**(-1)*x
                'C' or 'c'      Compute x = A'**(-1)*x
       DIAG   - Character specifing if the matrix has a unit diagonal.
                (Input)
                If DIAG is 'U' or 'u' then the elements in diagonal of A
                are assumed to be one and are not referenced.  If DIAG is
                'N' or 'n' then the actual diagonal elements of A are
                used
       N      - Order of the matrix A.  (Input)
       A      - Triangular matrix of order N.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                calling routine.  (Input)
       X      - Real vector of length (N-1)*IABS(INCX)+1.  (Input/Output)
       INCX   - Displacement between elements of X.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_strsv(Mchar *uplo, Mchar *trans, Mchar *diag, Mint n, Mfloat *a,
                Mint lda, Mfloat *x, Mint incx)
#else
void imsl_strsv(uplo, trans, diag, n, a, lda, x, incx)
	Mchar           *uplo;
	Mchar           *trans;
	Mchar           *diag;
	Mint             n;
	Mfloat           a[];
	Mint             lda;
	Mfloat           x[];
	Mint             incx;
#endif
{
	Mint	        nounit, tran, upper;
	Mint             i, ix;

	/*
	 * Quick return if possible.
	 */
	if (n == 0)
		return;
	nounit = diag[0]=='N' || diag[0]=='n';
	upper = uplo[0]=='U' || uplo[0]=='u';
	tran = trans[0]=='T' || trans[0]=='t' || trans[0]=='C' || trans[0]=='c';

	if (upper) {
		if (tran) {
			if (incx > 0) {
				ix = 1;
				for (i = 1; i <= n; i++) {
					x[ix - 1] -= imsl_sdot(i - 1, &a[lda * (i - 1)], 1,
							  x, incx);
					if (nounit)
						x[ix - 1] /= a[i + lda * (i - 1) - 1];
					ix += incx;
				}
			} else {
				ix = (-n + 1) * incx + 1;
				for (i = 1; i <= n; i++) {
					x[ix - 1] -= imsl_sdot(i - 1, &a[lda * (i - 1)], 1,
						   &x[ix - incx - 1], incx);
					if (nounit)
						x[ix - 1] /= a[i + lda * (i - 1) - 1];
					ix += incx;
				}
			}
		} else if (incx > 0) {
			ix = (n - 1) * incx + 1;
			for (i = n; i >= 1; i--) {
				if (i < n)
					x[ix - 1] -= imsl_sdot(n - i, &a[i + lda * i - 1], lda,
						   &x[ix + incx - 1], incx);
				if (nounit)
					x[ix - 1] /= a[i + lda * (i - 1) - 1];
				ix -= incx;
			}
		} else {
			ix = 1;
			for (i = n; i >= 1; i--) {
				if (i < n)
					x[ix - 1] -= imsl_sdot(n - i, &a[i + lda * i - 1], lda,
							  x, incx);
				if (nounit)
					x[ix - 1] /= a[i + lda * (i - 1) - 1];
				ix -= incx;
			}
		}
	} else if (tran) {
		if (incx > 0) {
			ix = (n - 1) * incx + 1;
			for (i = n; i >= 1; i--) {
				if (i < n)
					x[ix - 1] -= imsl_sdot(n - i, &a[i + lda * (i - 1)],
						1, &x[ix + incx - 1], incx);
				if (nounit)
					x[ix - 1] /= a[i + lda * (i - 1) - 1];
				ix -= incx;
			}
		} else {
			ix = 1;
			for (i = n; i >= 1; i--) {
				if (i < n)
					x[ix - 1] -= imsl_sdot(n - i, &a[i + lda * (i - 1)],
							  1, x, incx);
				if (nounit)
					x[ix - 1] /= a[i + lda * (i - 1) - 1];
				ix -= incx;
			}
		}
	} else if (incx > 0) {
		ix = 1;
		for (i = 1; i <= n; i++) {
			x[ix - 1] -= imsl_sdot(i - 1, &a[i - 1], lda, x, incx);
			if (nounit)
				x[ix - 1] /= a[i + lda * (i - 1) - 1];
			ix += incx;
		}
	} else {
		ix = (-n + 1) * incx + 1;
		for (i = 1; i <= n; i++) {
			x[ix - 1] -= imsl_sdot(i - 1, &a[i - 1], lda, &x[ix - incx - 1],
					  incx);
			if (nounit)
				x[ix - 1] /= a[i + lda * (i - 1) - 1];
			ix += incx;
		}
	}

	return;
}				/* end of function */
