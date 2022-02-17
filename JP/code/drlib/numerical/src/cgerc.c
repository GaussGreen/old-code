#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 17:34:50
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CGERC (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    July 18, 1989

    Purpose:    Perform the rank-one matrix update:
                    A = A + alpha*x*ctrans(y),
                where ctrans(y) is the conjugate transpose of the vector.

    Usage:      CALL CGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

    Arguments:
       M      - Number of rows in A.  (Input)
       N      - Number of columns in A.  (Input)
       ALPHA  - Complex scalar.  (Input)
       X      - Complex vector of length (M-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       Y      - Complex vector of length (N-1)*IABS(INCY)+1.  (Input)
       INCY   - Displacement between elements of Y.  (Input)
       A      - Complex array of size M by N.  (Input/Output)
                On input, A contains the matrix to be updated.
                On output, A contains the updated matrix.
       LDA    - Leading dimension of A exactly as specified in the
                calling routine.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_cgerc(Mint *m, Mint *n, Mf_complex *alpha, Mf_complex x[],
		Mint *incx, Mf_complex y[], Mint *incy, Mf_complex a[],
		Mint *lda)
#else
void imsl_cgerc(m, n, alpha, x, incx, y, incy, a, lda)
	Mint            *m, *n;
	Mf_complex      *alpha, x[];
	Mint            *incx;
	Mf_complex       y[];
	Mint            *incy;
	Mf_complex       a[];
	Mint            *lda;
#endif
{
	Mint             _l0, imsl_i1x, iy, j;
	Mf_complex       _cx0;

	/*
	 * Test the input parameters.
	 */
	if (*m < 0) {
		imsl_e1psh("CGERC ");
		imsl_e1sti(1, *m);

/*		imsl_ermes(5, 1, "M must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_M_GE_ZERO);
		imsl_e1pop("CGERC ");
		goto L_9000;
	} else if (*n < 0) {
		imsl_e1psh("CGERC ");
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 2, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("CGERC ");
		goto L_9000;
	} else if ((*lda < *m) || (*lda == 0)) {
		imsl_e1psh("CGERC ");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *m);

/*		imsl_ermes(5, 3, "LDA must be greater than or equal to M and greater than zero while LDA = %(i1) and M = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDA_GE_M);
		imsl_e1pop("CGERC ");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("CGERC ");
		imsl_e1sti(1, *incx);

/*		imsl_ermes(5, 4, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("CGERC ");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("CGERC ");
		imsl_e1sti(1, *incy);

/*		imsl_ermes(5, 5, "INCY must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("CGERC ");
		goto L_9000;
	}
	/*
	 * Quick return if possible
	 */
	if ((*m == 0 || *n == 0) || imsl_c_eq(*alpha, imsl_cf_convert(F_ZERO, F_ZERO)))
		goto L_9000;

	iy = 1;
	if (*incy < 0)
		iy = (-*n + 1) ** incy + 1;

	imsl_i1x = 1;
	_l0 = 1;
	for (j = 1; j <= *n; j++) {
		_cx0 = imsl_c_mul(*alpha, imsl_c_conjg(y[iy - 1]));
		imsl_caxpy(m, &_cx0, x, incx,
			   &a[imsl_i1x - 1], &_l0);
		iy += *incy;
		imsl_i1x += *lda;
	}

L_9000:
	return;
}				/* end of function */
