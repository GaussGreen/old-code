#define DOUBLE
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 08/02/90 at 15:34:57
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  DSYR2  (Double precision version)

    Computer:   FORC/DOUBLE

    Revised:    September 26, 1989

    Purpose:    Perform the rank-two symmetric matrix update
                A = A + alpha*x*trans(y) + alpha*y*trans(x), where A is
                a symmetric matrix, and trans( ) represents the
                transpose of the vectors.

    Usage:      CALL DSYR2 (UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA)

    Arguments:
       UPLO  -  Character indicating storage form of the matrix.  (Input)
                If UPLO is 'U' or 'u' then the matrix is stored in the
                upper half of A.  If UPLO is 'L' or 'l' then the matrix
                is stored in the lower half of A.
       N      - Order of the symmetric matrix A.  (Input)
       ALPHA  - Double precision scalar multiplier for the vector-vector
                product.  (Input)
       X      - Double precision vector of length (N-1)*IABS(INCX)+1.
                (Input)
       INCX   - Displacement between elements of X.  (Input)
       Y      - Double precision vector of length (N-1)*IABS(INCY)+1.
                (Input)
       INCY   - Displacement between elements of Y.  (Input)
       A      - Symmetric matrix of order N.  (Input/Output)
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
void imsl_dsyr2(Mchar *uplo, Mint uplo_s, Mint *n, Mdouble *alpha, Mdouble x[],
                Mint *incx, Mdouble y[], Mint *incy, Mdouble *a, Mint *lda)
#else
void imsl_dsyr2(uplo, uplo_s, n, alpha, x, incx, y, incy,
	   a, lda)
	Mchar           *uplo;
	unsigned        uplo_s;
	Mint            *n;
	Mdouble         *alpha, x[];
	Mint            *incx;
	Mdouble          y[];
	Mint            *incy;
	Mdouble         *a;
	Mint            *lda;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
/*	LOGICAL32       imsl_l1ame(), lower, upper; */
	Mint             lower, upper;
	Mint             ix, iy, j;


	upper = imsl_l1ame(uplo, uplo_s, "U", sizeof("U"));
	lower = imsl_l1ame(uplo, uplo_s, "L", sizeof("L"));
	/*
	 * Test the input parameters.
	 */
	if (*n < 0) {
		imsl_e1psh("DSYR2 ");
		imsl_e1sti(1, *n);

/*		(5, 1, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("DSYR2 ");
		goto L_9000;
	} else if ((*lda < *n) || (*lda == 0)) {
		imsl_e1psh("DSYR2 ");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		(5, 2, "LDA must be greater than or equal to N and greater than zero while LDA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_LDA_VALUE_GIVEN);
		imsl_e1pop("DSYR2 ");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("DSYR2 ");
		imsl_e1sti(1, *incx);

/*		(5, 3, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("DSYR2 ");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("DSYR2 ");
		imsl_e1sti(1, *incy);

/*		(5, 4, "INCY must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("DSYR2 ");
		goto L_9000;
	} else if ((!upper) && (!lower)) {
		imsl_e1psh("DSYR2 ");
		imsl_e1stl(1, uplo);

/*		(5, 5, "UPLO must be set equal to U or L while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
		imsl_e1pop("DSYR2 ");
		goto L_9000;
	}
	/* Quick return if possible */
	if (*n == 0 || *alpha == F_ZERO)
		goto L_9000;

	ix = 1;
	iy = 1;
	if (*incx < 0)
		ix = (-*n + 1) ** incx + 1;
	if (*incy < 0)
		iy = (-*n + 1) ** incy + 1;

	for (j = 1; j <= *n; j++) {
		if (upper) {
			if (*incx >= 0) {
				saxpy(j, *alpha * y[iy - 1], x, *incx, A(j - 1, 0),
					   1);
			} else {
				saxpy(j, *alpha * y[iy - 1], &x[ix - 1], *incx, A(j - 1, 0),
					   1);
			}
			if (*incy >= 0) {
				saxpy(j, *alpha * x[ix - 1], y, *incy, A(j - 1, 0),
					   1);
			} else {
				saxpy(j, *alpha * x[ix - 1], &y[iy - 1], *incy, A(j - 1, 0),
					   1);
			}
		} else {
			if (*incx >= 0) {
				saxpy(j, *alpha * y[iy - 1], x, *incx, A(0, j - 1),
					   *lda);
			} else {
				saxpy(j, *alpha * y[iy - 1], &x[ix - 1], *incx, A(0, j - 1),
					   *lda);
			}
			if (*incy >= 0) {
				saxpy(j, *alpha * x[ix - 1], y, *incy, A(0, j - 1),
					   *lda);
			} else {
				saxpy(j, *alpha * x[ix - 1], &y[iy - 1], *incy, A(0, j - 1),
					   *lda);
			}
		}
		ix += *incx;
		iy += *incy;
	}

L_9000:
	return;
}				/* end of function */
