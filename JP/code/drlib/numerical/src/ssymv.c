#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/26/90 at 14:42:27
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SSYMV (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 14, 1989

    Purpose:    Perform the matrix-vector multiplication:
                y = alpha*A*x + imsl_beta*y, where A is a symmetric matrix.

    Usage:      CALL SSYMV (UPLO, N, ALPHA, A, LDA, X, INCX, BETA,
                            Y, INCY)

    Arguments:
       UPLO   - Character specifing the storage structure.
                (Input)
                If UPLO is 'U' or 'u', then the upper part of A is used.
                If UPLO is 'L' or 'l', then the lower part of A is used.
       N      - Order of the matrix A.  (Input)
       ALPHA  - Scalar multiplier for the matrix-vector product.
                (Input)
       A      - Array of size LDA by N containing the symmetric matrix.
                (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling routine.  (Input)
       X      - Vector of length (N-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       BETA   - Scalar multiplier for Y.  (Input)
                When BETA is zero, Y is not referenced.  In that case,
                BETA*Y is defined as the zero vector.
       Y      - Vector of length (N-1)*IABS(INCY)+1.  (Input/Output)
                On output, Y contains the updated Y.
       INCY   - Displacement between elements of Y.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_ssymv(Mchar *uplo, Mint uplo_s, Mint *n, Mfloat *alpha, Mfloat *a,
                Mint *lda, Mfloat *x, Mint *incx, Mfloat *imsl_beta, Mfloat *y,
                Mint *incy)
#else
void imsl_ssymv(uplo, uplo_s, n, alpha, a, lda, x, incx, imsl_beta, y, incy)
	Mchar           *uplo;
	Mint            uplo_s;
	Mint            *n;
	Mfloat          *alpha, *a;
	Mint            *lda;
	Mfloat           x[];
	Mint            *incx;
	Mfloat          *imsl_beta, y[];
	Mint            *incy;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             lower, upper;
	Mint             ix, iy, j, ky;
	Mfloat           temp;


	upper = imsl_l1ame(uplo, uplo_s, "U", sizeof("U"));
	lower = imsl_l1ame(uplo, uplo_s, "L", sizeof("L"));
	/*
	 * Test the input parameters.
	 */
	if (*n < 0) {
		imsl_e1psh("SSYMV ");
		imsl_e1sti(1, *n);

/*		(5, 1, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("SSYMV ");
		goto L_9000;
	} else if ((*lda < *n) || (*lda == 0)) {
		imsl_e1psh("SSYMV ");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		(5, 2, "LDA must be greater than or equal to N and greater than zero while LDA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_LDA_VALUE_GIVEN);
		imsl_e1pop("SSYMV ");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("SSYMV ");
		imsl_e1sti(1, *incx);

/*		(5, 3, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("SSYMV ");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("SSYMV ");
		imsl_e1sti(1, *incy);

/*		(5, 4, "INCY must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("SSYMV ");
		goto L_9000;
	} else if ((!upper) && (!lower)) {
		imsl_e1psh("SSYMV ");
		imsl_e1stl(1, uplo);

/*		(5, 5, "UPLO must be set equal to U or L while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
		imsl_e1pop("SSYMV ");
		goto L_9000;
	}
	/* Quick return if possible */
	if (*n == 0 || (*alpha == F_ZERO && *imsl_beta == F_ONE))
		goto L_9000;

	ix = 1;
	iy = 1;
	if (*incx < 0)
		ix = (-*n + 1) ** incx + 1;
	if (*incy < 0)
		iy = (-*n + 1) ** incy + 1;

	if (*imsl_beta != F_ONE) {
		if (*imsl_beta == F_ZERO) {
			sset(*n, F_ZERO, y, abs(*incy));
		} else {
			sscal(*n, *imsl_beta, y, abs(*incy));
		}
	}
	if (*alpha == F_ZERO)
		goto L_9000;

	if (upper) {
		for (j = 1; j <= *n; j++) {
			temp = *alpha * x[ix - 1];
			ky = iy + (j - 2) * imsl_i_min(*incy, 0);
			saxpy(j - 1, temp, A(j - 1, 0), 1, &y[ky - 1], *incy);
			ky = iy + (j - 1) ** incy + (*n - j) * imsl_i_min(*incy, 0);
			saxpy(*n - j + 1, temp, A(j - 1, j - 1), *lda, &y[ky - 1],
				   *incy);
			ix += *incx;
		}
	} else {
		for (j = 1; j <= *n; j++) {
			temp = *alpha * x[ix - 1];
			ky = iy + (j - 2) * imsl_i_min(*incy, 0);
			saxpy(j - 1, temp, A(0, j - 1), *lda, &y[ky - 1], *incy);
			ky = iy + (j - 1) ** incy + (*n - j) * imsl_i_min(*incy, 0);
			saxpy(*n - j + 1, temp, A(j - 1, j - 1), 1, &y[ky - 1],
				   *incy);
			ix += *incx;
		}
	}

L_9000:
	return;
}				/* end of function */
