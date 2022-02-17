#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/30/90 at 09:58:41
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SGEMV  (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 26, 1989

    Purpose:    Perform one of the matrix-vector multiplications:
                    y = alpha*A*x + imsl_beta*y
                    y = alpha*trans(A)*x + imsl_beta*y
                Here trans(A) is the transpose of the matrix.

    Usage:      CALL SGEMV (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y,
                            INCY)

    Arguments:
       TRANS  - Character specifing the operation to be performed.
                (Input)
                   TRANS               Operation
                'N' or 'n'      y = alpha*A*x + imsl_beta*y
                'T' or 't'      y = alpha*trans(A)*x + imsl_beta*y
                'C' or 'c'      y = alpha*trans(A)*x + imsl_beta*y
       M      - Number of rows in A.  (Input)
       N      - Number of columns in A.  (Input)
       ALPHA  - Scalar multiplier for a matrix-vector product.  (Input)
       A      - Array of size M by N.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                calling routine.  (Input)
       X      - Vector of length (N-1)*IABS(INCX)+1 when TRANS is
                'N' or 'n' and of length (M-1)*IABS(INCX)+1 otherwise.
                (Input)
       INCX   - Displacement between elements of X.  (Input)
       BETA   - Scalar multiplier for Y.  (Input)
                When BETA is zero, Y need not be set on input.
                In that case, BETA*Y is defined as the zero vector.
       Y      - Vector of length (N-1)*IABS(INCY)+1 when TRANS is
                'M' or 'm' and of length (M-1)*IABS(INCY)+1 otherwise.
                (Input/Output)
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
void imsl_sgemv(Mchar *trans, Mint trans_s, Mint *m, Mint *n, Mfloat *alpha,
                Mfloat *a, Mint *lda, Mfloat *x, Mint *incx, Mfloat *imsl_beta,
                Mfloat *y, Mint *incy)
#else
void imsl_sgemv(trans, trans_s, m, n, alpha, a, lda, x, incx, imsl_beta, y, incy)
	Mchar           *trans;
	Mint             trans_s;
	Mint            *m, *n;
	Mfloat          *alpha, a[];
	Mint            *lda;
	Mfloat           x[];
	Mint            *incx;
	Mfloat          *imsl_beta, y[];
	Mint            *incy;
#endif
{
	Mint	        imsl_ctran, ntran, tran;
	Mint             i, ix, iy, kx, ky, lenx, leny;


	ntran = imsl_l1ame(trans, trans_s, "N", sizeof("N"));
	tran = imsl_l1ame(trans, trans_s, "T", sizeof("T"));
	imsl_ctran = imsl_l1ame(trans, trans_s, "C", sizeof("C"));
	/*
	 * Test the input parameters.
	 */
	if (*m < 0) {
		imsl_e1psh("SGEMV ");
		imsl_e1sti(1, *m);

/*		(5, 1, "M must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_M_GE_ZERO);
		imsl_e1pop("SGEMV ");
		goto L_9000;
	} else if (*n < 0) {
		imsl_e1psh("SGEMV ");
		imsl_e1sti(1, *n);

/*		(5, 2, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("SGEMV ");
		goto L_9000;
	} else if ((*lda < *m) || (*lda == 0)) {
		imsl_e1psh("SGEMV ");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *m);

/*		(5, 3, "LDA must be greater than or equal to M and greater than zero while LDA = %(i1) and M = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDA_GE_M);
		imsl_e1pop("SGEMV ");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("SGEMV ");
		imsl_e1sti(1, *incx);

/*		(5, 4, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("SGEMV ");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("SGEMV ");
		imsl_e1sti(1, *incy);

/*		(5, 5, "INCY must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("SGEMV ");
		goto L_9000;
	} else if (((!ntran) && (!tran)) && (!imsl_ctran)) {
		imsl_e1psh("SGEMV ");
		imsl_e1stl(1, trans);

/*		(5, 6, "TRANS must be set equal to N or T or C while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TRANS_MUST_EQUAL_N_T_OR_C);
		imsl_e1pop("SGEMV ");
		goto L_9000;
	}
	/* Quick return if possible */
	if ((*m == 0 || *n == 0) || ((*alpha == F_ZERO) && (*imsl_beta == F_ONE)))
		goto L_9000;

	if (ntran) {
		lenx = *n;
		leny = *m;
	} else {
		lenx = *m;
		leny = *n;
	}

	ix = 1;
	iy = 1;
	if (*incx < 0)
		ix = (-lenx + 1) ** incx + 1;
	if (*incy < 0)
		iy = (-leny + 1) ** incy + 1;

	if (*imsl_beta == F_ONE) {
	} else if (*incy == 0) {
		if (*imsl_beta == F_ZERO) {
			y[0] = F_ZERO;
		} else {
			y[0] *= powl(*imsl_beta, leny);
		}
	} else if (*imsl_beta == F_ZERO) {
		sset(leny, F_ZERO, y, abs(*incy));
	} else {
		sscal(leny, *imsl_beta, y, abs(*incy));
	}

	if (*alpha == F_ZERO)
		goto L_9000;

	/* Not transpose */
	if (ntran) {
		kx = ix;
		for (i = 1; i <= *n; i++) {
			saxpy(*m, *alpha * x[kx - 1], &a[*lda * (i - 1)], 1, y, *incy);
			kx += *incx;
		}
	} else {
		/* Transpose */
		ky = iy;
		for (i = 1; i <= *n; i++) {
			y[ky - 1] += *alpha * imsl_sdot(*m, &a[*lda * (i - 1)], 1, x,
							*incx);
			ky += *incy;
		}
	}

L_9000:
	return;
}				/* end of function */
