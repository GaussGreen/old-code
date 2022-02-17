#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static void PROTO(l_saxpy,(Mint n, Mfloat sa, Mfloat *sx, Mint incx, Mfloat *sy, Mint incy));
#ifdef INCL_PROTOTYPES
#include "/home/usr2/imsl1/clib/newclib/include/imsl_int.h"
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/14/90 at 13:56:29
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  STBSV (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 26, 1989

    Purpose:    Solve one of the triangular systems,
                    x = inv(A)*x,
                 or
                    x = inv(trans(A))*x,
                where A is a triangular matrix in band storage mode,
                and trans(A) is the transpose of the matrix.

    Usage:      CALL STBSV (UPLO, TRANS, DIAG, N, NCODA, A, LDA, X, INCX)

    Arguments:
       UPLO   - Character specifing the storage structure.
                (Input)
                   UPLO              Structure
                'U' or 'u'      Matrix is upper triangular
                'L' or 'l'      Matrix is lower triangular
       TRANS  - Character specifing if the transpose solution is to be
                computed.  (Input)
                   TRANS              Meaning
                'N' or 'n'      Compute x = inv(A)*x
                'T' or 't'      Compute x = inv(trans(A))*x
                'C' or 'c'      Compute x = inv(trans(A))*x
       DIAG   - Character specifing if A is unit triangular.
                (Input)
                   DIAG               Meaning
                'U' or 'u'      A is assumed to be unit triangular.
                'N' or 'n'      A is not assumed to be unit triangular.
                If DIAG is 'U' or 'u' then the diagonal elements
                of A are assumed to be one and are not referenced.  If
                DIAG is 'N' or 'n' then the actual diagonal elements of
                array A are used.  (Input)
       N      - Order of the matrix A.  (Input)
       NCODA  - Number of codiagonals in A.  (Input)
       A      - Band triangular matrix of order N in band storage mode.
                (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       X      - Vector of length (N-1)*IABS(INCX)+1.  (Input/Output)
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
void imsl_stbsv(Mchar *uplo, Mchar *trans, Mchar *diag,
	    Mint *n, Mint *ncoda, Mfloat *a, Mint *lda, Mfloat *x, Mint *incx)
#else
void imsl_stbsv(uplo, trans, diag,
	   n, ncoda, a, lda, x, incx)
	Mchar           *uplo;
	Mchar           *trans;
	Mchar           *diag;
	Mint            *n, *ncoda;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           x[];
	Mint            *incx;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
/*	LOGICAL32       imsl_ctran, imsl_l1ame(), lower, ndiag, ntran,
	                tran, udiag, upper;*/
	Mint            imsl_ctran, lower, ndiag, ntran,
	                tran, udiag, upper;
	Mint             ix, j, jx, kx, m;


	upper = imsl_l1ame(uplo, sizeof(uplo), "U", sizeof("U"));
	lower = imsl_l1ame(uplo, sizeof(uplo), "L", sizeof("L"));
	udiag = imsl_l1ame(diag, sizeof(diag), "U", sizeof("U"));
	ndiag = imsl_l1ame(diag, sizeof(diag), "N", sizeof("N"));
	ntran = imsl_l1ame(trans, sizeof(trans), "N", sizeof("N"));
	tran = imsl_l1ame(trans, sizeof(trans), "T", sizeof("T"));
	imsl_ctran = imsl_l1ame(trans, sizeof(trans), "C", sizeof("C"));

	/*
	 * Test the input parameters.
	 */
	if (*n < 0) {
		imsl_e1psh("STBSV ");
		imsl_e1sti(1, *n);

/*		(5, 1, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("STBSV ");
		goto L_9000;
	} else if ((*ncoda < 0) && (*n > 0)) {
		imsl_e1psh("STBSV ");
		imsl_e1sti(1, *ncoda);

/*		(5, 2, "NCODA must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NCODA_GE_ZERO);
		imsl_e1pop("STBSV ");
		goto L_9000;
	} else if (*lda < *ncoda + 1) {
		imsl_e1psh("STBSV ");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *ncoda);

/*		(5, 3, "LDA must be greater than NCODA while LDA = %(i1) and NCODA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDA_GT_NCODA);
		imsl_e1pop("STBSV ");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("STBSV ");
		imsl_e1sti(1, *incx);

/*		(5, 4, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("STBSV ");
		goto L_9000;
	} else if (((!ntran) && (!tran)) && (!imsl_ctran)) {
		imsl_e1psh("STBSV ");
		imsl_e1stl(1, trans);

/*		(5, 5, "TRANS must be set equal to N or T or C while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TRANS_MUST_EQUAL_N_T_OR_C);
		imsl_e1pop("STBSV ");
		goto L_9000;
	} else if ((!upper) && (!lower)) {
		imsl_e1psh("STBSV ");
		imsl_e1stl(1, uplo);

/*		(5, 6, "UPLO must be set equal to U or L while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
		imsl_e1pop("STBSV ");
		goto L_9000;
	} else if ((!udiag) && (!ndiag)) {
		imsl_e1psh("STBSV ");
		imsl_e1stl(1, diag);

/*		(5, 7, "DIAG must be set equal to U or N while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DIAG_MUST_EQUAL_U_OR_N);
		imsl_e1pop("STBSV ");
		goto L_9000;
	}
	/* Quick return if possible. */
	if (*n == 0)
		goto L_9000;
	/* Start point in X */
	kx = 1;
	if (*incx <= 0)
		kx = 1 - (*n - 1) ** incx;

	if (ntran) {
		/* Form  x = inv( A )*x */
		if (upper) {
			kx += (*n - 1) ** incx;
			jx = kx;
			for (j = *n; j >= 1; j--) {
				kx -= *incx;
				if (x[jx - 1] != F_ZERO) {
					if (ndiag)
						x[jx - 1] /= *A(j - 1, *ncoda);
					m = imsl_i_max(*ncoda + 1 - j, 0);
					ix = kx + (*ncoda - m - 1) * imsl_i_min(-*incx, 0);
					l_saxpy(*ncoda - m, -x[jx - 1], A(j - 1, m), -1,
						   &x[ix - 1], -*incx);
				}
				jx -= *incx;
			}
		} else {
			jx = kx;
			for (j = 1; j <= *n; j++) {
				kx += *incx;
				if (x[jx - 1] != F_ZERO) {
					if (ndiag)
						x[jx - 1] /= *A(j - 1, 0);
					m = imsl_i_min(*ncoda, *n - j);
					ix = kx + (m - 1) * imsl_i_min(*incx, 0);
					l_saxpy(m, -x[jx - 1], A(j - 1, 1), 1, &x[ix - 1],
						   *incx);
				}
				jx += *incx;
			}
		}
	} else {
		/* Form  x = inv(trans(A))*x */
		if (upper) {
			jx = kx;
			for (j = 1; j <= *n; j++) {
				m = imsl_i_max(*ncoda + 1 - j, 0);
				ix = kx + (*ncoda - m - 1) * imsl_i_min(*incx, 0);
				x[jx - 1] -= imsl_sdot(*ncoda - m, A(j - 1, m), 1, &x[ix - 1],
						       *incx);
				if (ndiag)
					x[jx - 1] /= *A(j - 1, *ncoda);
				jx += *incx;
				if (j > *ncoda)
					kx += *incx;
			}
		} else {
			kx += (*n - 1) ** incx;
			jx = kx;
			for (j = *n; j >= 1; j--) {
				m = imsl_i_min(*ncoda, *n - j);
				ix = kx + (m - 1) * imsl_i_min(-*incx, 0);
				x[jx - 1] -= imsl_sdot(m, A(j - 1, 1), -1, &x[ix - 1],
						       -*incx);
				if (ndiag)
					x[jx - 1] /= *A(j - 1, 0);
				jx -= *incx;
				if ((*n - j) >= *ncoda)
					kx -= *incx;
			}
		}
	}
L_9000:
	return;
}				/* end of function */




/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 18:50:57
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SAXPY (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Compute the scalar times a vector plus a vector,
                y = ax + y, all single precision.

    Usage:      CALL SAXPY (N, SA, SX, INCX, SY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SA     - Real scalar.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be
                   SX(1+(I-1)*INCX) if INCX.GE.0  or
                   SX(1+(I-N)*INCX) if INCX.LT.0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).
                (Input/Output)
                SAXPY replaces Y(I) with SA*X(I) + Y(I) for I=1,...,N.
                X(I) and Y(I) refer to specific elements of SX and SY.
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be
                   SY(1+(I-1)*INCY) if INCY.GE.0  or
                   SY(1+(I-N)*INCY) if INCY.LT.0.

    Keywords:   Level 1 BLAS; SAXPY

    GAMS:       D1a7

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_saxpy(Mint n, Mfloat sa, Mfloat *sx, Mint incx, Mfloat *sy, Mint incy)
#else
static void l_saxpy(n, sa, sx, incx, sy, incy)
	Mint             n;
	Mfloat           sa, sx[];
	Mint             incx;
	Mfloat           sy[];
	Mint             incy;
#endif
{
	Mint             i, ix, iy;


	if (n > 0) {
		if (sa != F_ZERO) {
			if (incx != 1 || incy != 1) {
				/*
				 * CODE FOR UNEQUAL INCREMENTS OR EQUAL
				 * INCREMENTS NOT EQUAL TO 1
				 */
				ix = 1;
				iy = 1;
				if (incx < 0)
					ix = (-n + 1) * incx + 1;
				if (incy < 0)
					iy = (-n + 1) * incy + 1;
				for (i = 1; i <= n; i++) {
					sy[iy - 1] += sa * sx[ix - 1];
					ix += incx;
					iy += incy;
				}
			} else {
				for (i = 1; i <= n; i++) {
					sy[i - 1] += sa * sx[i - 1];
				}
			}
		}
	}
	return;
}				/* end of function */
