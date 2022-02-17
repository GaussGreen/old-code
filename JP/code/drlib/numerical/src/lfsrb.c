#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  LFSRB/DLFSRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Solve a real system of linear equations given the LU
                factorization of the coefficient matrix in band storage
                mode.

    Usage:      CALL LFSRB (N, FAC, LDFAC, NLCA, NUCA, IPVT, B, IPATH, X)

    Arguments:
       N      - Number of equations.  (Input)
       FAC    - (2*NLCA+NUCA+1) by N array containing the LU
                factorization of the coefficient matrix A as output from
                subroutine LFCRB/DLFCRB or LFTRB/DLFTRB.  (Input)
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)
       NLCA   - Number of lower codiagonals of A.  (Input)
       NUCA   - Number of upper codiagonals of A.  (Input)
       IPVT   - Vector of length N containing the pivoting information
                for the LU factorization of A as output from subroutine
                LFCRB/DLFCRB or LFTRB/DLFTRB.  (Input)
       B      - Vector of length N containing the right-hand side of the
                linear system.  (Input)
       IPATH  - Path indicator.  (Input)
                IPATH = 1 means the system A*X = B is solved.
                IPATH = 2 means the system trans(A)*X = B is solved,
                          where trans(A) is the transpose of A.
       X      - Vector of length N containing the solution to the linear
                system.  (Output)
                If B is not needed, B and X can share the same storage
                locations.

    GAMS:       D2a2

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_lfsrb(Mint *n, Mfloat *imsl_fac, Mint *ldfac, Mint *nlca,
	   Mint *nuca, Mint ipvt[], Mfloat b[], Mint *ipath, Mfloat x[])
#else
void imsl_lfsrb(n, imsl_fac, ldfac, nlca, nuca, ipvt, b, ipath,
	   x)
	Mint            *n;
	Mfloat          *imsl_fac;
	Mint            *ldfac, *nlca, *nuca, ipvt[];
	Mfloat           b[];
	Mint            *ipath;
	Mfloat           x[];
#endif
{
#define FAC(I_,J_)	(imsl_fac+(I_)*(*ldfac)+(J_))
	Mint             _l0, _l1, k, l, lm, m, nrfac;
	Mfloat           big, small, t;


	imsl_e1psh("IMSL_LFSRB ");

	nrfac = 2 ** nlca + *nuca + 1;
	if (*n <= 0) {
		imsl_e1sti(1, *n);

/*		(5, 1, "The order of the matrix must be positive while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	} else if (*nlca < 0 || *nlca >= *n) {
		imsl_e1sti(1, *nlca);

/*		(5, 2, "The number of lower codiagonals must be greater than or equal to zero and less than N while NLCA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_NLCA_VALUE_SPECIFIED);

	} else if (*nuca < 0 || *nuca >= *n) {
		imsl_e1sti(1, *nuca);

/*		(5, 3, "The number of upper codiagonals must be greater than or equal to zero and less than N while NUCA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_NUCA_VALUE_SPECIFIED);

	} else if (nrfac > *ldfac) {
		imsl_e1sti(1, nrfac);
		imsl_e1sti(2, *ldfac);

/*		(5, 4, "The number of rows of matrix FAC must be less than or equal to its leading dimension while 2*NLCA + NUCA + 1 = %(i1) and LDFAC = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_ROWS_IN_FAC);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	/*
	 * COPY B INTO X AND USE X TO PRESERVE INPUT
	 */
	scopy(*n, b, 1, x, 1);

	m = *nuca + *nlca + 1;
	small = imsl_amach(1);
	big = imsl_amach(2);
	if (small * big < F_ONE)
		small = F_ONE / big;
	if (*ipath == 1) {
		/*
		 * IPATH = 1 , SOLVE A * X = B FIRST SOLVE L*Y = B
		 */
		if (*nlca != 0) {
			for (k = 1; k <= (*n - 1); k++) {
				lm = imsl_i_min(*nlca, *n - k);
				l = ipvt[k - 1];
				t = x[l - 1];
				if (l != k) {
					x[l - 1] = x[k - 1];
					x[k - 1] = t;
				}
				saxpy(lm, t, FAC(k - 1, m), 1, &x[k], 1);
			}
		}
		/* NOW SOLVE U*X = Y */
		for (k = *n; k >= 1; k--) {
			if (fabs(*FAC(k - 1, m - 1)) <= small) {

/*				(5, 6, "The input matrix is singular.  Some of the diagonal elements of the upper triangular matrix U of the LU factorization are close to zero.");
*/
                               imsl_ermes(IMSL_FATAL,
                               IMSL_SINGULAR_MATRIX);
				goto L_9000;
			}
		}
                _l0 =  *nuca + *nlca;
                _l1 = 1;
		imsl_stbsv("U", "N", "N", 
		n, &_l0, imsl_fac, ldfac, x, &_l1);

	} else if (*ipath == 2) {
		/*
		 * IPATH = 2, SOLVE TRANS(A) * X = B FIRST SOLVE TRANS(U)*Y =
		 * B
		 */
		for (k = 1; k <= *n; k++) {
			if (fabs(*FAC(k - 1, m - 1)) <= small) {

/*				(5, 7, "The input matrix is singular.  Some of the diagonal elements of the upper triangular matrix U of the LU factorization are close to zero.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_SINGULAR_MATRIX);
				goto L_9000;
			}
		}
                _l0 =  *nuca + *nlca;
                _l1 = 1;
		imsl_stbsv("U", "T",  "N",
		n, &_l0, imsl_fac, ldfac, x, &_l1);
		/* NOW SOLVE TRANS(L)*X = Y */
		if (*nlca != 0) {
			for (k = *n - 1; k >= 1; k--) {
				lm = imsl_i_min(*nlca, *n - k);
				x[k - 1] += imsl_sdot(lm, FAC(k - 1, m), 1, &x[k], 1);
				l = ipvt[k - 1];
				if (l != k) {
					t = x[l - 1];
					x[l - 1] = x[k - 1];
					x[k - 1] = t;
				}
			}
		}
	} else {
		imsl_e1sti(1, *ipath);

/*		(5, 5, "IPATH must be either 1 or 2 while a value of %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_IPATH_RANGE_3);
	}

L_9000:
	imsl_e1pop("IMSL_LFSRB ");
	return;
}				/* end of function */
