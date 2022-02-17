#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  L2LRB/DL2LRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Solve a real system of linear equations in band storage
                mode without iterative refinement.

    Usage:      CALL L2LRB (N, A, LDA, NLCA, NUCA, B, IPATH, X, FAC,
                            IPVT, WK)

    Arguments:  See LSLRB/DLSLRB.

    Remarks:    See LSLRB/DLSLRB.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_l2lrb(Mint *n, Mfloat *a, Mint *lda, Mint *nlca, Mint *nuca, Mfloat b[], Mint *ipath,
	   Mfloat x[], Mfloat imsl_fac[], Mint ipvt[], Mfloat wk[])
#else
void imsl_l2lrb(n, a, lda, nlca, nuca, b, ipath,
	   x, imsl_fac, ipvt, wk)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda, *nlca, *nuca;
	Mfloat           b[];
	Mint            *ipath;
	Mfloat           x[], imsl_fac[];
	Mint             ipvt[];
	Mfloat           wk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             ldfac, nra;
	Mfloat           rcond;


	imsl_e1psh("IMSL_L2LRB ");

	nra = *nlca + *nuca + 1;
	if (*n <= 0) {
		imsl_e1sti(1, *n);

/*		(5, 1, "The number of equations must be positive while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NUM_OF_EQUATIONS);
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

	} else if (nra > *lda) {
		imsl_e1sti(1, nra);
		imsl_e1sti(2, *lda);

/*		(5, 4, "The number of rows of matrix A must be less than or equal to its leading dimension while NLCA + NUCA + 1 = %(i1) and LDA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_ROWS_IN_MATRIX_A);

	} else if (*ipath != 1 && *ipath != 2) {
		imsl_e1sti(1, *ipath);

/*		(5, 5, "IPATH must be either 1 or 2 while a value of %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_IPATH_RANGE_3);
	} else {
		/*
		 * FACTOR AND ESTIMATE RECIPROCAL OF CONDITION NUMBER OF A
		 */
		ldfac = 2 ** nlca + *nuca + 1;
		imsl_l2crb(n, a, lda, nlca, nuca, imsl_fac, &ldfac, ipvt, &rcond, wk);
		/* SOLUTION STEP */
		if (imsl_n1rty(1) != 4) {
			imsl_lfsrb(n, imsl_fac, &ldfac, nlca, nuca, ipvt, b, ipath, x);

			if (rcond <= imsl_amach(4)) {
				imsl_e1str(1, rcond);

/*				(3, 1, "The matrix is too ill-conditioned.  An estimate of the reciprocal of its L1 condition number is RCOND = %(r1).  The solution might not be accurate.");
*/
                                imsl_ermes(IMSL_WARNING,
				IMSL_ILL_CONDITIONED);
			}
		}
	}

	imsl_e1pop("IMSL_L2LRB ");
	return;
}				/* end of function */
