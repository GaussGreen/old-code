#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*  -----------------------------------------------------------------------
    IMSL Name:  CRBRB/DCRBRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 12, 1986

    Purpose:    Copy a real band matrix stored in band storage mode.

    Usage:      CALL CRBRB (N, A, LDA, NLCA, NUCA, B, LDB, NLCB, NUCB)

    Arguments:
       N      - Order of the matrices A and B.  (Input)
       A      - Real band matrix of order N.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       NLCA   - Number of lower codiagonals in A.  (Input)
       NUCA   - Number of upper codiagonals in A.  (Input)
       B      - Real band matrix of order N containing a copy of A.
                (Output)
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)
       NLCB   - Number of lower codiagonals in B.  (Input)
                NLCB must be at least as large as NLCA.
       NUCB   - Number of upper codiagonals in B.  (Input)
                NUCB must be at least as large as NUCA.

    Keyword:    Basic matrix operation

    GAMS:       D1b8

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_crbrb(Mint *n, Mfloat *a, Mint *lda, Mint *nlca, Mint *nuca,
                Mfloat *b, Mint *ldb, Mint *nlcb, Mint *nucb)
#else
void imsl_crbrb(n, a, lda, nlca, nuca, b, ldb, nlcb, nucb)
        Mint            *n;
        Mfloat          *a;
        Mint            *lda, *nlca, *nuca;
        Mfloat          *b;
        Mint            *ldb, *nlcb, *nucb;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define B(I_,J_)	(b+(I_)*(*ldb)+(J_))
	Mint             j, k, ml, mu, ndiaga;


	imsl_e1psh("IMSL_CRBRB ");
	ndiaga = *nuca + *nlca + 1;
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  It must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
		goto L_9000;
	}
	/* Check LDA */
	if (*lda < ndiaga) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *nlca);
		imsl_e1sti(3, *nuca);
		imsl_e1sti(4, ndiaga);

/*		imsl_ermes(5, 2, "The argument LDA = %(i1).  It must be at least as large as NLCA+NUCA+1 = %(i2)+%(i3)+1 = %(i4).");
*/
                imsl_ermes(IMSL_TERMINAL,
		IMSL_LEADING_DIM_OF_A_TOO_SMALL);
		goto L_9000;
	}
	/* Check NLCA */
	if (*nlca < 0 || *nlca >= *n) {
		imsl_e1sti(1, *nlca);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The argument NLCA = %(i1).  It must be at least 0 and less than N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ZERO_LE_NLCA_LT_N);
		goto L_9000;
	}
	/* Check NUCA */
	if (*nuca < 0 || *nuca >= *n) {
		imsl_e1sti(1, *nuca);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 4, "The argument NUCA = %(i1).  It must be at least 0 and less than N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ZERO_LE_NUCA_LT_N);
		goto L_9000;
	}
	/* Check LDB */
	if (*ldb < *nlcb + *nucb + 1) {
		imsl_e1sti(1, *ldb);
		imsl_e1sti(2, *nlcb);
		imsl_e1sti(3, *nucb);
		imsl_e1sti(4, *nlcb + *nucb + 1);

/*		imsl_ermes(5, 5, "The argument LDB = %(i1).  It must be at least as large as NLCB+NUCB+1 = %(i2)+%(i3)+1 = %(i4).");
*/
                imsl_ermes(IMSL_TERMINAL,
		IMSL_LEADING_DIM_OF_B_TOO_SMALL);
		goto L_9000;
	}
	/* Check NLCB */
	if (*nlcb < *nlca || *nlcb >= *n) {
		imsl_e1sti(1, *nlcb);
		imsl_e1sti(2, *nlca);
		imsl_e1sti(3, *n);

/*		imsl_ermes(5, 6, "The argument NLCB = %(i1).  It must be at least as large as NLCA = %(i2) and less than N = %(i3).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NLCA_LE_NLCB_LT_N);
		goto L_9000;
	}
	/* Check NUCB */
	if (*nucb < *nuca || *nucb >= *n) {
		imsl_e1sti(1, *nucb);
		imsl_e1sti(2, *nuca);
		imsl_e1sti(3, *n);

/*		imsl_ermes(5, 7, "The argument NUCB = %(i1).  It must be at least as large as NUCA = %(i2) and less than N = %(i3).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NUCA_LE_NUCB_LT_N);
		goto L_9000;
	}
	if (*lda > *ldb) {
		/* Copy diagonals */
		for (j = 1; j <= *n; j++) {
			mu = imsl_i_min(j - 1, *nuca);
			ml = imsl_i_min(*n - j, *nlca);
			k = *nuca - mu + 1;
			imsl_scopy(mu + ml + 1, A(j - 1, k - 1), 1, B(j - 1, k - 1),
				   1);
			if (*nuca - mu != 0) {
				sset(*nuca - mu, F_ZERO, B(j - 1, 0), 1);
			}
			if (*nlca - ml != 0) {
				sset(*nlca - ml, F_ZERO, B(j - 1, *nuca + ml + 1), 1);
			}
		}
	} else {
		/* Copy backwards */
		for (j = *n; j >= 1; j--) {
			mu = imsl_i_min(j - 1, *nuca);
			ml = imsl_i_min(*n - j, *nlca);
			k = *nuca - mu + 1;
			imsl_scopy(mu + ml + 1, A(j - 1, k - 1), -1, B(j - 1, k - 1),
				   -1);
			if (*nuca - mu != 0) {
				sset(*nuca - mu, F_ZERO, B(j - 1, 0), 1);
			}
			if (*nlca - ml != 0) {
				sset(*nlca - ml, F_ZERO, B(j - 1, *nuca + ml + 1), 1);
			}
		}
	}
	/* Move diagonals into place */
	if (*nucb > *nuca) {
		for (j = ndiaga; j >= 1; j--) {
			imsl_scopy(*n, B(0, j - 1), *ldb, B(0, *nucb - *nuca + j - 1),
				   *ldb);
		}
	}
	/* Zero out extra upper B codiagonals */
	for (j = 1; j <= (*nucb - *nuca); j++) {
		sset(*n, F_ZERO, B(0, j - 1), *ldb);
	}
	/* Zero out extra lower B codiagonals */
	for (j = 1; j <= (*nlcb - *nlca); j++) {
		sset(*n, F_ZERO, B(0, *nucb - *nuca + ndiaga + j - 1), *ldb);
	}

L_9000:
	;
	imsl_e1pop("IMSL_CRBRB ");
	return;
}				/* end of function */
