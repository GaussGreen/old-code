#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/31/90 at 18:30:55
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  TRNRR/DTRNRR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 8, 1985

    Purpose:    Transpose a rectangular matrix.

    Usage:      CALL TRNRR (NRA, NCA, A, LDA, NRB, NCB, B, LDB)

    Arguments:
       NRA    - Number of rows of A.  (Input)
       NCA    - Number of columns of A.  (Input)
       A      - Real NRA by NCA matrix in full storage mode.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       NRB    - Number of rows of B.  (Input)
                NRB must be equal to NCA.
       NCB    - Number of columns of B.  (Input)
                NCB must be equal to NRA.
       B      - Real NRB by NCB matrix in full storage mode containing
                the transpose of A.  (Output)
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       If LDA=LDB and NRA=NCA then A and B may be the same matrix,
       otherwise A and B must be different matrices.

    Keyword:    Basic matrix operation

    GAMS:       D1b3

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_trnrr(Mint nra, Mint nca, Mfloat *a, Mint lda, Mint nrb, Mint ncb,
                Mfloat *b, Mint ldb)
#else
void imsl_trnrr(nra, nca, a, lda, nrb, ncb, b, ldb)
	Mint             nra, nca;
	Mfloat          *a;
	Mint             lda, nrb, ncb;
	Mfloat          *b;
	Mint             ldb;
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
#define B(I_,J_)	(b+(I_)*(ldb)+(J_))
	Mint             i;


	imsl_e1psh("TRNRR ");

	if (nra <= 0 || nca <= 0) {
		imsl_e1sti(1, nra);
		imsl_e1sti(2, nca);

/*		(5, 1, "Both the number of rows and the number of columns of a matrix have to be positive while NRA = %(i1) and NCA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NRA_AND_NCA_GT_ZERO);
	}
	if (nra > lda) {
		imsl_e1sti(1, nra);
		imsl_e1sti(2, lda);

/*		(5, 2, "The number of rows of the matrix must be less than or equal to its leading dimension while NRA = %(i1) and LDA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NRA_EXCEEDS_LDA_VALUE);
	}
	if (nrb <= 0 || ncb <= 0) {
		imsl_e1sti(1, nrb);
		imsl_e1sti(2, ncb);

/*		(5, 3, "Both the number of rows and the number of columns of a matrix have to be positive while NRB = %(i1) and NCB = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_POSITIVE_NRB_NCB_VALUES);
	}
	if (nrb > ldb) {
		imsl_e1sti(1, nrb);
		imsl_e1sti(2, ldb);

/*		(5, 4, "The number of rows of the matrix must be less than or equal to its leading dimension while NRB = %(i1) and LDB = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NRB_EXCEEDS_LDB_VALUE);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	if (nrb != nca || ncb != nra) {
		imsl_e1sti(1, nra);
		imsl_e1sti(2, nca);
		imsl_e1sti(3, nrb);
		imsl_e1sti(4, ncb);

/*		(5, 5, "The following must hold NRB = NCA and NCB = NRA while NRB = %(i3), NCA = %(i2), NCB = %(i4) and NRA = %(i1).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_DIMENSIONS_FOR_TRANS);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * If LDA=LDB and NRA=NCA then A and B may be the same, copy A to B
	 * and transpose B with SSWAP.
	 */
	if (nra == nca && lda == ldb) {
		for (i = 1; i <= nca; i++) {
			scopy(nra, A(i - 1, 0), 1, B(i - 1, 0), 1);
		}
		for (i = 1; i <= (nca - 1); i++) {
			sswap(nra - i, B(i - 1, i), 1, B(i, i - 1), ldb);
		}
	} else {
		for (i = 1; i <= nca; i++) {
			scopy(nra, A(i - 1, 0), 1, B(0, i - 1), ldb);
		}
	}
L_9000:
	imsl_e1pop("TRNRR ");
	return;
}				/* end of function */
