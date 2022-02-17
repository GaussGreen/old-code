#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  NR1RB/DNR1RB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 18, 1985

    Purpose:    Compute the 1-norm of a real band matrix in band storage
                mode.

    Usage:      CALL NR1RB (N, A, LDA, NLCA, NUCA, ANORM)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Real (NUCA + NLCA + 1) by N array containing the N by N
                band matrix in band storage mode.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       NLCA   - Number of lower codiagonals of A.  (Input)
       NUCA   - Number of upper codiagonals of A.  (Input)
       ANORM  - Real scalar containing the 1-norm of A.  (Output)

    Keyword:    Basic matrix operation

    GAMS:       D1b2

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_nr1rb(Mint *n, Mfloat *a, Mint *lda, Mint *nlca, Mint *nuca, Mfloat *anorm)
#else
void imsl_nr1rb(n, a, lda, nlca, nuca, anorm)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda, *nlca, *nuca;
	Mfloat          *anorm;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             is, j, l, nra;
	Mfloat           csum;


	imsl_e1psh("IMSL_NR1RB ");

	if (*n <= 0) {
		imsl_e1sti(1, *n);

/*		(5, 1, "The order of the matrix must be positive while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
		goto L_9000;
	}
	nra = *nlca + *nuca + 1;
	if (nra > *lda) {
		imsl_e1sti(1, nra);
		imsl_e1sti(2, *lda);

/*		(5, 2, "The number of upper codiagonals plus the number of lower codiagonals plus 1 must be less than or equal to its leading dimension while NUCA+NLCA+1 = %(i1) and LDA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_CODIAGONALS);
		goto L_9000;
	}
	if (*nlca < 0 || *nlca >= *n) {
		imsl_e1sti(1, *nlca);
		imsl_e1sti(2, *n);

/*		(5, 3, "The number of lower codiagonals must be greater than or equal to zero and less than N while NLCA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LOWER_CODIAGONALS);
	}
	if (*nuca < 0 || *nuca >= *n) {
		imsl_e1sti(1, *nuca);
		imsl_e1sti(2, *n);

/*		(5, 4, "The number of upper codiagonals must be greater than or equal to zero and less than N while NUCA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_UPPER_CODIAGONALS);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	*anorm = F_ZERO;
	l = *nlca + 1;
	is = *nuca + 1;
	for (j = 1; j <= *n; j++) {
		csum = imsl_sasum(l, A(j - 1, is - 1), 1);
		*anorm = imsl_f_max(*anorm, csum);
		if (is > 1)
			is -= 1;
		if (j <= *nuca)
			l += 1;
		if (j >= (*n - *nlca))
			l -= 1;
	}

L_9000:
	imsl_e1pop("IMSL_NR1RB ");
	return;
}				/* end of function */
