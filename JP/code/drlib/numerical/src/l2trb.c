#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  L2TRB/DL2TRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the LU factorization of a real matrix in band
                storage mode.

    Usage:      CALL L2TRB (N, A, LDA, NLCA, NUCA, FAC, LDFAC, IPVT,
                            SCALE)

    Arguments:  See LFTRB/DLFTRB.

    Remarks:    See LFTRB/DLFTRB.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_l2trb(Mint *n, Mfloat *a, Mint *lda, Mint *nlca, Mint *nuca, 
                Mfloat *imsl_fac, Mint *ldfac, Mint ipvt[], Mfloat scale[])
#else
void imsl_l2trb(n, a, lda, nlca, nuca, imsl_fac, ldfac, ipvt, scale)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda, *nlca, *nuca;
	Mfloat          *imsl_fac;
	Mint            *ldfac, ipvt[];
	Mfloat           scale[];
#endif
{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(*ldfac)+(J_))
	Mint             i, i0, ia, indj, info, iscale, j, j0, j1, ja, ju,
	                jz, k, klen, l, ldb, lm, m, nrfac;
	Mfloat           big, curmax, small, t, value;


	imsl_e1psh("IMSL_L2TRB ");

	m = *nlca + *nuca + 1;
	nrfac = 2 ** nlca + *nuca + 1;
	if (*n <= 0) {
		imsl_e1sti(1, *n);

/*		(5, 1, "The order of the matrix must be positive while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	} else if (*nlca < 0 || *nlca >= *n) {
		imsl_e1sti(1, *nlca);

/*	        (5, 2, "The number of lower codiagonals must be greater than or equal to zero and less than N while NLCA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_NLCA_VALUE_SPECIFIED);

	} else if (*nuca < 0 || *nuca >= *n) {
		imsl_e1sti(1, *nuca);

/*		(5, 3, "The number of upper codiagonals must be greater than or equal to zero and less than N while NUCA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_NUCA_VALUE_SPECIFIED);

	} else if (m > *lda) {
		imsl_e1sti(1, m);
		imsl_e1sti(2, *lda);

/*		(5, 4, "The number of rows of matrix A must be less than or equal to its leading dimension while NUCA + NLCA + 1 = %(i1) and LDA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_ROWS_IN_MATRIX_A);

	} else if (nrfac > *ldfac) {
		imsl_e1sti(1, nrfac);
		imsl_e1sti(2, *ldfac);

/*		(5, 5, "The number of rows of matrix FAC must be less than or equal to its leading dimension while 2*NLCA + NUCA + 1 = %(i1) and LDFAC = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_ROWS_IN_FAC);
	}
	if (imsl_n1rcd(0) == 0) {
		/*
		 * MAKE A COPY OF A IN FAC TO WORK WITH AND PRESERVE THE
		 * INPUT A. LOOP IS IN REVERSE ORDER TO ALLOW FOR STORAGE
		 * SHARING BETWEEN A AND FAC IF NEEDED.
		 */
		info = 0;
		imsl_crbrb(n, a, lda, nlca, nuca, imsl_fac, ldfac, nlca, nuca);
		if (*nlca > 0) {
			for (j = m; j >= 1; j--) {
				scopy(*n, FAC(0, j - 1), *ldfac, FAC(0, *nlca + j - 1),
					   *ldfac);
				sset(*n, F_ZERO, FAC(0, j - 1), *ldfac);
			}
		}

		small = imsl_amach(1);
		big = imsl_amach(2);
		if (small * big < F_ONE)
			small = F_ONE / big;
		/*
		 * COMPUTE THE INFINITY NORM OF EACH ROW OF A FOR SCALING
		 * PURPOSE
		 */
		ldb = *ldfac - 1;
		for (i = 1; i <= *n; i++) {
			ia = imsl_i_min(*nuca + i, m) + *nlca;
			ja = imsl_i_max(1, i - *nlca);
			klen = imsl_i_min(i - 1, *nlca) + imsl_i_min(*n - i, *nuca) + 1;
			indj = imsl_isamax(klen, FAC(ja - 1, ia - 1), ldb);
			scale[i - 1] = fabs(*FAC(ja + indj - 2, ia - indj));
		}
		/* ZERO INITIAL FILL-IN COLUMNS */
		j0 = *nuca + 2;
		j1 = imsl_i_min(*n, m) - 1;
		for (jz = j0; jz <= j1; jz++) {
			i0 = m + 1 - jz;
			sset(*nlca - i0 + 1, F_ZERO, FAC(jz - 1, i0 - 1), 1);
		}
		jz = j1;
		ju = 0;
		/*
		 * GAUSSIAN ELIMINATION WITH SCALED PARTIAL PIVOTING
		 */
		for (k = 1; k <= (*n - 1); k++) {
			/* ZERO NEXT FILL-IN COLUMN */
			jz += 1;
			if (jz <= *n)
				sset(*nlca, F_ZERO, FAC(jz - 1, 0), 1);
			/* FIND L = PIVOT INDEX */
			lm = imsl_i_min(*nlca, *n - k);
			l = m;
			iscale = k;
			curmax = F_ZERO;
			for (i = m; i <= (m + lm); i++) {
				if (scale[iscale - 1] > small) {
					value = fabs(*FAC(k - 1, i - 1)) / scale[iscale - 1];
				} else {
					value = fabs(*FAC(k - 1, i - 1));
				}
				iscale += 1;
				if (value > curmax) {
					curmax = value;
					l = i;
				}
			}
			ipvt[k - 1] = l + k - m;
			/*
			 * ZERO PIVOT IMPLIES THIS COLUMN ALREADY
			 * TRIANGULARIZED
			 */
			if (fabs(*FAC(k - 1, l - 1)) > small) {
				/* INTERCHANGE IF NECESSARY */
				if (l != m) {
					t = *FAC(k - 1, l - 1);
					*FAC(k - 1, l - 1) = *FAC(k - 1, m - 1);
					*FAC(k - 1, m - 1) = t;
				}

				/* COMPUTE MULTIPLIERS */
				t = -F_ONE / *FAC(k - 1, m - 1);
				if (lm != 0)
					sscal(lm, t, FAC(k - 1, m), 1);
				/* ROW ELIMINATION WITH COLUMN INDEXING */
				ju = imsl_i_min(imsl_i_max(ju, *nuca + ipvt[k - 1]), *n);
				if (m > 1 || k < ju)
					sswap(ju - k, FAC(k, m - 2), *ldfac - 1, FAC(k, l - 2),
						   *ldfac - 1);
				if (lm != 0 || k < ju)
					imsl_sger(lm, ju - k, F_ONE, FAC(k - 1, m), 1, FAC(k, m - 2),
						  *ldfac - 1, FAC(k, m - 1), *ldfac - 1);
			} else {
				info = k;

			}
		}
		ipvt[*n - 1] = *n;

		if (fabs(*FAC(*n - 1, m - 1)) <= small)
			info = *n;

		if (info != 0) {

/*		        (4, 2, "The input matrix is singular.  Some of the diagonal elements of the upper triangular matrix U of the LU factorization are equal to zero.");
*/
                        imsl_ermes(IMSL_FATAL, IMSL_SINGULAR_MATRIX);
		}
	}
	imsl_e1pop("IMSL_L2TRB ");
	return;
}				/* end of function */
