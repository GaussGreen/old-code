#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  L2CRB/DL2CRB (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the LU factorization of a real matrix in band
                storage mode and estimate its L1 condition number.

    Usage:      CALL L2CRB (N, A, LDA, NLCA, NUCA, FAC, LDFAC, IPVT,
                            RCOND, Z)

    Arguments:  See LFCRB/DLFCRB.

    Remarks:    See LFCRB/DLFCRB.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_l2crb(Mint *n, Mfloat *a, Mint *lda, Mint *nlca, Mint *nuca, Mfloat *imsl_fac,
	        Mint *ldfac, Mint ipvt[], Mfloat *rcond, Mfloat z[])
#else
void imsl_l2crb(n, a, lda, nlca, nuca, imsl_fac, ldfac, ipvt, rcond,
	   z)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda, *nlca, *nuca;
	Mfloat          *imsl_fac;
	Mint            *ldfac, ipvt[];
	Mfloat          *rcond, z[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(*ldfac)+(J_))
	Mint             j, ju, k, l, la, lm, lz, m, mm, nra,
	                nrfac;
	Mfloat           anorm, big, ek, s, sm, small, t,
	                wk, wkm, ynorm;


	imsl_e1psh("IMSL_L2CRB ");

	nra = *nlca + *nuca + 1;
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
	} else if (nra > *lda) {
		imsl_e1sti(1, nra);
		imsl_e1sti(2, *lda);

/*		(5, 4, "The number of rows of matrix A must be less than or equal to its leading dimension while NLCA + NUCA + 1 = %(i1) and LDA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_ROWS_IN_MATRIX_A);

	} else if (nrfac > *ldfac) {
		imsl_e1sti(1, nrfac);
		imsl_e1sti(2, *ldfac);

/*		(5, 5, "The number of rows of matrix FAC must be less than or equal to its leading dimension while 2*NLCA + NUCA + 1 = %(i1) and LDFAC = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_ROWS_IN_FAC);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	*rcond = F_ZERO;
	/* COMPUTE 1-NORM OF A */
	imsl_nr1rb(n, a, lda, nlca, nuca, &anorm);
	/* FACTORIZATION STEP */
	imsl_l2trb(n, a, lda, nlca, nuca, imsl_fac, ldfac, ipvt, z);
	if (imsl_n1rty(1) == 4)
		goto L_9000;
	/*
	 * RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . ESTIMATE =
	 * NORM(Z)/NORM(Y) WHERE A*Z = Y AND TRANS(A)*Y = E . TRANS(A) IS THE
	 * TRANSPOSE OF A . THE COMPONENTS OF E ARE CHOSEN TO CAUSE MAXIMUM
	 * LOCAL GROWTH IN THE ELEMENTS OF W WHERE TRANS(U)*W = E . THE
	 * VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. SOLVE
	 * TRANS(U)*W = E
	 */
	ek = F_ONE;
	sset(*n, F_ZERO, z, 1);
	m = *nlca + *nuca + 1;
	ju = 0;
	small = imsl_amach(1);
	big = imsl_amach(2);
	if (small * big < F_ONE)
		small = F_ONE / big;

	for (k = 1; k <= *n; k++) {
		if (z[k - 1] != F_ZERO)
			ek = sign(ek, -z[k - 1]);
		if (fabs(ek - z[k - 1]) > fabs(*FAC(k - 1, m - 1))) {
			s = fabs(*FAC(k - 1, m - 1)) / fabs(ek - z[k - 1]);
			sscal(*n, s, z, 1);
			ek *= s;
		}
		wk = ek - z[k - 1];
		wkm = -ek - z[k - 1];
		s = fabs(wk);
		sm = fabs(wkm);
		if (fabs(*FAC(k - 1, m - 1)) > small) {
			wk /= *FAC(k - 1, m - 1);
			wkm /= *FAC(k - 1, m - 1);
		} else {
			wk = F_ONE;
			wkm = F_ONE;
		}
		ju = imsl_i_min(imsl_i_max(ju, *nuca + ipvt[k - 1]), *n);
		mm = m;
		if (k + 1 <= ju) {
			for (j = k + 1; j <= ju; j++) {
				mm -= 1;
				sm += fabs(z[j - 1] + wkm ** FAC(j - 1, mm - 1));
			}
			if (m > 1)
				saxpy(ju - k, wk, FAC(k, m - 2), *ldfac - 1, &z[k],
					   1);
			s += imsl_sasum( ju - k, &z[k],  1);

			if (s < sm) {
				t = wkm - wk;
				wk = wkm;
				if (m > 1)
					saxpy(ju - k, t, FAC(k, m - 2), *ldfac - 1, &z[k],
						   1);
			}
		}
		z[k - 1] = wk;
	}

	s = F_ONE / imsl_sasum(*n, z, 1);
	sscal(*n, s, z, 1);
	/* SOLVE TRANS(L)*Y = W */
	for (k = *n; k >= 1; k--) {
		lm = imsl_i_min(*nlca, *n - k);
		if (k < *n && lm > 0)
			z[k - 1] += imsl_sdot(lm, FAC(k - 1, m), 1, &z[k], 1);
		if (fabs(z[k - 1]) > F_ONE) {
			s = F_ONE / fabs(z[k - 1]);
			sscal(*n, s, z, 1);
		}
		l = ipvt[k - 1];
		t = z[l - 1];
		z[l - 1] = z[k - 1];
		z[k - 1] = t;
	}
	s = F_ONE / imsl_sasum(*n, z, 1);
	sscal(*n, s, z, 1);
	ynorm = F_ONE;
	/* SOLVE L*V = Y */
	for (k = 1; k <= *n; k++) {
		l = ipvt[k - 1];
		t = z[l - 1];
		z[l - 1] = z[k - 1];
		z[k - 1] = t;
		lm = imsl_i_min(*nlca, *n - k);
		if (k < *n && lm > 0)
			saxpy(lm, t, FAC(k - 1, m), 1, &z[k], 1);
		if (fabs(z[k - 1]) > F_ONE) {
			s = F_ONE / fabs(z[k - 1]);
			sscal(*n, s, z, 1);
			ynorm *= s;
		}
	}
	s = F_ONE / imsl_sasum(*n, z, 1);
	sscal(*n, s, z, 1);
	ynorm *= s;
	/* SOLVE U*Z = V */
	for (k = *n; k >= 1; k--) {
		if (fabs(z[k - 1]) > fabs(*FAC(k - 1, m - 1))) {
			s = fabs(*FAC(k - 1, m - 1)) / fabs(z[k - 1]);
			sscal(*n, s, z, 1);
			ynorm *= s;
		}
		if (fabs(*FAC(k - 1, m - 1)) > small) {
			z[k - 1] /= *FAC(k - 1, m - 1);
		} else {
			z[k - 1] = F_ONE;
		}

		lm = imsl_i_min(k, m) - 1;
		la = m - lm;
		lz = k - lm;
		t = -z[k - 1];
		saxpy(lm, t, FAC(k - 1, la - 1), 1, &z[lz - 1], 1);
	}
	/* MAKE ZNORM = 1.0E0 */
	s = F_ONE / imsl_sasum(*n, z, 1);
	sscal(*n, s, z, 1);
	ynorm *= s;

	if (anorm > small)
		*rcond = ynorm / anorm;
	if (*rcond <= imsl_amach(4)) {
		imsl_e1str(1, *rcond);

/*		(3, 1, "The matrix is algorithmically singular.  An estimate of the reciprocal of its L1 condition number is RCOND = %(r1).");
*/
                imsl_ermes(IMSL_WARNING, IMSL_ILL_CONDITIONED);
	}
L_9000:
	imsl_e1pop("IMSL_L2CRB ");
	return;
}				/* end of function */



