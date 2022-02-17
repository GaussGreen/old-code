#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  G2RCF/DG2RCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 1, 1985

    Purpose:    Compute a Gauss, Gauss-Radau or Gauss-Lobatto quadrature
                rule given the recurrence coefficients for the monic
                polynomials orthogonal with respect to the weight
                function.

    Usage:      CALL G2RCF (N, B, C, NFIX, QXFIX, QX, QW, WK)

    Arguments:  (See GQRCF)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_g2rcf(Mint *n, Mfloat b[], Mfloat c[], Mint *nfix, Mfloat qxfix[],
                Mfloat qx[], Mfloat qw[], Mfloat wk[])
#else
void imsl_g2rcf(n, b, c, nfix, qxfix, qx, qw, wk)
	Mint            *n;
	Mfloat           b[], c[];
	Mint            *nfix;
	Mfloat           qxfix[], qx[], qw[], wk[];
#endif
{
	Mint             i;
	Mfloat           gam, t1;


	imsl_e1psh("G2RCF ");
	/* CHECK N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		(5, 1, "The number of quadrature points N = %(i1).  It must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NUM_QUADRATURE_POINTS);
		goto L_9000;
	}
	/* CHECK NFIX */
	if (*nfix < 0 || *nfix > 2) {
		imsl_e1sti(1, *nfix);

/*		(5, 3, "The number of fixed points NFIX = %(i1).  It must be 0, 1 or 2.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_NFIX_VALUE);
		goto L_9000;
	}
	/* N = 1 AND NFIX = 2 */
	if (*n <= *nfix) {
		imsl_e1sti(1, *n);
		imsl_e1sti(2, *nfix);

/*		(5, 4, "The number of quadrature points N = %(i1) and the number of fixed points NFIX = %(i2).  N must be greater than NFIX.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_MORE_QUAD_POINTS);
		goto L_9000;
	}
	/* CHECK THAT C(I) .GT. 0 FOR ALL I */
	for (i = 1; i <= *n; i++) {
		if (c[i - 1] <= F_ZERO) {
			imsl_e1sti(1, i);
			imsl_e1str(1, c[i - 1]);

/*		(5, 5, "The recurrence coefficient C(%(i1)) = %(r1).  Each element of the vector C must be greater than zero. ");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_RECURRENCE_COEFF_C);
			goto L_9000;
		}
	}
	/* COPY B ONTO QX, C ONTO QW */
	scopy(*n, b, 1, qx, 1);
	scopy(*n, c, 1, qw, 1);

	if (*nfix == 1) {
		qx[*n - 1] = imsl_g3rcf(&qxfix[0], n, qx, qw) * qw[*n - 1] + qxfix[0];
	} else if (*nfix == 2) {
		gam = imsl_g3rcf(&qxfix[0], n, qx, qw);
		t1 = (qxfix[0] - qxfix[1]) / (imsl_g3rcf(&qxfix[1], n, qx, qw) -
					      gam);
		qw[*n - 1] = t1;
		qx[*n - 1] = qxfix[0] + gam * t1;
	}
	/* CHECK C .GT. 0 WHEN NFIX IS 2 */
	if (*nfix == 2) {
		for (i = 1; i <= *n; i++) {
			if (qw[i - 1] <= F_ZERO) {
				imsl_e1sti(1, *n - 2);
				imsl_e1str(1, qxfix[0]);

/*				(5, 6, "Let I be the smallest interval containing the zeros of the Kth orthogonal polynomial, where K = N-2 = %(i1).  Either one of the fixed points in the vector QXFIX lies in I, or both lie on one side of I.  They must be on opposite sides of I.");
*/
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_FIXED_POINT_LOCATION);
				goto L_9000;
			}
		}
	}
	imsl_g4rcf(n, qx, qw, qx, qw, wk);

L_9000:
	;
	imsl_e1pop("G2RCF ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  G3RCF/DG3RCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 1, 1985

    Purpose:    Solve an equation to obtain the appropriate changes in
                the lower 2 by 2 submatrix of coefficients for
                orthogonal polynomials as used by G2RUL.

    Usage:      G3RCF(SHIFT, N, B, C)

    Arguments:
       SHIFT  - Constant in equation.  (Input)
                See remark.
       N      - Order of the tridiagonal matrix.  (Input)
       B      - Array of length N containing recurrence coefficients.
                (Input)
       C      - Array of length N containing recurrence coefficients.
                (Input)
       G3RCF  - N-th component of the solution DELTA.  (Output)
                See Remark.

    Remark:
       G3RCF performs elimination to solve for the N-th component
       of the solution DELTA to the equation
                   (JN - SHIFT*IDENTITY) * DELTA = EN,
       where EN is the vector of all zeros except for 1 in the N-th
       JN is the symmetric tridiagonal matrix with diagonal elements B
       and off-diagonal elements C.

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_g3rcf(Mfloat *shift, Mint *n, Mfloat b[], Mfloat c[])
#else
Mfloat imsl_g3rcf(shift,n,b,  c)
	Mfloat          *shift;
	Mint            *n;
	Mfloat           b[], c[];
#endif
{
	Mint             i;
	Mfloat           alpha, g3rcf_v;


	alpha = b[0] - *shift;
	for (i = 2; i <= (*n - 1); i++) {
		alpha = b[i - 1] - *shift - c[i - 1] / alpha;
	}

	g3rcf_v = F_ONE / alpha;

	return (g3rcf_v);
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  G4RCF/DG4RCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 1, 1985

    Purpose:    Derive the abscissae and weights of a Gauss-
                Christoffel quadrature formula given the recurrence
                coefficients for the polynomials orthogonal with
                respect to a user-supplied weight.

    Usage:      CALL G4RCF (N, B, C, X, W, CTMP)

    Arguments:
       N      - Number of abscissae and weights to be found.  (Input)
                The highest order polynomial defined above is N+1.
       B      - Array of length N containing recurrence coefficients.
                (Input)
       C      - Array of length N containing recurrence coefficients.
                (Input)
       X      - Array of length N containing Gauss-Christoffel abscissae.
                (Output)
                These are the N zeros of P(N+1).
       W      - Array of length N containing Gauss-Christoffel weights.
                (Output)
       CTMP   - Work array of length N.

    Remarks:
    1. The arrays X and B need not be distinct.  Likewise with W and C.

    2. The orthogonal polynomials are given by
       P(-1)=0, P(0)=1 AND P(M) = (X-B(M))*P(M-1) - C(M)*P(M-2).

    3. The abscissae and weights are the eigenvalues and first
       components of the eigenvectors of a symmetric tridiagonal
       matrix, which is solved using the implicit QL method.  This
       subroutine is based on the EISPACK routine IMTQL2 as modified
       by Gene Golub and further modified by Wayne Fullerton.

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_g4rcf(Mint *n, Mfloat b[], Mfloat c[], Mfloat x[], Mfloat w[], Mfloat ctmp[])
#else
void imsl_g4rcf(n,  b, c, x, w, ctmp)
	Mint            *n;
	Mfloat           b[], c[], x[], w[], ctmp[];
#endif
{
	Mint             i, ii, iter, j, k, m, mmk;
	Mfloat           bb, cc, eps, f, g, p, r, s, wtsum;


	eps = imsl_amach(4);

	wtsum = c[0];
	for (i = 1; i <= *n; i++) {
		if (i > 1)
			ctmp[i - 2] = sqrt(c[i - 1]);
		x[i - 1] = b[i - 1];
		w[i - 1] = F_ZERO;
	}
	w[0] = F_ONE;
	ctmp[*n - 1] = F_ZERO;

	for (k = 1; k <= *n; k++) {
		for (iter = 1; iter <= 100; iter++) {
			/* LOOK FOR SMALL SUB-DIAGONAL ELEMENT */
			for (m = k; m <= *n; m++) {
				if (m == *n)
					goto L_30;
				if (fabs(ctmp[m - 1]) <= eps * (fabs(x[m - 1]) +
								fabs(x[m])))
					goto L_30;
			}

	L_30:
			p = x[k - 1];
			if (m == k)
				goto L_80;
			/* FORM SHIFT */
			g = (x[k] - p) / (F_TWO * ctmp[k - 1]);
			r = sqrt(g * g + F_ONE);
			s = r;
			if (g < F_ZERO)
				s = -r;
			g = x[m - 1] - p + ctmp[k - 1] / (g + s);
			s = F_ONE;
			cc = F_ONE;
			p = F_ZERO;

			mmk = m - k;
			for (ii = 1; ii <= mmk; ii++) {
				i = m - ii;
				f = s * ctmp[i - 1];
				bb = cc * ctmp[i - 1];
				if (fabs(f) < fabs(g))
					goto L_40;
				cc = g / f;
				r = sqrt(cc * cc + F_ONE);
				ctmp[i] = f * r;
				s = F_ONE / r;
				cc *= s;
				goto L_50;

		L_40:
				s = f / g;
				r = sqrt(s * s + F_ONE);
				ctmp[i] = g * r;
				cc = F_ONE / r;
				s *= cc;

		L_50:
				g = x[i] - p;
				r = (x[i - 1] - g) * s + (F_TWO * cc) * bb;
				p = s * r;
				x[i] = g + p;
				g = cc * r - bb;
				/* FORM FIRST COMPONENT OF VECTOR */
				f = w[i];
				w[i] = s * w[i - 1] + cc * f;
				w[i - 1] = cc * w[i - 1] - s * f;
			}

			x[k - 1] -= p;
			ctmp[k - 1] = g;
			ctmp[m - 1] = F_ZERO;
		}

/*		(4, 1, "No convergence in 100 iterations.");  */
                imsl_ermes(IMSL_FATAL, IMSL_NO_CONVERGE_100_ITERATIONS);
		goto L_120;

L_80:
		w[k - 1] *= wtsum * w[k - 1];
	}
	/*
	 * ORDER THE EIGENVALUES AND EIGENVECTORS
	 */
	for (ii = 2; ii <= *n; ii++) {
		i = ii - 1;
		k = i;
		p = x[i - 1];

		for (j = ii; j <= *n; j++) {
			if (x[j - 1] >= p)
				goto L_100;
			k = j;
			p = x[j - 1];
	L_100:
			;
		}

		if (k == i)
			goto L_110;
		x[k - 1] = x[i - 1];
		x[i - 1] = p;
		p = w[i - 1];
		w[i - 1] = w[k - 1];
		w[k - 1] = p;
L_110:
		;
	}

L_120:
	;
	return;
}				/* end of function */
