    


























#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  RECCF/DRECCF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 1, 1985

    Purpose:    Compute recurrence coefficients for various monic
                polynomials.

    Usage:      CALL RECCF (N, IWEIGH, ALPHA, BETA, B, C)

    Arguments:
       N      - Number of recurrence coefficients.  (Input)
       IWEIGH - Index of the weight function.  (Input)
                IWEIGH    WT(X)           Interval    Name
                1           1              (-1,+1)   Legendre
                2       1/SQRT(1-X**2)     (-1,+1)   Chebyshev 1st kind
                3       SQRT(1-X**2)       (-1,+1)   Chebyshev 2nd kind
                4       EXP(-X**2)       (-inf,+inf) Hermite
                5 (1-X)**ALPHA*(1+X)**BETA (-1,+1)   Jacobi
                6     EXP(-X)*X**ALPHA      (0,+inf) Generalized Laguerre
                7       1/COSH(X)        (-inf,+inf) COSH
       ALPHA  - Parameter used in the weight function with some values
                of IWEIGH, otherwise it is ignored.  (Input)
       BETA   - Parameter used in the weight function with some values
                of IWEIGH, otherwise it is ignored.  (Input)
       B      - Array of length N containing recurrence coefficients.
                (Output)
       C      - Array of length N containing recurrence coefficients.
                (Output)

    Remark:
       The recurrence coefficients B(I) and C(I) define the monic
       polynomials via the relation
          P(I) = (X - B(I+1))*P(I-1) - C(I+1)*P(I-2)
       The zero-th moment (INTEGRAL WT(X) dX) of the weight function
       is returned in C(1).

    Keywords:   Univariate quadrature; Three-term recurrence; Orthogonal
                polynomials; Legendre; Chebyshev; Hermite; Jacobi;
                Laguerre; Numerical integration

    GAMS:       H2c

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_reccf(Mint *n, Mint *iweigh, Mfloat *alpha, Mfloat *imsl_beta, Mfloat b[], Mfloat c[])
#else
void imsl_reccf(n, iweigh, alpha, imsl_beta, b, c)
	Mint            *n, *iweigh;
	Mfloat          *alpha, *imsl_beta, b[], c[];
#endif
{
	Mint             i;
	Mfloat           a2b2, ab, abi, xi;
	static Mfloat    pi = 3.1415926535897932384626433831;



	imsl_e1psh("RECCF ");

	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		(5, 1, "The number of recurrence coefficients N = %(i1).  It must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_RECURRENCE_COEFF);
		goto L_130;
	}
	if (*iweigh < 1 || *iweigh > 7) {
		imsl_e1sti(1, *iweigh);

/*		(5, 2, "The index of the weight function IWEIGH = %(i1).  It must be in the range 1 to 7.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_IWEIGHT_OUT_OF_RANGE);
		goto L_130;
	}
	if (*iweigh == 2)
		goto L_30;
	if (*iweigh == 3)
		goto L_40;
	if (*iweigh == 4)
		goto L_50;
	if (*iweigh == 5)
		goto L_70;
	if (*iweigh == 6)
		goto L_90;
	if (*iweigh == 7)
		goto L_110;
	/*
	 * IWEIGH = 1.  Legendre polynomials P(X) on (-1,+1) with W(X) = 1.
	 */
	c[0] = F_TWO;
	for (i = 1; i <= (*n - 1); i++) {
		b[i - 1] = F_ZERO;
		c[i] = (float) (i * i) / (imsl_fi_power((float) (2 * i), 2) - F_ONE);
	}
	b[*n - 1] = F_ZERO;
	goto L_130;
	/*
	 * IWEIGH = 2.  Chebyshev polynomials of the first kind T(X) on
	 * (-1,+1) with W(X) = 1/SQRT(1-X*X)
	 */
L_30:
	c[0] = pi;
	sset(*n - 1, F_ZERO, b, 1);
	sset(*n - 1, 0.25, &c[1], 1);
	if (*n > 1)
		c[1] = F_HALF;
	b[*n - 1] = F_ZERO;
	goto L_130;
	/*
	 * IWEIGH = 3.  Chevyshev polynomials of the second kind U(X) on
	 * (-1,+1) with W(X) = SQRT(1-X*X).
	 */
L_40:
	c[0] = pi / F_TWO;
	sset(*n - 1, F_ZERO, b, 1);
	sset(*n - 1, 0.25, &c[1], 1);
	b[*n - 1] = F_ZERO;
	goto L_130;
	/*
	 * IWEIGH = 4.  Hermite polynomials H(X) on (-INF,+INF) with W(X) =
	 * EXP(-X**2)
	 */
L_50:
	c[0] = sqrt(pi);
	sset(*n - 1, F_ZERO, b, 1);
	for (i = 1; i <= (*n - 1); i++) {
		c[i] = (float) (i) / F_TWO;
	}
	b[*n - 1] = F_ZERO;
	goto L_130;
	/*
	 * IWEIGH = 5.  Jacobi polynomials P(ALPHA, BETA)(X) on (-1,+1) with
	 * W(X) = (1-X)**ALPHA * (1+X)**BETA, ALPHA and BETA .GT. -1.0 . Note
	 * that when ALPHA = BETA, Gegenbauer polynomials result.
	 */
L_70:
	if (*alpha <= -F_ONE) {
		imsl_e1str(1, *alpha);
                imsl_e1stl(1, "alpha");
/*		(5, 3, "The parameter in weight function ALPHA = %(r1).  It must be greater than -1.0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WEIGHT_FCN_PARAMETER);
		goto L_130;
	}
	if (*imsl_beta <= -F_ONE) {
		imsl_e1str(1, *imsl_beta);
                imsl_e1stl(1, "beta");
/*		(5, 4, "The parameter in weight function BETA = %(r1).  It must be greater than -1.0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WEIGHT_FCN_PARAMETER);
		goto L_130;
	}
	ab = *alpha + *imsl_beta;
	abi = F_TWO + ab;
	c[0] = pow(F_TWO, ab + F_ONE) * imsl_f_gamma(*alpha + F_ONE) * imsl_f_gamma(*imsl_beta +
						   F_ONE) / imsl_f_gamma(abi);
	b[0] = (*imsl_beta - *alpha) / abi;
	c[1] = F_FOUR * (F_ONE + *alpha) * (F_ONE + *imsl_beta) / ((abi + F_ONE) * abi *
								  abi);
	a2b2 = *imsl_beta ** imsl_beta - *alpha ** alpha;
	for (i = 2; i <= (*n - 1); i++) {
		xi = i;
		abi = F_TWO * xi + ab;
		b[i - 1] = a2b2 / ((abi - F_TWO) * abi);
		c[i] = F_FOUR * xi * (xi + *alpha) * (xi + *imsl_beta) * (xi + ab) / ((abi *
						  abi - F_ONE) * abi * abi);
	}
	abi = F_TWO * (float) (*n) + ab;
	b[*n - 1] = a2b2 / ((abi - F_TWO) * abi);
	goto L_130;
	/*
	 * IWEIGH = 6.  Laguerre polynomials L(ALPHA)(X) on (0,+INF) with
	 * W(X) = X**ALPHA * EXP(-X), ALPHA .GT. -1.
	 */
L_90:
	if (*alpha <= -F_ONE) {
		imsl_e1str(1, *alpha);
                imsl_e1stl(1, "alpha");
/*		(5, 3, "The parameter in weight function ALPHA = %(r1).  It must be greater than -1.0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WEIGHT_FCN_PARAMETER);
		goto L_130;
	}
	c[0] = imsl_f_gamma(*alpha + F_ONE);
	for (i = 1; i <= (*n - 1); i++) {
		xi = i;
		b[i - 1] = F_TWO * xi - F_ONE + *alpha;
		c[i] = xi * (xi + *alpha);
	}
	b[*n - 1] = F_TWO * (float) (*n) - F_ONE + *alpha;
	goto L_130;
	/*
	 * IWEIGH = 7.  Polynomials orthogonal on (-INF,+INF) with weight
	 * W(X) = 1/COSH(X).
	 */
L_110:
	c[0] = pi;
	for (i = 1; i <= (*n - 1); i++) {
		b[i - 1] = F_ZERO;
		c[i] = imsl_fi_power((float) (i) * pi / F_TWO, 2);
	}
	b[*n - 1] = F_ZERO;

L_130:
	;
	imsl_e1pop("RECCF ");
	return;
}				/* end of function */
