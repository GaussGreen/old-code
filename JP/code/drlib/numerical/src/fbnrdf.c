#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef TRUE
#undef TRUE
#define TRUE 1
#else
#define TRUE 1
#endif

#ifdef FALSE
#undef FALSE
#define FALSE 0
#else
#define FALSE 0
#endif


#ifdef ANSI
static Mfloat l_bnrdf(Mfloat *x, Mfloat *y, Mfloat *rho);
static Mfloat l_b1nnt(Mfloat *y, Mfloat *z);
#else
static Mfloat l_bnrdf();
static Mfloat l_b1nnt();
#endif

/****************************** USER INTERFACE ************************/

static Mfloat lv_answer;

#ifdef ANSI
Mfloat imsl_f_bivariate_normal_cdf(Mfloat x, Mfloat y, Mfloat rho)
#else
Mfloat imsl_f_bivariate_normal_cdf(x, y, rho)
    Mfloat      x, y, rho;
#endif
{
    Mfloat x2, y2, rho2;
    E1PSH("imsl_f_bivariate_normal_cdf","imsl_d_bivariate_normal_cdf");
    lv_answer = imsl_amach(6);
    x2 = (Mfloat) x; y2 = (Mfloat) y; rho2 = (Mfloat) rho;
    lv_answer = l_bnrdf(&x2, &y2, &rho2);
    E1POP("imsl_f_bivariate_normal_cdf","imsl_d_bivariate_normal_cdf");
    return lv_answer;
}


/*Translated by FOR_C++, v0.1, on 09/09/91 at 09:48:40 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 09/09/91 at 09:48:39
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BNRDF/DBNRDF (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    April 29, 1991

    Purpose:    Evaluate the bivariate normal distribution function.

    Usage:      BNRDF(X, Y, RHO)

    Arguments:
       X      - One argument for which the bivariate normal
                distribution function is to be evaluated.  (Input)
       Y      - The other argument for which the bivariate normal
                distribution function is to be evaluated.  (Input)
       RHO    - Correlation coefficient.  (Input)
       BNRDF  - Function value, the probability that a bivariate normal
                random variable with correlation RHO takes a value less
                than or equal to X and less than or equal to Y.
                (Output)

    Keywords:   P-value; Probability integral; Probability distribution;
                Continuous random variables; Cumulative distribution
                function; CDF

    GAMS:       L5b1n

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                Inverses
                MATH/LIBRARY SPECIAL FUNCTIONS Statistical Distribution
                Functions and Inverses

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_bnrdf(Mfloat *x, Mfloat *y, Mfloat *rho)
#else
static Mfloat l_bnrdf(x, y, rho)
	Mfloat          *x, *y, *rho;
#endif
{
	Mint             iax, iay, ind, ner;
	Mfloat           ax, ay, bnrdf_v,f1, p, qx, qy, tx,
	                ty, xy, z;


	imsl_e1psh("BNRDF");
	bnrdf_v = imsl_amach(6);
	ner = 1;
	/* Check RHO */
	if (fabs(*rho) > 1.0e0) {
		imsl_e1str(1, *rho);
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ABS_RHO_LE_1);
/*		imsl_e1mes(5, ner, "The absolute value of the correlation coefficient, RHO = %(r1), must be less than or equal to
1.");*/
		goto L_9000;
	}
	/* ABS(RHO) = 1 */
	if (fabs(*rho) == 1.0e0) {
        imsl_ermes (IMSL_WARNING_IMMEDIATE, IMSL_ABS_RHO_EQ_1); 
/*		imsl_e1mes(6, ner + 1, "Since the absolute value of the correlation coefficient, RHO, is equal to 1.0, the distribution is singular and the function is computed based on the univariate normal
distribution.");*/
		if (*rho > 0.0e0) {
			z = imsl_f_min(*x, *y);
			bnrdf_v = imsl_f_normal_cdf(z);
		} else {
			bnrdf_v = imsl_f_max(imsl_f_normal_cdf(*x) + imsl_f_normal_cdf(*y) - 1.0e0, 0.0e0);
		}
		goto L_9000;
	}
	f1 = 1.0e0 / sqrt(1.0e0 - imsl_fi_power(*rho, 2));
	xy = *x ** y;
	iax = 0;
	iay = 0;
	ind = 0;
	if (xy != 0.0e0) {
		ax = f1 * (*y / *x - *rho);
		ay = f1 * (*x / *y - *rho);
	} else {
		if (*x == 0.0e0 && *y == 0.0e0) {
			/*
			 * For X=Y=0 and AX=AY=(1-RHO)/SQRT(1-RHO**2)
			 */
			ax = f1 * (1.0e0 - *rho);
			ay = ax;
		} else {
			if (*x != 0.0e0) {
				/*
				 * For Y=0,X less than 0     TY = -1/4 For
				 * Y=0,X greater than 0  TY =  1/4
				 */
				ty = 0.25e0;
				if (*x < 0.0e0)
					ty = -ty;
				ax = -f1 ** rho;
				ind = 1;
			} else {
				/*
				 * For X=0,Y less than 0     TX = -1/4 For
				 * X=0,Y greater than 0  TX =  1/4
				 */
				tx = 0.25e0;
				if (*y < 0.0e0)
					tx = -tx;
				ay = -f1 ** rho;
				if (ay < 0.0e0) {
					iay = 1;
					ay = -ay;
				}
				ty = l_b1nnt(y, &ay);
				if (iay != 0)
					ty = -ty;
				goto L_10;
			}
		}
	}
	if (ax < 0.0e0) {
		iax = 1;
		ax = -ax;
	}
	tx = l_b1nnt(x, &ax);
	if (iax != 0)
		tx = -tx;
	if (ind == 0) {
		if (ay < 0.0e0) {
			iay = 1;
			ay = -ay;
		}
		ty = l_b1nnt(y, &ay);
		if (iay != 0)
			ty = -ty;
	}
L_10:
	if (*x <= 0.0e0) {
		qx = imsl_f_normal_cdf(*x);
	} else {
		qx = imsl_f_normal_cdf(-*x);
		qx = 1.0e0 - qx;
	}
	if (*y <= 0.0e0) {
		qy = imsl_f_normal_cdf(*y);
	} else {
		qy = imsl_f_normal_cdf(-*y);
		qy = 1.0e0 - qy;
	}
	/* Now evaluate P */
	p = 0.5e0 * (qx + qy) - tx - ty;
	if (xy <= 0.0e0 && (xy != 0.0e0 || *x + *y < 0.0e0))
		p -= 0.5e0;
	p = imsl_f_min(imsl_f_max(0.0e0, p), 1.0e0);

	bnrdf_v = p;

L_9000:
	imsl_e1pop("BNRDF" );
	return (bnrdf_v);
}				/* end of function */


/*Translated by FOR_C++, v0.1, on 09/09/91 at 13:49:44 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 09/09/91 at 13:49:43
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  B1NNT/DB1NNT (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    February 8, 1991

    Purpose:    Integral related to calculation of non-central T and
                bivariate normal probability distribution functions

    Usage:      B1NNT(Y, Z)

    Arguments:
       Y      - Real variable to be used in expression to be integrated.
                (Input)
       Z      - Real variable containing the upper limit of integration.
                (Input)
       B1NNT  - The value of the integral, from 0 to Z, of
                (EXP((-Y**2/2)*(1+X**2))/(2*PI*(1+X**2)))DX.  (Output)

    Keywords:   ???????

    GAMS:       ?????

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                             Inverses
                SFUN/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_b1nnt(Mfloat *y, Mfloat *z)
#else
static Mfloat l_b1nnt(y, z)
	Mfloat          *y, *z;
#endif

{
	Mfloat           a, a4, a4b4, ab4, aeps, ahsqb, asq, b, b1nnt_v,
	                b4, ber, bexp, big1, big2, d, d1, d2, eps, expov,
	                f, g, g1, hsqb, sum, t, ta, ter;
	static Mfloat    oned2p, twopi;
#ifdef COMPUTER_DECOSF
	static int first = TRUE;
#else
	static long first = TRUE;
#endif
/*        Mint            true = 1;
        Mint            false = 0;*/




	if (first) {
		twopi = 4.0e0 * asin(1.0e0);
		oned2p = 1.0e0 / twopi;
		first = FALSE;
	}
	expov = -log(imsl_amach(1)) - imsl_amach(4);
	eps = imsl_amach(4);
	b = fabs(*y);
	a = fabs(*z);

	if (a == 0.0e0) {
		t = 0.0e0;
		goto L_9000;
	}
	big1 = pow(imsl_amach(2), 0.25e0);
	big2 = big1 / log(big1);
	if (a > big2) {
		t = 0.5e0 * (1.0e0 - imsl_f_normal_cdf(b));
		if (*z < 0.0e0)
			t = -t;
		goto L_9000;
	}
	ta = atan(a);
	if (b == 0.0e0) {
		t = oned2p * ta;
		if (*z < 0.0e0)
			t = -t;
		goto L_9000;
	}
	if (a * b > 4.0e0) {
		t = 0.25e0 - 0.5e0 * (imsl_f_normal_cdf(b) - 0.5e0);
		if (*z < 0.0e0)
			t = -t;
		goto L_9000;
	}
	hsqb = 0.5e0 * b * b;
	if (hsqb <= expov) {
		bexp = exp(-hsqb);
		asq = a * a;
		a4 = asq * asq;
		b4 = hsqb * hsqb;
		a4b4 = a4 * b4;
		ahsqb = a * hsqb;
		ab4 = a * b4 * 0.5e0;
		f = 1.0e0;
		sum = 0.0e0;
		g = 3.0e0;
		/* Begin series expansion */
L_10:
		g1 = g;
		ber = 0.0e0;
		ter = ab4;
L_20:
		ber += ter;
		if (ter > ber * eps) {
			/* Develop coefficient series */
			ter *= hsqb / g1;
			g1 += 1.0e0;
			goto L_20;
		}
		d1 = (ber + ahsqb) / f;
		d2 = ber * asq / (f + 2.0e0);
		d = d1 - d2;
		sum += d;
		t = ta - sum * bexp;
		if (t > 0.0e0) {
			aeps = eps * t;
		} else {
			aeps = eps;
		}
		ahsqb = ahsqb * a4b4 / ((g - 1.0e0) * g);
		ab4 = ab4 * a4b4 / ((g + 1.0e0) * g);
		f += 4.0e0;
		g += 2.0e0;
		/* Should series expansion be terminated */
		if (d2 * bexp >= aeps)
			goto L_10;
		t *= oned2p;
		if (*z < 0.0e0)
			t = -t;
	}
L_9000:
	b1nnt_v = t;
	return (b1nnt_v);
}				/* end of function */




