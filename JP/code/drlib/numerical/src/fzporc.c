#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK  l_zeros_poly(Mint ndeg, Mfloat *coeff,  va_list argptr);
static void     l_zporc(Mint *ndeg, Mfloat coeff[], Mf_complex root[]);
static void     l_z3lrc(Mfloat *a, Mfloat *b, Mfloat *c, 
                        Mf_complex *zsm, Mf_complex *zlg);
static void     l_z3orc(Mint *l2, Mint *nz);
static void     l_z4orc(Mfloat *uu, Mfloat *vv, Mint *nz);
static void     l_z5orc(Mfloat *sss, Mint *nz, Mint *iflag);
static void     l_z6orc(Mint *itype);
static void     l_z7orc(Mint *itype);
static void     l_z8orc(Mint *itype, Mfloat *uu, Mfloat *vv);
static void     l_z9orc(Mint *nn, Mfloat *u, Mfloat *v, Mfloat p[],
                        Mfloat q[], Mfloat *ra, Mfloat *rb);
static void     l_z10rc(Mfloat *ra, Mfloat *b1, Mfloat *c,
                        Mf_complex *s, Mf_complex *rl);
static void     l_sdadd(Mfloat *a, Mdouble acc[]);
static void     l_sdini(Mfloat *s, Mdouble acc[]);
static void     l_sdmul(Mfloat *a, Mfloat *b, Mdouble acc[]);
static void     l_sdsto(Mdouble acc[], Mfloat *s);
static void     l_sswap(Mint n, Mfloat sx[], Mint incx, Mfloat sy[],
                        Mint incy);
void            imsl_cset(Mint *n, Mf_complex *ca, Mf_complex *cx,
                          Mint *incx);
#else
static VA_LIST_HACK l_zeros_poly();
static void     l_zporc();
static void     l_z3lrc();
static void     l_z3orc();
static void     l_z4orc();
static void     l_z5orc();
static void     l_z6orc();
static void     l_z7orc();
static void     l_z8orc();
static void     l_z9orc();
static void     l_z10rc();
static void     l_sdadd();
static void     l_sdini();
static void     l_sdmul();
static void     l_sdsto();
static void     l_sswap();
void            imsl_cset();
#endif


static Mf_complex   *lv_roots;
 
#ifdef ANSI
Mf_complex *imsl_f_zeros_poly(Mint ndeg, Mfloat *coeff, ...)
#else
Mf_complex *imsl_f_zeros_poly(ndeg, coeff, va_alist)
    Mint        ndeg;
    Mfloat      *coeff;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, coeff);
    E1PSH("imsl_f_zeros_poly", "imsl_d_zeros_poly");
    lv_roots = NULL;
    IMSL_CALL(l_zeros_poly(ndeg, coeff, argptr));
    va_end(argptr);
    E1POP("imsl_f_zeros_poly", "imsl_d_zeros_poly"); 
    return lv_roots;
}
 
 
#ifdef ANSI
static VA_LIST_HACK l_zeros_poly(Mint ndeg, Mfloat *coef, va_list argptr)
#else
static VA_LIST_HACK l_zeros_poly(ndeg, coef, argptr)
    Mint        ndeg;
    Mfloat      *coef;
    va_list     argptr;
#endif
{
    Mint            code;
    Mint            arg_number  = 2;
    Mint            return_user = 0;
    Mint            use_companion = 0;
    Mfloat         *companion = NULL;
    Mint            last_row_offset;
    Mint            i;

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch (code) {
                 /* 
                  * In order to use the  companion matrix method of finding 
                  * the roots the follwing case was added.  This change was 
                  * for IMSL/IDL.
                  */
            case IMSL_COMPANION_METHOD:
                use_companion = 1;
                arg_number++;
                break;
            case IMSL_RETURN_USER:
                lv_roots = va_arg(argptr, Mf_complex*);
                arg_number++;
                return_user = 1;
                break;
            case 0:
                break;
            default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
                break;
        }
    } 
    if (imsl_n1rty(0)) goto RETURN;
    if (use_companion) {
      /* 
       * NOTE: for IMSL/IDL, the ndeg is always positive, thus
       * there is no need to check it.  C/Math should check it though. 
       */
      companion = (Mfloat*) calloc(ndeg*ndeg, sizeof(*companion));
      if (companion == NULL) {
        imsl_e1sti(1, ndeg);
        imsl_e1stl(1, "ndeg");
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto RETURN;
      }
      
      if (lv_roots == NULL) {
        lv_roots = (Mf_complex *)imsl_malloc(ndeg*sizeof(*lv_roots));
        if (lv_roots == NULL) {
	  imsl_e1sti(1, ndeg);
	  imsl_e1stl(1,"ndeg");
	  imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	  goto RETURN;
        }
      }
      
      /* build companion matrix */
      
      /* fill super-diagonal */
      
      for (i=0; i<ndeg-1; i++) {
        *(companion+(ndeg+1)*i+1) = 1.0;
      }
      
      /* fill last row with -coef[i]/coef[ndeg], i = 0 ,..., n-1 */
      
      last_row_offset = (ndeg-1) * ndeg;
      for (i=0; i<ndeg; i++) 
        *(companion+last_row_offset+i) = -coef[i]/coef[ndeg];
      
      /* roots of p(z) = e-values of companion matrix */
      
      imsl_f_eig_gen(ndeg, companion, IMSL_RETURN_USER, lv_roots, 0);

    } else {
      /* CHECK FOR DEGREE > 100 OR < = ZERO.  */
      if (ndeg > 100 || ndeg <= 0) {
        imsl_e1sti(1, ndeg);
	/*      (5, 1, "The degree of the polynomial must be less */
	/*              than 100 and greater than zero.  NDEG is  */
	/*              given as %(i1)."); */
        imsl_ermes(IMSL_TERMINAL, IMSL_POLYNOMIAL_DEGREE);
        goto RETURN;
      }
      
      if (imsl_n1rty(0)) goto RETURN;
      
      if (lv_roots == NULL) {
        lv_roots = (Mf_complex *)imsl_malloc(ndeg*sizeof(*lv_roots));
        if (lv_roots == NULL) {
	  imsl_e1sti(1, ndeg);
	  imsl_e1stl(1,"ndeg");
	  imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	  goto RETURN;
        }
      }
      
      l_zporc(&ndeg, coef, lv_roots);
    }
RETURN:
    if (imsl_n1rty(0) > 3) {
        if (!return_user && lv_roots != NULL)   imsl_free(lv_roots);
        lv_roots = NULL;
    }
    return (argptr);
}


/* -----------------------------------------------------------------------
    IMSL Name:  ZPORC/DZPORC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Find the zeros of a polynomial with real coefficients
                using the Jenkins-Traub three-stage algorithm.

    Usage:      CALL ZPORC (NDEG, COEFF, ROOT)

    Arguments:
       NDEG   - Degree of the polynomial.  (Input)
       COEFF  - Vector of length NDEG+1 containing the coefficients of
                the polynomial in increasing order by degree.  (Input)
                The polynomial is COEFF(NDEG+1)*Z**NDEG +
                COEFF(NDEG)*Z**(NDEG-1) + ... + COEFF(1).
       ROOT   - Complex vector of length NDEG containing the zeros of
                the polynomial.  (Output)

    Remark:
       Informational errors
       Type Code
         3   1  The first several coefficients of the polynomial are
                equal to zero.  Several of the last roots will be set to
                machine infinity to compensate for this problem.
         3   2  Fewer than NDEG zeros were found.  The ROOT vector will
                contain the value for machine infinity in the locations
                which do not contain zeros.

    Keywords:   Polynomials; Roots; Deflation; Nonlinear equations

    GAMS:       F1a1a

    Chapter:    MATH/LIBRARY Nonlinear Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
static struct t_z11rc {
	Mfloat           p[101], qp[101], rk[101], qk[101], svk[101], sr,
	                imsl_si, u, v, ra, rb, c, d, a1, a2, a3, a6, a7,
	                e, f, g, h, szr, szi, rlzr, rlzi, eta, are, rmre;
	Mint             n, nn;
}               lv_z11rc;

#ifdef ANSI
static void l_zporc(Mint *ndeg, Mfloat coeff[], Mf_complex root[])
#else
static void l_zporc(ndeg, coeff, root)
	Mint            *ndeg;
	Mfloat           coeff[];
	Mf_complex       root[];
#endif
{
#ifdef COMPUTER_DECOSF
	int             zerok;
#else
	long            zerok;
#endif
	Mint             _l0, i, icnt, j, jj, k, l, n2, newdeg,
	                nm1, nz;
	Mfloat           aa, bb, bnd, cc, cosr, df, dx, factor, ff, finity,
	                fn, pt[101], radix, repsp, repsr1, rinfp, rlo,
	                rmax, rmin, sc, sinr, t, temp[101], x, xm, xx,
	                xxx, yy;
	Mf_complex       _cx0;



	imsl_e1psh("l_zporc");

	/* REVERSE THE ORDER OF THE COEFFICIENT VECTOR.  */

	lv_z11rc.n = *ndeg;
	lv_z11rc.nn = *ndeg + 1;
	k = mod(lv_z11rc.nn, 2) + 1 + lv_z11rc.nn / 2;
	n2 = *ndeg;
	l_sswap(lv_z11rc.nn / 2, &coeff[0], 1, &coeff[k - 1], -1);
	/*
	 * CALL Z3LRC IF POLYNOMIAL IS QUADRATIC
	 */
	if (coeff[0] != F_ZERO && *ndeg == 2) {
		l_z3lrc(&coeff[0], &coeff[1], &coeff[2], &root[0], &root[1]);
		goto L_160;
	}

	/* CHECK FOR LEADING COEFFICIENTS EQUAL TO ZERO.  */

	l = 0;
	rinfp = imsl_amach(2);
	finity = imsl_amach(7);
	if (coeff[0] == F_ZERO) {
L_10:
		l += 1;
		if (l <= *ndeg) {
			root[n2 - 1] = imsl_cf_convert(finity, F_ZERO);
			n2 -= 1;
			if (coeff[l] == F_ZERO)
				goto L_10;
		}
		if (l <= *ndeg) {
			imsl_e1sti(1, l);
			imsl_e1sti(2, l);
		} else {
			imsl_e1sti(1, l);
			imsl_e1sti(2, l - 1);
		}
/*		imsl_ermes(3, 2, "The %(i1) leading coefficients of the    */
/*                                polynomial are equal to zero.  The last  */
/*                                %(i2) roots will be set to               */
/*                                ("infinity",0.0) where "infinity" is the */
/*                                largest machine constant.");             */
		imsl_ermes(IMSL_WARNING, IMSL_ZERO_COEFF);
	}
	if (l >= *ndeg)
		goto L_160;
	newdeg = *ndeg - l;
	lv_z11rc.nn = newdeg + 1;
	/* INITIALIZE REMAINING CONSTANTS */
	radix = imsl_i_machine(6);
	repsp = imsl_amach(1);
	repsr1 = imsl_amach(4);
	lv_z11rc.eta = repsr1;
	lv_z11rc.are = lv_z11rc.eta;
	lv_z11rc.rmre = lv_z11rc.eta;
	rlo = repsp / lv_z11rc.eta;
	xx = sqrt(F_HALF);
	yy = -xx;
	cosr = -.06975647;
	sinr = .9975641;
	/*
	 * REMOVE THE ZEROS AT THE ORIGIN IF ANY EXIST.
	 */
L_20:
	if (coeff[lv_z11rc.nn + l - 1] != F_ZERO)
		goto L_30;
	j = newdeg - lv_z11rc.nn + 2;
	root[j - 1] = imsl_cf_convert(F_ZERO, F_ZERO);
	lv_z11rc.nn -= 1;
	if (lv_z11rc.nn == 1)
		goto L_160;
	goto L_20;
	/* MAKE A COPY OF THE COEFFICIENTS */
L_30:
	lv_z11rc.n = lv_z11rc.nn - 1;
	scopy(lv_z11rc.nn, &coeff[l], 1, lv_z11rc.p, 1);
	/* START THE ALGORITHM FOR ONE ZERO */
L_40:
	if (lv_z11rc.n <= 2) {
		/*
		 * CALCULATE THE FINAL ZERO OR PAIR OF ZEROS
		 */
		if (lv_z11rc.n == 1) {
			root[newdeg - 1] = imsl_cf_convert(-lv_z11rc.p[1] / lv_z11rc.p[0], F_ZERO);
		} else if (lv_z11rc.n == 2) {
			l_z10rc(&lv_z11rc.p[0], &lv_z11rc.p[1], &lv_z11rc.p[2], &root[newdeg - 2],
				   &root[newdeg - 1]);
		}
		if (lv_z11rc.n < 1)
			goto L_160;
		goto L_150;
	}
	/*
	 * FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
	 */
	rmax = F_ZERO;
	rmin = rinfp;
	for (i = 1; i <= lv_z11rc.nn; i++) {
		x = fabs(lv_z11rc.p[i - 1]);
		if (x > rmax)
			rmax = x;
		if (x != F_ZERO && x < rmin)
			rmin = x;
	}
	/*
	 * SCALE IF THERE ARE LARGE OR VERY SMALL COEFFICIENTS COMPUTES A
	 * SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL. THE
	 * SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID UNDETECTED
	 * UNDERFLOW INTERFERING WITH THE CONVERGENCE CRITERION. THE FACTOR
	 * IS A POWER OF THE BASE.
	 */
	sc = rlo / rmin;

	if (sc <= F_ONE && rmax >= F_TEN) {
		if (sc == F_ZERO)
			sc = repsp * radix * radix;
		l = log(sc) / log(radix) + F_HALF;
		if (l != 0) {
			factor = imsl_fi_power(radix, l);
			sscal(lv_z11rc.nn, factor, lv_z11rc.p, 1);
		}
	} else if (sc > F_ONE) {
		if (rinfp / sc >= rmax) {
			l = log(sc) / log(radix) + F_HALF;
			if (l != 0) {
				factor = imsl_fi_power(radix, l);
				sscal(lv_z11rc.nn, factor, lv_z11rc.p, 1);
			}
		}
	}
	/*
	 * COMPUTE LOWER BOUND ON MODULI OF ZEROS.
	 */
	for (i = 1; i <= lv_z11rc.nn; i++) {
		pt[i - 1] = fabs(lv_z11rc.p[i - 1]);
	}
	pt[lv_z11rc.nn - 1] = -pt[lv_z11rc.nn - 1];
	/* COMPUTE UPPER ESTIMATE OF BOUND */
	x = exp((log(-pt[lv_z11rc.nn - 1]) - log(pt[0])) / lv_z11rc.n);
	if (pt[lv_z11rc.n - 1] != F_ZERO) {
		/*
		 * IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
		 */
		xm = -pt[lv_z11rc.nn - 1] / pt[lv_z11rc.n - 1];
		if (xm < x)
			x = xm;
	}
	/*
	 * CHOP THE INTERVAL (0,X) UNTIL FF.LE.0
	 */
L_70:
	xm = x * 0.1;
	ff = pt[0];
	for (i = 2; i <= lv_z11rc.nn; i++) {
		ff = ff * xm + pt[i - 1];
	}
	if (ff > F_ZERO) {
		x = xm;
		goto L_70;
	}
	dx = x;
	/*
	 * DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES
	 * DOWHILE (ABS(DX/X).GT..005)
	 */
L_10034:
	if (!(fabs(dx / x) > .005))
		goto L_10036;
	ff = pt[0];
	df = ff;
	for (i = 2; i <= lv_z11rc.n; i++) {
		ff = ff * x + pt[i - 1];
		df = x * df + ff;
	}
	ff = ff * x + pt[lv_z11rc.nn - 1];
	dx = ff / df;
	x -= dx;
	/* ENDWHILE */

	goto L_10034;
L_10036:
	;
	bnd = x;
	/*
	 * COMPUTE THE DERIVATIVE AS THE INTIAL K POLYNOMIAL AND DO 5 STEPS
	 * WITH NO SHIFT
	 */
	nm1 = lv_z11rc.n - 1;
	fn = F_ONE / lv_z11rc.n;
	for (i = 2; i <= lv_z11rc.n; i++) {
		lv_z11rc.rk[i - 1] = (lv_z11rc.nn - i) * lv_z11rc.p[i - 1] * fn;
	}
	lv_z11rc.rk[0] = lv_z11rc.p[0];
	aa = lv_z11rc.p[lv_z11rc.nn - 1];
	bb = lv_z11rc.p[lv_z11rc.n - 1];
	zerok = lv_z11rc.rk[lv_z11rc.n - 1] == F_ZERO;
	for (jj = 1; jj <= 5; jj++) {
		cc = lv_z11rc.rk[lv_z11rc.n - 1];
		if (!zerok) {
			/*
			 * USE SCALED FORM OF RECURRENCE IF VALUE OF K AT 0
			 * IS NONZERO
			 */
			t = -aa / cc;
			for (i = 1; i <= nm1; i++) {
				j = lv_z11rc.nn - i;
				lv_z11rc.rk[j - 1] = t * lv_z11rc.rk[j - 2] + lv_z11rc.p[j - 1];
			}
			lv_z11rc.rk[0] = lv_z11rc.p[0];
			zerok = fabs(lv_z11rc.rk[lv_z11rc.n - 1]) <= fabs(bb) * lv_z11rc.eta *
				F_TEN;
		} else {
			/* USE UNSCALED FORM OF RECURRENCE */
			scopy(nm1, &lv_z11rc.rk[0], -1, &lv_z11rc.rk[1], -1);
			lv_z11rc.rk[0] = F_ZERO;
			zerok = lv_z11rc.rk[lv_z11rc.n - 1] == F_ZERO;
		}
	}
	/* SAVE K FOR RESTARTS WITH NEW SHIFTS */
	scopy(lv_z11rc.n, lv_z11rc.rk, 1, temp, 1);
	/*
	 * LOOP TO SELECT THE QUADRATIC CORRESPONDING TO EACH NEW SHIFT
	 */
	for (icnt = 1; icnt <= 20; icnt++) {
		/*
		 * QUADRATIC CORRESPONDS TO A DOUBLE SHIFT TO A NON-REAL
		 * POINT AND ITS COMPLEX CONJUGATE. THE POINT HAS MODULUS BND
		 * AND AMPLITUDE ROTATED BY 94 DEGREES FROM THE PREVIOUS
		 * SHIFT.
		 */
		xxx = cosr * xx - sinr * yy;
		yy = sinr * xx + cosr * yy;
		xx = xxx;
		lv_z11rc.sr = bnd * xx;
		lv_z11rc.imsl_si = bnd * yy;
		lv_z11rc.u = -lv_z11rc.sr - lv_z11rc.sr;
		lv_z11rc.v = bnd * bnd;
		/*
		 * SECOND STAGE CALCULATION, FIXED QUADRATIC
		 */
                _l0 = 20 * icnt;
		l_z3orc(&_l0, &nz);
		if (nz == 0)
			goto L_130;
		/*
		 * THE SECOND STAGE JUMPS DIRECTLY TO ONE OF THE THIRD STAGE
		 * ITERATIONS AND RETURNS HERE IF SUCCESSFUL. DEFLATE THE
		 * POLYNOMIAL, STORE THE ZERO OR ZEROS AND RETURN TO THE MAIN
		 * ALGORITHM.
		 */
		j = newdeg - lv_z11rc.n + 1;
		root[j - 1] = imsl_cf_convert(lv_z11rc.szr, lv_z11rc.szi);
		lv_z11rc.nn -= nz;
		lv_z11rc.n = lv_z11rc.nn - 1;
		scopy(lv_z11rc.nn, lv_z11rc.qp, 1, lv_z11rc.p, 1);
		if (nz == 1)
			goto L_40;
		root[j] = imsl_cf_convert(lv_z11rc.rlzr, lv_z11rc.rlzi);
		goto L_40;
		/*
		 * IF THE ITERATION IS UNSUCCESSFUL ANOTHER QUADRATIC IS
		 * CHOSEN AFTER RESTORING K
		 */
L_130:
		scopy(lv_z11rc.n, temp, 1, lv_z11rc.rk, 1);
	}
	/*
	 * RETURN WITH FAILURE IF NO CONVERGENCE WITH 20 SHIFTS
	 */
	imsl_e1sti(1, newdeg - lv_z11rc.n);
	imsl_e1sti(2, lv_z11rc.n);

/*	(3, 3,"Only %(i1) roots were found.  The "root" vector will contain   */
/*             the value for machine infinity in the last %(i2) locations."); */
        imsl_ermes(IMSL_WARNING, IMSL_FEWER_ZEROS_FOUND);
	/*
	 * SET UNFOUND ROOTS TO MACHINE INFINITY
	 */
L_150:
	if (imsl_n1rcd(0) == (Mint)IMSL_FEWER_ZEROS_FOUND ) {
		n2 = newdeg - lv_z11rc.n + 1;
                _cx0 = imsl_cf_convert(finity, F_ZERO);
                _l0 = 1;
		imsl_cset(&lv_z11rc.n, &_cx0, &root[n2 - 1], &_l0);
	}
L_160:
	lv_z11rc.nn = *ndeg + 1;
	k = mod(lv_z11rc.nn, 2) + 1 + lv_z11rc.nn / 2;
	l_sswap(lv_z11rc.nn / 2, &coeff[0], 1, &coeff[k - 1], -1);

	imsl_e1pop("l_zporc");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z3LRC/DZ3LRC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus called by ZPLRC to find the zeros of a quadratic
                with real coefficients.

    Usage:      CALL Z3LRC (A, B, C, ZSM, ZLG)

    Arguments:
       A      - Coefficient of the quadratic Ax**2 + Bx + C.  (Input)
       B      - Coefficient of the quadratic Ax**2 + Bx + C.  (Input)
       C      - Coefficient of the quadratic Ax**2 + Bx + C.  (Input)
       ZSM    - The smallest root of the quadratic in absolute value.
                  (Output)
       ZLG    - The largest root of the quadratic in absolute value.
                  (Output)

    Remark:
       Informational errors
       Type Code
         3   1  Implies A and B equal zero.  On output ZLG will be set to
                CMPLX(FINITY,0.0) and ZSM will be set to -ZLG.  FINITY is
                largest machine constant.
         3   2  Implies A equals zero.  On output ZLG will be set to
                CMPLX(FINITY,0.0).  FINITY is the largest machine con-
                stant.

    Chapter:    MATH/LIBRARY Nonlinear Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z3lrc(Mfloat *a, Mfloat *b, Mfloat *c, Mf_complex *zsm,
                    Mf_complex *zlg)
#else
static void l_z3lrc(a, b, c, zsm, zlg)
	Mfloat          *a, *b, *c;
	Mf_complex      *zsm, *zlg;
#endif
{
	Mint             is;
	Mfloat           _f0, a0, b0, b1, c0, d, d1, dd, finity, radix,
	                 rnlgrx, s, scale;
	Mdouble          acc[2];
	Mf_complex       zl, zs;


	imsl_e1psh("l_z3lrc");
	/* INITIALIZE MACHINE CONSTANTS */
	radix = imsl_i_machine(6);
	rnlgrx = log(radix);
	finity = imsl_amach(7);
	/*
	 * PUT THE COEFFICIENTS IN TEMPORARY STORAGE TO SAVE EXECUTION TIME.
	 */
	a0 = *a;
	b1 = -*b;
	c0 = *c;
	/*
	 * CHECK FOR A AND/OR B EQUAL TO ZERO.
	 */
	if (a0 == F_ZERO) {
		if (b1 == F_ZERO) {
                        imsl_e1sti(1, 2);
                        imsl_e1sti(2, 2);
/*                      (3, 2, "The %(i1) leading coefficients of the       */
/*                              polynomial are equal to zero.  The last     */
/*                              %(i2) roots will be set to ("infinity",0.0) */
/*                              where "infinity" is the largest machine     */
/*                              constant.");                                */
                        imsl_ermes(IMSL_WARNING, IMSL_ZERO_COEFF);

			*zlg = imsl_cf_convert(finity, F_ZERO);
			*zsm = imsl_c_neg(*zlg);
		} else {
                        imsl_e1sti(1, 1); 
                        imsl_e1sti(2, 1); 
                        imsl_ermes(IMSL_WARNING, IMSL_ZERO_COEFF); 

			*zlg = imsl_cf_convert(finity, F_ZERO);
			*zsm = imsl_cf_convert(c0 / b1, F_ZERO);
		}
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	/* CHECK FOR C EQUAL TO ZERO */
	if (c0 != F_ZERO) {

		/*
		 * SCALING TO AVOID OVERFLOW OR UNDERFLOW. SCALE THE
		 * COEFFICIENTS SO THAT A*C IS APPROXIMATELY ONE. THE SCALE
		 * FACTOR CSQRT(A*C) FITS THIS REQUIREMENT BUT MAY CAUSE
		 * OVERFLOW OR UNDERFLOW IN THE SCALING PROCEDURE. LET
		 * AMAX1(ABS(AR),ABS(AI)) BE REPRESENTED BY RADIX**IA AND LET
		 * AMAX1(ABS(CR),ABS(CI) BE REPRESENTED BY RADIX**IC. THE
		 * SCALE FACTOR, SCALE, IS DEFINED BY THE FOLLOWING FORMULA:
		 * SCALE=RADIX**IS, WHERE IS=ENTIER((IA+IC+1)/2) AND ENTIER
		 * IS THE MATHEMATICAL GREATEST INTEGER FUNCTION.
		 */
		is = (Mint) ((log(fabs(a0)) + log(fabs(c0)) + rnlgrx) /
			    (rnlgrx + rnlgrx));
		scale = imsl_fi_power(radix, is);
		/*
		 * IF THE SCALE FACTOR .LE. DEPS*ABS(B1) DO NOT SCALE THE
		 * COEFFICIENTS.
		 */
		d1 = fabs(b1);
		l_sdini(&d1, acc);
		l_sdadd(&scale, acc);
                _f0 = -d1;
		l_sdadd(&_f0, acc);
		l_sdsto(acc, &d);

		if (d != F_ZERO) {
			/*
			 * IF ABS(B1) .GE. EPS*SCALE FACTOR THEN SCALE B0.
			 * OTHERWISE SET B0 = ZERO.
			 */
			b0 = F_ZERO;
			l_sdini(&d1, acc);
			l_sdadd(&scale, acc);
                        _f0 = -scale;
			l_sdadd(&_f0, acc);
			l_sdsto(acc, &d);
			if (d != F_ZERO)
				b0 = (b1 / scale) * F_HALF;
			a0 /= scale;
			c0 /= scale;
			/* SOLVE A0*Z**2-2.0*B0*Z+C0=ZERO */
                        _f0 = F_ZERO;
			l_sdini(&_f0, acc);
			l_sdmul(&b0, &b0, acc);
                        _f0 = -a0;
			l_sdmul(&_f0, &c0, acc);
			l_sdsto(acc, &dd);
			s = sqrt(fabs(dd));

			if (dd <= F_ZERO) {
				/*
				 * COINCIDENT OR COMPLEX ROOTS (D .LE. ZERO).
				 */
				*zlg = imsl_cf_convert(b0 / a0, fabs(s / a0));
				*zsm = imsl_c_conjg(*zlg);
				goto L_9000;
			}
			/* DISTINCT REAL ROOTS (D .GT. ZERO). */
			b1 = sign(s, b0) + b0;
		}
		zs = imsl_cf_convert(c0 / b1, F_ZERO);
	} else {
		zs = imsl_cf_convert(F_ZERO, F_ZERO);
	}
	zl = imsl_cf_convert(b1 / a0, F_ZERO);
	if (fabs(imsl_fc_convert(zl)) < fabs(imsl_fc_convert(zs)))
		zs = imsl_c_neg(zl);
	*zlg = zl;
	*zsm = zs;

L_9000:
	imsl_e1pop("l_z3lrc");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z3ORC/DZ3ORC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus called by subroutine ZPORC to compute up to L2
                fixed shift k-polynomials and tests for convergence.

    Usage:      CALL Z3ORC (L2, NZ)

    Arguments:
       L2     - Limit of fixed shift steps.  (Input)
       NZ     - Number of zeros found.  (Output)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z3orc(Mint *l2, Mint *nz)
#else
static void l_z3orc(l2, nz)
	Mint            *l2, *nz;
#endif
{
#ifdef COMPUTER_DECOSF
	int              spass, stry, vpass, vtry;
#else
	long             spass, stry, vpass, vtry;
#endif
	Mint             iflag, itype, j;
	Mfloat           betas, betav, oss, ots, otv, ovv, s, ss, svu, svv,
	                ts, tss, tv, tvv, ui, vi, vv;
#ifdef COMPUTER_DECOSF
        int              IMSLFALSE = 0, IMSLTRUE = 1;
#else
        long             IMSLFALSE = 0, IMSLTRUE = 1;
#endif

	*nz = 0;
	betav = .25;
	betas = .25;
	oss = lv_z11rc.sr;
	ovv = lv_z11rc.v;
	/*
	 * EVALUATE POLYNOMIAL BY SYNTHETIC DIVISION
	 */
	l_z9orc(&lv_z11rc.nn, &lv_z11rc.u, &lv_z11rc.v, lv_z11rc.p, lv_z11rc.qp, &lv_z11rc.ra,
		   &lv_z11rc.rb);
	l_z6orc(&itype);
	for (j = 1; j <= *l2; j++) {
		/*
		 * CALCULATE NEXT K POLYNOMIAL AND ESTIMATE V
		 */
		l_z7orc(&itype);
		l_z6orc(&itype);
		l_z8orc(&itype, &ui, &vi);
		vv = vi;
		/* ESTIMATE S */
		ss = F_ZERO;
		if (lv_z11rc.rk[lv_z11rc.n - 1] != F_ZERO)
			ss = -lv_z11rc.p[lv_z11rc.nn - 1] / lv_z11rc.rk[lv_z11rc.n - 1];
		tv = F_ONE;
		ts = F_ONE;
		if (j == 1 || itype == 3)
			goto L_40;
		/*
		 * COMPUTE RELATIVE MEASURES OF CONVERGENCE OF S AND V
		 * SEQUENCES
		 */
		if (vv != F_ZERO)
			tv = fabs((vv - ovv) / vv);
		if (ss != F_ZERO)
			ts = fabs((ss - oss) / ss);
		/*
		 * IF DECREASING, MULTIPLY TWO MOST RECENT CONVERGENCE
		 * MEASURES
		 */
		tvv = F_ONE;
		if (tv < otv)
			tvv = tv * otv;
		tss = F_ONE;
		if (ts < ots)
			tss = ts * ots;
		/* COMPARE WITH CONVERGENCE CRITERIA */
		vpass = tvv < betav;
		spass = tss < betas;
		if (!(spass || vpass))
			goto L_40;
		/*
		 * AT LEAST ONE SEQUENCE HAS PASSED THE CONVERGENCE TEST.
		 * STORE VARIABLES BEFORE ITERATING
		 */
		svu = lv_z11rc.u;
		svv = lv_z11rc.v;
		scopy(lv_z11rc.n, lv_z11rc.rk, 1, lv_z11rc.svk, 1);
		s = ss;
		/*
		 * CHOOSE ITERATION ACCORDING TO THE FASTEST CONVERGING
		 * SEQUENCE
		 */
		vtry = IMSLFALSE;
		stry = IMSLFALSE;
		if (spass && ((!vpass) || tss < tvv))
			goto L_20;
L_10:
		l_z4orc(&ui, &vi, nz);
		if (*nz > 0)
			return;
		/*
		 * QUADRATIC ITERATION HAS FAILED. FLAG THAT IT HAS BEEN
		 * TRIED AND DECREASE THE CONVERGENCE CRITERION.
		 */
		vtry = IMSLTRUE;
		betav *= .25;
		/*
		 * TRY LINEAR ITERATION IF IT HAS NOT BEEN TRIED AND THE S
		 * SEQUENCE IS CONVERGING
		 */
		if (stry || (!spass))
			goto L_30;
		scopy(lv_z11rc.n, lv_z11rc.svk, 1, lv_z11rc.rk, 1);
L_20:
		l_z5orc(&s, nz, &iflag);
		if (*nz > 0)
			return;
		/*
		 * LINEAR ITERATION HAS FAILED. FLAG THAT IT HAS BEEN TRIED
		 * AND DECREASE THE CONVERGENCE CRITERION
		 */
		stry = IMSLTRUE;
		betas *= .25;
		if (iflag == 0)
			goto L_30;
		/*
		 * IF LINEAR ITERATION SIGNALS AN ALMOST DOUBLE REAL ZERO
		 * ATTEMPT QUADRATIC INTERATION
		 */
		ui = -(s + s);
		vi = s * s;
		goto L_10;
		/* RESTORE VARIABLES */
L_30:
		lv_z11rc.u = svu;
		lv_z11rc.v = svv;
		scopy(lv_z11rc.n, lv_z11rc.svk, 1, lv_z11rc.rk, 1);
		/*
		 * TRY QUADRATIC ITERATION IF IT HAS NOT BEEN TRIED AND THE V
		 * SEQUENCE IS CONVERGING
		 */
		if (vpass && (!vtry))
			goto L_10;
		/*
		 * RECOMPUTE QP AND SCALAR VALUES TO CONTINUE THE SECOND
		 * STAGE
		 */
		l_z9orc(&lv_z11rc.nn, &lv_z11rc.u, &lv_z11rc.v, lv_z11rc.p, lv_z11rc.qp, &lv_z11rc.ra,
			   &lv_z11rc.rb);
		l_z6orc(&itype);
L_40:
		ovv = vv;
		oss = ss;
		otv = tv;
		ots = ts;
	}

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z4ORC/DZ4ORC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus subroutine called by ZPORC to do a variable-shift
                k-polynomial iteration.

    Usage:      CALL Z4ORC (UU, VV, NZ)

    Arguments:
       UU     - Coefficient of the starting quadratic.  (Input)
       VV     - Coefficient of the starting quadratic.  (Input)
       NZ     - Number of zeros found.  (Output)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z4orc(Mfloat *uu, Mfloat *vv, Mint *nz)
#else
static void l_z4orc(uu, vv, nz)
	Mfloat          *uu, *vv;
	Mint            *nz;
#endif
{
#ifdef COMPUTER_DECOSF
	int              tried;
#else
	long             tried;
#endif
	Mint             i, itype, j;
	Mfloat           _f0, ee, omp, relstp, rmp, t, ui, vi, zm;
	Mdouble          acc[2];
	Mf_complex       rlz, sz;
#ifdef COMPUTER_DECOSF
        int              IMSLFALSE = 0, IMSLTRUE = 1;
#else
        long             IMSLFALSE = 0, IMSLTRUE = 1;
#endif

	*nz = 0;
	tried = IMSLFALSE;
	lv_z11rc.u = *uu;
	lv_z11rc.v = *vv;
	j = 0;
	/* MAIN LOOP */
L_10:
        _f0 = F_ONE;
	l_z10rc(&_f0, &lv_z11rc.u, &lv_z11rc.v, &sz, &rlz);
	lv_z11rc.szr = imsl_fc_convert(sz);
	lv_z11rc.szi = imsl_c_aimag(sz);
	lv_z11rc.rlzr = imsl_fc_convert(rlz);
	lv_z11rc.rlzi = imsl_c_aimag(rlz);
	/*
	 * RETURN IF ROOTS OF THE QUADRATIC ARE REAL AND NOT CLOSE TO
	 * MULTIPLE OR NEARLY EQUAL AND OF OPPOSITE SIGN
	 */
	if (fabs(fabs(lv_z11rc.szr) - fabs(lv_z11rc.rlzr)) > 0.01 * fabs(lv_z11rc.rlzr))
		return;
	/*
	 * EVALUATE POLYNOMIAL BY QUADRATIC SYNTHETIC DIVISION
	 */
	l_z9orc(&lv_z11rc.nn, &lv_z11rc.u, &lv_z11rc.v, lv_z11rc.p, lv_z11rc.qp, &lv_z11rc.ra,
		   &lv_z11rc.rb);
	rmp = fabs(lv_z11rc.ra - lv_z11rc.szr * lv_z11rc.rb) + fabs(lv_z11rc.szi *
							   lv_z11rc.rb);
	/*
	 * COMPUTE A RIGOROUS BOUND ON THE ROUNDING ERROR IN EVALUTING P
	 */
	zm = sqrt(fabs(lv_z11rc.v));
	ee = F_TWO * fabs(lv_z11rc.qp[0]);
	t = -lv_z11rc.szr * lv_z11rc.rb;
	for (i = 2; i <= lv_z11rc.n; i++) {
		ee = ee * zm + fabs(lv_z11rc.qp[i - 1]);
	}
	ee = ee * zm + fabs(lv_z11rc.ra + t);
	ee = (F_FIVE * lv_z11rc.rmre + F_FOUR * lv_z11rc.are) * ee - (F_FIVE * lv_z11rc.rmre + F_TWO *
	     lv_z11rc.are) * (fabs(lv_z11rc.ra + t) + fabs(lv_z11rc.rb) * zm) + F_TWO *
		lv_z11rc.are * fabs(t);
	/*
	 * ITERATION HAS CONVERGED SUFFICIENTLY IF THE POLYNOMIAL VALUE IS
	 * LESS THAN 20 TIMES THIS BOUND
	 */
	if (rmp > 20.0 * ee)
		goto L_30;
	*nz = 2;
	return;
L_30:
	j += 1;
	/* STOP ITERATION AFTER 20 STEPS */
	if (j > 20)
		return;
	if (j >= 2) {
		if ((relstp <= 0.01 && rmp >= omp) && !tried) {
			/*
			 * A CLUSTER APPEARS TO BE STALLING THE CONVERGENCE.
			 * FIVE FIXED SHIFT STEPS ARE TAKEN WITH A U,V CLOSE
			 * TO THE CLUSTER
			 */
			if (relstp < lv_z11rc.eta)
				relstp = lv_z11rc.eta;
			relstp = sqrt(relstp);
			/* U = U - U*RELSTP */
			l_sdini(&lv_z11rc.u, acc);
                        _f0 = -lv_z11rc.u;
			l_sdmul(&_f0, &relstp, acc);
			l_sdsto(acc, &lv_z11rc.u);
			/* V = V + V*RELSTP */
			l_sdini(&lv_z11rc.v, acc);
			l_sdmul(&lv_z11rc.v, &relstp, acc);
			l_sdsto(acc, &lv_z11rc.v);
			l_z9orc(&lv_z11rc.nn, &lv_z11rc.u, &lv_z11rc.v, lv_z11rc.p, lv_z11rc.qp,
				   &lv_z11rc.ra, &lv_z11rc.rb);
			for (i = 1; i <= 5; i++) {
				l_z6orc(&itype);
				l_z7orc(&itype);
			}
			tried = IMSLTRUE;
			j = 0;
		}
	}
	omp = rmp;
	/*
	 * CALCULATE NEXT K POLYNOMIAL AND NEW U AND V
	 */
	l_z6orc(&itype);
	l_z7orc(&itype);
	l_z6orc(&itype);
	l_z8orc(&itype, &ui, &vi);
	/*
	 * IF VI IS ZERO THE ITERATION IS NOT CONVERGING
	 */
	if (vi == F_ZERO)
		return;
	relstp = fabs((vi - lv_z11rc.v) / vi);
	lv_z11rc.u = ui;
	lv_z11rc.v = vi;
	goto L_10;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z5ORC/DZ5ORC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus subroutine called by ZPORC to do a variable-shift
                h-polynomial iteration.

    Usage:      CALL Z5ORC (SSS, NZ, IFLAG)

    Arguments:
       SSS    - Starting iterate.  (Input)
       NZ     - Number of zeros found.  (Output)
       IFLAG  - Flag to indicate a pair of zeros near the real axis.
                   (Output)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z5orc(Mfloat *sss, Mint *nz, Mint *iflag)
#else
static void l_z5orc(sss, nz, iflag)
	Mfloat          *sss;
	Mint            *nz, *iflag;
#endif
{
	Mint             i, j;
	Mfloat           ee, omp, pv, rkv, rmp, rms, s, t;
	Mdouble          acc[2];


	*nz = 0;
	s = *sss;
	*iflag = 0;
	j = 0;
	/* MAIN LOOP */
L_10:
	pv = lv_z11rc.p[0];
	/* EVALUATE P AT S */
	lv_z11rc.qp[0] = pv;
	/* PV = PV*S + P(I), FOR I = 2, NN */
	for (i = 2; i <= lv_z11rc.nn; i++) {
		l_sdini(&lv_z11rc.p[i - 1], acc);
		l_sdmul(&pv, &s, acc);
		l_sdsto(acc, &pv);
		lv_z11rc.qp[i - 1] = pv;
	}
	rmp = fabs(pv);
	/*
	 * COMPUTE A RIGOROUS BOUND ON THE ERROR IN EVALUATING P
	 */
	rms = fabs(s);
	ee = (lv_z11rc.rmre / (lv_z11rc.are + lv_z11rc.rmre)) * fabs(lv_z11rc.qp[0]);
	for (i = 2; i <= lv_z11rc.nn; i++) {
		ee = ee * rms + fabs(lv_z11rc.qp[i - 1]);
	}
	/*
	 * ITERATION HAS CONVERGED SUFFICIENTLY IF THE POLYNOMIAL VALUE IS
	 * LESS THAN 20 TIMES THIS BOUND
	 */
	if (rmp <= 20.0 * ((lv_z11rc.are + lv_z11rc.rmre) * ee - lv_z11rc.rmre * rmp)) {
		*nz = 1;
		lv_z11rc.szr = s;
		lv_z11rc.szi = F_ZERO;
		goto L_9000;
	}
	j += 1;
	/* STOP ITERATION AFTER 10 STEPS */
	if (j > 10)
		goto L_9000;
	if (j >= 2) {
		if (fabs(t) <= 0.001 * fabs(s - t) && rmp > omp) {
			/*
			 * A CLUSTER OF ZEROS NEAR THE REAL AXIS HAS BEEN
			 * ENCOUNTERED RETURN WITH IFLAG SET TO INITIATE A
			 * QUADRATIC ITERATION
			 */
			*iflag = 1;
			*sss = s;
			goto L_9000;
		}
	}
	/*
	 * RETURN IF THE POLYNOMIAL VALUE HAS INCREASED SIGNIFICANTLY
	 */
	omp = rmp;
	/*
	 * COMPUTE T, THE NEXT POLYNOMIAL, AND THE NEW ITERATE
	 */
	rkv = lv_z11rc.rk[0];
	lv_z11rc.qk[0] = rkv;
	/* RKV = RKV*S + RK(I), FOR I = 2, N */
	for (i = 2; i <= lv_z11rc.n; i++) {
		l_sdini(&lv_z11rc.rk[i - 1], acc);
		l_sdmul(&rkv, &s, acc);
		l_sdsto(acc, &rkv);
		lv_z11rc.qk[i - 1] = rkv;
	}
	if (fabs(rkv) > fabs(lv_z11rc.rk[lv_z11rc.n - 1]) * F_TEN * lv_z11rc.eta) {
		/*
		 * USE THE SCALED FORM OF THE RECURRENCE IF THE VALUE OF K AT
		 * S IS NONZERO
		 */
		t = -pv / rkv;
		lv_z11rc.rk[0] = lv_z11rc.qp[0];
		/* RK(I) = T*QK(I-1) + QP(I), FOR I=2,N */
		for (i = 2; i <= lv_z11rc.n; i++) {
			l_sdini(&lv_z11rc.qp[i - 1], acc);
			l_sdmul(&t, &lv_z11rc.qk[i - 2], acc);
			l_sdsto(acc, &lv_z11rc.rk[i - 1]);
		}
	} else {
		/* USE UNSCALED FORM */
		lv_z11rc.rk[0] = F_ZERO;
		scopy(lv_z11rc.n - 1, lv_z11rc.qk, 1, &lv_z11rc.rk[1], 1);
	}

	rkv = lv_z11rc.rk[0];
	/* RKV = RKV*S + RK(I), FOR I = 2, N */
	for (i = 2; i <= lv_z11rc.n; i++) {
		l_sdini(&lv_z11rc.rk[i - 1], acc);
		l_sdmul(&rkv, &s, acc);
		l_sdsto(acc, &rkv);
	}
	t = F_ZERO;
	if (fabs(rkv) > fabs(lv_z11rc.rk[lv_z11rc.n - 1]) * F_TEN * lv_z11rc.eta)
		t = -pv / rkv;
	/* S = S + T */
	l_sdini(&t, acc);
	l_sdadd(&s, acc);
	l_sdsto(acc, &s);
	goto L_10;

L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z6ORC/DZ6ORC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus subroutine called by ZPORC to calculate scalar
                quantities for the next k-polynomial and new estimates
                of the quadratic coefficients.

    Usage:      CALL Z6ORC (ITYPE)

    Arguments:
       ITYPE  - Indicator to determine how the calculations are normal-
                ized to avoid overflow.  (Output)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
static void l_z6orc(Mint *itype)
#else
static void l_z6orc(itype)
	Mint            *itype;
#endif
{
	Mfloat           _f0;
	Mdouble          acc[2];


	/*
	 * SYNTHETIC DIVISION OF K BY THE QUADRATIC 1,U,V
	 */
	l_z9orc(&lv_z11rc.n, &lv_z11rc.u, &lv_z11rc.v, lv_z11rc.rk, lv_z11rc.qk, &lv_z11rc.c,
		   &lv_z11rc.d);
	if (fabs(lv_z11rc.c) <= fabs(lv_z11rc.rk[lv_z11rc.n - 1]) * 100.0 * lv_z11rc.eta) {
		if (fabs(lv_z11rc.d) <= fabs(lv_z11rc.rk[lv_z11rc.n - 2]) * 100.0 *
		    lv_z11rc.eta) {
			*itype = 3;
			/*
			 * TYPE=3 INDICATES THE QUADRATIC IS ALMOST A FACTOR
			 * OF K
			 */
			goto L_9000;
		}
	}
	if (fabs(lv_z11rc.d) >= fabs(lv_z11rc.c)) {
		*itype = 2;
		/*
		 * TYPE=2 INDICATES THAT ALL FORMULAS ARE DIVIDED BY D
		 */
		lv_z11rc.e = lv_z11rc.ra / lv_z11rc.d;
		lv_z11rc.f = lv_z11rc.c / lv_z11rc.d;
		lv_z11rc.g = lv_z11rc.u * lv_z11rc.rb;
		lv_z11rc.h = lv_z11rc.v * lv_z11rc.rb;
		/* A3 = (RA+G)*E + H*(RB/D) */
                _f0 = F_ZERO;
		l_sdini(&_f0, acc);
		lv_z11rc.a3 = lv_z11rc.rb / lv_z11rc.d;
		l_sdmul(&lv_z11rc.h, &lv_z11rc.a3, acc);
		l_sdmul(&lv_z11rc.g, &lv_z11rc.e, acc);
		l_sdmul(&lv_z11rc.ra, &lv_z11rc.e, acc);
		l_sdsto(acc, &lv_z11rc.a3);
		/* A1 = RB*F - RA */
                _f0 = -lv_z11rc.ra;
		l_sdini(&_f0, acc);
		l_sdmul(&lv_z11rc.rb, &lv_z11rc.f, acc);
		l_sdsto(acc, &lv_z11rc.a1);
		/* A7 = (F+U)*RA + H */
		l_sdini(&lv_z11rc.h, acc);
		l_sdmul(&lv_z11rc.f, &lv_z11rc.ra, acc);
		l_sdmul(&lv_z11rc.u, &lv_z11rc.ra, acc);
		l_sdsto(acc, &lv_z11rc.a7);
		goto L_9000;
	}
	*itype = 1;
	/*
	 * TYPE=1 INDICATES THAT ALL FORMULAS ARE DIVIDED BY C
	 */
	lv_z11rc.e = lv_z11rc.ra / lv_z11rc.c;
	lv_z11rc.f = lv_z11rc.d / lv_z11rc.c;
	lv_z11rc.g = lv_z11rc.u * lv_z11rc.e;
	lv_z11rc.h = lv_z11rc.v * lv_z11rc.rb;
	/* A3 = RA*E + (H/C+G)*RB */
	lv_z11rc.a3 = lv_z11rc.h / lv_z11rc.c;
        _f0 = F_ZERO;
	l_sdini(&_f0, acc);
	l_sdmul(&lv_z11rc.rb, &lv_z11rc.a3, acc);
	l_sdmul(&lv_z11rc.rb, &lv_z11rc.g, acc);
	l_sdmul(&lv_z11rc.ra, &lv_z11rc.e, acc);
	l_sdsto(acc, &lv_z11rc.a3);
	/* A1 = RB - RA*(D/C) */
	lv_z11rc.a1 = lv_z11rc.d / lv_z11rc.c;
	l_sdini(&lv_z11rc.rb, acc);
        _f0 = -lv_z11rc.ra;
	l_sdmul(&_f0, &lv_z11rc.a1, acc);
	l_sdsto(acc, &lv_z11rc.a1);
	/* A7 = RA + G*D + H*F */
	l_sdini(&lv_z11rc.ra, acc);
	l_sdmul(&lv_z11rc.g, &lv_z11rc.d, acc);
	l_sdmul(&lv_z11rc.h, &lv_z11rc.f, acc);
	l_sdsto(acc, &lv_z11rc.a7);

L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z7ORC/DZ7ORC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus subroutine called by ZPORC to compute the next
                k-polynomials using scalars from Z5ORC.

    Usage:      CALL Z7ORC (ITYPE)

    Arguments:
       ITYPE  - Indicator to determine how the calculations are normal-
                ized to avoid overflow.  (Input)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z7orc(Mint *itype)
#else
static void l_z7orc(itype)
	Mint            *itype;
#endif
{
	Mint             i;
	Mfloat           _f0, temp;
	Mdouble          acc[2];


	if (*itype != 3) {
		temp = lv_z11rc.ra;
		if (*itype == 1)
			temp = lv_z11rc.rb;
		/*
		 * IF A1 IS NEARLY ZERO THEN USE A SPECIAL FORM OF THE
		 * RECURRENCE
		 */
		lv_z11rc.rk[0] = F_ZERO;
		if (fabs(lv_z11rc.a1) <= fabs(temp) * lv_z11rc.eta * F_TEN) {
			lv_z11rc.rk[1] = -lv_z11rc.a7 * lv_z11rc.qp[0];
			/*
			 * RK(I) = A3*QK(I-2) - A7*QP(I-1) FOR I = 3, N
			 */
			for (i = 3; i <= lv_z11rc.n; i++) {
				_f0 = F_ZERO;
				l_sdini(&_f0, acc);
				l_sdmul(&lv_z11rc.a3, &lv_z11rc.qk[i - 3], acc);
				_f0 = -lv_z11rc.a7;
				l_sdmul(&_f0, &lv_z11rc.qp[i - 2], acc);
				l_sdsto(acc, &lv_z11rc.rk[i - 1]);
			}
			goto L_9000;
		}
		/* USE SCALED FORM OF THE RECURRENCE */
		lv_z11rc.a7 /= lv_z11rc.a1;
		lv_z11rc.a3 /= lv_z11rc.a1;
		lv_z11rc.rk[0] = lv_z11rc.qp[0];
		/* RK(2) = QP(2) - A7*QP(1) */
		l_sdini(&lv_z11rc.qp[1], acc);
		_f0 = -lv_z11rc.a7;
		l_sdmul(&_f0, &lv_z11rc.qp[0], acc);
		l_sdsto(acc, &lv_z11rc.rk[1]);
		/*
		 * RK(I) = A3*QK(I-2) - A7*QP(I-1) + QP(I) FOR I = 3, N
		 */
		for (i = 3; i <= lv_z11rc.n; i++) {
			l_sdini(&lv_z11rc.qp[i - 1], acc);
			_f0 = -lv_z11rc.a7;
			l_sdmul(&_f0, &lv_z11rc.qp[i - 2], acc);
			l_sdmul(&lv_z11rc.a3, &lv_z11rc.qk[i - 3], acc);
			l_sdsto(acc, &lv_z11rc.rk[i - 1]);
		}

		goto L_9000;
	}
	/*
	 * USE UNSCALED FORM OF THE RECURRENCE IF TYPE IS 3
	 */
	lv_z11rc.rk[0] = F_ZERO;
	lv_z11rc.rk[1] = F_ZERO;
	scopy(lv_z11rc.n - 2, lv_z11rc.qk, 1, &lv_z11rc.rk[2], 1);

L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z8ORC/DZ8ORC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus subroutine called by ZPORC to compute new
                estimates of the quadratic coefficients.

    Usage:      CALL Z8ORC (ITYPE, UU, VV)

    Arguments:
       ITYPE  - Indicator to determine how the calculations are normal-
                ized to avoid overflow.  (Input)
       UU     - Coefficient of the starting quadratic.  (Output)
       VV     - Coefficient of the starting quadratic.  (Output)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z8orc(Mint *itype, Mfloat *uu, Mfloat *vv)
#else
static void l_z8orc(itype, uu, vv)
	Mint            *itype;
	Mfloat          *uu, *vv;
#endif
{
	Mfloat           _f0, a4, a5, atemp, b1, b2, c1, c2, c3, c4, temp;
	Mdouble          acc[2];


	if (*itype != 3) {
		if (*itype != 2) {
			/* A4 = RA + U*RB + H*F */
			l_sdini(&lv_z11rc.ra, acc);
			l_sdmul(&lv_z11rc.u, &lv_z11rc.rb, acc);
			l_sdmul(&lv_z11rc.h, &lv_z11rc.f, acc);
			l_sdsto(acc, &a4);
			/* A5 = C + (U+V*F)*D */
			l_sdini(&lv_z11rc.u, acc);
			l_sdmul(&lv_z11rc.v, &lv_z11rc.f, acc);
			l_sdsto(acc, &atemp);
			l_sdini(&lv_z11rc.c, acc);
			l_sdmul(&lv_z11rc.d, &atemp, acc);
			l_sdsto(acc, &a5);
		} else {
			/* A4 = (RA+G)*F + H */
			l_sdini(&lv_z11rc.g, acc);
			l_sdadd(&lv_z11rc.ra, acc);
			l_sdsto(acc, &atemp);
			l_sdini(&lv_z11rc.h, acc);
			l_sdmul(&lv_z11rc.f, &atemp, acc);
			l_sdsto(acc, &a4);
			/* A5 = (F+U)*C + V*D */
			l_sdini(&lv_z11rc.u, acc);
			l_sdadd(&lv_z11rc.f, acc);
			l_sdsto(acc, &atemp);
                        _f0 = F_ZERO;
			l_sdini(&_f0, acc);
			l_sdmul(&lv_z11rc.c, &atemp, acc);
			l_sdmul(&lv_z11rc.v, &lv_z11rc.d, acc);
			l_sdsto(acc, &a5);
		}
		/* EVALUATE NEW QUADRATIC COEFFICIENTS. */
		b1 = -lv_z11rc.rk[lv_z11rc.n - 1] / lv_z11rc.p[lv_z11rc.nn - 1];
		/* B2 = -(RK(N-1)+B1*P(N))/P(NN) */
		l_sdini(&lv_z11rc.rk[lv_z11rc.n - 2], acc);
		l_sdmul(&b1, &lv_z11rc.p[lv_z11rc.n - 1], acc);
		l_sdsto(acc, &b2);
		b2 = -b2 / lv_z11rc.p[lv_z11rc.nn - 1];
		c1 = lv_z11rc.v * b2 * lv_z11rc.a1;
		c2 = b1 * lv_z11rc.a7;
		c3 = b1 * b1 * lv_z11rc.a3;
		/* C4 = C1-C2-C3 */
		l_sdini(&c1, acc);
                _f0 = -c2;
		l_sdadd(&_f0, acc);
                _f0 = -c3;
		l_sdadd(&_f0, acc);
		l_sdsto(acc, &c4);
		/* TEMP = A5 + B1*A4 - C4 */
		l_sdini(&a5, acc);
		l_sdmul(&b1, &a4, acc);
                _f0 = -c4;
		l_sdadd(&_f0, acc);
		l_sdsto(acc, &temp);

		if (F_ONE * temp != F_ZERO) {
			/*
			 * UU = U - (U*(C3+C2)+V*(B1*A1+B2*A7)) /TEMP
			 */
			l_sdini(&c2, acc);
			l_sdadd(&c3, acc);
			l_sdsto(acc, &atemp);
                        _f0 = F_ZERO;
			l_sdini(&_f0, acc);
			l_sdmul(&lv_z11rc.u, &atemp, acc);
			l_sdsto(acc, &atemp);
                        _f0 = F_ZERO;
			l_sdini(&_f0, acc);
			l_sdmul(&b2, &lv_z11rc.a7, acc);
			l_sdmul(&b1, &lv_z11rc.a1, acc);
			l_sdsto(acc, uu);
			l_sdini(&atemp, acc);
			l_sdmul(uu, &lv_z11rc.v, acc);
			l_sdsto(acc, uu);
			*uu /= temp;
                        _f0 = -*uu;
			l_sdini(&_f0, acc);
			l_sdadd(&lv_z11rc.u, acc);
			l_sdsto(acc, uu);
			/* VV = V*(1+C4/TEMP) */
			*vv = c4 / temp;
			l_sdini(&lv_z11rc.v, acc);
			l_sdmul(&lv_z11rc.v, vv, acc);
			l_sdsto(acc, vv);
			goto L_9000;
		}
	}
	/* IF TYPE=3 THE QUADRATIC IS ZEROED */
	*uu = F_ZERO;
	*vv = F_ZERO;

L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z9ORC/DZ9ORC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus subroutine called by ZPORC to divide the poly-
                nomial by the quadratic 1, U, V.

    Usage:      CALL Z9ORC (NN, U, V, P, Q, RA, RB)

    Arguments:
       NN     - Number of coefficients in the polynomial.  (Input)
       U      - Coefficient of the starting quadratic.  (Input)
       V      - Coefficient of the starting quadratic.  (Input)
       P      - Vector of length NN containing the coefficients of the
                polynomial.  (Input)
       Q      - Vector of length NN containing the quotient.  (Output)
       RA     - Scalar containing the first coefficient of the remainder.
                   (Output)
       RB     - Scalar containing the second coefficient of the remain-
                der.  (Output)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z9orc(Mint *nn, Mfloat *u, Mfloat *v, Mfloat p[], Mfloat q[],
                    Mfloat *ra, Mfloat *rb)
#else
static void l_z9orc(nn, u, v, p, q, ra, rb)
	Mint            *nn;
	Mfloat          *u, *v, p[], q[], *ra, *rb;
#endif
{
	Mint             i;
	Mfloat           _f0, c;
	Mdouble          acc[2];


	*rb = p[0];
	q[0] = *rb;
	/* RA = P(2) - U*RB */
	l_sdini(&p[1], acc);
        _f0 = -*u;
	l_sdmul(&_f0, rb, acc);
	l_sdsto(acc, ra);
	q[1] = *ra;
	/* C = P(I) - U*RA - V*RB FOR I = 3, NN */
	for (i = 3; i <= *nn; i++) {
		l_sdini(&p[i - 1], acc);
		_f0 = -*u;
		l_sdmul(&_f0, ra, acc);
	        _f0 = -*v;
		l_sdmul(&_f0, rb, acc);
		l_sdsto(acc, &c);
		q[i - 1] = c;
		*rb = *ra;
		*ra = c;
	}

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  Z10RC/DZ10RC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Nucleus subroutine called by ZPORC to calculate the zeros
                of a quadratic.

    Usage:      CALL Z10RC (RA, B1, C, S, RL)

    Arguments:
       RA     - Coefficient of the starting quadratic.  (Input)
       B1     - Coefficient of the starting quadratic.  (Input)
       C      - Coefficient of the starting quadratic.  (Input)
       S      - Small zero. (Output)
       RL     - Large zero. (Output)

    Copyright:  1985 by IMSL, Inc. All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_z10rc(Mfloat *ra, Mfloat *b1, Mfloat *c, Mf_complex *s,
                    Mf_complex *rl)
#else
static void l_z10rc(ra, b1, c, s, rl)
	Mfloat          *ra, *b1, *c;
	Mf_complex      *s, *rl;
#endif
{
	Mfloat           d, e, rb, rlr, imsl_si, sr;


	if (*ra == F_ZERO) {
		*s = imsl_cf_convert(F_ZERO, F_ZERO);
		if (*b1 != F_ZERO)
			*s = imsl_cf_convert(-*c / *b1, F_ZERO);
		*rl = imsl_cf_convert(F_ZERO, F_ZERO);
		goto L_9000;
	}
	if (*c == F_ZERO) {
		*s = imsl_cf_convert(F_ZERO, F_ZERO);
		*rl = imsl_cf_convert(-*b1 / *ra, F_ZERO);
		goto L_9000;
	}
	/*
	 * COMPUTE DISCRIMINANT AVOIDING OVERFLOW
	 */
	rb = *b1 / F_TWO;
	if (fabs(rb) >= fabs(*c)) {
		e = F_ONE - (*ra / rb) * (*c / rb);
		d = sqrt(fabs(e)) * fabs(rb);
	} else {
		e = *ra;
		if (*c < F_ZERO)
			e = -*ra;
		e = rb * (rb / fabs(*c)) - e;
		d = sqrt(fabs(e)) * sqrt(fabs(*c));
	}

	if (e >= F_ZERO) {
		/* REAL ZEROS */
		if (rb >= F_ZERO)
			d = -d;
		rlr = (-rb + d) / *ra;
		sr = F_ZERO;
		if (rlr != F_ZERO)
			sr = (*c / rlr) / *ra;
		*s = imsl_cf_convert(sr, F_ZERO);
		*rl = imsl_cf_convert(rlr, F_ZERO);
		goto L_9000;
	}
	/* COMPLEX CONJUGATE ZEROS */
	sr = -rb / *ra;
	imsl_si = fabs(d / *ra);
	*s = imsl_cf_convert(sr, imsl_si);
	*rl = imsl_cf_convert(sr, -imsl_si);

L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SDADD (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Double precision add.

    Usage:      CALL SDADD (A, ACC)

    Arguments:
       A      - Single precision number to be added to the accumulator.
                (Input)
       ACC    - Accumulator. (Input/Output)
                ACC is a double precision vector of length 2. On output,
                ACC contains the sum of C input ACC and A.

    Remark:
       SDADD adds the single precision number A to the double
       precision accumulator, ACC. The subroutine assumes that a
       double precision number is already in the accumulator.

    GAMS:       A3b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_sdadd(Mfloat *a, Mdouble acc[])
#else
static void l_sdadd(a, acc)
	Mfloat          *a;
	Mdouble          acc[];
#endif
{


	acc[0] += (double)*a;
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SDINI (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Single precision initialization.

    Usage:      CALL SDINI (S, ACC)

    Arguments:
       S      - Single precision scalar. (Input)
                On Input, S contains a single precision value that the
                accumulator is to be initialized to.
       ACC    - Accumulator. (Output)
                ACC is a double precision vector of length 2.  ACC is
                initialized to S.

    GAMS:       A3b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_sdini(Mfloat *s, Mdouble acc[])
#else
static void l_sdini(s, acc)
	Mfloat          *s;
	Mdouble          acc[];
#endif
{


	acc[0] = (double)*s;
	acc[1] = F_ZERO;

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SDMUL (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Double precision multiply.

    Usage:      CALL SDMUL (A, B, ACC)

    Arguments:
       A      - Single precision multiplier. (Input)
       B      - Single precision multiplicand. (Input)
       ACC    - Accumulator. (Input/Output)
                ACC is a double precision vector of length 2.
                On output, ACC contains the sum of input ACC and
                A*B.

    Remark:
       SDMUL adds the product A*B to the double precision
       accumulator, ACC. The subroutine assumes that a double
       precision number is already in the accumulator.

    GAMS:       A3b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_sdmul(Mfloat *a, Mfloat *b, Mdouble acc[])
#else
static void l_sdmul(a, b, acc)
	Mfloat          *a, *b;
	Mdouble          acc[];
#endif
{


	acc[0] += ((double)*a) * ((double)*b);
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  SDSTO (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Single precision store.

    Usage:      CALL SDSTO (ACC, S)

    Arguments:
       ACC    - Accumulator. (Input)
                ACC is a double precision vector of length 2. ACC is
                assumed to be the result computed by calling IMSL
                double precision routines.
       S      - Single precision scalar. (Output)
                On output, S contains a single precision approximation
                to the value of the double precision accumulator.

    GAMS:       A3b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_sdsto(Mdouble acc[], Mfloat *s)
#else
static void l_sdsto(acc, s)
	Mdouble          acc[];
	Mfloat          *s;
#endif
{


	*s = acc[0];

	return;
}				/* end of function */
/* -----------------------------------------------------------------------
    IMSL Name:  SSWAP (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Interchange vectors X and Y, both single precision.

    Usage:      CALL SSWAP (N, SX, INCX, SY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).
                (Input/Output)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be
                   SX(1+(I-1)*INCX) if INCX.GE.0  or
                   SX(1+(I-N)*INCX) if INCX.LT.0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).
                (Input/Output)
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be
                   SY(1+(I-1)*INCY) if INCY.GE.0  or
                   SY(1+(I-N)*INCY) if INCY.LT.0.

    Keywords:   Level 1 BLAS; SSWAP; Swap; Exchange

    GAMS:       D1a5

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_sswap(Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy)
#else
static void l_sswap(n, sx, incx, sy, incy)
	Mint             n;
	Mfloat           sx[];
	Mint             incx;
	Mfloat           sy[];
	Mint             incy;
#endif
{
	Mint             i, ix, iy;
	Mfloat           stemp;


	if (n > 0) {
		if (incx != 1 || incy != 1) {
			/*
			 * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
			 * NOT EQUAL TO 1
			 */
			ix = 1;
			iy = 1;
			if (incx < 0)
				ix = (-n + 1) * incx + 1;
			if (incy < 0)
				iy = (-n + 1) * incy + 1;
			for (i = 1; i <= n; i++) {
				stemp = sx[ix - 1];
				sx[ix - 1] = sy[iy - 1];
				sy[iy - 1] = stemp;
				ix += incx;
				iy += incy;
			}
		} else {
			for (i = 1; i <= n; i++) {
				stemp = sx[i - 1];
				sx[i - 1] = sy[i - 1];
				sy[i - 1] = stemp;
			}
		}
	}
	return;
}				/* end of function */
