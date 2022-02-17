#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#undef	TRUE

#undef 	FALSE



#ifdef ANSI

static VA_LIST_HACK l_bessel_Ix(Mfloat xnu, Mf_complex z,

                Mint n, va_list argptr);

static VA_LIST_HACK l_bessel_Ix_adr(Mfloat *xnu, Mf_complex *z,

                Mint n, va_list argptr);

static void l_cbis(Mfloat *, Mf_complex *, Mint *, Mf_complex[]);

void imsl_c3is(Mf_complex*, Mfloat*, Mint*, Mf_complex[],

		Mf_complex[], Mf_complex[], Mf_complex[],

		Mint*);

static Mf_complex l_c4is(Mf_complex*, Mfloat*, Mfloat*, Mfloat*,

		Mint*, Mint*, Mf_complex*, Mint*);

#else

static VA_LIST_HACK l_bessel_Ix();

static VA_LIST_HACK l_bessel_Ix_adr();

static void l_cbis();

void imsl_c3is();

static Mf_complex l_c4is();

#endif



static Mf_complex *lv_cbs = NULL;

#ifdef ANSI

Mf_complex *imsl_c_bessel_Ix(Mfloat xnu, Mf_complex z,

		Mint n, ...)

#else

Mf_complex *imsl_c_bessel_Ix(xnu, z, n, va_alist)

    Mfloat 	xnu;

    Mf_complex 	z;

    Mint        n;

    va_dcl

#endif

{

    va_list     argptr;



    VA_START(argptr, n);



    E1PSH("imsl_c_bessel_Ix", "imsl_z_bessel_Ix");



    lv_cbs = NULL;

    IMSL_CALL(l_bessel_Ix(xnu, z, n, argptr));

    va_end(argptr);



    E1POP("imsl_c_bessel_Ix", "imsl_z_bessel_Ix");



    return(lv_cbs);

}



#ifdef ANSI

Mf_complex *imsl_c_bessel_Ix_adr(Mfloat *xnu, Mf_complex *z,

		Mint n, ...)

#else

Mf_complex *imsl_c_bessel_Ix_adr(xnu, z, n, va_alist)

    Mfloat 	*xnu;

    Mf_complex 	*z;

    Mint        n;

    va_dcl

#endif

{

    va_list     argptr;



    VA_START(argptr, n);



    E1PSH("imsl_c_bessel_Ix_adr", "imsl_z_bessel_Ix_adr");



    lv_cbs = NULL;

    IMSL_CALL(l_bessel_Ix_adr(xnu, z, n, argptr));

    va_end(argptr);



    E1POP("imsl_c_bessel_Ix_adr", "imsl_z_bessel_Ix_adr");



    return(lv_cbs);

}



#ifdef ANSI

static VA_LIST_HACK l_bessel_Ix(Mfloat xnu, Mf_complex z,

		Mint n, va_list argptr)

#else

static VA_LIST_HACK l_bessel_Ix(xnu, z, n, argptr)

    Mfloat 	xnu;

    Mf_complex 	z;

    Mint        n;

    va_list     argptr;

#endif

{

    Mfloat float_xnu;

    Mint        code, user_cbs = 0, arg_number = 3;



    code = va_arg(argptr, Mint);

    if (code == (Mint)IMSL_RETURN_USER) {

        lv_cbs = va_arg(argptr, Mf_complex*);

        user_cbs = 1;

    }

    else if (code != 0){

        imsl_e1sti (1, code);

        imsl_e1sti (2, arg_number);

        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);

        goto RETURN;

    }



    if (n <= 0) {

        imsl_e1sti (1, n);



        imsl_e1mes (5, 1, "The number of function values to be calculated, N = %(i1), must be at least 1.");

	goto RETURN;

    }

    if (!user_cbs) {

        lv_cbs = (Mf_complex *) imsl_malloc (n*sizeof(*lv_cbs));

        if (!lv_cbs){

            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);

            goto RETURN;



        }

    }

    float_xnu = xnu;

    l_cbis (&float_xnu, &z, &n, lv_cbs);

    if (imsl_n1rty(0)>3 && !user_cbs) {

        imsl_free(lv_cbs);

        lv_cbs = NULL;

    }



RETURN:



    return(argptr);

}



#ifdef ANSI

static VA_LIST_HACK l_bessel_Ix_adr(Mfloat *xnu, Mf_complex *z,

		Mint n, va_list argptr)

#else

static VA_LIST_HACK l_bessel_Ix_adr(xnu, z, n, argptr)

    Mfloat 	*xnu;

    Mf_complex 	*z;

    Mint        n;

    va_list     argptr;

#endif

{

    Mfloat float_xnu;

    Mint        code, user_cbs = 0, arg_number = 3;



    code = va_arg(argptr, Mint);

    if (code == (Mint)IMSL_RETURN_USER) {

        lv_cbs = va_arg(argptr, Mf_complex*);

        user_cbs = 1;

    }

    else if (code != 0){

        imsl_e1sti (1, code);

        imsl_e1sti (2, arg_number);

        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);

        goto RETURN;

    }



    if (n <= 0) {

        imsl_e1sti (1, n);



        imsl_e1mes (5, 1, "The number of function values to be calculated, N = %(i1), must be at least 1.");

	goto RETURN;

    }

    if (!user_cbs) {

        lv_cbs = (Mf_complex *) imsl_malloc (n*sizeof(*lv_cbs));

        if (!lv_cbs){

            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);

            goto RETURN;



        }

    }

    float_xnu = *xnu;

    l_cbis (&float_xnu, z, &n, lv_cbs);

    if (imsl_n1rty(0)>3 && !user_cbs) {

        imsl_free(lv_cbs);

        lv_cbs = NULL;

    }



RETURN:



    return(argptr);

}









/*Translated by FOR_C++, v0.1, on 06/06/91 at 15:21:52 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 06/06/91 at 15:21:50

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  CBIS/DCBIS (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    May 1, 1991



    Purpose:    Evaluate a sequence of Modified Bessel functions of the

                first kind with real order and d_complex arguments.



    Usage:      CALL CBIS (XNU, Z, N, CBS)



    Arguments:

       XNU    - Real argument which is the lowest order desired.

                (Input)

                XNU must be greater than -1/2.

       Z      - Complex argument for which the sequence of Bessel

                functions is to be evaluated.  (Input)

       N      - Number of elements in the sequence.  (Input)

       CBS    - Vector of length N containing the values of the

                function through the series.  (Output)

                CBS(I) contains the value of the Bessel function of

                order XNU+I-1 at Z for I = 1 to N.



    Remark:

       Informational errors

       Type Code

         3   1  One of the continued fractions failed.

         4   2  Only the first several entries in CBS are valid.



    GAMS:       C10b4



    Chapter:    MATH/LIBRARY SPECIAL FUNCTIONS Bessel Functions



    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_cbis(Mfloat*xnu, Mf_complex *z, Mint *n, Mf_complex cbs[])

#else

static void l_cbis (xnu, z, n, cbs)

    Mfloat      *xnu;

    Mf_complex  *z;

    Mint        *n;

    Mf_complex   cbs[];

#endif

{

    Mint         mode;

    Mf_complex   fip[1], fk[1], fkp[1];





    imsl_e1psh ("CBIS ");



    mode = 4;

    imsl_c3is (z, xnu, n, cbs, fk, fip, fkp, &mode);



    imsl_e1pop ("CBIS ");



    return;

}				/* end of function */







/*Translated by FOR_C++, v0.1, on 06/06/91 at 15:39:04 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 06/06/91 at 15:39:00

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  C3IS/DC3IS (Single/Double precision version)



    Computer:   SUN4/SINGLE



    Revised:    January 8, 1990



    Purpose:    Evaluate a sequence of I, K Bessel functions of

                with real order and d_complex arguments using Steed's

                method.



    Usage:      CALL C3IS (ZZ, XNU, NL, FI, FK, FIP, FKP, MODE1)



    Arguments:

       ZZ     - Complex argument for which the sequence of Bessel

                functions is to be evaluated.  (Input)

       XNU    - Real argument which is the lowest order desired.

                (Input)

                The first order XNU must be .GT. -0.5, and

                XNU .LT. 2/ln(abs(Z)) if abs(Z) is much less than 1

                (i.e. not too negative if abs(Z) is much less than 1)

       NL     - Number of elements in the sequence.  (Input)

       FI     - Complex array containing values of the Bessel function

                I(z) starting with order XNU, if MODE1 indicates that

                I is returned.  (Output)

       FK     - Complex array containing values of the Bessel function

                K(z) starting with order XNU, if MODE1 indicates that

                K is returned.  (Output)

       FIP    - Complex array containing values of the derivative of the

                Bessel function I(z) starting with order XNU, if MODE1

                indicates that I' is returned.  (Output)

       FKP    - Complex array containing values of the derivative of the

                Bessel function K(z) starting with order XNU, if MODE1

                indicates that K' is returned.  (Output)

       MODE1  - Flag indicating which series is to be computed.  (Input)

                if abs(MODE1) = 1  get I,K,I',K'

                              = 2      I,K

                              = 3      I,  I'

                              = 4      I

                if MODE1 .LT. 0 then the values returned are scaled by an

                exponential factor (dependent only on Z) to bring nearer

                unity the functions for large abs(Z), small

                abs(XN) .LT. abs(Z)

                So define SCALE = (  0        if MODE1 .GT. 0

                                  (  REAL(Z)  if MODE1 .LT. 0

                        then FI = EXP(-ABS(SCALE)) * I

                            FIP = EXP(-ABS(SCALE)) * I'

                         and FK = EXP(SCALE)       * K

                            FKP = EXP(SCALE)       * K'



    Remarks:

    1. The algorithm is based on the code BESSCC.

       I. J. Thompson          Bristol     JULY    1986



       Original program  RCWFN       in    CPC  8 (1974) 377-395

                      +  RCWFF       in    CPC 11 (1976) 141-142

                      +  COULFG      in    CPC 27 (1982) 147-166

                      +  COULCC      in    CPC 36 (1985) 363-372

       Description of real algorithm in    CPC 21 (1981) 297-314

       Description of d_complex algorithm    JCP 64 (1986) 490-509

       This version written up       in    CPC 47 (1987) 245-257



    2. Results to within 1-2 decimals of 'machine accuracy'.



    3. LIMIT is the imsl_i_max. no. iterations for CF1, CF2 continued fractions

       If XMU .GT. 0.35*abs(Z),  then abs(Z) is limited to LIMIT.



    4. GAM, CSC are the coefficients of the continued fraction form

       of the diagonal Pade approximants for

             ln(Gamma(1+nu))/nu and  1/sin(pi.nu) - 1/(pi.nu) resp.

       The number of terms required is a linear function of the number

       of digits accuracy, i.e. of log(AMACH(4)).

       The given GAM & CSC parameters are sufficient for ACC .GT. 1D-24

       For fixed accuracy worse than this, expressions in the code

       involving ACC may be pre-evaluated, and NGAM & NCSC reduced.



    Chapter:    SFUN/LIBRARY Bessel Functions



    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#define	NCSC	12

#define	NGAM	33



#ifdef ANSI

void imsl_c3is (Mf_complex*zz, Mfloat *xnu, Mint *nl, Mf_complex fi[],

		Mf_complex fk[], Mf_complex fip[], Mf_complex fkp[],

		Mint *mode1)

#else

void imsl_c3is (zz, xnu, nl, fi, fk, fip, fkp, mode1)

    Mf_complex  *zz;

    Mfloat      *xnu;

    Mint        *nl;

    Mf_complex   fi[], fk[], fip[], fkp[];

    Mint        *mode1;

#endif

{

    Mint   ifip;

    Mint         _l0, i, ifail, ist, j, kase, l, l0, m, mode, nfp, nl1,

                npq[2];

    Mfloat       _f0, _f1, absz, c1, c2, ek, g, g1, pi, rk, scale, x2m,

                xll, xm1, xmu, xn0, xx, yl0;

    static Mfloat acclog, csc[NCSC], fphmin, fplmin, fpmax, fpmin, gam[NGAM], hugh,

                pi1;

    Mf_complex   alpha, d0, f1, f2, fil, fipl, fk0, fkp0,

                g2, pq1, pq2, rat, sl, sl1, sum[2], x, z, zi, ztemp;

    static Mint  limit = 20000;

    static Mfloat accur = -1.0;

    static Mfloat pi2 = -8.90891020676153735661672e-06;

    static Mint  _aini = 1;

    if (_aini) {

	gam[0] = -5.772156649015328606065121e-01;

	gam[1] = 1.42488688965920173739608e00;

	gam[2] = -9.377115767249094049120554e-01;

	gam[3] = -9.773477272948880533901762e-02;

	gam[4] = 6.130615858357687237828236e-01;

	gam[5] = 1.035548353597846452535594e-01;

	gam[6] = 3.833453676357374951071144e-01;

	gam[7] = 1.676098533180169272436803e-01;

	gam[8] = 3.423153070879739741315024e-01;

	gam[9] = 1.872096384664854899310796e-01;

	gam[10] = 3.044176452057531884181363e-01;

	gam[11] = 2.038057743675276230967228e-01;

	gam[12] = 3.036922239952957869745095e-01;

	gam[13] = 2.111473436553448271828676e-01;

	gam[14] = 2.815492532986102124308247e-01;

	gam[15] = 2.186620399336151796598028e-01;

	gam[16] = 2.887599926903399843284945e-01;

	gam[17] = 2.211365686356764978140098e-01;

	gam[18] = 2.711931229606902430385447e-01;

	gam[19] = 2.278516401771275692711102e-01;

	gam[20] = 2.79828901683609140009025e-01;

	gam[21] = 2.257719734937412432938853e-01;

	gam[22] = 2.670653349467782427924255e-01;

	gam[23] = 2.338192475772592332514537e-01;

	gam[24] = 2.72079712938743098157715e-01;

	gam[25] = 2.297304315464212589631377e-01;

	gam[26] = 2.662328902179656243145202e-01;

	gam[27] = 2.355772901124245270936971e-01;

	gam[28] = 2.665129963814811071666658e-01;

	gam[29] = 2.352423215988759927243035e-01;

	gam[30] = 2.639101271705042366932099e-01;

	gam[31] = 2.349813957292258269500878e-01;

	gam[32] = 2.659470304082729651745135e-01;

	csc[0] = 1.666666666666666666666667e-01;

	csc[1] = -1.166666666666666666666667e-01;

	csc[2] = 1.122448979591836734693878e-02;

	csc[3] = -2.839620696763553906411049e-02;

	csc[4] = 3.881831014317402702157693e-03;

	csc[5] = -1.25669016828591456240659e-02;

	csc[6] = 1.959015840039126317490664e-03;

	csc[7] = -7.058674448331350123193777e-03;

	csc[8] = 1.179797436463512584131754e-03;

	csc[9] = -4.514563690718488905367012e-03;

	csc[10] = 7.88000856804251553256682e-04;

	csc[11] = -3.133992819253743972358894e-03;

	_aini = 0;

    }







#define TRUE	1

#define FALSE	0



#define TIMESI(x)	(imsl_cf_convert( -imsl_c_aimag( (x) ), imsl_fc_convert( (x) ) ))

#define CTR(x,g)	(imsl_cf_convert( imsl_fc_convert( (x) )*(g), imsl_c_aimag( (x) )* (g) ))

#define ABSC(x)	(Mfloat)(fabs( imsl_fc_convert( (x) ) ) + fabs( imsl_c_aimag( (x) ) ))

#define CSINH(x)	(imsl_c_neg(TIMESI( imsl_c_sin( TIMESI( (x) ) ) )))



    imsl_e1psh ("C3IS ");



    if (accur < 0.0) {

	accur = imsl_amach (3);

	fpmin = imsl_amach (1) / accur;

	hugh = imsl_amach (2);

	fpmax = hugh * accur;

	acclog = log (accur);

	fphmin = sqrt (fpmin);

	fplmin = log (fpmin);

	pi1 = 3217.0e0 / 1024.0e0;

    }

    /* Check N */

    if (*nl <= 0) {

	imsl_e1sti (1, *nl);



	imsl_e1mes (5, 1, "The number of function values to be calculated, N = %(i1), must be at least 1.");

	/* Check XNU */

    }

    else if (*xnu <= -0.5e0) {

	imsl_e1str (1, *xnu);



	imsl_e1mes (5, 2, "The value of the input argument XNU = %(r1), must be greater than -1/2.");

	ztemp = imsl_cf_convert (hugh, 0.0e0);

	_l0 = 1;

	imsl_cset (nl, &ztemp, fi, &_l0);

    }

    if (imsl_n1rty (0) > 0)

	goto L_9000;



    mode = abs (*mode1);

    ifip = mod (mode, 2) == 1;

    nfp = 0;

    kase = 0;

    npq[0] = 0;

    npq[1] = 0;

    pi = pi1 + pi2;



    rk = 0;

    ek = 1.0e0;

    if (imsl_fc_convert (*zz) < 0.0e0) {

	ek = -1.0e0;

	rk = sign (1.0e0, imsl_c_aimag (*zz));

    }

    z = CTR (*zz, ek);

    x = TIMESI (z);

    scale = 0.0e0;

    if (*mode1 < 0)

	scale = imsl_fc_convert (z);

    absz = imsl_c_abs (z);



    if (absz > 1.0e0 / fphmin) {

	imsl_e1stc (1, z);

	imsl_e1str (1, 1.0e0 / fphmin);



	imsl_e1mes (5, 2, "The complex absolute value of the input argument, Z = %(c1), must not be greater than %(r1).");

	ztemp = imsl_cf_convert (hugh, 0.0e0);

	_l0 = 1;

	imsl_cset (nl, &ztemp, fi, &_l0);

	goto L_9000;

    }



    zi = imsl_c_div (imsl_cf_convert (1.0e0, 0.), z);

    ifail = -1;



    /*

     * XNU  is initial nu value, XMU  is final nu value, XN0  is nu value

     * near zero, and XLL  is much greater that XMU, such that XLL-XN0 = 2m

     * is even and .GT. 0. (The Y.. values are the X.. values + HALF, ie

     * Coulomb Ls)

     */

    xmu = *xnu + *nl - 1;

    ist = -imsl_i_max ((int) (*xnu + 0.4999), 0);

    if (absz < 1.0e-3 && (*xnu + ist) * log (absz) > 2.0)

	ist += 1;

    if (ist > 0)

	goto L_120;

    xn0 = *xnu + ist;

    yl0 = xn0 - 0.5e0;

    l0 = 1 + ist;



    /*

     * Choose KASE of solution required :

     */

    if (absz < 3.0 || (fabs (imsl_fc_convert (z)) < 2.0 && fabs (imsl_c_aimag (z)) <

	    20.0)) {

	kase = 3;

    }

    else if (absz > 1.5 * xmu && absz > 25.0) {

	kase = 1;

    }

    else {

	kase = 2;

    }



    /*

     * Normalisations for modified cylindrical Bessel functions

     */

    alpha = CTR (zi, 0.5e0);

    fil = imsl_c_sqrt (CTR (alpha, pi));

    if (imsl_fc_convert (fil) < 0.0e0)

	fil = imsl_c_neg (fil);



    /*

     * Start Calculation with K(xn0) (for xn0 near zero) :

     */

    if (kase <= 2) {



	_f0 = 1.0;

	_l0 = TRUE;

	pq1 = l_c4is (&x, &yl0, &_f0, &accur, &limit, &_l0,

	    &f1, &npq[0]);

	c1 = scale - imsl_fc_convert (z);

/*







	g = 2.0e0 * trunc (imsl_c_aimag (z) / (pi + pi));







*/

	g = 2.0e0 * (int) (imsl_c_aimag (z) / (pi + pi));



	c2 = (imsl_c_aimag (z) - g * pi1) - g * pi2;

	if (c1 < fplmin)

	    goto L_120;

	g2 = CTR (imsl_c_mul (fil, imsl_cf_convert (cos (c2), -sin (c2))), exp (c1));

	fk0 = imsl_c_mul (g2, f1);

	fkp0 = imsl_c_mul ((imsl_c_sub (TIMESI (pq1), alpha)), fk0);



	if (kase == 1) {

	    /*

	     * Find I(xn0) if OK to recur it upwards to XMU

	     */

	    _f0 = -1.0;

	    _l0 = FALSE;

	    pq2 = l_c4is (&x, &yl0, &_f0, &accur, &limit,

		&_l0, &f2, &npq[1]);

	    f2 = imsl_c_div (imsl_cf_convert (0.0e0, 2.0e0), (imsl_c_mul (f1, (imsl_c_sub (pq1,

				pq2)))));

	    sl = imsl_cf_convert (-2.0e0 * imsl_fc_convert (z), yl0 * pi - c2 - c2);

	    rat = imsl_cf_convert (0.0e0, 0.);

	    if (imsl_fc_convert (sl) > acclog - 2.0e0)

		rat = imsl_c_div (imsl_c_mul (imsl_c_exp (sl), f1), f2);

	    d0 = imsl_c_div (imsl_c_mul (alpha, f2), g2);

	    sl1 = imsl_c_sub (imsl_cf_convert (1.0e0, 0.), rat);

	    fil = imsl_c_mul (d0, sl1);

	    fipl = imsl_c_mul ((imsl_c_sub (TIMESI (imsl_c_sub (pq2, imsl_c_mul (rat, pq1))),

			imsl_c_mul (sl1, alpha))), d0);

	    goto L_80;

	}

    }

    /*

     * For KASE = 2 or 3, XMU too large to recur I up to it, so find XLL much

     * greater than XMU, and recur downwards -stable NB. in these KASES until

     * stmt 140, FIPL = I'/I

     */

L_10:

    i = imsl_i_vmax (3, *nl + 1, 3 + ist, (int) (absz));

    fil = imsl_cf_convert (0.0e0, 0.);

    fipl = imsl_cf_convert (1.0e0, 0.);

    sl = CTR (zi, 2.0e0);

    g2 = CTR (sl, *xnu + i - 1);

    g = 10.0 / accur;

    if (kase == 2)

	g = sqrt (g);

    for (nl1 = i; nl1 <= limit; nl1++) {

	f1 = imsl_c_sub (fil, imsl_c_mul (fipl, g2));

	if (ABSC (f1) > g)

	    goto L_30;

	fil = fipl;

	fipl = f1;

	g2 = imsl_c_add (g2, sl);

    }

    goto L_120;

L_30:

    nfp = nl1 - *nl;

    m = (nl1 - ist) / 2;

    if (nl1 - ist == m * 2)

	nl1 += 1;

    xll = *xnu + nl1 - 1;

    fil = imsl_cf_convert (1.0e0, 0.);



    fipl = imsl_c_add (CTR (zi, xll), CTR (z, 0.5e0 / (xll + 1.0e0)));



    sum[0] = imsl_cf_convert (1.0e0, 0.);

    sum[1] = imsl_cf_convert (1.0e0, 0.);

    c1 = 1.0e0;

    c2 = 1.0e0;

    yl0 = fpmax * imsl_f_min (1.0e0, pow (absz, xn0 + 1.0e0));



    /*

     * downward recurrence  of I and I'/I to nu = XN0

     */

    for (l = nl1 - 1; l >= l0; l--) {

	sl = CTR (zi, *xnu + l);

	sl1 = imsl_c_sub (sl, zi);

	if (ABSC (fil) > yl0) {

	    /* renormalise here & previously */

	    fil = imsl_c_mul (fil, imsl_cf_convert (fphmin, 0.));

	    if (kase == 3) {

		sum[0] = imsl_c_mul (sum[0], imsl_cf_convert (fphmin, 0.));

		sum[1] = imsl_c_mul (sum[1], imsl_cf_convert (fphmin, 0.));

	    }

	    for (i = imsl_i_max (l + 1, 1); i <= *nl; i++) {

		if (ABSC (fi[i - 1]) < fphmin)

		    fi[i - 1] = imsl_cf_convert (0.0e0, 0.);

		fi[i - 1] = CTR (fi[i - 1], fphmin);

	    }

	}

	d0 = imsl_c_add (sl, fipl);

	fil = imsl_c_mul (fil, d0);

	fipl = imsl_c_add (sl1, imsl_c_div (imsl_cf_convert (1.0e0, 0.), d0));

	if (l >= 1 && l <= *nl) {

	    fi[l - 1] = fil;

	    if (ifip)

		fip[l - 1] = fipl;

	}

	if (kase == 3 && mod (l - l0, 2) == 0) {

	    xm1 = xn0 + (m - 1);

	    x2m = xn0 + (m + m);

	    xx = x2m;

	    if (m > 1)

		xx = xm1 * x2m / (m * (x2m - 2.0e0));

	    c1 = -c1 / xx;

	    sum[0] = imsl_c_add (sum[0], CTR (fil, c1));

	    if (l >= l0 + 2 && mode <= 2) {

		c2 = c2 * (m - xn0) / (xx * (xn0 + xm1));

		sum[1] = imsl_c_add (sum[1], CTR (fil, c2));

	    }

	    m -= 1;

	}

    }

    if (imsl_fc_convert (fil) == 0.0e0)

	fil = imsl_cf_convert (accur, 0.0);



    if (kase == 2) {

	/*

	 * Normalise I(nu=XN0) by Wronskian with K(nu) by determining RAT =

	 * I(XN0) / FIL

	 */

	rat = imsl_c_div (zi, (imsl_c_mul ((imsl_c_sub (imsl_c_mul (fipl, fk0), fkp0)), fil)));



    }

    else {

	/*

	 * Normalise I and K by I sums & Von Neumann series: Begin by

	 * calculating G(xnu) = log(Gamma(1+xnu))/xnu  at XN0

	 */

	xm1 = -acclog * 0.4343;

	l = imsl_i_min (NGAM, (int) (1.1 + 1.3 * xm1) + 1);

	g = gam[l - 1];

	for (i = l - 1; i >= 1; i--) {

	    g = gam[i - 1] / (1.0e0 + xn0 * g);

	}



	/*

	 * Next calculate the first coefficient for the I sum :

	 */

#if 0

	g2 = imsl_c_add (imsl_cf_convert (g, 0.), imsl_c_log (imsl_c_add (zi, zi)));

#else

	{

	  Mf_complex logtmp;

	  Mf_complex addtmp;



	  addtmp = imsl_c_add (zi, zi);

	  logtmp.re = log(hypot(addtmp.re, addtmp.im));

	  if (addtmp.re) {

	    logtmp.im = atan2(addtmp.im, addtmp.re);

	  } else if (addtmp.im) {

	    logtmp.im = imsl_f_constant("pi", 0)/2.0;

	    if (addtmp.im < 0) logtmp.im *= -1;

	  } else {

	    logtmp.im = 0;

	  }

	  g2 = imsl_c_add (imsl_cf_convert (g, 0.), logtmp);

	}

#endif

	f2 = imsl_c_mul (g2, imsl_cf_convert (xn0, 0.));

	if (imsl_fc_convert (f2) < fplmin)

	    goto L_120;

	sl = imsl_c_exp(f2);

	rat = imsl_c_div (imsl_cf_convert (c1 * exp (-fabs (scale)), 0.), (imsl_c_mul (sl, sum[0])));



	if (mode <= 2) {

	    /* second coefficient for the K sum : */

	    sl1 = imsl_c_mul (imsl_c_mul (sl, sl), imsl_cf_convert ((xn0 + 2.0e0) / (1.0e0 -

			xn0), 0.));

	    /*

	     * first coefficient for K sum G1 = - pi*0.5e0 * (1/sin(pi*xn0) -

	     * 1/(pi*xn0))

	     */

	    xx = pi * xn0;

	    l = imsl_i_min (NCSC, (int) (2 + 0.45 * xm1) + 1);

	    g1 = csc[l - 1];

	    for (i = l - 1; i >= 1; i--) {

		g1 = csc[i - 1] / (1.0e0 + (xx * xx) * g1);

	    }

	    g1 *= -pi * 0.5e0 * xx;



	    d0 = g2;

	    if (xn0 != 0.0e0)

		d0 = imsl_c_add (imsl_cf_convert (g1, 0.), imsl_c_div (imsl_c_mul (sl, CSINH (f2)),

			imsl_cf_convert (xn0, 0.)));

	    fk0 = imsl_c_mul ((imsl_c_add (imsl_c_mul (d0, (imsl_c_mul (rat, fil))), imsl_c_mul (imsl_c_mul ((imsl_c_mul (rat,

					sl1)), sum[1]), imsl_cf_convert (1.0e0 / c2, 0.)))), imsl_cf_convert (exp (scale +

			fabs (scale)), 0.));

	    if (ABSC (fipl) < 30.0 || absz < 0.5e0) {

		fkp0 = imsl_c_sub (imsl_c_mul (fipl, fk0), imsl_c_div (zi, (imsl_c_mul (rat, fil))));

	    }

	    else {

		/*

		 * Calculate K' independently of Wronskian I'K - K'I=1/z, as

		 * I(xn0) is near a zero.   Use CF2 to give K'/K

		 */

		_f0 = xn0 - .5;

		_f1 = 1.0;

		_l0 = FALSE;

		pq1 = l_c4is (&x, &_f0, &_f1, 

		    &accur, &limit, &_l0, &f1, &npq[0]);

		fkp0 = imsl_c_mul ((imsl_c_sub (TIMESI (pq1), alpha)), fk0);

		kase = 4;

	    }

	}

    }

    yl0 = ABSC (rat);

    rat = imsl_c_div (rat, imsl_cf_convert (yl0, 0.));



    /*

     * Upward recurrence of K from FK0,FKP0 Upward recurrence of I from

     * FIL,FIPL if KASE=1 or   Renormalise FI,FIP at each nu if KASE=2,3

     */

L_80:

    ifail = *nl;

    g1 = xn0 * rk * pi;

    f2 = imsl_cf_convert (cos (g1), sin (g1));

    g2 = imsl_cf_convert (0.0e0, rk * pi);

    g = exp (imsl_f_max (-scale, fplmin));

    c2 = fpmin / g;

    c1 = ABSC (fil) * 0.2;

    j = 0;

    sl = imsl_c_mul (imsl_cf_convert (xn0 - 1.0e0, 0.), zi);

    for (l = l0; l <= *nl; l++) {

	sl1 = imsl_c_add (sl, zi);

	if (kase >= 2) {

	    if (l < 1)

		goto L_90;

	    fil = CTR (imsl_c_mul (rat, fi[l - 1]), yl0);

	    if (ifip)

		fipl = imsl_c_mul (fil, fip[l - 1]);

	    /*

	     * KASE=1 : I & I' at XN0 given, so recur upwards:

	     */

	}

	else if (l > l0) {

	    f1 = imsl_c_sub (fipl, imsl_c_mul (sl, fil));

	    fipl = imsl_c_sub (fil, imsl_c_mul (sl1, f1));

	    fil = f1;

	    if (ABSC (fil) < c1) {

		j += 1;

	    }

	    else {

		j = 0;

	    }

	    if (j > 3) {

		/*

		 * LOSING ACCURACY in KASE = 1 !!!!! Redo KASE=2 from large

		 * nu=XMU !

		 */

		l0 = l - 1;

		ist = l0 - 1;

		xn0 = *xnu + ist;

		kase = 2;

		goto L_10;

	    }

	}

	if (ABSC (fil) < fpmin)

	    goto L_120;

	if (l >= 1) {

	    fi[l - 1] = imsl_c_mul (fil, f2);

	    if (ifip)

		fip[l - 1] = CTR (imsl_c_mul (fipl, f2), ek);

	}

L_90:

	if (mode >= 3)

	    goto L_100;

	if (l > l0) {

	    f1 = imsl_c_sub (imsl_c_mul (sl, fk0), fkp0);

	    if (ABSC (f1) > fpmax * imsl_f_min (absz, 1.0e0))

		goto L_120;

	    fkp0 = imsl_c_neg (imsl_c_add (fk0, imsl_c_mul (sl1, f1)));

	    fk0 = f1;

	}

	if (l >= 1 && ek > 0.0e0) {

	    fk[l - 1] = fk0;

	    if (mode == 1)

		fkp[l - 1] = fkp0;

	}

	else if (l >= 1 && ek < 0.0e0) {

	    d0 = imsl_cf_convert (0.0e0, 0.);

	    if (g * ABSC (fk0) > c2)

		d0 = imsl_cf_convert (g * imsl_fc_convert (f2), -g * imsl_c_aimag (f2));

	    fk[l - 1] = imsl_c_sub (imsl_c_mul (d0, (imsl_c_mul (imsl_cf_convert (g, 0.), fk0))), imsl_c_mul (g2,

		    fil));

	    if (mode == 1)

		fkp[l - 1] = imsl_c_neg (imsl_c_sub (imsl_c_mul (d0, (imsl_c_mul (imsl_cf_convert (g, 0.),

				    fkp0))), imsl_c_mul (g2, fipl)));

	}

L_100:

	ifail = imsl_i_min (ifail, *nl - l);

	f2 = CTR (f2, ek);

	sl = sl1;

    }



L_120:

    if (ifail >= *nl)

	ifail = -1;

    /* Check ifail */

    if (ifail == -1) {



	imsl_e1mes (3, 1, "One of the continued fractions failed.");

    }

    else if (ifail > 0) {

	imsl_e1sti (1, *nl - ifail);



	imsl_e1mes (4, 2, "Only the first %(i1) entries in CBS are valid.");

    }



L_9000:

    imsl_e1pop ("C3IS ");

    return;

#undef	TRUE

#undef 	FALSE

#undef	CSINH

#undef	ABSC

#undef	CTR

#undef	TIMESI

}				/* end of function */

/* Structured by FOR_STRUCT, v0.2, on 06/06/91 at 18:28:44

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  C4IS/DC4IS (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    May 22, 1989



    Purpose:    Evaluate a continued fraction using the Steed version of

                Temme's algorithm.



    Usage:      C4IS(X, XL, OMEGA, EPS, LIMIT, NORM, F20, N)



    Arguments:

       X      - Complex point at which CF2 is to be evaluated.  (Input)

       XL     - Real parameter in CF2.  (Input)

       OMEGA  - Real parameter in CF2.  (Input)

       EPS    - Tolerance for error control.  (Input)

       LIMIT  - Maximum number of iterations for the continued fraction.

                (Input)

       NORM   - Logical for normalized result.  (Input)

       F20    - Complex hypergeometic function value.  (Output)

       N      - Number of iterations used.  (Output)

       C4IS   - Complex result.  (Output)



    Remark:

       C4IS evaluates



                                      (omega)        (omega)

                 CF2  = p + PM.q  =  H   (0.0,X)' / H   (0.0,X)

                                      XL             XL

       where PM = omega*i



       and evaluates H  (0.0,X) itself using the Steed version of Temme's

                      XL                     algorithm.



    Chapter:    SFUN/LIBRARY Bessel Functions



    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static Mf_complex l_c4is (Mf_complex *x, Mfloat *xl, Mfloat *omega,

		Mfloat *eps, Mint *limit, Mint *norm, Mf_complex *f20,

		Mint *n)

#else

static Mf_complex l_c4is (x, xl, omega, eps, limit, norm, f20, n)

    Mf_complex  *x;

    Mfloat      *xl, *omega, *eps;

    Mint        *limit;

    Mint  *norm;

    Mf_complex  *f20;

    Mint        *n;

#endif

{

    Mfloat       aa, asn, err;

    Mf_complex   bb, c4is_v, dd, dl, ds, em1, em2, en, pm, pq, qq, sn;





#define ABSC(w)	(float)(fabs( imsl_fc_convert( (w) ) ) + fabs( imsl_c_aimag( (w) ) ))



    pq = *x;

    aa = -*xl * (*xl + 1.0e0);

    pm = imsl_cf_convert (0.0e0, *omega);

    bb = imsl_c_mul (imsl_cf_convert (2.0e0, 0.), (imsl_c_add (pq, pm)));

    dd = imsl_c_div (imsl_cf_convert (1.0e0, 0.), bb);

    dl = imsl_c_mul (imsl_cf_convert (aa, 0.), dd);

    em1 = imsl_cf_convert (0.0e0, 0.);

    en = pm;

    *n = 1;

    qq = en;

    sn = imsl_c_add (imsl_cf_convert (1.0e0, 0.), imsl_c_mul (qq, dl));

    /* ASM = ABSC(SN) */

L_10:

    pq = imsl_c_add (pq, dl);

    *n += 1;

    em2 = em1;

    em1 = en;

    en = imsl_c_sub (imsl_c_mul (imsl_c_mul (imsl_cf_convert (0.0e0, -imsl_c_aimag (pm) / *n), bb), em1),

	imsl_c_mul (em2, imsl_cf_convert (aa / ((*n - 1) ** n), 0.)));

    aa += *n + *n - 2;

    bb = imsl_c_add (bb, (imsl_c_add (pm, pm)));

    dd = imsl_c_div (imsl_cf_convert (1.0e0, 0.), (imsl_c_add (imsl_c_mul (imsl_cf_convert (aa, 0.), dd), bb)));

    dl = imsl_c_mul (dl, (imsl_c_sub (imsl_c_mul (bb, dd), imsl_cf_convert (1.0e0, 0.))));

    err = ABSC (dl) / ABSC (pq);

    if (*norm) {

	qq = imsl_c_add (qq, en);

	ds = imsl_c_mul (qq, dl);

	sn = imsl_c_add (sn, ds);

	asn = ABSC (sn);

	err = imsl_f_max (err, ABSC (ds) / asn);

	/* ASM = MAX(ASM,ASN) */

    }

    if (err >= *eps && *n <= *limit)

	goto L_10;



    c4is_v = imsl_c_div (imsl_c_mul ((imsl_c_add (pq, dl)), pm), *x);

    if (!*norm)

	return (c4is_v);

    *f20 = imsl_c_div (imsl_cf_convert (1.0e0, 0.), sn);

    /* ERR = ACCUR* ASM / ASN */

    return (c4is_v);

#undef	ABSC

}				/* end of function */



