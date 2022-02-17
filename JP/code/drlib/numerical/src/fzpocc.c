#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif



#ifdef ANSI

static VA_LIST_HACK l_zeros_poly(Mint ndeg, Mf_complex *coeff, va_list argptr);

static void l_zpocc(Mint *ndeg, Mf_complex coeff[], Mf_complex root[]);

static void l_z3occ(Mint *l1);

static void l_z4occ (Mint *l2, Mf_complex *z, Mint *conv);

static void l_z5occ (Mint *l3, Mf_complex *z, Mint *conv);

static void l_z6occ (Mint *bowl);

static void l_z7occ (Mint *bowl);

static void l_z8occ (Mint *nn, Mf_complex *s, Mf_complex p[],

                Mf_complex q[], Mf_complex *pv);

static Mfloat l_z9occ (Mint *nn, Mf_complex q[], Mfloat *rms, Mfloat *rmp,

                Mfloat *are, Mfloat *rmre);

static Mfloat l_z10cc (Mint *nn, Mf_complex p[], Mfloat q[]);

static Mfloat l_z11cc (Mint *nn, Mf_complex p[], Mfloat *repsr1,

                Mfloat *rinfp, Mfloat *repsp, Mfloat *radix);

static void l_z4lrc (Mf_complex *a, Mf_complex *b, Mf_complex *c,

                Mf_complex *zsm, Mf_complex *zlg);

static void l_sdsto (Mdouble acc[], Mfloat *s);

static void l_sdini (Mfloat *s, Mdouble acc[]);

static void l_sdadd (Mfloat *a, Mdouble acc[]);

static void l_czsto (Mdouble acc[], Mf_complex *z);

static void l_czini (Mf_complex *z, Mdouble acc[]);

static void l_czmul (Mf_complex *a, Mf_complex *b, Mdouble acc[]);

static void l_sdmul (Mfloat *a, Mfloat *b, Mdouble acc[]);

#else

static VA_LIST_HACK l_zeros_poly();

static void l_zpocc();

static void l_z3occ();

static void l_z4occ();

static void l_z5occ();

static void l_z6occ();

static void l_z7occ();

static void l_z8occ();

static Mfloat l_z9occ();

static Mfloat l_z10cc();

static Mfloat l_z11cc();

static void l_z4lrc();

static void l_sdsto();

static void l_sdini();

static void l_sdadd();

static void l_czsto();

static void l_czini();

static void l_czmul();

static void l_sdmul();

#endif



#define ADR(t,x)    ( t = x, &t )

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



static Mf_complex   *lv_roots;



#ifdef ANSI

Mf_complex *imsl_c_zeros_poly(Mint ndeg, Mf_complex *coeff, ...)

#else

Mf_complex *imsl_c_zeros_poly(ndeg, coeff, va_alist)

    Mint        ndeg;

    Mf_complex  *coeff;

    va_dcl

#endif

{

    va_list     argptr;



    VA_START(argptr, coeff);

    E1PSH("imsl_c_zeros_poly", "imsl_z_zeros_poly");

    lv_roots = NULL;

    IMSL_CALL(l_zeros_poly(ndeg, coeff, argptr));

    va_end(argptr);

    E1POP("imsl_c_zeros_poly", "imsl_z_zeros_poly");

    return lv_roots;

}





#ifdef ANSI

static VA_LIST_HACK l_zeros_poly(Mint ndeg, Mf_complex *coef, va_list argptr)

#else

static VA_LIST_HACK l_zeros_poly(ndeg, coef, argptr)

    Mint        ndeg;

    Mf_complex      *coef;

    va_list     argptr;

#endif

{

    Mint            code;

    Mint            arg_number  = 2;

    Mint            return_user = 0;

    Mint            use_companion = 0;

    Mf_complex     *companion = NULL;

    Mint            last_row_offset;

    Mint            i;



    code = 1;

    while (code > 0) {

        code = va_arg(argptr, Mint);

        arg_number++;

        switch (code) {

            case IMSL_RETURN_USER:

                lv_roots = va_arg(argptr, Mf_complex*);

                arg_number++;

                return_user = 1;

                break;

                 /*

                  * In order to use the  companion matrix method of finding

                  * the roots the follwing case was added.  This change was

                  * for IMSL/IDL.

                  */

            case IMSL_COMPANION_METHOD:

                use_companion = 1;

                arg_number++;

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

      companion = (Mf_complex*) calloc(ndeg*ndeg, sizeof(*companion));

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

        (companion+(ndeg+1)*i+1)->re = 1.0;

        (companion+(ndeg+1)*i+1)->im = 0.0;

      }



      /* fill last row with -coef[i]/coef[ndeg], i = 0 ,..., n-1 */



      last_row_offset = (ndeg-1) * ndeg;

      for (i=0; i<ndeg; i++)

        *(companion+last_row_offset+i) =

	  imsl_c_div(imsl_c_neg(coef[i]), coef[ndeg]);

      /* roots of p(z) = e-values of companion matrix */



      imsl_c_eig_gen(ndeg, companion, IMSL_RETURN_USER, lv_roots, 0);

      if (companion != NULL) imsl_free(companion);





    } else {



    /* CHECK FOR DEGREE > 50 OR < = ZERO.  */

    if (ndeg > 50 || ndeg <= 0) {

        imsl_e1sti(1, ndeg);

/*      (5, 1, "The degree of the polynomial must be less */

/*              than 50 and greater than zero.  NDEG is  */

/*              given as %(i1)."); */

        imsl_ermes(IMSL_TERMINAL, IMSL_POLYNOMIAL_DEGREE_50);

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



    l_zpocc(&ndeg, coef, lv_roots);

  }

RETURN:

    if (imsl_n1rty(0) > 3) {

        if (!return_user && lv_roots != NULL)   imsl_free(lv_roots);

        lv_roots = NULL;

    }

    return (argptr);

}





/*Translated by FOR_C++, v0.1, on 06/07/91 at 11:09:55 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 11:09:51

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  ZPOCC/DZPOCC (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    March 25, 1991



    Purpose:    Find the zeros of a polynomial with d_complex

                coefficients using the Jenkins-Traub three-stage

                algorithm.



    Usage:      CALL ZPOCC (NDEG, COEFF, ROOT)



    Arguments:

       NDEG   - Degree of the polynomial.  (Input)

       COEFF  - Complex vector of length NDEG+1 containing the

                coefficients of the polynomial in increasing order by

                degree.  (Input)

                The polynomial is COEFF(NDEG+1)*Z**NDEG+

                COEFF(NDEG)*Z**(NDEG-1)+...+COEFF(1).

       ROOT   - Complex vector of length NDEG containing the zeros of

                the polynomial.  (Output)



    Remark:

       Informational errors

       Type Code

         3   1  The first several coefficients of the polynomial are

                equal to zero.  Several of the last roots will be set

                to machine infinity to compensate for this problem.

         3   2  Fewer than NDEG zeros were found.  The ROOT vector will

                contain the value for machine infinity in the locations

                which do not contain zeros.



    Keywords:   Jenkins-Traub; Polynomials; Roots; Increasing order of

                modulus; Deflation; Nonlinear equations



    GAMS:       F1a1b



    Chapter:    MATH/LIBRARY Nonlinear Equations



    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

static struct t_z13cc {

    Mf_complex   p[50], h[50], qp[50], qh[50], s, t, pv;

    Mfloat       are, rmre, repsr1, rinfp;

    Mf_complex   sh[50];

    Mfloat       shr[50];

    Mint         nn;

}           z13cc;



#ifdef ANSI

static void l_zpocc(Mint *ndeg, Mf_complex coeff[], Mf_complex root[])

#else

static void l_zpocc (ndeg, coeff, root)

    Mint        *ndeg;

    Mf_complex   coeff[], root[];

#endif

{

    Mint 	conv;

    Mint         _l0, _l1, _l2, icnt1, icnt2, inx, k, l,

                n1, n2, newdeg;

    Mfloat       bnd, cosr, finity, radix, repsp, sinr, sqr2, xx, xxx, yy;

    Mf_complex   _cx0, z;





    imsl_e1psh ("ZPOCC ");

    /*

     * REVERSE THE ORDER OF THE COEFFICIENT VECTOR.

     */

    z13cc.nn = *ndeg + 1;

    k = mod (z13cc.nn, 2) + 1 + z13cc.nn / 2;

    n2 = *ndeg;

    imsl_cswap (ADR (_l0, z13cc.nn / 2), &coeff[0], ADR (_l1, 1), &coeff[k - 1],

	ADR (_l2, -1));

    /* CALL Z4LRC IF POLYNOMIAL IS QUADRATIC */

    if ((imsl_fc_convert (coeff[0]) != 0.0 || imsl_c_aimag (coeff[0]) != 0.0) &&

	*ndeg == 2) {

	l_z4lrc (&coeff[0], &coeff[1], &coeff[2], &root[0], &root[1]);

	goto L_110;

    }

    /*

     * CHECK FOR LEADING COEFFICIENTS EQUAL TO ZERO.

     */

    l = 0;

    z13cc.rinfp = imsl_amach (2);

    finity = imsl_amach (7);

    if (imsl_fc_convert (coeff[0]) == 0.0 && imsl_c_aimag (coeff[0]) == 0.0) {

L_20:

	l += 1;

	if (l <= *ndeg) {

	    root[n2 - 1] = imsl_cf_convert (finity, 0.0);

	    n2 -= 1;

	    if (imsl_fc_convert (coeff[l]) == 0.0 && imsl_c_aimag (coeff[l]) == 0.0)

		goto L_20;

	}



	if (l == 1) {



/*

	    imsl_e1mes (3, 1, "The leading coefficient of the polynomial is equal to zero.  The last root will be set to CMPLX(FINITY,0.0) where FINITY is the largest machine constant.");

*/

	    imsl_ermes(IMSL_WARNING, IMSL_ZERO_COEFF_1);

	}

	else {

	    imsl_e1sti (1, l);



/*

	    imsl_e1mes (3, 1, "The %(i1) leading coefficients of the polynomial are equal to zero.  The last %(i1) roots will be set to CMPLX(FINITY,0.0) where FINITY is the largest machine constant.");

*/

	    imsl_ermes(IMSL_WARNING, IMSL_ZERO_COEFF);

	}

    }

    if (l == *ndeg)

	goto L_110;

    newdeg = *ndeg - l;

    z13cc.nn = newdeg + 1;

    /* INITIALIZE REMAINING CONSTANTS */

    radix = imsl_imach (6);

    repsp = imsl_amach (1);

    z13cc.repsr1 = imsl_amach (4);

    sqr2 = sqrt (2.0e0);

    z13cc.are = z13cc.repsr1;

    z13cc.rmre = 2.0e0 * sqr2 * z13cc.repsr1;

    xx = sqrt (0.5e0);

    yy = -xx;

    cosr = -.06975647;

    sinr = .9975641;

    /*

     * REMOVE THE ZEROS AT THE ORIGIN IF ANY EXIST.

     */

L_30:

    if (imsl_fc_convert (coeff[z13cc.nn - 1]) != 0.0 || imsl_c_aimag (coeff[z13cc.nn - 1]) !=

	0.0)

	goto L_40;

    inx = newdeg - z13cc.nn + 2;

    root[inx - 1] = imsl_cf_convert (0.0, 0.0);

    z13cc.nn -= 1;

    if (z13cc.nn == 1)

	goto L_110;

    goto L_30;

    /* MAKE A COPY OF THE COEFFICIENTS */

L_40:

#ifdef COMPUTER_MIPNT

    imsl_ccopy (&z13cc.nn, &coeff[l], ADR (_l0, 1), &(z13cc.p)[0], ADR (_l1, 1));

#else

    imsl_ccopy (&z13cc.nn, &coeff[l], ADR (_l0, 1), z13cc.p, ADR (_l1, 1));

#endif

    /* SCALE THE POLYNOMIAL */

    bnd = l_z11cc (&z13cc.nn, z13cc.p, &z13cc.repsr1, &z13cc.rinfp,

	&repsp, &radix);

    if (bnd != 1.0e0)

	imsl_cscal (&z13cc.nn, ADR (_cx0, imsl_cf_convert (bnd, 0.0)), z13cc.p, ADR (_l0, 1));

    /* START THE ALGORITHM FOR ONE ZERO */

L_50:

    if (z13cc.nn > 2)

	goto L_60;

    /* CALCULATE THE FINAL ZERO AND RETURN */

    root[newdeg - 1] = imsl_c_neg (imsl_c_div (z13cc.p[1], z13cc.p[0]));

    goto L_100;

    /*

     * CALCULATE BND, A LOWER BOUND ON THE MODULUS OF THE ZEROS

     */

L_60:

#ifdef COMPUTER_MIPNT

    bnd = l_z10cc (&z13cc.nn, &(z13cc.p)[0], &(z13cc.shr)[0]);

#else

    bnd = l_z10cc (&z13cc.nn, z13cc.p, z13cc.shr);

#endif

    /*

     * OUTER LOOP TO CONTROL 2 MAJOR PASSES WITH DIFFERENT SEQUENCES OF

     * SHIFTS.

     */

    for (icnt1 = 1; icnt1 <= 2; icnt1++) {

	/* FIRST STAGE CALCULATION, NO SHIFT */

	l_z3occ (ADR (_l0, 5));

	/* INNER LOOP TO SELECT A SHIFT */

	for (icnt2 = 1; icnt2 <= 9; icnt2++) {

	    /*

	     * SHIFT IS CHOSEN WITH MODULUS BND AND AMPLITUDE ROTATED BY 94

	     * DEGREES FROM THE PREVIOUS SHIFT

	     */

	    xxx = cosr * xx - sinr * yy;

	    yy = sinr * xx + cosr * yy;

	    xx = xxx;

	    z13cc.s = imsl_c_mul (imsl_cf_convert (bnd, 0.), imsl_cf_convert (xx, yy));

	    /*

	     * SECOND STAGE CALCULATION, FIXED SHIFT.

	     */

	    l_z4occ (ADR (_l0, 10 * icnt2), &z, &conv);

	    if (!conv)

		goto L_70;

	    /*

	     * THE SECOND STAGE JUMPS DIRECTLY TO THE THIRD STAGE ITERATION.

	     * IF SUCCESSFUL THE ZERO IS STORED AND THE POLYNOMIAL DEFLATED.

	     */

	    inx = newdeg + 2 - z13cc.nn;

	    root[inx - 1] = z;

	    z13cc.nn -= 1;



#ifdef COMPUTER_MIPNT

	    imsl_ccopy (&z13cc.nn, &(z13cc.qp)[0], ADR (_l0, 1), &(z13cc.p)[0], ADR (_l1, 1));

#else

	    imsl_ccopy (&z13cc.nn, z13cc.qp, ADR (_l0, 1), z13cc.p, ADR (_l1, 1));

#endif



	    goto L_50;

    L_70:

	    ;

	    /*

	     * IF THE ITERATION IS UNSUCCESSFUL ANOTHER SHIFT IS CHOSEN.

	     */

	}

	/*

	 * IF 9 SHIFTS FAIL, THE OUTER LOOP IS REPEATED WITH ANOTHER SEQUENCE

	 * OF SHIFTS.

	 */

    }

    /*

     * THE ZEROFINDER HAS FAILED ON TWO MAJOR PASSES. RETURN EMPTY HANDED.

     */

    n1 = z13cc.nn - 1;

    imsl_e1sti (1, newdeg - n1);

    imsl_e1sti (2, n1);

/*

    imsl_e1mes (3, 2, "Only %(i1) zeros were found.  The last %(i2) locations of the ROOT vector will be set to CMPLX(FINITY,0.0), where FINITY is the largest machine constant.");

*/

    imsl_ermes(IMSL_WARNING, IMSL_FEWER_ZEROS_FOUND);

L_100:

    if (imsl_n1rcd (0) != 0) {

	/* SET UNFOUND ZEROS TO MACHINE INFINITY */

	n1 = z13cc.nn - 1;

	n2 = newdeg - z13cc.nn + 2;

	imsl_cset (&n1, ADR (_cx0, imsl_cf_convert (finity, 0.0)), &root[n2 - 1],

	    ADR (_l0, 1));

    }

    /* RESET COEFFICIENTS TO ORIGINAL ORDER */

L_110:

    z13cc.nn = *ndeg + 1;

    k = mod (z13cc.nn, 2) + 1 + z13cc.nn / 2;

    imsl_cswap (ADR (_l0, z13cc.nn / 2), &coeff[0], ADR (_l1, 1), &coeff[k - 1],

	ADR (_l2, -1));



    imsl_e1pop ("ZPOCC ");

    return;

}				/* end of function */







/*Translated by FOR_C++, v0.1, on 06/07/91 at 11:21:25 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 11:21:23

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z3OCC/DZ3OCC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Compute the derivative polynomial and the no-shift H

                polynomials.



    Usage:      CALL Z3OCC (L1)



    Arguments:

       L1     - The number of H no-shift polynomials.  (Input)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------



 */

#ifdef ANSI

static void l_z3occ(Mint *l1)

#else

static void l_z3occ (l1)

    Mint        *l1;

#endif

{

    Mint         _l0, _l1, _l2, i, j, jj, n;

    Mfloat       onedn, xni;



    n = z13cc.nn - 1;

    onedn = 1.0e0 / n;

    /*

     * COMPUTES THE DERIVATIVE POLYNOMIAL AS THE INITIAL H POLYNOMIAL AND

     * COMPUTES L1 NO-SHIFT H POLYNOMIALS.

     */

    for (i = 1; i <= n; i++) {

	xni = z13cc.nn - i;

	z13cc.h[i - 1] = imsl_c_mul (imsl_c_mul (imsl_cf_convert (xni, 0.), z13cc.p[i - 1]),

	    imsl_cf_convert (onedn, 0.));

    }



    for (jj = 1; jj <= *l1; jj++) {

	if (imsl_c_abs (z13cc.h[n - 1]) <= z13cc.repsr1 * 10.0 * imsl_c_abs (z13cc.p[n - 1])) {

	    /*

	     * IF THE CONSTANT TERM IS ESSENTIALLY ZERO, SHIFT H COEFFICIENTS

	     */

	    imsl_ccopy (ADR (_l0, n - 1), &z13cc.h[0], ADR (_l1, -1), &z13cc.h[1],

		ADR (_l2, -1));

	    z13cc.h[0] = imsl_cf_convert (0.0, 0.0);

	}

	else {

	    z13cc.t = imsl_c_neg (imsl_c_div (z13cc.p[z13cc.nn - 1], z13cc.h[n - 1]));



	    for (j = n; j >= 2; j--) {

		z13cc.h[j - 1] = imsl_c_add (imsl_c_mul (z13cc.t, z13cc.h[j - 2]),

		    z13cc.p[j - 1]);

	    }



	    z13cc.h[0] = z13cc.p[0];

	}

    }



    return;

}				/* end of function */











/*Translated by FOR_C++, v0.1, on 06/07/91 at 11:22:03 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 11:22:02

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z4OCC/DZ4OCC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPOCC to compute the fixed-shift H

                polynomials and test for convergence.



    Usage:      CALL Z4OCC (L2, Z, CONV)



    Arguments:

       L2     - Limit of fixed shift steps.  (Input)

       Z      - Approximate zero if CONV is true.  (Output)

       CONV   - Logical scalar indicating if convergence has occurred.

                   (Output)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */



#ifdef ANSI

static void l_z4occ (Mint *l2, Mf_complex *z, Mint *conv)

#else

static void l_z4occ (l2, z, conv)

    Mint        *l2;

    Mf_complex  *z;

    Mint  *conv;

#endif

{

    Mint   bowl, pasd, test;

    Mint         _l0, _l1, j, n;

    Mfloat       temp1, temp2;

    Mf_complex   ot, svs;





    n = z13cc.nn - 1;

    /*

     * COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TEST FOR CONVERGENCE.

     * INITIATES A VARIABLE-SHIFT ITERATION AND RETURN WITH THE APPROXIMATE

     * ZERO IF SUCCESSFUL. EVALUATE P AT S

     */



#ifdef COMPUTER_MIPNT

    l_z8occ (&z13cc.nn, &z13cc.s, &(z13cc.p)[0], &(z13cc.qp)[0], &z13cc.pv);

#else

    l_z8occ (&z13cc.nn, &z13cc.s, z13cc.p, z13cc.qp, &z13cc.pv);

#endif



    test = TRUE;

    pasd = FALSE;

    /* CALCULATE FIRST T = -P(S)/H(S) */

    l_z6occ (&bowl);

    /* MAIN LOOP FOR ONE SECOND STAGE STEP */

    for (j = 1; j <= *l2; j++) {

	ot = z13cc.t;

	/* COMPUTE NEXT H POLYNOMIAL AND NEW T */

	l_z7occ (&bowl);

	l_z6occ (&bowl);

	*z = imsl_c_add (z13cc.s, z13cc.t);

	/*

	 * TEST FOR CONVERGENCE UNLESS STAGE 3 HAS FAILED ONCE OR THIS IS THE

	 * LAST H POLYNOMIAL

	 */

	if ((!bowl && test) && j != *l2) {

	    temp1 = imsl_c_abs (imsl_c_sub (z13cc.t, ot));

	    temp2 = 0.5e0 * imsl_c_abs (*z);

	    if (temp1 >= temp2) {

		pasd = FALSE;

	    }

	    else if (!pasd) {

		pasd = TRUE;

		/*

		 * THE WEAK CONVERGENCE TEST HAS BEEN PASSED TWICE, START THE

		 * THIRD STAGE ITERATION, AFTER SAVING THE CURRENT H

		 * POLYNOMIAL AND SHIFT.

		 */

	    }

	    else {



#ifdef COMPUTER_MIPNT

		imsl_ccopy (&n, &(z13cc.h)[0], ADR (_l0, 1), &(z13cc.sh)[0], ADR (_l1, 1));

#else

		imsl_ccopy (&n, z13cc.h, ADR (_l0, 1), z13cc.sh, ADR (_l1, 1));

#endif



		svs = z13cc.s;

		l_z5occ (ADR (_l0, 10), z, conv);

		if (*conv)

		    goto L_9000;

		/*

		 * THE ITERATION FAILED TO CONVERGE. TURN OFF TESTING AND

		 * RESTORE H,S,PV AND T.

		 */

		test = FALSE;

		imsl_ccopy (&n, z13cc.sh, ADR (_l0, 1), z13cc.h, ADR (_l1, 1));

		z13cc.s = svs;

		l_z8occ (&z13cc.nn, &z13cc.s, z13cc.p, z13cc.qp, &z13cc.pv);

		l_z6occ (&bowl);

	    }

	}

    }

    /*

     * ATTEMPT AN ITERATION WITH FINAL H POLYNOMIAL FROM SECOND STAGE

     */

    l_z5occ (ADR (_l0, 10), z, conv);

L_9000:

    return;

}				/* end of function */











/*Translated by FOR_C++, v0.1, on 06/07/91 at 11:26:31 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 11:26:29

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z5OCC/DZ5OCC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPOCC to carry out the stage three

                iteration.



    Usage:      CALL Z5OCC (L3, Z, CONV)



    Arguments:

       L3     - Limit of steps in stage 3.  (Input)

       Z      - Contains the initial iterate on input and the final

                iterate if convergence is obtained on output.

                (Input/Output)

       CONV   - Logical scalar indicating if convergence has occurred.

                   (Output)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */



#ifdef ANSI

static void l_z5occ (Mint *l3, Mf_complex *z, Mint *conv)

#else

static void l_z5occ (l3, z, conv)

    Mint        *l3;

    Mf_complex  *z;

    Mint  *conv;

#endif

{

    Mint   b, bowl;

    Mint         i, j;

    Mfloat       omp, r1, relstp, rmp, rms, tp;





    *conv = FALSE;

    b = FALSE;

    z13cc.s = *z;

    /*

     * CARRY OUT THE THIRD STAGE ITERATION. MAIN LOOP FOR STAGE THREE

     */

    for (i = 1; i <= *l3; i++) {

	/*

	 * EVALUATE P AT S AND TEST FOR CONVERGENCE

	 */



#ifdef 	COMPUTER_MIPNT

	l_z8occ (&z13cc.nn, &z13cc.s, &(z13cc.p)[0], &(z13cc.qp)[0], &z13cc.pv);

#else

	l_z8occ (&z13cc.nn, &z13cc.s, z13cc.p, z13cc.qp, &z13cc.pv);

#endif

	rmp = imsl_c_abs (z13cc.pv);

	rms = imsl_c_abs (z13cc.s);



#ifdef COMPUTER_MIPNT

	if (rmp > 20.0e0 * l_z9occ (&z13cc.nn, &(z13cc.qp)[0], &rms, &rmp, &z13cc.are,

		&z13cc.rmre))

	    goto L_10;

#else

	if (rmp > 20.0e0 * l_z9occ (&z13cc.nn, z13cc.qp, &rms, &rmp, &z13cc.are,

		&z13cc.rmre))

	    goto L_10;

#endif

	/*

	 * POLYNOMIAL VALUE IS SMALLER IN VALUE THAN A BOUND ON THE ERROR IN

	 * EVALUATING P, TERMINATE THE ITERATION

	 */

	*conv = TRUE;

	*z = z13cc.s;

	return;

L_10:

	if (i == 1)

	    goto L_40;

	if ((b || rmp < omp) || relstp >= 0.05e0)

	    goto L_30;

	/*

	 * ITERATION HAS STALLED. PROBABLY A CLUSTER OF ZEROS. DO 5 FIXED

	 * SHIFT STEPS INTO THE CLUSTER TO FORCE ONE ZERO TO DOMINATE.

	 */

	tp = relstp;

	b = TRUE;

	if (relstp < z13cc.repsr1)

	    tp = z13cc.repsr1;

	/* 1       R1 = DSQRT(TP) */

	r1 = sqrt (tp);

	z13cc.s = imsl_c_mul (imsl_cf_convert (1.0 + r1, r1), z13cc.s);

	l_z8occ (&z13cc.nn, &z13cc.s, z13cc.p, z13cc.qp, &z13cc.pv);

	for (j = 1; j <= 5; j++) {

	    l_z6occ (&bowl);

	    l_z7occ (&bowl);

	}

	omp = z13cc.rinfp;

	goto L_50;

	/*

	 * EXIT IF POLYNOMIAL VALUE INCREASES SIGNIFICANTLY

	 */

L_30:

	if (0.1e0 * rmp > omp)

	    return;

L_40:

	omp = rmp;

	/* CALCULATE NEXT ITERATE */

L_50:

	l_z6occ (&bowl);

	l_z7occ (&bowl);

	l_z6occ (&bowl);

	if (bowl)

	    goto L_60;

	relstp = imsl_c_abs (z13cc.t) / imsl_c_abs (z13cc.s);

	z13cc.s = imsl_c_add (z13cc.s, z13cc.t);

L_60:

	;

    }

    return;

}				/* end of function */





















/*Translated by FOR_C++, v0.1, on 06/07/91 at 11:27:00 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 11:26:59

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z6OCC/DZ6OCC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPOCC to compute T = -P(S)/H(S).



    Usage:      CALL Z6OCC (BOWL)



    Arguments:

       BOWL   - Logical scalar set to TRUE if H(S) is zero.  (Output)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */



#ifdef ANSI

static void l_z6occ (Mint *bowl)

#else

static void l_z6occ (bowl)

    Mint	  *bowl;

#endif

{

    Mint         n;

    Mf_complex   th;





    n = z13cc.nn - 1;

    /* EVALUATE H(S) */



#ifdef COMPUTER_MIPNT

    l_z8occ (&n, &z13cc.s, &(z13cc.h)[0], &(z13cc.qh)[0], &th);

#else

    l_z8occ (&n, &z13cc.s, z13cc.h, z13cc.qh, &th);

#endif



    *bowl = imsl_c_abs (th) <= z13cc.are * 10.0 * imsl_c_abs (z13cc.h[n - 1]);



    if (*bowl) {

	z13cc.t = imsl_cf_convert (0.0, 0.0);

    }

    else {

	/* COMPUTE -P(S)/H(S) */

	z13cc.t = imsl_c_neg (imsl_c_div (z13cc.pv, th));

    }



    return;

}				/* end of function */









/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 11:27:28

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z7OCC/DZ7OCC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPOCC to calculate the next shifted H

                polynomial.



    Usage:      CALL Z7OCC (BOWL)



    Arguments:

       BOWL   - Logical scalar set to TRUE if H(S) is zero.  (Output)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */



#ifdef ANSI

static void l_z7occ (Mint *bowl)

#else

static void l_z7occ (bowl)

    Mint  *bowl;

#endif

{

    Mint         _l0, _l1, _l2, j, n;





    n = z13cc.nn - 1;

    if (!*bowl) {

	for (j = 2; j <= n; j++) {

	    z13cc.h[j - 1] = imsl_c_add (imsl_c_mul (z13cc.t, z13cc.qh[j - 2]),

		z13cc.qp[j - 1]);

	}

	z13cc.h[0] = z13cc.qp[0];

    }

    else {

	/* IF H(S) IS ZERO REPLACE H WITH QH */

	imsl_ccopy (ADR (_l0, n - 1), &z13cc.qh[0], ADR (_l1, 1), &z13cc.h[1],

	    ADR (_l2, 1));

	z13cc.h[0] = imsl_cf_convert (0.0, 0.0);

    }



    return;

}				/* end of function */















/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 11:27:54

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z8OCC/DZ8OCC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPOCC to evaluate a polynomial P at S

                by the Horner Recurrence.



    Usage:      CALL Z8OCC (NN, S, P, Q, PV)



    Arguments:

       NN     - Length of vectors P, and Q.  (Input)

       S      - Scalar containing the point the polynomial is to be

                evaluated at.  (Input)

       P      - Vector of length NN containing the coefficients of the

                polynomial.  (Input)

       Q      - Vector of length NN containing the partial sums.

                (Output)

       PV     - The computed polynomial evaluated at S.  (Output)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_z8occ (Mint *nn, Mf_complex *s, Mf_complex p[],

		Mf_complex q[], Mf_complex *pv)

#else

static void l_z8occ (nn, s, p, q, pv)

    Mint        *nn;

    Mf_complex  *s, p[], q[], *pv;

#endif

{

    Mint         i;





    q[0] = p[0];

    *pv = q[0];



    for (i = 2; i <= *nn; i++) {

	*pv = imsl_c_add (imsl_c_mul (*pv, *s), p[i - 1]);

	q[i - 1] = *pv;

    }



    return;

}				/* end of function */













/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 14:59:38

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z9OCC/DZ9OCC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPOCC to bound the error in evaluating

                the polynomial by the Horner Recurrence.



    Usage:      Z9OCC (NN, Q, RMS, RMP, ARE, RMRE)



    Arguments:

       NN     - Length of vector Q.  (Input)

       Q      - Vector containing the partial sums.  (Input)

       RMS    - Modulus of the point.  (Input)

       RMP    - Modulus of the polynomial value.  (Input)

       ARE    - Error bound on d_complex addition and multiplication.

                   (Input)

       RMRE   - Error bound on d_complex addition and multiplication.

                   (Input)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static Mfloat l_z9occ (Mint *nn, Mf_complex q[], Mfloat *rms, Mfloat *rmp,

		Mfloat *are, Mfloat *rmre)

#else

static Mfloat l_z9occ (nn, q, rms, rmp, are, rmre)

    Mint        *nn;

    Mf_complex   q[];

    Mfloat      *rms, *rmp, *are, *rmre;

#endif

{

    Mint         i;

    Mfloat       e, z9occ_v;





    e = imsl_c_abs (q[0]) ** rmre / (*are + *rmre);



    for (i = 1; i <= *nn; i++) {

	e = e ** rms + imsl_c_abs (q[i - 1]);

    }



    z9occ_v = e * (*are + *rmre) - *rmp ** rmre;



    return (z9occ_v);

}				/* end of function */







/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:03:19

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z10CC/DZ10CC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPOCC to compute a lower bound on

                the moduli of the zeros of a polynomial.



    Usage:      Z10CC (NN, P, Q)



    Arguments:

       NN     - Length of vectors PT and Q.  (Input)

       P      - Vector containing the coefficients.  (Input)

       Q      - Real work vector.  (Output)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static Mfloat l_z10cc (Mint *nn, Mf_complex p[], Mfloat q[])

#else

static Mfloat l_z10cc (nn, p, q)

    Mint        *nn;

    Mf_complex   p[];

    Mfloat       q[];

#endif

{

    Mint         i, n;

    Mfloat       df, dx, f, pt[50], x, xm, z10cc_v;





    for (i = 1; i <= *nn; i++) {

	pt[i - 1] = imsl_c_abs (p[i - 1]);

    }

    n = *nn - 1;

    pt[*nn - 1] = -pt[*nn - 1];

    /* COMPUTE UPPER ESTIMATE OF BOUND */

    x = exp ((log (-pt[*nn - 1]) - log (pt[0])) / n);

    if (pt[n - 1] == 0.0e0)

	goto L_20;

    /*

     * IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.

     */

    xm = -pt[*nn - 1] / pt[n - 1];

    if (xm < x)

	x = xm;

    /* CHOP THE INTERVAL (0,X) UNITL F.LE.0 */

L_20:

    xm = 0.1e0 * x;

    f = pt[0];

    for (i = 2; i <= *nn; i++) {

	f = f * xm + pt[i - 1];

    }

    if (f <= 0.0e0)

	goto L_40;

    x = xm;

    goto L_20;

L_40:

    dx = x;

    /*

     * DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES

     */

L_50:

    if (x == 0.0)

	goto L_80;

    if (fabs (dx / x) <= 0.005e0)

	goto L_80;

    q[0] = pt[0];

    for (i = 2; i <= *nn; i++) {

	q[i - 1] = q[i - 2] * x + pt[i - 1];

    }

    f = q[*nn - 1];

    df = q[0];

    for (i = 2; i <= n; i++) {

	df = df * x + q[i - 1];

    }

    dx = f / df;

    x -= dx;

    goto L_50;

L_80:

    z10cc_v = x;



    return (z10cc_v);

}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:03:45

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z11CC/DZ11CC  (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nuclues called by ZPOCC that returns a scale factor to

                multiply the coefficients of the polynomial.



    Usage:      Z11CC (NN, P, REPSR1, RINFP, REPSP, RADIX)



    Arguments:

       NN     - Length of vector P.  (Input)

       P      - Vector containing the coefficients of the polynomial.

                (Input)

       REPSR1 - Machine constant.  (Input)

       RINFP  - Machine constant.  (Input)

       REPSP  - Machine constant.  (Input)

       RADIX  - Machine constant.  (Input)



    Copyright:  1985 by IMSL, Inc. All Rights Reserved



    Warranty:   IMSL warrants only that IMSL testing has been applied to

                this code.  No other warranty, expressed or implied, is

                applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static Mfloat l_z11cc (Mint *nn, Mf_complex p[], Mfloat *repsr1,

		Mfloat *rinfp, Mfloat *repsp, Mfloat *radix)

#else

static Mfloat l_z11cc (nn, p, repsr1, rinfp, repsp, radix)

    Mint        *nn;

    Mf_complex   p[];

    Mfloat      *repsr1, *rinfp, *repsp, *radix;

#endif

{

    Mint         i, l;

    Mfloat       rhi, rlo, rmax, rmin, sc, x, z11cc_v;





#if defined(COMPUTER_ALFAC_IEEE) && defined(DOUBLE)

    rhi = 1.340781e+154;

#else

    rhi = sqrt (*rinfp);

#endif

    rlo = *repsp / *repsr1;

    rmax = 0.0e0;

    rmin = *rinfp;

    /*

     * FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.

     */

    for (i = 1; i <= *nn; i++) {

	x = imsl_c_abs (p[i - 1]);

	if (x > rmax)

	    rmax = x;

	if (x != 0.0e0 && x < rmin)

	    rmin = x;

    }

    /*

     * SCALE ONLY IF THERE ARE VERY LARGE OR VERY SMALL COMPONENTS

     */

    z11cc_v = 1.0e0;

    if (rmin < rlo || rmax > rhi) {

	x = rlo / rmin;



	if (x <= 1.0e0) {

	    sc = 1.0e0 / (sqrt (rmax) * sqrt (rmin));

	}

	else {

	    sc = x;

	    if (*rinfp / sc < rmax)

		sc = 1.0e0;

	}



	l = log (sc) / log (*radix) + 0.5e0;

	z11cc_v = imsl_fi_power (*radix, l);

    }



    return (z11cc_v);

}				/* end of function */







/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:08:01

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Z4LRC/DZ4LRC (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Nucleus called by ZPLRC to find the zeros of a quadratic

                with d_complex coefficients.



    Usage:      CALL Z4LRC (A, B, C, ZSM, ZLG)



    Arguments:

       A      - Complex coefficient of the quadritic Ax**2 + Bx + C.

                   (Input)

       B      - Complex coefficient of the quadritic Ax**2 + Bx + C.

                   (Input)

       C      - Complex coefficient of the quadritic Ax**2 + Bx + C.

                   (Input)

       ZSM    - The smallest root of the quadratic, in absolute value.

                   (Output)

       ZLG    - The largest root of the quadratic, in absolute value.

                   (Output)



    Remark:

       Informational errors

       Type Code

         3   1  Implies A and B equal zero.  On output ZLG will be set to

                CMPLX (FINITY,0.0) and ZSM will be set to -ZLG.  FINITY

                is the largest machine constant.

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

static void l_z4lrc (Mf_complex *a, Mf_complex *b, Mf_complex *c,

		Mf_complex *zsm, Mf_complex *zlg)

#else

static void l_z4lrc (a, b, c, zsm, zlg)

    Mf_complex  *a, *b, *c, *zsm, *zlg;

#endif

{

    Mint         is;

    Mfloat       _f0, imsl_ai, ar, imsl_bi, br, imsl_ci, cr, d, d1, finity,

                radix, rnlgrx, scale;

    Mdouble      acc[2], cacc[4];

    Mf_complex   _cx0, a0, b0, c0, czeros, zl, zs;





    imsl_e1psh ("Z4LRC ");

    /* INITIALIZE MACHINE CONSTANTS */

    radix = imsl_imach (6);

    rnlgrx = log (radix);

    finity = imsl_amach (7);

    /*

     * PUT THE COEFFICIENTS IN TEMPORARY STORAGE TO SAVE EXECUTION TIME.

     */

    a0 = *a;

    b0 = imsl_c_neg (*b);

    c0 = *c;

    ar = imsl_fc_convert (a0);

    imsl_ai = imsl_c_aimag (a0);

    br = imsl_fc_convert (b0);

    imsl_bi = imsl_c_aimag (b0);

    cr = imsl_fc_convert (c0);

    imsl_ci = imsl_c_aimag (c0);

    /*

     * CHECK FOR A AND/OR B EQUAL TO ZERO.

     */

    if (imsl_fc_convert (a0) == 0.0 && imsl_c_aimag (a0) == 0.0) {

	if (imsl_fc_convert (b0) == 0.0 && imsl_c_aimag (b0) == 0.0) {



                        imsl_e1sti(1, 2);

                        imsl_e1sti(2, 2);

                        imsl_ermes(IMSL_WARNING, IMSL_ZERO_COEFF);

/*

	    imsl_e1mes (3, 1, " The two leading coefficients are equal to zero.  On output ZLG will be set to CMPLX(FINITY,0.0) and ZSM will be set to -ZLG.  FINITY is the largest machine constant.");

*/

	    *zlg = imsl_cf_convert (finity, 0.0);

	    *zsm = imsl_c_neg (*zlg);

	}

	else {



                        imsl_e1sti(1, 1);

                        imsl_e1sti(2, 1);

                        imsl_ermes(IMSL_WARNING, IMSL_ZERO_COEFF);

/*

	    imsl_e1mes (3, 1, " The leading coefficient is equal to zero.  On output ZLG will be set to CMPLX(FINITY,0.0).  FINITY is the largest machine constant.");

*/

	    *zlg = imsl_cf_convert (finity, 0.0);

	    *zsm = imsl_c_div (c0, b0);

	}

    }

    if (imsl_n1rcd (0) != 0)

	goto L_9000;



    /* CHECK IF C IS EQUAL TO ZERO. */

    if (imsl_fc_convert (c0) != 0.0 || imsl_c_aimag (c0) != 0.0) {



	/*

	 * SCALING TO AVOID OVERFLOW OR UNDERFLOW. SCALE THE COEFFICIENTS SO

	 * THAT A*C IS APPROXIMATELY ONE. THE SCALE FACTOR CSQRT(A*C) FITS

	 * THIS REQUIREMENT BUT MAY CAUSE OVERFLOW OR UNDERFLOW IN THE

	 * SCALING PROCEDURE. LET AMAX1(ABS(AR),ABS(AI)) BE REPRESENTED BY

	 * RADIX**IA AND LET AMAX1(ABS(CR),ABS(CI) BE REPRESENTED BY

	 * RADIX**IC. THE SCALE FACTOR, SCALE, IS DEFINED BY THE FOLLOWING

	 * FORMULA:  SCALE=RADIX**IS, WHERE IS=ENTIER((IA+IC+1)/2) AND ENTIER

	 * IS THE MATHEMATICAL GREATEST INTEGER FUNCTION.

	 */

	is = (log (imsl_f_max (fabs (ar), fabs (imsl_ai))) + log (imsl_f_max (fabs (cr),

		    fabs (imsl_ci))) + rnlgrx) / (rnlgrx + rnlgrx);

	scale = imsl_fi_power (radix, is);

	/*

	 * IF THE SCALE FACTOR .LE. EPS*MAX(ABS(BR),ABS(BI)) DO NOT SCALE THE

	 * COEFFICIENTS.

	 */

	d1 = imsl_f_max (fabs (br), fabs (imsl_bi));

	l_sdini (&d1, acc);

	l_sdadd (&scale, acc);

	l_sdadd (ADR (_f0, -d1), acc);

	l_sdsto (acc, &d);



	if (d != 0.0) {

	    /*

	     * IF MAX(ABS(BR),ABS(BI)) .GE. DEPS*SCALE FACTOR THEN SCALE B0.

	     * OTHERWISE SET B0 = ZERO.

	     */

	    l_sdini (&d1, acc);

	    l_sdadd (&scale, acc);

	    l_sdadd (ADR (_f0, -scale), acc);

	    l_sdsto (acc, &d);



	    if (d == 0.0) {

		b0 = imsl_cf_convert (0.0, 0.0);

	    }

	    else {

		b0 = imsl_c_mul ((imsl_c_div (b0, imsl_cf_convert (scale, 0.))), imsl_cf_convert (0.5, 0.));

	    }



	    a0 = imsl_c_div (a0, imsl_cf_convert (scale, 0.));

	    c0 = imsl_c_div (c0, imsl_cf_convert (scale, 0.));

	    /* SOLVE A0*Z**2-2.0*B0*Z+C0=ZERO */

	    czeros = imsl_cf_convert (0.0, 0.0);

	    l_czini (&czeros, cacc);

	    l_czmul (&b0, &b0, cacc);

	    l_czmul (ADR (_cx0, imsl_c_neg (a0)), &c0, cacc);

	    l_czsto (cacc, &zs);

	    zs = imsl_c_sqrt (zs);

	    /*

	     * CHOOSE THE SIGN OF ZS SUCH THAT

	     * CABS(B)=AMAX1(CABS(B+ZS),CABS(B-ZS)).

	     */

	    if (imsl_c_abs (imsl_c_add (b0, zs)) < imsl_c_abs (imsl_c_sub (b0, zs)))

		zs = imsl_c_neg (zs);

	    b0 = imsl_c_add (b0, zs);

	}



	zs = imsl_c_div (c0, b0);

    }

    else {

	zs = imsl_cf_convert (0.0, 0.0);

    }



    zl = imsl_c_div (b0, a0);

    *zlg = zl;

    *zsm = zs;



L_9000:

    imsl_e1pop ("Z4LRC ");

    return;

}				/* end of function */







/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:15:25

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  SDSTO (Single precision version)



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

static void l_sdsto (Mdouble acc[], Mfloat *s)

#else

static void l_sdsto (acc, s)

    Mdouble      acc[];

    Mfloat      *s;

#endif

{





    *s = acc[0];



    return;

}				/* end of function */









/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:15:55

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  SDINI (Single precision version)



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

static void l_sdini (Mfloat *s, Mdouble acc[])

#else

static void l_sdini (s, acc)

    Mfloat      *s;

    Mdouble      acc[];

#endif

{





    acc[0] = (Mdouble)*s;

    acc[1] = 0.0e0;



    return;

}				/* end of function */











/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:16:20

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  SDADD (Single precision version)



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

static void l_sdadd (Mfloat *a, Mdouble acc[])

#else

static void l_sdadd (a, acc)

    Mfloat      *a;

    Mdouble      acc[];

#endif

{





    acc[0] += (Mdouble) *a;

    return;

}				/* end of function */







/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:16:44

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  CZSTO (Single precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Complex store.



    Usage:      CALL CZSTO (ACC, Z)



    Arguments:

       ACC    - Accumulator. (Input)

                ACC is a double precision vector of length 4. ACC is

                assumed to be the result computed by calling IMSL

                extended precision routines.

       Z      - Complex scalar. (Output)

                On output, Z contains a d_complex approximation

                to the value of the extended precision accumulator.



    GAMS:       A4a



    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_czsto (Mdouble acc[], Mf_complex *z)

#else

static void l_czsto (acc, z)

    Mdouble      acc[];

    Mf_complex  *z;

#endif

{





    *z = imsl_cf_convert ((Mfloat) acc[0], (Mfloat) acc[2]);



    return;

}				/* end of function */















/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:17:08

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  CZINI (Single precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Double precision initialization.



    Usage:      CALL CZINI (Z, ACC)



    Arguments:

       Z      - Complex scalar. (Input)

                On Input, Z contains a d_complex value that the

                accumulator is to be initialized to.

       ACC    - Accumulator. (Output)

                ACC is a double precision vector of length 4. ACC is

                initialized to Z.



    GAMS:       A4a



    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_czini (Mf_complex *z, Mdouble acc[])

#else

static void l_czini (z, acc)

    Mf_complex  *z;

    Mdouble      acc[];

#endif

{





    acc[0] = imsl_fc_convert (*z);

    acc[1] = 0.0e0;

    acc[2] = imsl_c_aimag (*z);

    acc[3] = 0.0e0;



    return;

}				/* end of function */











/* Structured by FOR_STRUCT, v0.2, on 06/07/91 at 15:17:32

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  CZMUL (Single precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Double precision multiply.



    Usage:      CALL CZMUL (A, B, ACC)



    Arguments:

       A      - Complex multiplier. (Input)

       B      - Complex multiplicand. (Input)

       ACC    - Accumulator. (Input/Output)

                ACC is a double precision vector of length 4.

                On output, ACC contains the sum of ACC and A*B.



    Remark:

       CZMUL adds the product A*B to the extended precision

       accumulator, ACC. The subroutine assumes that an double

       precision number is already in the accumulator.



    GAMS:       A4a



    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_czmul (Mf_complex *a, Mf_complex *b, Mdouble acc[])

#else

static void l_czmul (a, b, acc)

    Mf_complex  *a, *b;

    Mdouble      acc[];

#endif

{

    Mfloat       _f0, cy, cz, rw, rx;

    /*

     * SPLIT A AND B INTO REAL AND IMAGINARY COMPONENTS

     */

    rx = imsl_fc_convert (*a);

    rw = imsl_fc_convert (*b);

    cy = imsl_c_aimag (*a);

    cz = imsl_c_aimag (*b);

    /*

     * MULTIPLY CORRESPONDING COMPONENTS AND STORE IN THE ACCUMULATOR

     */

    l_sdmul (&rx, &rw, &acc[0]);

    l_sdmul (&cy, ADR (_f0, -cz), &acc[0]);

    l_sdmul (&rx, &cz, &acc[2]);

    l_sdmul (&rw, &cy, &acc[2]);



    return;

}				/* end of function */



#ifdef ANSI

static void l_sdmul (Mfloat *a, Mfloat *b, Mdouble acc[])

#else

static void l_sdmul(a, b, acc)

        Mfloat          *a, *b;

        Mdouble          acc[];

#endif

{



        acc[0] += ((Mdouble)*a) * ((Mdouble)*b);

        return;

}

