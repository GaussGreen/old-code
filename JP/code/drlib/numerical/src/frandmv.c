#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static void PROTO (l_c1trg, (Mint *n, Mfloat *a, Mint *lda));
static Mfloat   *lv_random;

static void PROTO (l_chfac, (Mint *n, Mfloat *a, Mint *lda, Mfloat *tol,
                    Mint *irank, Mfloat *r, Mint *ldr));
static VA_LIST_HACK PROTO(l_random_normal_multivariate,(Mint, Mint,
	Mfloat*, va_list));
static void PROTO(l_rnmvn, (Mint*, Mint*, Mfloat*, Mint*, Mfloat*,
        Mint*));
static void PROTO(l_rnnoa,(Mint*, Mfloat[]));


#ifdef ANSI
Mfloat *imsl_f_random_normal_multivariate(Mint n_vectors, Mint l_vectors,
	Mfloat *covariances, ...)
#else
Mfloat *imsl_f_random_normal_multivariate(n_vectors, l_vectors, covariances,
	va_alist)
    Mint        n_vectors;
    Mint        l_vectors;
    Mfloat	*covariances;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, covariances);

    E1PSH("imsl_f_random_normal_multivariate", "imsl_d_random_normal_multivariate");

    lv_random = NULL;
    IMSL_CALL(l_random_normal_multivariate(n_vectors, l_vectors, covariances,
		argptr));
    va_end(argptr);

    E1POP("imsl_f_random_normal_multivariate", "imsl_d_random_normal_multivariate");

    return(lv_random);
}

#ifdef ANSI
static VA_LIST_HACK l_random_normal_multivariate(Mint n_vectors, Mint l_vectors,
	Mfloat *covariances, va_list argptr)
#else
static VA_LIST_HACK l_random_normal_multivariate(n_vectors, l_vectors, covariances,
	argptr)
    Mint        n_vectors;
    Mint        l_vectors;
    Mfloat 	*covariances;
    va_list     argptr;
#endif
{
    Mint        code, return_user = 0, ner = 1, arg_number = 4;
    Mfloat	*rsig = NULL;
    Mfloat	tol = .00001;
    Mint	rank;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_random = va_arg(argptr, Mfloat*);
        return_user = 1;
    }
    else if (code != 0){
        imsl_e1sti (1, code);
        imsl_e1sti (2, arg_number);
        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
        goto RETURN;
    }

    imsl_c1iarg(n_vectors, "n_vectors", 1, 0, &ner);
    imsl_c1iarg(l_vectors, "l_vectors", 1, 0, &ner);
    if (imsl_n1rty(0) > 0) goto RETURN;

    if (!return_user) {
        lv_random = (Mfloat *) imsl_malloc (n_vectors*l_vectors*sizeof(Mfloat));
        if (!lv_random){
	    imsl_e1sti(1, n_vectors);
	    imsl_e1stl(1, "n_vectors");
	    imsl_e1sti(2, l_vectors);
	    imsl_e1stl(2, "l_vectors");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto RETURN;

        }
    }

    rsig = (Mfloat*) imsl_malloc (l_vectors*l_vectors*sizeof(Mfloat));
    if (!rsig) {
        imsl_e1sti(1, l_vectors);
	imsl_e1stl(1, "l_vectors");
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto RETURN;
    }
	
    l_chfac(&l_vectors, covariances, &l_vectors, &tol, &rank,
	rsig, &l_vectors);

    l_rnmvn (&n_vectors, &l_vectors, rsig, &l_vectors, lv_random, &n_vectors);

    imsl_f_m1ran (l_vectors, n_vectors, lv_random, lv_random);

RETURN:
    if (imsl_n1rty(0)>3 && !return_user) {
        imsl_free(lv_random);
        lv_random = NULL;
    }
    if (rsig != NULL) imsl_free (rsig);
    return (argptr);
}

/*Translated by FOR_C++, v0.1, on 01/27/92 at 13:46:13 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 01/27/92 at 13:46:11
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  RNMVN/DRNMVN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 30, 1991

    Purpose:    Generate pseudorandom numbers from a multivariate
                normal distribution.

    Usage:      CALL RNMVN (NR, K, RSIG, LDRSIG, R, LDR)

    Arguments:
       NR     - Number of random multivariate normal vectors to
                generate.  (Input)
       K      - Length of the multivariate normal vectors.  (Input)
       RSIG   - Upper triangular matrix, K by K, containing the
                Cholesky factor of the variance-covariance matrix.
                (Input)
                The variance-covariance matrix is equal to the product
                of the transpose of RSIG and RSIG.  RSIG can be
                obtained from the variance-covariance matrix using
                routine CHFAC.
       LDRSIG - Leading dimension of RSIG exactly as specified in the
                dimension statement in the calling program.  (Input)
       R      - NR by K matrix containing the random multivariate
                normal vectors in its rows.  (Output)
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       The routine RNSET  can be used to initialize the seed
       of the random number generator.  The routine RNOPT  can
       be used to select the form of the generator.

    Keywords:   Monte Carlo; Simulation; Bivariate normal distribution;
                Variance-covariance matrix; Cholesky factor; Random
                numbers

    GAMS:       L6b14

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rnmvn (Mint *nr, Mint *k, Mfloat *rsig, Mint *ldrsig,
	Mfloat *r, Mint *ldr)
#else
static void l_rnmvn (nr, k, rsig, ldrsig, r, ldr)
    Mint        *nr, *k;
    Mfloat      *rsig;
    Mint        *ldrsig;
    Mfloat      *r;
    Mint        *ldr;
#endif
{
#define RSIG(I_,J_)	(rsig+(I_)*(*ldrsig)+(J_))
#define R(I_,J_)	(r+(I_)*(*ldr)+(J_))
    Mint         _l0, i, j, ner;


    imsl_e1psh ("RNMVN ");
    ner = 1;

    imsl_c1dim (1, *nr, "NR", *ldr, "LDR", &ner);

    imsl_c1dim (1, *k, "K", *ldrsig, "LDRSIG", &ner);
    imsl_c1r (*k, rsig, *ldrsig, &ner);
    if (imsl_n1rty (0) > 0)
	goto L_9000;
    /* Generate NR*K normal random deviates. */
    if (*nr == *ldr) {
	_l0 = *nr ** k;
	l_rnnoa (&_l0, r);
    }
    else {
	for (j = 1; j <= *k; j++) {
	    l_rnnoa (nr, R (j - 1, 0));
	}
    }
    /* Form R*RSIG in place. */
    for (i = 1; i <= *nr; i++) {
	for (j = *k; j >= 1; j--) {
	    *R (j - 1, i - 1) = imsl_sdot (j, R (0, i - 1), *ldr, RSIG (j - 1, 0),
		1);
	}
    }
L_9000:
    imsl_e1pop ("RNMVN ");
    return;
}				/* end of function */
#undef RSIG
#undef R


/* Structured by FOR_STRUCT, v0.2, on 01/27/92 at 13:55:34
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  RNNOA/DRNNOA (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 30, 1991

    Purpose:    Generate pseudorandom numbers from a standard normal
                distribution using an acceptance/rejection method.

    Usage:      CALL RNNOA (NR, R)

    Arguments:
       NR     - Number of random numbers to generate.  (Input)
       R      - Vector of length NR containing the random standard
                normal deviates.  (Output)

    Remark:
       The routine RNSET  can be used to initialize the seed
       of the random number generator.  The routine RNOPT  can
       be used to select the form of the generator.

    Keywords:   Monte Carlo; Simulation; Gaussian distribution;
                Univariate continuous

    GAMS:       L6a14

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rnnoa (Mint *nr, Mfloat r[])
#else
static void l_rnnoa (nr, r)
    Mint        *nr;
    Mfloat       r[];
#endif
{
    Mint         _l0, _l1, i, ner;
    Mfloat       f, t, u, v, w, z;
    static Mfloat rsqr2p;
    static Mfloat b1 = 1.13113163544418e0;
    static Mfloat b2 = 0.63083480192196e0;
    static Mfloat b3 = 0.034240503750111e0;
    static Mfloat b4 = 0.180025191068563e0;
    static Mfloat b5 = 0.479727404222441e0;
    static Mfloat b6 = 1.10547366102207e0;
    static Mfloat b7 = 0.049264496373128e0;
    static Mfloat b8 = 0.59550713801594e0;
    static Mfloat b9 = 0.053377549506886e0;
    static Mfloat p1 = 0.884070402298758e0;
    static Mfloat p2 = 0.973310954173898e0;
    static Mfloat p3 = 0.986655477086949e0;
    static Mfloat p4 = 0.958720824790463e0;
    static Mfloat p5 = 0.755591531667601e0;
    static Mfloat p6 = 0.911312780288703e0;
    static Mfloat p7 = 0.87283497667179e0;
    static Mfloat p8 = 0.805577924423817e0;
    static Mfloat xi = 2.216035867166471e0;
    static Mfloat xi2 = 2.45540748228412e0;
    static Mint first = 1;

    if (first) {
	rsqr2p = 1.0e0 / (sqrt (4.0e0 * asin (1.0e0)));
	first = 0;
    }

    if (*nr <= 0) {
	imsl_e1psh ("RNNOA ");
	ner = 1;
	imsl_c1iarg (*nr, "NR", 1, 0, &ner);
	imsl_e1pop ("RNNOA ");
	goto L_9000;
    }
    else {
	for (i = 1; i <= *nr; i++) {
	    /* Step 1 */
	    imsl_f_random_uniform(1, IMSL_RETURN_USER, &u, 0);
	    if (u < p1) {
		imsl_f_random_uniform(1, IMSL_RETURN_USER, &v, 0);
		r[i - 1] = xi * (b1 * u + v - 1.0e0);
		goto L_60;
		/* Step 2 */
	    }
	    else if (u >= p2) {
		/* Step 3 */
	L_10:
		imsl_f_random_uniform(1, IMSL_RETURN_USER, &v, 0);
		imsl_f_random_uniform(1, IMSL_RETURN_USER, &w, 0);
		t = xi2 - log (w);
		if (v * v * t > xi2)
		    goto L_10;
		r[i - 1] = sqrt (2.0e0 * t);
		if (u >= p3)
		    r[i - 1] = -r[i - 1];
		goto L_60;
		/* Step 4 */
	    }
	    else if (u >= p4) {
		/* Step 5 */
	L_20:
		imsl_f_random_uniform(1, IMSL_RETURN_USER, &v, 0);
                imsl_f_random_uniform(1, IMSL_RETURN_USER, &w, 0);
		z = v - w;
		t = xi - b2 * imsl_f_min (v, w);
		if (imsl_f_max (v, w) <= p5)
		    goto L_50;
		f = rsqr2p * exp (-0.5e0 * t * t) - b4 * (xi - fabs (t));
		if (b3 * fabs (z) <= f)
		    goto L_50;
		goto L_20;
		/* Step 6 */
	    }
	    else if (u >= p6) {
		/* Step 7 */
	L_30:
		imsl_f_random_uniform(1, IMSL_RETURN_USER, &v, 0);
                imsl_f_random_uniform(1, IMSL_RETURN_USER, &w, 0);
		z = v - w;
		t = b5 + b6 * imsl_f_min (v, w);
		if (imsl_f_max (v, w) <= p7)
		    goto L_50;
		f = rsqr2p * exp (-0.5e0 * t * t) - b4 * (xi - fabs (t));
		if (b7 * fabs (z) <= f)
		    goto L_50;
		goto L_30;
		/* Step 8 */
	    }
	    else {
	L_40:
		imsl_f_random_uniform(1, IMSL_RETURN_USER, &v, 0);
                imsl_f_random_uniform(1, IMSL_RETURN_USER, &w, 0);
		z = v - w;
		t = b5 - b8 * imsl_f_min (v, w);
		if (t < 0.0e0)
		    goto L_40;
		if (imsl_f_max (v, w) <= p8)
		    goto L_50;
		f = rsqr2p * exp (-0.5e0 * t * t) - b4 * (xi - fabs (t));
		if (b9 * fabs (z) <= f)
		    goto L_50;
		goto L_40;
		/* Step 9 */
	    }
    L_50:
	    r[i - 1] = t;
	    if (z >= 0.0e0)
		r[i - 1] = -r[i - 1];
    L_60:
	    ;
	}
    }

L_9000:
    return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 08/28/90 at 16:28:49
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CHFAC/DCHFAC (Single/Double precision version)
 
    Computer:   FORC/SINGLE
 
    Revised:    April 9, 1985
 
    Purpose:    Compute an upper triangular factorization of a real
                symmetric nonnegative definite matrix.
 
    Usage:      CALL CHFAC (N, A, LDA, TOL, IRANK, R, LRD)
 
    Arguments:
       N      - Order of the matrix.  (Input)
       A      - N by N symmetric nonnegative definite matrix for which an
                upper triangular factorization is desired.  (Input)
                Only elements in the upper triangle of A are referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement in the calling program.  (Input)
       TOL    - Tolerance used in determining linear dependence.  (Input)
                For CHFAC, TOL = 100.*AMACH(4) is a common choice.   For
                DCHFAC, TOL = 100.D0*DMACH(4) is a common choice.  See
                documentation for IMSL routine AMACH/DMACH.
       IRANK  - Rank of A.  (Output)
                N - IRANK is the number of effective zero pivots.
       R      - N by N upper triangular matrix containing the R matrix
                from a Cholesky decomposition trans(R)*R of A.  (Output)
                The elements of the appropriate rows of R are set to 0.0
                if linear dependence of the columns of A is declared.
                (There are N-IRANK rows of R whose elements are set
                to 0.0.)  If A is not needed, then R and A can share
                the same storage locations.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement in the calling program.  (Input)
 
    Remarks:
    1. Informational error
       Type Code
         3   1  The input matrix is not nonnegative definite within the
                tolerance defined by TOL.
 
    2. Elements of row I of R are set to 0.0 if a linear dependence is
       declared.  Linear dependence is declared if
                      I-1
         ABS(A(I,I) - SUM R(J,I)**2) .LE. TOL*ABS(A(I,I))
                      J=1
 
    Keywords:   Cholesky decomposition; Square root method
 
    GAMS:       D9; D2b1b
 
    Chapter:    STAT/LIBRARY Mathematical Support
 
    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
 
    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.
 
  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_chfac (Mint *n, Mfloat *a, Mint *lda, Mfloat *tol, Mint *irank,
                Mfloat *r, Mint *ldr)
#else
static void l_chfac (n, a, lda, tol, irank, r, ldr)
    Mint       *n;
    Mfloat      a[];
    Mint       *lda;
    Mfloat     *tol;
    Mint       *irank;
    Mfloat      r[];
    Mint       *ldr;
#endif
{
    Mint        info, j, k, ner;
    Mfloat      s, t, x;
 
 
    imsl_e1psh ("l_chfac");
    /* CHECK FOR TERMINAL ERRORS */
    ner = 1;
 
    imsl_c1dim (1, *n, "n", *lda, "lda", &ner);
 
    imsl_c1dim (1, *n, "*n", *ldr, "ldr", &ner);
    if (*tol < 0.0 || *tol > 1.0) {
        imsl_e1str (1, *tol);
        imsl_ermes (IMSL_TERMINAL, IMSL_TOL_OUT_OF_RANGE);
    }
    if (imsl_n1rty (0) != 0)
        goto L_9000;
    info = 0;
    /* PRESERVE A COPY OF THE INPUT MATRIX */
    for (j = 1; j <= *n; j++) {
        scopy (j, &a[*lda * (j - 1)], 1, &r[*ldr * (j - 1)], 1);
    }
 
    *irank = 0;
    for (j = 1; j <= *n; j++) {
        s = 0.0e0;
        x = *tol * sqrt (fabs (r[j + *ldr * (j - 1) - 1]));
        for (k = 1; k <= (j - 1); k++) {
            t = r[k + *ldr * (j - 1) - 1] - imsl_sdot (k - 1, &r[*ldr * (k - 1)],
                1, &r[*ldr * (j - 1)], 1);
            if (r[k + *ldr * (k - 1) - 1] != 0.0) {
                t /= r[k + *ldr * (k - 1) - 1];
                r[k + *ldr * (j - 1) - 1] = t;
                s += t * t;
            }
            else {
                if (info == 0) {
                    if (fabs (t) > x * imsl_snrm2 (k - 1, &r[*ldr * (k - 1)],
                            1)) {
                        info = j;
                    }
                }
                r[k + *ldr * (j - 1) - 1] = 0.0;
            }
        }
        s = r[j + *ldr * (j - 1) - 1] - s;
        if (fabs (s) <= *tol * fabs (r[j + *ldr * (j - 1) - 1])) {
            s = 0.0e0;
        }
        else if (s < 0.0e0) {
            s = 0.0e0;
            if (info == 0)
                info = j;
        }
        else {
            *irank += 1;
        }
        r[j + *ldr * (j - 1) - 1] = sqrt (s);
    }
 
    if (info != 0) {
        imsl_e1sti (1, info);
        imsl_e1str (1, *tol);
        imsl_ermes (IMSL_WARNING, IMSL_NOT_NONNEG_DEFINITE);
    }
    /* FILL LOWER TRIANGLE OF R WITH ZEROS */
    l_c1trg (n, r, ldr);
    /* EXIT SECTION */
L_9000:
    imsl_e1pop ("l_chfac");
    return;
}                               /* end of function */
/*---------------------------------------------------------------------- */
 
/*  IMSL Name:  C1TRG/DC1TRG  (Single/Double precision version)
 
    Computer:   FORC/SINGLE
 
    Revised:    June 7, 1985
 
    Purpose:    Fill out the lower triangular portion of an upper
                triangular matrix with zeros.
 
    Usage:      C1TRG (N, A, LDA)
 
    Arguments:
       N      - Order of the matrix A.  (Input)
       A      - N x N upper triangular matrix to be filled out.
                (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
 
    GAMS:       D1b9;
 
    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Utilities
 
    Copyright:  1985 by IMSL, Inc.  All Rights Reserved
 
    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.
 
  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c1trg (Mint *n, Mfloat *a, Mint *lda)
#else
static void l_c1trg (n, a, lda)
    Mint       *n;
    Mfloat      a[];
    Mint       *lda;
#endif
{
    Mint        i;
 
 
    imsl_e1psh ("l_c1trg");
 
    if (*n <= 0) {
        imsl_e1sti (1, *n);
 
/*              imsl_ermes(5, 1, "N = %(i1).  The order of A, N, must be greater than 0.");
*/
        imsl_ermes (IMSL_TERMINAL, IMSL_NEED_A_AND_N_GT_ZERO);
    }
    if (*n > *lda) {
        imsl_e1sti (1, *n);
        imsl_e1sti (2, *lda);
        imsl_ermes (IMSL_TERMINAL, IMSL_NEED_N_LE_LDA);
    }
    if (imsl_n1rcd (0) != 0)
        goto L_9000;
    /* Fill lower triangle with zeros */
    for (i = 1; i <= (*n - 1); i++) {
        sset (*n - i, 0.0, &a[(i + 1) + *lda * (i - 1) - 1], 1);
    }
    /* Exit section */
L_9000:
    imsl_e1pop ("l_c1trg");
    return;
}                               /* end of function */
