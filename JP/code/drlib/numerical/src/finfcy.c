#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn_cauchy (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Mfloat c, va_list argptr);
static void l_q2awc (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *c,
                Mfloat *errabs, Mfloat *errrel, Mfloat *result,
                Mfloat *errest, Mint *maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[]);
static void l_q3awc (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *c,
                Mfloat *epsabs, Mfloat *epsrel, Mint *limit,
                Mfloat *result, Mfloat *abserr, Mint *neval,
                Mint *ier, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint *last);
static void l_q4awc (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *c,
                Mfloat *result, Mfloat *abserr, Mint *krul,
                Mint *neval);
static Mfloat l_q5awc (Mfloat *x, Mfloat *c, Mfloat *p2, Mfloat *p3,
                Mfloat *p4, Mint *kp);
#else
static VA_LIST_HACK l_int_fcn_cauchy ();
static void l_q2awc ();
static void l_q3awc ();
static void l_q4awc ();
static Mfloat l_q5awc ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn_cauchy (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Mfloat c,...)
#else
Mfloat      imsl_f_int_fcn_cauchy (fcn, a, b, c, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    Mfloat      c;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, c);
    E1PSH ("imsl_f_int_fcn_cauchy", "imsl_d_int_fcn_cauchy");
    lv_value = F_ZERO;
    IMSL_CALL (l_int_fcn_cauchy (fcn, a, b, c, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn_cauchy", "imsl_d_int_fcn_cauchy");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn_cauchy (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Mfloat c, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn_cauchy (fcn, a, b, c, argptr)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    Mfloat      c;
    va_list     argptr;
#endif
{
    Mfloat      a_float;
    Mfloat      b_float;
    Mfloat      c_float;
    Mint        code;
    Mint        arg_number = 4;
    Mfloat      err_abs;
    Mfloat      err_rel;
    Mfloat     *err_est = NULL;
    Mfloat     *alist = NULL;
    Mfloat     *blist = NULL;
    Mfloat     *rlist = NULL;
    Mfloat     *elist = NULL;
    Mint        max_subinter = 500;
    Mint       *n_subinter = NULL;
    Mint       *n_evals = NULL;
    Mint       *iord = NULL;
    Mint        rule = 2;
    Mint        user_err_list = 0;
    Mint        user_err_order = 0;
    Mfloat      temp_err_est;
    Mint        temp_n_subinter;
    Mint        temp_n_evals;
    err_abs = sqrt( (double)imsl_amach (4));
    err_rel = sqrt( (double)imsl_amach (4));

    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RULE:
	    arg_number++;
	    rule = (Mint) va_arg (argptr, Mint);
	    break;
        case IMSL_ERR_ABS_ADR:
            arg_number++;
            err_abs = *(va_arg (argptr, Mfloat *));
            break;
        case IMSL_ERR_REL_ADR:
            arg_number++;
            err_rel = *(va_arg (argptr, Mfloat *));
            break;
	case IMSL_ERR_ABS:
	    arg_number++;
	    err_abs = (Mfloat) va_arg (argptr, double);
	    break;
	case IMSL_ERR_REL:
	    arg_number++;
	    err_rel = (Mfloat) va_arg (argptr, double);
	    break;
	case IMSL_ERR_EST:
	    arg_number++;
	    err_est = va_arg (argptr, Mfloat *);
	    break;
	case IMSL_MAX_SUBINTER:
	    arg_number++;
	    max_subinter = va_arg (argptr, Mint);
	    break;
	case IMSL_N_SUBINTER:
	    arg_number++;
	    n_subinter = va_arg (argptr, Mint *);
	    break;
	case IMSL_N_EVALS:
	    arg_number++;
	    n_evals = va_arg (argptr, Mint *);
	    break;
	case IMSL_ERR_LIST:
	    user_err_list = 1;
	    arg_number++;
	    elist = va_arg (argptr, Mfloat *);
	    break;
	case IMSL_ERR_ORDER:
	    user_err_order = 1;
	    arg_number++;
	    iord = va_arg (argptr, Mint *);
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
    /* CHECK MAXSUB */
    if (max_subinter < 1) {
	imsl_e1sti (1, max_subinter);
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_SUBINTER_SMALL);
    }
    if (imsl_n1rty (0))
	goto RETURN;
    if (*fcn == NULL) {
	imsl_e1stl (1, "fcn");
	/* (5, 1, "The required argument %(L1) is NULL."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
    }
    if (imsl_n1rty (0))
	goto RETURN;


    if (!user_err_list)
	elist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*elist));

    if (!user_err_order)
	iord = (Mint *) imsl_malloc (max_subinter * sizeof (*iord));

    alist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*alist));
    blist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*blist));
    rlist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*rlist));

    if (elist == NULL || iord == NULL || alist == NULL || blist == NULL ||
	rlist == NULL) {
	/* Not enough memory, with %(L1) = %(I1). */
	imsl_e1sti (1, max_subinter);
	imsl_e1stl (1, "max_subinter");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }

    if (err_est == NULL)
	err_est = &temp_err_est;
    if (n_subinter == NULL)
	n_subinter = &temp_n_subinter;
    if (n_evals == NULL)
	n_evals = &temp_n_evals;

    a_float = a;
    b_float = b;
    c_float = c;

    l_q2awc (fcn, &a_float, &b_float, &c_float, &err_abs, &err_rel,
	&lv_value, err_est, &max_subinter, n_evals, n_subinter,
	alist, blist, rlist, elist, iord);
FREE_SPACE:
    ;
    if (!user_err_list && elist != NULL)
	imsl_free (elist);
    if (!user_err_order && iord != NULL)
	imsl_free (iord);
    if (alist != NULL)
	imsl_free (alist);
    if (blist != NULL)
	imsl_free (blist);
    if (rlist != NULL)
	imsl_free (rlist);
RETURN:
    ;
    if (imsl_n1rty(0) > 3)
	lv_value = imsl_amach(6);
    return (argptr);
}



























/*Translated by FOR_C++, v0.1, on 08/16/90 at 17:23:54 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 17:23:52
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q2AWC/DQ2AWC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function F(X)/(X-C) in the Cauchy Principle
                Value sense.

    Usage:      CALL Q2AWC (F, A, B, C, ERRABS, ERRREL, RESULT, ERREST,
                            MAXSUB, NEVAL, NSUBIN, ALIST, BLIST, RLIST,
                            ELIST, IORD)

    Arguments:  (See QDAWC)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q2awc (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *c,
                Mfloat *errabs, Mfloat *errrel, Mfloat *result,
                Mfloat *errest, Mint *maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[])
#else
static void l_q2awc (f, a, b, c, errabs, errrel, result, errest,
                maxsub, neval, nsubin, alist, blist, rlist, elist, iord)
    Mfloat      (*f) (), *a, *b, *c, *errabs, *errrel, *result, *errest;
    Mint       *maxsub, *neval, *nsubin;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[];
#endif
{
    Mint        ier;


    imsl_e1psh ("l_q2awc");
    /* CHECK C */
    if (*c == *a) {
	imsl_e1str (1, *c);
	/*
	 * (5, 2, "Argument C and argument A both equal to %(r1).  They must
	 * be different.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_C_AND_A_DIFFERENT);
    }
    else if (*c == *b) {
	imsl_e1str (1, *c);
	/*
	 * (5, 3, "Argument C and argument B both equal to %(r1).  They must
	 * be different.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_C_AND_B_DIFFERENT);
    }
    /* CHECK ERRABS */
    if (*errabs < F_ZERO) {
	imsl_e1str (1, *errabs);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_ABS_SMALL);
    }
    /* CHECK ERRREL */
    if (*errrel < F_ZERO) {
	imsl_e1str (1, *errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_SMALL);
    }
    /*
     * CHECK ERRABS AND ERRREL
     */
    if (*errabs == F_ZERO && *errrel == F_ZERO) {
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_TOL_ZERO);
    }
    /* CHECK ERRREL .GE. 1 */
    if (*errrel >= F_ONE) {
	imsl_e1str (1, *errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_BIG);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;

    l_q3awc (f, a, b, c, errabs, errrel, maxsub, result, errest, neval,
	&ier, alist, blist, rlist, elist, iord, nsubin);

    if (ier == 1) {
	imsl_e1sti (1, *maxsub);
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_SUBINTERVALS);
    }
    else if (ier == 2) {
	imsl_e1str (1, *errabs);
	imsl_e1str (2, *errrel);
	imsl_ermes (IMSL_WARNING, IMSL_ROUNDOFF_CONTAMINATION);
    }
    else if (ier == 3) {
	imsl_e1str (1, alist[iord[0] - 1]);
	imsl_e1str (2, blist[iord[0] - 1]);
	imsl_ermes (IMSL_WARNING, IMSL_PRECISION_DEGRADATION);
    }
L_9000:
    ;
    imsl_e1pop ("l_q2awc");
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/16/90 at 17:24:56 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 17:24:54
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3AWC/DQ3AWC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q3AWC (F, A, B, C, EPSABS, EPSREL, LIMIT, RESULT,
                            ABSERR, NEVAL, IER, ALIST, BLIST, RLIST,
                            ELIST, IORD, LAST)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


  .......................................................................

   1.     Q3AWC
          COMPUTATION OF A CAUCHY PRINCIPAL VALUE
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A
             CAUCHY PRINCIPAL VALUE I = INTEGRAL OF F*W OVER (A,B)
             (W(X) = 1/(X-C), (C.NE.A, C.NE.B), HOPEFULLY SATISFYING
             FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

   3.     CALLING SEQUENCE
             CALL Q3AWC(F,A,B,C,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
                         NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)

          PARAMETERS
           ON ENTRY
              F      - REAL
                       FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                       FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                       DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

              A      - REAL
                       LOWER LIMIT OF INTEGRATION

              B      - REAL
                       UPPER LIMIT OF INTEGRATION

              C      - REAL
                       PARAMETER IN THE WEIGHT FUNCTION, C.NE.A, C.NE.B
                       IF C = A OR C = B, THE ROUTINE WILL END WITH
                       IER = 6.

              EPSABS - REAL
                       ABSOLUTE ACCURACY REQUESTED
              EPSREL - REAL
                       RELATIVE ACCURACY REQUESTED
                       IF  EPSABS.LT.0 AND EPSREL.LT.0,
                       THE ROUTINE WILL END WITH IER = 6.

              LIMIT  - INTEGER
                       GIVES AN UPPER BOUND ON THE NUMBER OF SUBINTERVALS
                       IN THE PARTITION OF (A,B), LIMIT.GE.1.

           ON RETURN
              RESULT - REAL
                       APPROXIMATION TO THE INTEGRAL

              ABSERR - REAL
                       ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                       WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

              NEVAL  - INTEGER
                       NUMBER OF INTEGRAND EVALUATIONS

              IER    - INTEGER
                       IER = 0 NORMAL AND RELIABLE TERMINATION OF THE
                               ROUTINE. IT IS ASSUMED THAT THE REQUESTED
                               ACCURACY HAS BEEN ACHIEVED.
                       IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE
                               THE ESTIMATES FOR INTEGRAL AND ERROR ARE
                               LESS RELIABLE. IT IS ASSUMED THAT THE
                               REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
                       IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
                               HAS BEEN ACHIEVED. ONE CAN ALLOW MORE SUB-
                               DIVISIONS BY INCREASING THE VALUE OF
                               LIMIT.
                               HOWEVER, IF THIS YIELDS NO IMPROVEMENT
                               IT IS ADVISED TO ANALYZE THE INTEGRAND,
                               IN ORDER TO DETERMINE THE INTEGRATION
                               DIFFICULTIES.  IF THE POSITION OF A
                               LOCAL DIFFICULTY CAN BE DETERMINED
                               (E.G. SINGULARITY, DISCONTINUITY WITHIN
                               THE INTERVAL) ONE WILL PROBABLY GAIN
                               FROM SPLITTING UP THE INTERVAL AT THIS
                               POINT AND CALLING APPROPRIATE INTEGRATORS
                               ON THE SUBRANGES.
                           = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS DETEC-
                               TED, WHICH PREVENTS THE REQUESTED
                               TOLERANCE FROM BEING ACHIEVED.
                           = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR
                               OCCURS AT SOME INTERIOR POINTS OF
                               THE INTEGRATION INTERVAL.
                           = 6 THE INPUT IS INVALID, BECAUSE
                               C = A OR C = B OR
                               EPSABS.LT.0 AND EPSREL.LT.0,
                               OR LIMIT.LT.1.
                               RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
                               IORD(1) AND LAST ARE SET TO ZERO.
                               ALIST(1) AND BLIST(1) ARE SET TO
                               A AND B RESPECTIVELY.

              ALIST   - REAL
                        VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                         LAST  ELEMENTS OF WHICH ARE THE LEFT END POINTS
                        OF THE SUBINTERVALS IN THE PARTITION OF THE
                        GIVEN INTEGRATION RANGE (A,B)

              BLIST   - REAL
                        VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                         LAST  ELEMENTS OF WHICH ARE THE RIGHT END POINTS
                        OF THE SUBINTERVALS IN THE PARTITION OF THE
                        GIVEN INTEGRATION RANGE (A,B)

              RLIST   - REAL
                        VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                         LAST  ELEMENTS OF WHICH ARE THE INTEGRAL
                        APPROXIMATIONS ON THE SUBINTERVALS

              ELIST   - REAL
                        VECTOR OF DIMENSION LIMIT, THE FIRST
                         LAST  ELEMENTS OF WHICH ARE THE MODULI OF THE
                        ABSOLUTE ERROR ESTIMATES ON THE SUBINTERVALS

              IORD    - INTEGER
                        VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST K
                        ELEMENTS OF WHICH ARE POINTERS TO THE ERROR
                        ESTIMATES OVER THE SUBINTERVALS, SO THAT
                        ELIST(IORD(1)), ..., ELIST(IORD(K)) WITH
                        K = LAST IF LAST.LE.(LIMIT/2+2), AND
                        K = LIMIT+1-LAST OTHERWISE, FORM A
                        DECREASING SEQUENCE

              LAST    - INTEGER
                        NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN
                        THE SUBDIVISION PROCESS

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q4AWC
                - Q10G
                - Q8AWO
                - Q7AWO
                - Q5AWC
                - Q4NG
                - F (USER-PROVIDED FUNCTION)
                - FORTRAN ABS, AMAX1, AMIN1

   ......................................................................




              LIST OF MAJOR VARIABLES
              -----------------------

             ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
                         CONSIDERED UP TO NOW
             BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
                         CONSIDERED UP TO NOW
             RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
                         (ALIST(I),BLIST(I))
             ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
             MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
                         ERROR ESTIMATE
             ERRMAX    - ELIST(MAXERR)
             AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
             ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
             ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
                         ABS(RESULT))
             *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
             *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
             LAST      - INDEX FOR SUBDIVISION


              MACHINE DEPENDENT CONSTANTS
              ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q3awc (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *c,
                Mfloat *epsabs, Mfloat *epsrel, Mint *limit,
                Mfloat *result, Mfloat *abserr, Mint *neval,
                Mint *ier, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint *last)
#else
static void l_q3awc (f, a, b, c, epsabs, epsrel, limit, result,
                abserr, neval, ier, alist, blist, rlist, elist, iord, last)
    Mfloat      (*f) (), *a, *b, *c, *epsabs, *epsrel;
    Mint       *limit;
    Mfloat     *result, *abserr;
    Mint       *neval, *ier;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[], *last;
#endif
{
    Mint        iroff1, iroff2, k, krule, maxerr, nev, nrmax;
    Mfloat      a1, a2, aa, area, area1, area12, area2, b1, b2, bb, epmach,
                errbnd, errmax, erro12, error1, error2, errsum, oflow,
                uflow;


    imsl_q4ng (&epmach, &uflow, &oflow);
    /* TEST ON VALIDITY OF PARAMETERS */
    *ier = 6;
    *neval = 0;
    *last = 0;
    alist[0] = *a;
    blist[0] = *b;
    rlist[0] = F_ZERO;
    elist[0] = F_ZERO;
    iord[0] = 0;
    *result = F_ZERO;
    *abserr = F_ZERO;
    if ((*c == *a || *c == *b) || (*epsabs < F_ZERO && *epsrel < F_ZERO))
	goto L_90;
    /* FIRST APPROXIMATION TO THE INTEGRAL */
    aa = *a;
    bb = *b;
    if (*a <= *b)
	goto L_10;
    aa = *b;
    bb = *a;
L_10:
    *ier = 0;
    krule = 1;
    l_q4awc (f, &aa, &bb, c, result, abserr, &krule, neval);
    *last = 1;
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 1;
    alist[0] = *a;
    blist[0] = *b;
    /* TEST ON ACCURACY */
    errbnd = imsl_f_max (*epsabs, *epsrel * fabs (*result));
    if (*limit == 1)
	*ier = 1;
    if (*abserr < imsl_f_min (1.0e-02 * fabs (*result), errbnd) || *ier ==
	1)
	goto L_80;
    /* INITIALIZATION */
    alist[0] = aa;
    blist[0] = bb;
    rlist[0] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    nrmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    /* MAIN DO-LOOP */
    for (*last = 2; *last <= *limit; (*last)++) {
	/*
	 * BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE.
	 */
	a1 = alist[maxerr - 1];
	b1 = F_HALF * (alist[maxerr - 1] + blist[maxerr - 1]);
	b2 = blist[maxerr - 1];
	if (*c <= b1 && *c > a1)
	    b1 = F_HALF * (*c + b2);
	if (*c > b1 && *c < b2)
	    b1 = F_HALF * (a1 + *c);
	a2 = b1;
	krule = 2;
	l_q4awc (f, &a1, &b1, c, &area1, &error1, &krule, &nev);
	*neval += nev;
	l_q4awc (f, &a2, &b2, c, &area2, &error2, &krule, &nev);
	*neval += nev;
	/*
	 * IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR AND TEST FOR
	 * ACCURACY.
	 */
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum += erro12 - errmax;
	area += area12 - rlist[maxerr - 1];
	if ((fabs (rlist[maxerr - 1] - area12) < 1.0e-05 * fabs (area12) &&
		erro12 >= 9.9e-01 * errmax) && krule == 0)
	    iroff1 += 1;
	if ((*last > 10 && erro12 > errmax) && krule == 0)
	    iroff2 += 1;
	rlist[maxerr - 1] = area1;
	rlist[*last - 1] = area2;
	errbnd = imsl_f_max (*epsabs, *epsrel * fabs (area));
	if (errsum <= errbnd)
	    goto L_20;
	/*
	 * TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
	 */
	if (iroff1 >= 6 && iroff2 > 20)
	    *ier = 2;
	/*
	 * SET ERROR FLAG IN THE CASE THAT NUMBER OF INTERVAL BISECTIONS
	 * EXCEEDS LIMIT.
	 */
	if (*last == *limit)
	    *ier = 1;
	/*
	 * SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR AT A POINT
	 * OF THE INTEGRATION RANGE.
	 */
	if (imsl_f_max (fabs (a1), fabs (b2)) <= (F_ONE + 1.0e03 * epmach) *
	    (fabs (a2) + 1.0e03 * uflow))
	    *ier = 3;
	/*
	 * APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
	 */
L_20:
	if (error2 > error1)
	    goto L_30;
	alist[*last - 1] = a2;
	blist[maxerr - 1] = b1;
	blist[*last - 1] = b2;
	elist[maxerr - 1] = error1;
	elist[*last - 1] = error2;
	goto L_40;
L_30:
	alist[maxerr - 1] = a2;
	alist[*last - 1] = a1;
	blist[*last - 1] = b1;
	rlist[maxerr - 1] = area2;
	rlist[*last - 1] = area1;
	elist[maxerr - 1] = error2;
	elist[*last - 1] = error1;
	/*
	 * CALL SUBROUTINE Q10G TO MAINTAIN THE DESCENDING ORDERING IN THE
	 * LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL WITH NRMAX-TH
	 * LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
	 */
L_40:
	imsl_q10g (limit, last, &maxerr, &errmax, elist, iord, &nrmax);

	/* JUMP OUT OF DO-LOOP */
	if (*ier != 0 || errsum <= errbnd)
	    goto L_60;
    }
    /* COMPUTE FINAL RESULT. */
L_60:
    *result = F_ZERO;
    for (k = 1; k <= *last; k++) {
	*result += rlist[k - 1];
    }
    *abserr = errsum;
L_80:
    if (aa == *b)
	*result = -*result;
L_90:
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/16/90 at 17:25:58 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 17:25:56
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q4AWC/DQ4AWC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q4AWC (F, A, B, C, RESULT, ABSERR, KRUL, NEVAL)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


  .......................................................................

   1.     Q4AWC
          INTEGRATION RULES FOR THE COMPUTATION OF CAUCHY
          PRINCIPAL VALUE INTEGRALS
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             TO COMPUTE I = INTEGRAL OF F*W OVER (A,B) WITH
             ERROR ESTIMATE, WHERE W(X) = 1/(X-C)

   3.     CALLING SEQUENCE
             CALL Q4AWC(F,A,B,C,RESULT,ABSERR,KRUL,NEVAL)

          PARAMETERS
             F      - REAL
                      FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                      F(X). THE ACTUAL NAME FOR F NEEDS TO BE DECLARED
                      E X T E R N A L  IN THE DRIVER PROGRAM.

             A      - REAL
                      LEFT END POINT OF THE INTEGRATION INTERVAL

             B      - REAL
                      RIGHT END POINT OF THE INTEGRATION INTERVAL, B.GT.A

             C      - REAL
                      PARAMETER IN THE WEIGHT FUNCTION

             RESULT - REAL
                      APPROXIMATION TO THE INTEGRAL
                      RESULT IS COMPUTED BY USING A GENERALIZED
                      CLENSHAW-CURTIS METHOD IF C LIES WITHIN TEN PERCENT
                      OF THE INTEGRATION INTERVAL.  IN THE OTHER CASE THE
                      15-POINT KRONROD RULE OBTAINED BY OPTIMAL ADDITION
                      OF ABSCISSAE TO THE 7-POINT GAUSS RULE, IS APPLIED.

             ABSERR - REAL
                      ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                      WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

             KRUL   - INTEGER
                      KEY WHICH IS DECREASED BY 1 IF THE 15-POINT
                      GAUSS-KRONROD SCHEME HAS BEEN USED

             NEVAL  - INTEGER
                      NUMBER OF INTEGRAND EVALUATIONS

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q7AWO
                - Q8AWO
                - F (USER-PROVIDED FUNCTION)
                - Q5AWC
                - FORTRAN ABS, ALOG, AMAX1, AMIN1

  .......................................................................




             THE VECTOR X CONTAINS THE VALUES COS(K*PI/24),
             K = 1, ..., 11, TO BE USED FOR THE CHEBYSHEV SERIES
             EXPANSION OF F

             LIST OF MAJOR VARIABLES
             ----------------------
             FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
                      COS(K*PI/24),  K = 0, ..., 24
             CHEB12 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
                      FOR THE FUNCTION F, OF DEGREE 12
             CHEB24 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
                      FOR THE FUNCTION F, OF DEGREE 24
             RES12  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
                      TO THE USE OF CHEB12
             RES24  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
                      TO THE USE OF CHEB24
             Q5AWC - EXTERNAL FUNCTION SUBPROGRAM DEFINING
                      THE WEIGHT FUNCTION
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             CENTR  - MID POINT OF THE INTERVAL


             CHECK THE POSITION OF C.

  -----------------------------------------------------------------------
 */
static Mfloat lv_x[] = {
    0.991444861373810411144557526929e0,
    0.965925826289068286749743199729e0,
    0.923879532511286756128183189397e0,
    0.866025403784438646763723170753e0,
    0.793353340291235164579776961501e0,
    0.707106781186547524400844362105e0,
    0.608761429008720639416097542898e0,
    0.5e0,
    0.38268343236508977172845998403e0,
    0.258819045102520762348898837624e0,
    0.130526192220051591548406227896e0
};
#ifdef ANSI
static void l_q4awc (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *c,
                Mfloat *result, Mfloat *abserr, Mint *krul,
                Mint *neval)
#else
static void l_q4awc (f, a, b, c, result, abserr, krul, neval)
    Mfloat      (*f) (), *a, *b, *c, *result, *abserr;
    Mint       *krul, *neval;
#endif
{
    Mint        i, isym, k, kp;
    Mfloat      ak22, amom0, amom1, amom2, cc, centr, cheb12[13], cheb24[25],
                fval[25], hlgth, p2, p3, p4, res12, res24, resabs, resasc,
                u;
    cc = (F_TWO ** c - *b - *a) / (*b - *a);
    if (fabs (cc) < 1.1e00)
	goto L_10;
    /*
     * APPLY THE 15-POINT GAUSS-KRONROD SCHEME.
     */
    *krul -= 1;
    imsl_q8awo (f, l_q5awc, c, &p2, &p3, &p4, &kp, a, b, result, abserr,
	&resabs, &resasc);
    *neval = 15;
    if (resasc == *abserr)
	*krul += 1;
    goto L_50;
    /*
     * USE THE GENERALIZED CLENSHAW-CURTIS METHOD.
     */
L_10:
    hlgth = F_HALF * (*b - *a);
    centr = F_HALF * (*b + *a);
    *neval = 25;
    imsl_e1usr ("ON");
    fval[0] = F_HALF * (*f) (hlgth + centr);
    fval[12] = (*f) (centr);
    fval[24] = F_HALF * (*f) (centr - hlgth);
    for (i = 2; i <= 12; i++) {
	u = hlgth * lv_x[i - 2];
	isym = 26 - i;
	fval[i - 1] = (*f) (u + centr);
	fval[isym - 1] = (*f) (centr - u);
    }
    imsl_e1usr ("OFF");
    /*
     * COMPUTE THE CHEBYSHEV SERIES EXPANSION.
     */
    imsl_q7awo (lv_x, fval, cheb12, cheb24);
    /*
     * THE MODIFIED CHEBYSHEV MOMENTS ARE COMPUTED BY FORWARD RECURSION,
     * USING AMOM0 AND AMOM1 AS STARTING VALUES.
     */
    amom0 = log (fabs ((F_ONE - cc) / (F_ONE + cc)));
    amom1 = F_TWO + cc * amom0;
    res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
    res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
    for (k = 3; k <= 13; k++) {
	amom2 = F_TWO * cc * amom1 - amom0;
	ak22 = (k - 2) * (k - 2);
	if ((k / 2) * 2 == k)
	    amom2 += -F_FOUR / (ak22 - F_ONE);
	res12 += cheb12[k - 1] * amom2;
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
    }
    for (k = 14; k <= 25; k++) {
	amom2 = F_TWO * cc * amom1 - amom0;
	ak22 = (k - 2) * (k - 2);
	if ((k / 2) * 2 == k)
	    amom2 += -F_FOUR / (ak22 - F_ONE);
	res24 += cheb24[k - 1] * amom2;
	amom0 = amom1;
	amom1 = amom2;
    }
    *result = res24;
    *abserr = fabs (res24 - res12);
L_50:
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/16/90 at 17:26:49 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 17:26:48
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q5AWC/DQ5AWC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      Q5AWC(X, C, P2, P3, P4, KP)

    Arguments:
       X      - Independent variable of integration.  (Input)
       C      - Variable of the weight function 1.0/(X-C).  (Input)
       P2     - Parameter that is passed through, with no function in
                this option.
       P3     - Parameter that is passed through, with no function in
                this option.
       P4     - Parameter that is passed through, with no function in
                this option.
       KP     - Parameter that is passed through, with no function in
                this option.
       Q5AWC  - Weight function.  (Output)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* P2, P3, P4, and KP are not used here, but leave the calling sequence 
   intact. */
#ifdef ANSI
static Mfloat l_q5awc (Mfloat *x, Mfloat *c, Mfloat *p2, Mfloat *p3,
                Mfloat *p4, Mint *kp)
#else
static Mfloat l_q5awc (x, c, p2, p3, p4, kp)
    Mfloat     *x, *c, *p2, *p3, *p4;
    Mint       *kp;
#endif
{
    Mfloat      q5awc_v;


    q5awc_v = F_ONE / (*x - *c);
    return (q5awc_v);
}				/* end of function */
