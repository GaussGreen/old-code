#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn_alg_log (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Imsl_quad weight, Mfloat alpha, Mfloat beta, va_list argptr);
static void l_q2aws (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mint *iweigh,
                Mfloat *alpha, Mfloat *beta, Mfloat *errabs,
                Mfloat *errrel, Mfloat *result, Mfloat *errest,
                Mint *maxsub, Mint *neval, Mint *nsubin,
                Mfloat alist[], Mfloat blist[], Mfloat rlist[],
                Mfloat elist[], Mint iord[]);
static void l_q3aws (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *alfa,
                Mfloat *beta, Mint *integr, Mfloat *epsabs,
                Mfloat *epsrel, Mint *limit, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mint *ier,
                Mfloat alist[], Mfloat blist[], Mfloat rlist[],
                Mfloat elist[], Mint iord[], Mint *last);
static void l_q4aws (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *bl,
                Mfloat *br, Mfloat *alfa, Mfloat *beta, Mfloat ri[],
                Mfloat rj[], Mfloat rg[], Mfloat rh[], Mfloat *result,
                Mfloat *abserr, Mfloat *resasc, Mint *integr,
                Mint *nev);
static void l_q5aws (Mfloat *alfa, Mfloat *beta, Mfloat ri[], Mfloat rj[],
                Mfloat rg[], Mfloat rh[], Mint *integr);
    static Mfloat l_q6aws (Mfloat *x, Mfloat *a, Mfloat *b, Mfloat *alfa,
                Mfloat *beta, Mint *integr);
#else
static VA_LIST_HACK l_int_fcn_alg_log ();
static void l_q2aws ();
static void l_q3aws ();
static void l_q4aws ();
static void l_q5aws ();
static Mfloat l_q6aws ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn_alg_log (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Imsl_quad weight, Mfloat alpha, Mfloat beta,...)
#else
Mfloat      imsl_f_int_fcn_alg_log (fcn, a, b, weight, alpha, beta, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    Imsl_quad   weight;
    Mfloat      alpha;
    Mfloat      beta;
va_dcl
#endif
{
    va_list argptr;
    VA_START (argptr, beta);
    E1PSH ("imsl_f_int_fcn_alg_log", "imsl_d_int_fcn_alg_log");
    lv_value = F_ZERO;
    IMSL_CALL (l_int_fcn_alg_log (fcn, a, b, weight, alpha, beta, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn_alg_log", "imsl_d_int_fcn_alg_log");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn_alg_log (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Imsl_quad weight, Mfloat alpha, Mfloat beta, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn_alg_log (fcn, a, b, weight, alpha, beta, argptr)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    Imsl_quad   weight;
    Mfloat      alpha;
    Mfloat      beta;
    va_list     argptr;
#endif
{
    Mint        l_weight=0;
    Mfloat      a_float;
    Mfloat      b_float;
    Mfloat      alpha_float;
    Mfloat      beta_float;
    Mint        code;
    Mint        arg_number = 6;
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
    Mint        user_err_list = 0;
    Mint        user_err_order = 0;
    Mint        i;
    Mfloat      temp_err_est;
    Mint        temp_n_subinter;
    Mint        temp_n_evals;


    err_abs = sqrt((double) imsl_amach (4));
    err_rel = sqrt((double) imsl_amach (4));


    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
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
	    err_abs = (Mfloat) va_arg (argptr, Mdouble);
	    break;
	case IMSL_ERR_REL:
	    arg_number++;
	    err_rel = (Mfloat) va_arg (argptr, Mdouble);
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
    alpha_float = alpha;
    beta_float = beta;

    if (weight == IMSL_ALG)
	l_weight = 1;
    if (weight == IMSL_ALG_LEFT_LOG)
	l_weight = 2;
    if (weight == IMSL_ALG_RIGHT_LOG)
	l_weight = 3;
    if (weight == IMSL_ALG_LOG)
	l_weight = 4;

    l_q2aws (fcn, &a_float, &b_float, &l_weight, &alpha_float, &beta_float,
	&err_abs, &err_rel, &lv_value, err_est, &max_subinter, n_evals,
	n_subinter, alist, blist, rlist, elist, iord);
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
    if (imsl_n1rty (0) > 3)
	lv_value = imsl_amach(6);
    return (argptr);
}





/*Translated by FOR_C++, v0.1, on 08/10/90 at 16:30:57 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 16:30:53
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q2AWS/DQ2AWS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function with algebraic-logarithmic
                singularities.

    Usage:      CALL Q2AWS (F, A, B, IWEIGH, ALPHA, BETA, ERRABS, ERRREL,
                            RESULT, ERREST, MAXSUB, NEVAL, NSUBIN, ALIST,
                            BLIST, RLIST, ELIST, IORD)

    Arguments:  (See QDAWS)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q2aws (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mint *iweigh,
                Mfloat *alpha, Mfloat *beta, Mfloat *errabs,
                Mfloat *errrel, Mfloat *result, Mfloat *errest,
                Mint *maxsub, Mint *neval, Mint *nsubin,
                Mfloat alist[], Mfloat blist[], Mfloat rlist[],
                Mfloat elist[], Mint iord[])
#else
static void l_q2aws (f, a, b, iweigh, alpha, beta, errabs, errrel,
                result, errest, maxsub, neval, nsubin, alist, blist, rlist, elist,
                iord)
    Mfloat      (*f) (), *a, *b;
    Mint       *iweigh;
    Mfloat     *alpha, *beta, *errabs, *errrel, *result, *errest;
    Mint       *maxsub, *neval, *nsubin;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[];
#endif
{
    Mint        ier;
    /* FIRST EXECUTABLE STATEMENT */
    imsl_e1psh ("Q2AWS  ");
    /* CHECK MAXSUB */
    if (*maxsub < 1) {
	imsl_e1sti (1, *maxsub);
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_SUBINTER_SMALL);
    }
    /* CHECK IWEIGH */
    if (*iweigh < 1 || *iweigh > 4) {
	imsl_e1sti (1, *iweigh);

	/*
	 * (5, 2, "The weight choice weight = %(i1).  It must be in the range
	 * 1 to 4.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_BAD_WEIGHT_CHOICE);
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
    /* CHECK ALPHA */
    if (*alpha <= -F_ONE) {
	imsl_e1str (1, *alpha);

	/*
	 * (5, 6, "The argument alpha = %(r1).  It must greater than -1.0.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_ALPHA_ARGUMENT);
    }
    /* CHECK BETA */
    if (*beta <= -F_ONE) {
	imsl_e1str (1, *beta);

	/*
	 * (5, 7, "The argument beta = %(r1).  It must greater than -1.0.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_BETA_ARGUMENT);
    }
    /* CHECK A AND B */
    if (*a >= *b) {
	imsl_e1str (1, *a);
	imsl_e1str (2, *b);

	/*
	 * (5, 8, "The limits of integration are a = %(r1) and b = %(r2).
	 * The value of the lower limit, a,  must be less than b.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_INTEGRATION_LIMITS);
    }
    /* CHECK ERRREL .GE. 1 */
    if (*errrel >= F_ONE) {
	imsl_e1str (1, *errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_BIG);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;

    l_q3aws (f, a, b, alpha, beta, iweigh, errabs, errrel, maxsub, result,
	errest, neval, &ier, alist, blist, rlist, elist, iord, nsubin);

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
    imsl_e1pop ("Q2AWS  ");
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 16:32:13 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 16:32:08
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3AWS/DQ3AWS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q3AWS (F, A, B, ALFA, BETA, INTEGR, EPSABS, EPSREL,
                            LIMIT, RESULT, ABSERR, NEVAL, IER, ALIST,
                            BLIST, RLIST, ELIST, IORD, LAST)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


   ......................................................................

   1.     Q3AWS
          INTEGRATION OF FUNCTIONS HAVING ALGEBRAICO-LOGARITHMIC
          END POINT SINGULARITIES
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
             DEFINITE INTEGRAL   I = INTEGRAL OF F*W OVER (A,B), (WHERE W
             SHOWS A SINGULAR BEHAVIOUR AT THE END POINTS, SEE PARAMETER
             INTEGR), HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

   3.     CALLING SEQUENCE
             CALL Q3AWS(F,A,B,ALFA,BETA,INTEGR,EPSABS,EPSREL,LIMIT,
                         RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,
                         IORD,LAST)

          PARAMETERS
           ON ENTRY
              F      - REAL
                       FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                       FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                       DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

              A      - REAL
                       LOWER LIMIT OF INTEGRATION

              B      - REAL
                       UPPER LIMIT OF INTEGRATION, B.GT.A
                       IF B.LE.A, THE ROUTINE WILL END WITH IER = 6.

              ALFA   - REAL
                       PARAMETER IN THE WEIGHT FUNCTION, ALFA.GT.(-1)
                       IF ALFA.LE.(-1), THE ROUTINE WILL END WITH
                       IER = 6.

              BETA   - REAL
                       PARAMETER IN THE WEIGHT FUNCTION, BETA.GT.(-1)
                       IF BETA.LE.(-1), THE ROUTINE WILL END WITH
                       IER = 6.

              INTEGR - INTEGER
                       INDICATES WHICH WEIGHT FUNCTION IS TO BE USED
                       = 1  (X-A)**ALFA*(B-X)**BETA
                       = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
                       = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
                       = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X)
                       IF INTEGR.LT.1 OR INTEGR.GT.4, THE ROUTINE
                       WILL END WITH IER = 6.

              EPSABS - REAL
                       ABSOLUTE ACCURACY REQUESTED
              EPSREL - REAL
                       RELATIVE ACCURACY REQUESTED
                       IF  EPSABS.LT.0 AND EPSREL.LT.0,
                       THE ROUTINE WILL END WITH IER = 6.

              LIMIT  - INTEGER
                       GIVES AN UPPER BOUND ON THE NUMBER OF SUBINTERVALS
                       IN THE PARTITION OF (A,B), LIMIT.GE.2
                       IF LIMIT.LT.2, THE ROUTINE WILL END WITH IER = 6.

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
                               THE ESTIMATES FOR THE INTEGRAL AND ERROR
                               ARE LESS RELIABLE. IT IS ASSUMED THAT THE
                               REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
                           = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
                               HAS BEEN ACHIEVED. ONE CAN ALLOW MORE
                               SUBDIVISIONS BY INCREASING THE
                               VALUE OF LIMIT.
                               HOWEVER, IF THIS YIELDS NO IMPROVEMENT
                               IT IS ADVISED TO ANALYZE THE INTEGRAND,
                               IN ORDER TO DETERMINE THE INTEGRATION
                               DIFFICULTIES WHICH PREVENT THE REQUESTED
                               TOLERANCE FROM BEING ACHIEVED.
                               IN CASE OF A JUMP DISCONTINUITY OR A LOCAL
                               SINGULARITY OF ALGEBRAICO-LOGARITHMIC TYPE
                               AT ONE OR MORE INTERIOR POINTS OF THE
                               INTEGRATION RANGE, ONE SHOULD PROCEED BY
                               SPLITTING UP THE INTERVAL AT THESE POINTS
                               AND CALLING THE INTEGRATOR ON THE
                               SUBRANGES.
                           = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS
                               DETECTED, WHICH PREVENTS THE REQUESTED
                               TOLERANCE FROM BEING ACHIEVED.
                           = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS
                               AT SOME POINTS OF THE INTEGRATION
                               INTERVAL.
                           = 6 THE INPUT IS INVALID, BECAUSE
                               B.LE.A OR ALFA.LE.(-1) OR BETA.LE.(-1), OR
                               INTEGR.LT.1 OR INTEGR.GT.4, OR
                               EPSABS.LT.0 AND EPSREL.LT.0,
                               OR LIMIT.LT.2.
                               RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
                               IORD(1) AND LAST ARE SET TO ZERO.
                               ALIST(1) AND BLIST(1) ARE SET TO
                               A AND B RESPECTIVELY.

              ALIST  - REAL
                       VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                        LAST  ELEMENTS OF WHICH ARE THE LEFT END POINTS
                       OF THE SUBINTERVALS IN THE PARTITION OF THE
                       GIVEN INTEGRATION RANGE (A,B)

              BLIST  - REAL
                       VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                        LAST  ELEMENTS OF WHICH ARE THE RIGHT END POINTS
                       OF THE SUBINTERVALS IN THE PARTITION OF THE
                       GIVEN INTEGRATION RANGE (A,B)

              RLIST  - REAL
                       VECTOR OF DIMENSION AT LEAST LIMIT,THE FIRST
                        LAST  ELEMENTS OF WHICH ARE THE INTEGRAL
                       APPROXIMATIONS ON THE SUBINTERVALS

              ELIST  - REAL
                       VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                        LAST  ELEMENTS OF WHICH ARE THE MODULI OF THE
                       ABSOLUTE ERROR ESTIMATES ON THE SUBINTERVALS

              IORD   - INTEGER
                       VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST K
                       OF WHICH ARE POINTERS TO THE ERROR ESTIMATES
                       OVER THE SUBINTERVALS, SO THAT
                       ELIST(IORD(1)), ..., ELIST(IORD(K)) WITH
                       K = LAST IF LAST.LE.(LIMIT/2+2), AND
                       K = LIMIT+1-LAST OTHERWISE FORM A
                       DECREASING SEQUENCE

              LAST   - INTEGER
                       NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN
                       THE SUBDIVISION PROCESS

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q5AWS
                - Q4AWS
                - Q10G
                - Q8AWO
                - Q7AWO
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
static void l_q3aws (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *alfa,
                Mfloat *beta, Mint *integr, Mfloat *epsabs,
                Mfloat *epsrel, Mint *limit, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mint *ier,
                Mfloat alist[], Mfloat blist[], Mfloat rlist[],
                Mfloat elist[], Mint iord[], Mint *last)
#else
static void l_q3aws (f, a, b, alfa, beta, integr, epsabs, epsrel,
                limit, result, abserr, neval, ier, alist, blist, rlist, elist,
                iord, last)
    Mfloat      (*f) (), *a, *b, *alfa, *beta;
    Mint       *integr;
    Mfloat     *epsabs, *epsrel;
    Mint       *limit;
    Mfloat     *result, *abserr;
    Mint       *neval, *ier;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[], *last;
#endif
{
    Mint        iroff1, iroff2, k, maxerr, nev, nrmax;
    Mfloat      a1, a2, area, area1, area12, area2, b1, b2, centre, epmach,
                errbnd, errmax, erro12, error1, error2, errsum, oflow,
                resas1, resas2, rg[25], rh[25], ri[25], rj[25], uflow;


    imsl_q4ng (&epmach, &uflow, &oflow);
    /* TEST ON VALIDITY OF PARAMETERS */
    *ier = 6;
    *neval = 0;
    *last = 0;
    rlist[0] = F_ZERO;
    elist[0] = F_ZERO;
    iord[0] = 0;
    *result = F_ZERO;
    *abserr = F_ZERO;
    if ((((((*b <= *a || (*epsabs < F_ZERO && *epsrel < F_ZERO)) ||
			*alfa <= (-F_ONE)) || *beta <= (-F_ONE)) || *integr < 1) ||
	    *integr > 4) || *limit < 2)
	goto L_100;
    *ier = 0;
    /*
     * COMPUTE THE MODIFIED CHEBYSHEV MOMENTS.
     */
    l_q5aws (alfa, beta, ri, rj, rg, rh, integr);
    /*
     * INTEGRATE OVER THE INTERVALS (A,(A+B)/2) AND ((A+B)/2,B).
     */
    centre = F_HALF * (*b + *a);
    l_q4aws (f, a, b, a, &centre, alfa, beta, ri, rj, rg, rh, &area1,
	&error1, &resas1, integr, &nev);
    *neval = nev;
    l_q4aws (f, a, b, &centre, b, alfa, beta, ri, rj, rg, rh, &area2,
	&error2, &resas2, integr, &nev);
    *last = 2;
    *neval += nev;
    *result = area1 + area2;
    *abserr = error1 + error2;
    /* TEST ON ACCURACY. */
    errbnd = imsl_f_max (*epsabs, *epsrel * fabs (*result));
    /* INITIALIZATION */
    if (error2 > error1)
	goto L_10;
    alist[0] = *a;
    alist[1] = centre;
    blist[0] = centre;
    blist[1] = *b;
    rlist[0] = area1;
    rlist[1] = area2;
    elist[0] = error1;
    elist[1] = error2;
    goto L_20;
L_10:
    alist[0] = centre;
    alist[1] = *a;
    blist[0] = *b;
    blist[1] = centre;
    rlist[0] = area2;
    rlist[1] = area1;
    elist[0] = error2;
    elist[1] = error1;
L_20:
    iord[0] = 1;
    iord[1] = 2;
    if (*limit == 2)
	*ier = 1;
    if (*abserr <= errbnd || *ier == 1)
	goto L_100;
    errmax = elist[0];
    maxerr = 1;
    nrmax = 1;
    area = *result;
    errsum = *abserr;
    iroff1 = 0;
    iroff2 = 0;
    /* MAIN DO-LOOP */
    for (*last = 3; *last <= *limit; (*last)++) {
	/*
	 * BISECT THE SUBINTERVAL WITH LARGEST ERROR ESTIMATE.
	 */
	a1 = alist[maxerr - 1];
	b1 = F_HALF * (alist[maxerr - 1] + blist[maxerr - 1]);
	a2 = b1;
	b2 = blist[maxerr - 1];

	l_q4aws (f, a, b, &a1, &b1, alfa, beta, ri, rj, rg, rh, &area1,
	    &error1, &resas1, integr, &nev);
	*neval += nev;
	l_q4aws (f, a, b, &a2, &b2, alfa, beta, ri, rj, rg, rh, &area2,
	    &error2, &resas2, integr, &nev);
	*neval += nev;
	/*
	 * IMPROVE PREVIOUS APPROXIMATIONS INTEGRAL AND ERROR AND TEST FOR
	 * ACCURACY.
	 */
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum += erro12 - errmax;
	area += area12 - rlist[maxerr - 1];
	if (*a == a1 || *b == b2)
	    goto L_30;
	if (resas1 == error1 || resas2 == error2)
	    goto L_30;

	/*
	 * TEST FOR ROUNDOFF ERROR.
	 */
	if (fabs (rlist[maxerr - 1] - area12) < 1.0e-05 * fabs (area12) &&
	    erro12 >= 9.9e-01 * errmax)
	    iroff1 += 1;
	if (*last > 10 && erro12 > errmax)
	    iroff2 += 1;
L_30:
	rlist[maxerr - 1] = area1;
	rlist[*last - 1] = area2;
	/* TEST ON ACCURACY. */
	errbnd = imsl_f_max (*epsabs, *epsrel * fabs (area));
	if (errsum <= errbnd)
	    goto L_40;
	/*
	 * SET ERROR FLAG IN THE CASE THAT THE NUMBER OF INTERVAL BISECTIONS
	 * EXCEEDS LIMIT.
	 */
	if (*last == *limit)
	    *ier = 1;
	/*
	 * SET ERROR FLAG IN THE CASE OF ROUNDOFF ERROR.
	 */
	if (iroff1 >= 6 || iroff2 >= 20)
	    *ier = 2;
	/*
	 * SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR AT INTERIOR
	 * POINTS OF INTEGRATION RANGE.
	 */
	if (imsl_f_max (fabs (a1), fabs (b2)) <= (F_ONE + 1.0e03 * epmach) *
	    (fabs (a2) + 1.0e03 * uflow))
	    *ier = 3;
	/*
	 * APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
	 */
L_40:
	if (error2 > error1)
	    goto L_50;
	alist[*last - 1] = a2;
	blist[maxerr - 1] = b1;
	blist[*last - 1] = b2;
	elist[maxerr - 1] = error1;
	elist[*last - 1] = error2;
	goto L_60;
L_50:
	alist[maxerr - 1] = a2;
	alist[*last - 1] = a1;
	blist[*last - 1] = b1;
	rlist[maxerr - 1] = area2;
	rlist[*last - 1] = area1;
	elist[maxerr - 1] = error2;
	elist[*last - 1] = error1;
	/*
	 * CALL SUBROUTINE Q10G TO MAINTAIN THE DESCENDING ORDERING IN THE
	 * LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL WITH LARGEST
	 * ERROR ESTIMATE (TO BE BISECTED NEXT).
	 */
L_60:
	imsl_q10g (limit, last, &maxerr, &errmax, elist, iord, &nrmax);

	/* JUMP OUT OF DO-LOOP */
	if (*ier != 0 || errsum <= errbnd)
	    goto L_80;
    }
    /* COMPUTE FINAL RESULT. */
L_80:
    *result = F_ZERO;
    for (k = 1; k <= *last; k++) {
	*result += rlist[k - 1];
    }
    *abserr = errsum;
L_100:
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 16:35:20 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 16:35:16
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q4AWS/DQ4AWS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q4AWS (F, A, B, BL, BR, ALFA, BETA, RI, RJ, RG, RH,
                            RESULT, ABSERR, RESASC, INTEGR, NEV)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


  .......................................................................

   1.     Q4AWS
          INTEGRATION RULES FOR INTEGRANDS HAVING ALGEBRAICO-LOGARITHMIC
          END POINT SINGULARITIES
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             TO COMPUTE I = INTEGRAL OF F*W OVER (BL,BR), WITH ERROR
             ESTIMATE, WHERE THE WEIGHT FUNCTION W HAS A SINGULAR
             BEHAVIOUR OF ALGEBRAICO-LOGARITHMIC TYPE AT THE POINTS
             A AND/OR B. (BL,BR) IS A PART OF (A,B).

   3.     CALLING SEQUENCE
             CALL Q4AWS(F,A,B,BL,BR,ALFA,BETA,RI,RJ,RG,RH,RESULT,
                         ABSERR,RESASC,INTEGR,NEV)

          PARAMETERS
             F      - REAL
                      FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                      F(X). THE ACTUAL NAME FOR F NEEDS TO BE DECLARED
                      E X T E R N A L  IN THE DRIVER PROGRAM.

             A      - REAL
                      LEFT END POINT OF THE ORIGINAL INTERVAL

             B      - REAL
                      RIGHT END POINT OF THE ORIGINAL INTERVAL, B.GT.A.

             BL     - REAL
                      LOWER LIMIT OF INTEGRATION, BL.GE.A

             BR     - REAL
                      UPPER LIMIT OF INTEGRATION, BR.LE.B

             ALFA   - REAL
                      PARAMETER IN THE WEIGHT FUNCTION

             BETA   - REAL
                      PARAMETER IN THE WEIGHT FUNCTION

             RI,RJ,RG,RH - REAL
                      MODIFIED CHEBYSHEV MOMENTS FOR THE APPLICATION
                      OF THE GENERALIZED CLENSHAW-CURTIS METHOD
                      (COMPUTED IN SUBROUTINE Q6AWS)

             RESULT - REAL
                      APPROXIMATION TO THE INTEGRAL
                      RESULT IS COMPUTED BY USING A GENERALIZED
                      CLENSHAW-CURTIS METHOD IF B1 = A OR BR = B.
                      IN ALL OTHER CASES THE 15-POINT KRONROD
                      RULE IS APPLIED, OBTAINED BY OPTIMAL ADDITION OF
                      ABSCISSAE TO THE 7-POINT GAUSS RULE.

             ABSERR - REAL
                      ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                      WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

             RESASC - REAL
                      APPROXIMATION TO THE INTEGRAL OF ABS(F*W-I/(B-A))

             INTEGR - INTEGER
                      WHICH DETERMINES THE WEIGHT FUNCTION
                      = 1   W(X) = (X-A)**ALFA*(B-X)**BETA
                      = 2   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
                      = 3   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
                      = 4   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*
                                   LOG(B-X)

             NEV    - INTEGER
                      NUMBER OF INTEGRAND EVALUATIONS

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q7AWO
                - Q8AWO
                - F (USER-PROVIDED FUNCTION)
                - Q6AWS
                - FORTRAN ABS, ALOG, AMAX1

  .......................................................................




             THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
             K = 1, ..., 11, TO BE USED FOR THE COMPUTATION OF THE
             CHEBYSHEV SERIES EXPANSION OF F.


             LIST OF MAJOR VARIABLES
             -----------------------

             FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
                      (BR-BL)*0.5*COS(K*PI/24)+(BR+BL)*0.5
                      K = 0, ..., 24
             CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
                      OF DEGREE 12, FOR THE FUNCTION F, IN THE
                      INTERVAL (BL,BR)
             CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
                      OF DEGREE 24, FOR THE FUNCTION F, IN THE
                      INTERVAL (BL,BR)
             RES12  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB12
             RES24  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB24
             Q6AWS - EXTERNAL FUNCTION SUBPROGRAM DEFINING
                      THE FOUR POSSIBLE WEIGHT FUNCTIONS
             HLGTH  - HALF-LENGTH OF THE INTERVAL (BL,BR)
             CENTR  - MID POINT OF THE INTERVAL (BL,BR)

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
static void l_q4aws (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *bl,
                Mfloat *br, Mfloat *alfa, Mfloat *beta, Mfloat ri[],
                Mfloat rj[], Mfloat rg[], Mfloat rh[], Mfloat *result,
                Mfloat *abserr, Mfloat *resasc, Mint *integr,
                Mint *nev)
#else
static void l_q4aws (f, a, b, bl, br, alfa, beta, ri, rj, rg, rh,
                result, abserr, resasc, integr, nev)
    Mfloat      (*f) (), *a, *b, *bl, *br, *alfa, *beta, ri[], rj[], rg[], rh[],
           *result, *abserr, *resasc;
    Mint       *integr, *nev;
#endif
{
    Mint        i, isym, _i, _r;
    Mfloat      _f0, centr, cheb12[13], cheb24[25], dc, factor, fix, fval[25],
                hlgth, res12, res24, resabs, u;
    *nev = 25;
    if (*bl == *a && ((*alfa != F_ZERO || *integr == 2) || *integr ==
	    4))
	goto L_10;
    if (*br == *b && ((*beta != F_ZERO || *integr == 3) || *integr ==
	    4))
	goto L_140;
    /*
     * IF A.GT.BL AND B.LT.BR, APPLY THE 15-POINT GAUSS-KRONROD SCHEME.
     */
    imsl_q8awo (f, l_q6aws, a, b, alfa, beta, integr, bl, br, result, abserr,
	&resabs, resasc);
    *nev = 15;
    goto L_270;
    /*
     * THIS PART OF THE PROGRAM IS EXECUTED ONLY IF A = BL. COMPUTE THE
     * CHEBYSHEV SERIES EXPANSION OF THE FOLLOWING FUNCTION F1 =
     * (0.5*(B+B-BR-A)-0.5*(BR-A)*X)**BETA A F(0.5*(BR-A)*X+0.5*(BR+A))
     */
L_10:
    hlgth = F_HALF * (*br - *bl);
    centr = F_HALF * (*br + *bl);
    fix = *b - centr;
    imsl_e1usr ("ON");
    _f0 = hlgth + centr;
    fval[0] = F_HALF * (*f) (_f0) * pow (fix - hlgth, *beta);
    fval[12] = (*f) (centr) * (pow (fix, *beta));
    fval[24] = F_HALF * (*f) (centr - hlgth) * pow (fix + hlgth, *beta);
    for (i = 2; i <= 12; i++) {
	u = hlgth * lv_x[i - 2];
	isym = 26 - i;
	fval[i - 1] = (*f) (u + centr) * pow (fix - u, *beta);
	fval[isym - 1] = (*f) (centr - u) * pow (fix + u, *beta);
    }
    imsl_e1usr ("OFF");
    factor = pow (hlgth, *alfa + F_ONE);
    *result = F_ZERO;
    *abserr = F_ZERO;
    res12 = F_ZERO;
    res24 = F_ZERO;
    if (*integr > 2)
	goto L_70;
    imsl_q7awo (lv_x, fval, cheb12, cheb24);
    /* INTEGR = 1 (OR 2) */
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * ri[i - 1];
	res24 += cheb24[i - 1] * ri[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * ri[i - 1];
    }
    if (*integr == 1)
	goto L_130;
    /* INTEGR = 2 */
    dc = log (*br - *bl);
    *result = res24 * dc;
    *abserr = fabs ((res24 - res12) * dc);
    res12 = F_ZERO;
    res24 = F_ZERO;
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * rg[i - 1];
	res24 += cheb24[i - 1] * rg[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * rg[i - 1];
    }
    goto L_130;
    /*
     * COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE FOLLOWING FUNCTION F4 =
     * F1*LOG(0.5*(B+B-BR-A)-0.5*(BR-A)*X))
     */
L_70:
    fval[0] *= log (fix - hlgth);
    fval[12] *= log (fix);
    fval[24] *= log (fix + hlgth);
    for (i = 2; i <= 12; i++) {
	u = hlgth * lv_x[i - 2];
	isym = 26 - i;
	fval[i - 1] *= log (fix - u);
	fval[isym - 1] *= log (fix + u);
    }
    imsl_q7awo (lv_x, fval, cheb12, cheb24);
    /* INTEGR = 3 (OR 4) */
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * ri[i - 1];
	res24 += cheb24[i - 1] * ri[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * ri[i - 1];
    }
    if (*integr == 3)
	goto L_130;
    /* INTEGR = 4 */
    dc = log (*br - *bl);
    *result = res24 * dc;
    *abserr = fabs ((res24 - res12) * dc);
    res12 = F_ZERO;
    res24 = F_ZERO;
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * rg[i - 1];
	res24 += cheb24[i - 1] * rg[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * rg[i - 1];
    }
L_130:
    *result = (*result + res24) * factor;
    *abserr = (*abserr + fabs (res24 - res12)) * factor;
    goto L_270;
    /*
     * THIS PART OF THE PROGRAM IS EXECUTED ONLY IF B = BR. COMPUTE THE
     * CHEBYSHEV SERIES EXPANSION OF THE FOLLOWING FUNCTION F2 =
     * (0.5*(B+BL-A-A)+0.5*(B-BL)*X)**ALFA A F(0.5*(B-BL)*X+0.5*(B+BL))
     */
L_140:
    hlgth = F_HALF * (*br - *bl);
    centr = F_HALF * (*br + *bl);
    fix = centr - *a;
    imsl_e1usr ("ON");
    fval[0] = F_HALF * (*f) (hlgth + centr) * pow (fix + hlgth, *alfa);
    fval[12] = (*f) (centr) * (pow (fix, *alfa));
    fval[24] = F_HALF * (*f) (centr - hlgth) * pow (fix - hlgth, *alfa);
    for (i = 2; i <= 12; i++) {
	u = hlgth * lv_x[i - 2];
	isym = 26 - i;
	fval[i - 1] = (*f) (u + centr) * pow (fix + u, *alfa);
	fval[isym - 1] = (*f) (centr - u) * pow (fix - u, *alfa);
    }
    imsl_e1usr ("OFF");
    factor = pow (hlgth, *beta + F_ONE);
    *result = F_ZERO;
    *abserr = F_ZERO;
    res12 = F_ZERO;
    res24 = F_ZERO;
    if (*integr == 2 || *integr == 4)
	goto L_200;
    /* INTEGR = 1 (OR 3) */
    imsl_q7awo (lv_x, fval, cheb12, cheb24);
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * rj[i - 1];
	res24 += cheb24[i - 1] * rj[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * rj[i - 1];
    }
    if (*integr == 1)
	goto L_260;
    /* INTEGR = 3 */
    dc = log (*br - *bl);
    *result = res24 * dc;
    *abserr = fabs ((res24 - res12) * dc);
    res12 = F_ZERO;
    res24 = F_ZERO;
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * rh[i - 1];
	res24 += cheb24[i - 1] * rh[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * rh[i - 1];
    }
    goto L_260;
    /*
     * COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE FOLLOWING FUNCTION F3 =
     * F2*LOG(0.5*(B-BL)*X+0.5*(B+BL-A-A)))
     */
L_200:
    fval[0] *= log (hlgth + fix);
    fval[12] *= log (fix);
    fval[24] *= log (fix - hlgth);
    for (i = 2; i <= 12; i++) {
	u = hlgth * lv_x[i - 2];
	isym = 26 - i;
	fval[i - 1] *= log (u + fix);
	fval[isym - 1] *= log (fix - u);
    }
    imsl_q7awo (lv_x, fval, cheb12, cheb24);
    /* INTEGR = 2 (OR 4) */
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * rj[i - 1];
	res24 += cheb24[i - 1] * rj[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * rj[i - 1];
    }
    if (*integr == 2)
	goto L_260;
    dc = log (*br - *bl);
    *result = res24 * dc;
    *abserr = fabs ((res24 - res12) * dc);
    res12 = F_ZERO;
    res24 = F_ZERO;
    /* INTEGR = 4 */
    for (i = 1; i <= 13; i++) {
	res12 += cheb12[i - 1] * rh[i - 1];
	res24 += cheb24[i - 1] * rh[i - 1];
    }
    for (i = 14; i <= 25; i++) {
	res24 += cheb24[i - 1] * rh[i - 1];
    }
L_260:
    *result = (*result + res24) * factor;
    *abserr = (*abserr + fabs (res24 - res12)) * factor;
L_270:
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 16:33:16 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 16:33:13
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q5AWS/DQ5AWS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q5AWS (ALFA, BETA, RI, RJ, RG, RH, INTEGR)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


  .......................................................................

   1.     Q5AWS
          MODIFIED CHEBYSHEV MOMENTS
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THIS ROUTINE COMPUTES MODIFIED CHEBYSHEV MOMENTS.
             THE K-TH MODIFIED CHEBYSHEV MOMENT IS DEFINED AS THE
             INTEGRAL OVER (-1,1) OF W(X)*T(K,X), WHERE T(K,X)
             IS THE CHEBYSHEV POLYNOMIAL OF DEGREE K.

   3.     CALLING SEQUENCE
             CALL Q5AWS(ALFA,BETA,RI,RJ,RG,RH,INTEGR)

          PARAMETERS
             ALFA   - REAL
                      PARAMETER IN THE WEIGHT FUNCTION W(X),
                      ALFA.GT.(-1)

             BETA   - REAL
                      PARAMETER IN THE WEIGHT FUNCTION W(X),
                      BETA.GT.(-1)

             RI     - REAL
                      VECTOR OF DIMENSION 25
                      RI(K) IS THE INTEGRAL OVER (-1,1) OF
                      (1+X)**ALFA*T(K-1,X), K = 1, ..., 25.

             RJ     - REAL
                      VECTOR OF DIMENSION 25
                      RJ(K) IS THE INTEGRAL OVER (-1,1) OF
                      (1-X)**BETA*T(K-1,X), K = 1, ..., 25.

             RG     - REAL
                      VECTOR OF DIMENSION 25
                      RG(K) IS THE INTEGRAL OVER (-1,1) OF
                      (1+X)**ALFA*LOG((1+X)/2)*T(K-1,X), K = 1, ..., 25.

             RH     - REAL
                      VECTOR OF DIMENSION 25
                      RH(K) IS THE INTEGRAL OVER (-1,1) OF
                      (1-X)**BETA*LOG((1-X)/2)*T(K-1,X), K = 1, ..., 25.

             INTEGR - INTEGER
                      INPUT PARAMETER INDICATING THE MODIFIED
                      MOMENTS TO BE COMPUTED
                      INTEGR = 1 COMPUTE RI, RJ
                             = 2 COMPUTE RI, RJ, RG
                             = 3 COMPUTE RI, RJ, RH
                             = 4 COMPUTE RI, RJ, RG, RH

   4.     NO SUBROUTINES OR FUNCTIONS NEEDED

  .......................................................................


  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q5aws (Mfloat *alfa, Mfloat *beta, Mfloat ri[], Mfloat rj[],
                Mfloat rg[], Mfloat rh[], Mint *integr)
#else
static void l_q5aws (alfa, beta, ri, rj, rg, rh, integr)
    Mfloat     *alfa, *beta, ri[], rj[], rg[], rh[];
    Mint       *integr;
#endif
{
    Mint        i, im1;
    Mfloat      alfp1, alfp2, an, anm1, betp1, betp2, ralf, rbet;


    alfp1 = *alfa + F_ONE;
    betp1 = *beta + F_ONE;
    alfp2 = *alfa + F_TWO;
    betp2 = *beta + F_TWO;
    ralf = pow (F_TWO, alfp1);
    rbet = pow (F_TWO, betp1);
    /*
     * COMPUTE RI, RJ USING A FORWARD RECURRENCE RELATION.
     */
    ri[0] = ralf / alfp1;
    rj[0] = rbet / betp1;
    ri[1] = ri[0] ** alfa / alfp2;
    rj[1] = rj[0] ** beta / betp2;
    an = F_TWO;
    anm1 = F_ONE;
    for (i = 3; i <= 25; i++) {
	ri[i - 1] = -(ralf + an * (an - alfp2) * ri[i - 2]) / (anm1 * (an +
		alfp1));
	rj[i - 1] = -(rbet + an * (an - betp2) * rj[i - 2]) / (anm1 * (an +
		betp1));
	anm1 = an;
	an += F_ONE;
    }
    if (*integr == 1)
	goto L_60;
    if (*integr == 3)
	goto L_30;
    /*
     * COMPUTE RG USING A FORWARD RECURRENCE RELATION.
     */
    rg[0] = -ri[0] / alfp1;
    rg[1] = -(ralf + ralf) / (alfp2 * alfp2) - rg[0];
    an = F_TWO;
    anm1 = F_ONE;
    im1 = 2;
    for (i = 3; i <= 25; i++) {
	rg[i - 1] = -(an * (an - alfp2) * rg[im1 - 1] - an * ri[im1 - 1] +
	    anm1 * ri[i - 1]) / (anm1 * (an + alfp1));
	anm1 = an;
	an += F_ONE;
	im1 = i;
    }
    if (*integr == 2)
	goto L_60;
    /*
     * COMPUTE RH USING A FORWARD RECURRENCE RELATION.
     */
L_30:
    rh[0] = -rj[0] / betp1;
    rh[1] = -(rbet + rbet) / (betp2 * betp2) - rh[0];
    an = F_TWO;
    anm1 = F_ONE;
    im1 = 2;
    for (i = 3; i <= 25; i++) {
	rh[i - 1] = -(an * (an - betp2) * rh[im1 - 1] - an * rj[im1 - 1] +
	    anm1 * rj[i - 1]) / (anm1 * (an + betp1));
	anm1 = an;
	an += F_ONE;
	im1 = i;
    }
    for (i = 2; i <= 25; i += 2) {
	rh[i - 1] = -rh[i - 1];
    }
L_60:
    for (i = 2; i <= 25; i += 2) {
	rj[i - 1] = -rj[i - 1];
    }
L_80:
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 16:33:58 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 16:33:56
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q6AWS/DQ6AWS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      Q6AWS(X, A, B, ALFA, BETA, INTEGR)

    Arguments:
       X      - Scalar.  (Input)
       A      - Scalar.  (Input)
       B      - Scalar.  (Input)
       ALFA   - Scalar.  (Input)
       BETA   - Scalar.  (Input)
       INTEGR - Scalar.  (Input)
       Q6AWS  - Weight function of QDAWS.  (Output)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_q6aws (Mfloat *x, Mfloat *a, Mfloat *b, Mfloat *alfa,
                Mfloat *beta, Mint *integr)
#else
static Mfloat l_q6aws (x, a, b, alfa, beta, integr)
    Mfloat     *x, *a, *b, *alfa, *beta;
    Mint       *integr;
#endif
{
    Mfloat      bmx, q6aws_v, xma;


    xma = *x - *a;
    bmx = *b - *x;
    q6aws_v = pow (xma, *alfa) * pow (bmx, *beta);
    if (*integr == 1)
	goto L_40;
    if (*integr == 3)
	goto L_20;
    if (*integr == 4)
	goto L_30;
L_10:
    q6aws_v *= log (xma);
    goto L_40;
L_20:
    q6aws_v *= log (bmx);
    goto L_40;
L_30:
    q6aws_v *= log (xma) * log (bmx);
L_40:
    return (q6aws_v);
}				/* end of function */
