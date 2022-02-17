#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                va_list argptr);
static void l_q2ag (Mfloat (*f) (Mfloat), Mfloat a, Mfloat b, Mfloat errabs,
                Mfloat errrel, int irule, Mfloat *result,
                Mfloat *errest, Mint maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[]);
static void l_q3ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *epsabs,
                Mfloat *epsrel, Mint *key, Mint *limit,
                Mfloat *result, Mfloat *abserr, Mint *neval,
                Mint *ier, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint *last);
static void l_q4ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc);
static void l_q5ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc);
static void l_q6ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc);
static void l_q7ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc);
static void l_q8ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc);
#else
static VA_LIST_HACK l_int_fcn ();
static void l_q2ag ();
static void l_q3ag ();
static void l_q4ag ();
static void l_q5ag ();
static void l_q6ag ();
static void l_q7ag ();
static void l_q8ag ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,...)
#else
Mfloat      imsl_f_int_fcn (fcn, a, b, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, b);
    E1PSH ("imsl_f_int_fcn", "imsl_d_int_fcn");
    lv_value = F_ZERO;
    IMSL_CALL (l_int_fcn (fcn, a, b, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn", "imsl_d_int_fcn");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn (fcn, a, b, argptr)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        arg_number = 3;
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
    err_abs = sqrt( (double) imsl_amach (4));
    err_rel = sqrt( (double) imsl_amach (4));

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
	    /* Argument number %(I2) is an unknown */
	    /* optional argument %(I1). */
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    break;
	}
    }
    /* CHECK MAXSUB */
    if (max_subinter < 1) {
	/* The maximum number of subintervals MAXSUB = %(i1). */
	/* It must be at least 1. */
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

    l_q2ag (fcn, a, b, err_abs, err_rel, rule, &lv_value, err_est,
	max_subinter, n_evals, n_subinter, alist, blist, rlist,
	elist, iord);
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



/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:10:04 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:10:01
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q2AG/DQ2AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q2AG (F, A, B, ERRABS, ERRREL, IRULE, RESULT,
                           ERREST, MAXSUB, NEVAL, NSUBIN, ALIST, BLIST,
                           RLIST, ELIST, IORD)

    Arguments:  (See QDAG)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q2ag (Mfloat (*f) (Mfloat), Mfloat a, Mfloat b, Mfloat errabs,
                Mfloat errrel, int l_irule, Mfloat *result,
                Mfloat *errest, Mint maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[])
#else
static void l_q2ag (f, a, b, errabs, errrel, l_irule, result, errest,
                maxsub, neval, nsubin, alist, blist, rlist, elist, iord)
    Mfloat      (*f) (), a, b, errabs, errrel;
    Mint        l_irule;
    Mfloat     *result, *errest;
    Mint        maxsub, *neval, *nsubin;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[];
#endif
{
    Mfloat      temp_a = a;
    Mfloat      temp_b = b;
    Mfloat      temp_errabs = errabs;
    Mfloat      temp_errrel = errrel;
    Mint        ier;
    Mint        irule;
    imsl_e1psh ("l_q2ag");
    irule = l_irule;
    /* CHECK MAXSUB */
    if (maxsub < 1) {
	/* The maximum number of subintervals MAXSUB = %(i1). */
	/* It must be at least 1. */
	imsl_e1sti (1, maxsub);
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_SUBINTER_SMALL);
    }
    /* CHECK IRULE */
    if (irule < 1 || irule > 6) {
	/* The quadrature rule IRULE = %(i1).  */
	/* It must be in the range 1 to 6. */
	imsl_e1sti (1, irule);
	imsl_ermes (IMSL_TERMINAL, IMSL_RULE_SMALL);
    }
    /* CHECK ERRABS */
    if (errabs < F_ZERO) {
	/* The absolute error desired err_abs = %(R1). */
	/* It must be at least zero. */
	imsl_e1str (1, errabs);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_ABS_SMALL);
    }
    /* CHECK ERRREL */
    if (errrel < F_ZERO) {
	/* The absolute error desired err_rel = %(R1). */
	/* It must be at least zero. */
	imsl_e1str (1, errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_SMALL);
    }
    /*
     * CHECK ERRABS AND ERRREL
     */
    if (errabs == F_ZERO && errrel == F_ZERO) {
	/* The error tolerance arguments err_abs and */
	/* and err_rel are both equal to zero. */
	/* At least one of them must be greater */
	/* than zero. */
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_TOL_ZERO);
    }
    /* CHECK ERRREL .GE. 1 */
    if (errrel >= F_ONE) {
	/* The relative error desired err_rel = %(R1). */
	/* When err_rel is greater than or equal to */
	/* one zero can always be returned as an */
	/* answer. */
	imsl_e1str (1, errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_BIG);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;

    l_q3ag (f, &temp_a, &temp_b, &temp_errabs, &temp_errrel, &irule,
	&maxsub, result, errest, neval, &ier, alist, blist, rlist, elist,
	iord, nsubin);

    if (ier == 1) {
	/* The maximum number of subintervals allowed */
	/* maxsub = %(i1) has been reached. */
	/* Increase maxsub. */
	imsl_e1sti (1, maxsub);
	imsl_ermes (IMSL_FATAL, IMSL_MAX_SUBINTERVALS);
    }
    else if (ier == 2) {
	/* Roundoff error has been detected.  The */
	/* requested tolerances, err_abs = %(r1) */
	/* and err_rel = %(r2) cannot be reached. */
	imsl_e1str (1, errabs);
	imsl_e1str (2, errrel);
	imsl_ermes (IMSL_WARNING, IMSL_ROUNDOFF_CONTAMINATION);
    }
    else if (ier == 3) {
	/* Precision is degraded due to too fine a */
	/* subdivision relative to the requested */
	/* tolerance.  This may be due to bad integrand */
	/* behavior in the interval (%(r1),%(r2)). */
	/* Higher precision may alleviate this problem. */
	imsl_e1str (1, alist[iord[0] - 1]);
	imsl_e1str (2, blist[iord[0] - 1]);
	imsl_ermes (IMSL_WARNING, IMSL_PRECISION_DEGRADATION);
    }
L_9000:
    ;
    imsl_e1pop ("l_q2ag");
    return;
}				/* end of function */


/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:13:15 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:13:12
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3AG/DQ3AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q3AG (F, A, B, EPSABS, EPSREL, KEY, LIMIT, RESULT,
                           ABSERR, NEVAL, IER, ALIST, BLIST, RLIST,
                           ELIST, IORD, LAST)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


   ......................................................................

   1.     Q3AG
          COMPUTATION OF A DEFINITE INTEGRAL
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
             DEFINITE INTEGRAL   I = INTEGRAL OF  F  OVER (A,B),
             HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

   3.     CALLING SEQUENCE
             CALL Q3AG(F,A,B,EPSABS,EPSREL,KEY,LIMIT,RESULT,ABSERR,
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

              EPSABS - REAL
                       ABSOLUTE ACCURACY REQUESTED
              EPSREL - REAL
                       RELATIVE ACCURACY REQUESTED
                       IF  EPSABS.LT.0 AND EPSREL.LT.0,
                       THE ROUTINE WILL END WITH IER = 6.

              KEY    - INTEGER
                       KEY FOR CHOICE OF LOCAL INTEGRATION RULE
                       A GAUSS-KRONROD PAIR IS USED WITH
                            7 - 15 POINTS IF KEY.LT.2,
                           10 - 21 POINTS IF KEY = 2,
                           15 - 31 POINTS IF KEY = 3,
                           20 - 41 POINTS IF KEY = 4,
                           25 - 51 POINTS IF KEY = 5,
                           30 - 61 POINTS IF KEY.GT.5.

              LIMIT  - INTEGER
                       GIVES AN UPPERBOUND ON THE NUMBER OF SUBINTERVALS
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
                               THE ESTIMATES FOR RESULT AND ERROR ARE
                               LESS RELIABLE. IT IS ASSUMED THAT THE
                               REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
                       IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
                               HAS BEEN ACHIEVED. ONE CAN ALLOW MORE
                               SUBDIVISIONS BY INCREASING THE VALUE
                               OF LIMIT.
                               HOWEVER, IF THIS YIELDS NO IMPROVEMENT IT
                               IS RATHER ADVISED TO ANALYZE THE INTEGRAND
                               IN ORDER TO DETERMINE THE INTEGRATION
                               DIFFICULTIES. IF THE POSITION OF A LOCAL
                               DIFFICULTY CAN BE DETERMINED(E.G.
                               SINGULARITY, DISCONTINUITY WITHIN THE
                               INTERVAL) ONE WILL PROBABLY GAIN FROM
                               SPLITTING UP THE INTERVAL AT THIS POINT
                               AND CALLING THE INTEGRATOR ON THE
                               SUBRANGES. IF POSSIBLE, AN APPROPRIATE
                               SPECIAL-PURPOSE INTEGRATOR SHOULD BE USED
                               WHICH IS DESIGNED FOR HANDLING THE TYPE OF
                               DIFFICULTY INVOLVED.
                           = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS
                               DETECTED, WHICH PREVENTS THE REQUESTED
                               TOLERANCE FROM BEING ACHIEVED.
                           = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS
                               AT SOME POINTS OF THE INTEGRATION
                               INTERVAL.
                           = 6 THE INPUT IS INVALID, BECAUSE
                               EPSABS.LT.0 AND EPSREL.LT.0,
                               RESULT, ABSERR, NEVAL, LAST, RLIST(1),
                               ELIST(1) AND IORD(1) ARE SET TO ZERO.
                               ALIST(1) AND BLIST(1) ARE SET TO A AND B
                               RESPECTIVELY.

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
                        VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                         LAST  ELEMENTS OF WHICH ARE THE MODULI OF THE
                        ABSOLUTE ERROR ESTIMATES ON THE SUBINTERVALS

              IORD    - INTEGER
                        VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST K
                        ELEMENTS OF WHICH ARE POINTERS TO THE ERROR
                        ESTIMATES OVER THE SUBINTERVALS, SUCH THAT
                        ELIST(IORD(1)), ..., ELIST(IORD(K)) FORM A
                        DECREASING SEQUENCE, WITH
                        K = LAST IF LAST.LE.(LIMIT/2+2), AND
                        K = LIMIT+1-LAST OTHERWISE

              LAST    - INTEGER
                        NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN THE
                        SUBDIVISION PROCESS

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q4AG, Q9AG, Q5AG, Q6AG, Q7AG, Q8AG
                - Q10G
                - F (USER-PROVIDED FUNCTION)
                - Q4NG
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

             EPMACH  IS THE LARGEST RELATIVE SPACING.
             UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW  IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q3ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *epsabs,
                Mfloat *epsrel, Mint *key, Mint *limit,
                Mfloat *result, Mfloat *abserr, Mint *neval,
                Mint *ier, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint *last)
#else
static void l_q3ag (f, a, b, epsabs, epsrel, key, limit, result,
                abserr, neval, ier, alist, blist, rlist, elist, iord, last)
    Mfloat      (*f) (), *a, *b, *epsabs, *epsrel;
    Mint       *key, *limit;
    Mfloat     *result, *abserr;
    Mint       *neval, *ier;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[], *last;
#endif
{
    Mint        iroff1, iroff2, k, keyf, maxerr, nrmax;
    Mfloat      a1, a2, area, area1, area12, area2, b1, b2, c;
    Mfloat      defab1, defab2, defabs, epmach, errbnd, errmax;
    Mfloat      erro12, error1, error2, errsum, oflow, resabs, uflow;
    imsl_q4ng (&epmach, &uflow, &oflow);

    /*
     * TEST ON VALIDITY OF PARAMETERS
     */
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = F_ZERO;
    *abserr = F_ZERO;
    alist[0] = *a;
    blist[0] = *b;
    rlist[0] = F_ZERO;
    elist[0] = F_ZERO;
    iord[0] = 0;
    if (*epsabs < F_ZERO && *epsrel < F_ZERO)
	*ier = 6;
    if (*ier == 6)
	goto L_90;

    /*
     * FIRST APPROXIMATION TO THE INTEGRAL
     */
    keyf = *key;
    if (*key <= 0)
	keyf = 1;
    if (*key >= 7)
	keyf = 6;
    c = keyf;
    *neval = 0;
    if (keyf == 1)
	l_q4ag (f, a, b, result, abserr, &defabs, &resabs);
    if (keyf == 2)
	imsl_q9ag (f, a, b, result, abserr, &defabs, &resabs);
    if (keyf == 3)
	l_q5ag (f, a, b, result, abserr, &defabs, &resabs);
    if (keyf == 4)
	l_q6ag (f, a, b, result, abserr, &defabs, &resabs);
    if (keyf == 5)
	l_q7ag (f, a, b, result, abserr, &defabs, &resabs);
    if (keyf == 6)
	l_q8ag (f, a, b, result, abserr, &defabs, &resabs);
    *last = 1;
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 1;

    /*
     * TEST ON ACCURACY.
     */
    errbnd = imsl_f_max (*epsabs, *epsrel * fabs (*result));
    if (*abserr <= 5.0e01 * epmach * defabs && *abserr > errbnd)
	*ier = 2;
    if (*limit == 1)
	*ier = 1;
    if ((*ier != 0 || (*abserr <= errbnd && *abserr != resabs)) ||
	*abserr == F_ZERO)
	goto L_80;

    /*
     * INITIALIZATION
     */
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    nrmax = 1;
    iroff1 = 0;
    iroff2 = 0;

    /*
     * MAIN DO-LOOP
     */
    for (*last = 2; *last <= *limit; (*last)++) {

	/*
	 * BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE
	 */
	a1 = alist[maxerr - 1];
	b1 = F_HALF * (alist[maxerr - 1] + blist[maxerr - 1]);
	a2 = b1;
	b2 = blist[maxerr - 1];
	if (keyf == 1)
	    l_q4ag (f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	if (keyf == 2)
	    imsl_q9ag (f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	if (keyf == 3)
	    l_q5ag (f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	if (keyf == 4)
	    l_q6ag (f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	if (keyf == 5)
	    l_q7ag (f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	if (keyf == 6)
	    l_q8ag (f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	if (keyf == 1)
	    l_q4ag (f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	if (keyf == 2)
	    imsl_q9ag (f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	if (keyf == 3)
	    l_q5ag (f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	if (keyf == 4)
	    l_q6ag (f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	if (keyf == 5)
	    l_q7ag (f, &a2, &b2, &area2, &error2, &resabs, &defab2);
	if (keyf == 6)
	    l_q8ag (f, &a2, &b2, &area2, &error2, &resabs, &defab2);

	/*
	 * IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR; TEST FOR
	 * ACCURACY
	 */
	*neval += 1;
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum += erro12 - errmax;
	area += area12 - rlist[maxerr - 1];
	if (defab1 == error1 || defab2 == error2)
	    goto L_10;
	if (fabs (rlist[maxerr - 1] - area12) <= 1.0e-05 * fabs (area12) &&
	    erro12 >= 9.9e-01 * errmax)
	    iroff1 += 1;
	if (*last > 10 && erro12 > errmax)
	    iroff2 += 1;
L_10:
	rlist[maxerr - 1] = area1;
	rlist[*last - 1] = area2;
	errbnd = imsl_f_max (*epsabs, *epsrel * fabs (area));
	if (errsum <= errbnd)
	    goto L_20;

	/*
	 * TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG
	 */
	if (iroff1 >= 6 || iroff2 >= 20)
	    *ier = 2;

	/*
	 * SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS EQUALS
	 * LIMIT
	 */
	if (*last == *limit)
	    *ier = 1;

	/*
	 * SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOR AT A POINT OF
	 * THE INTEGRATION RANGE
	 */
	if (imsl_f_max (fabs (a1), fabs (b2)) <= (F_ONE + c * 1.0e03 * epmach) *
	    (fabs (a2) + 1.0e04 * uflow))
	    *ier = 3;

	/*
	 * APPEND THE NEWLY-CREATED INTERVALS TO THE LIST
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
	 * LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL WITH THE
	 * LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT)
	 */
L_40:
	imsl_q10g (limit, last, &maxerr, &errmax, elist, iord, &nrmax);
	/* JUMP OUT OF DO-LOOP */
	if (*ier != 0 || errsum <= errbnd)
	    goto L_60;
    }

    /*
     * COMPUTE FINAL RESULT.
     */
L_60:
    *result = F_ZERO;
    for (k = 1; k <= *last; k++) {
	*result += rlist[k - 1];
    }
    *abserr = errsum;
L_80:
    if (keyf != 1)
	*neval = (10 * keyf + 1) * (2 ** neval + 1);
    if (keyf == 1)
	*neval = 30 ** neval + 15;
L_90:
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:15:34 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:15:28
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q4AG/DQ4AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q4AG (F, A, B, RESULT, ABSERR, RESABS, RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q4AG
             INTEGRATION RULES
                STANDARD FORTRAN SUBROUTINE
                REAL VERSION

   2.        PURPOSE
                TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR
                               ESTIMATE
                           J = INTEGRAL OF ABS(F) OVER (A,B)

   3.        CALLING SEQUENCE
                CALL Q4AG(F,A,B,RESULT,ABSERR,RESABS,RESASC)

             PARAMETERS
              ON ENTRY
                F      - REAL
                         FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                         FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                         DECLARED E X T E R N A L IN THE CALLING PROGRAM.

                A      - REAL
                         LOWER LIMIT OF INTEGRATION

                B      - REAL
                         UPPER LIMIT OF INTEGRATION

              ON RETURN
                RESULT - REAL
                         APPROXIMATION TO THE INTEGRAL I
                         RESULT IS COMPUTED BY APPLYING THE 15-POINT
                         KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
                         OF ABSCISSAE TO THE 7-POINT GAUSS RULE(RESG).

                ABSERR - REAL
                         ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                         WHICH SHOULD NOT EXCEED ABS(I-RESULT)

                RESABS - REAL
                         APPROXIMATION TO THE INTEGRAL J

                RESASC - REAL
                         APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
                         OVER (A,B)

   4.        SUBROUTINES OR FUNCTIONS NEEDED
                   - F (USER-PROVIDED FUNCTION)
                   - Q4NG
                   - FORTRAN ABS, AMAX1, AMIN1

       ..................................................................



             THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
             BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
             CORRESPONDING WEIGHTS ARE GIVEN.

             XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
                      XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
                      GAUSS RULE
                      XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
                      ADDED TO THE 7-POINT GAUSS RULE

             WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE

             WG     - WEIGHTS OF THE 7-POINT GAUSS RULE


             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC   - ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
             RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
             RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
                      I.E. TO I/(B-A)

             MACHINE DEPENDENT CONSTANTS
             ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q4ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc)
#else
static void /* FUNCTION */ l_q4ag (f, a, b, result, abserr, resabs, resasc)
    Mfloat      (*f) (), *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    Mint        j, jtw, jtwm1;
    Mfloat      _f0, absc, centr, dhlgth, epmach, fc, fsum, fv1[7], fv2[7],
                fval1, fval2, hlgth, oflow, resg, resk, reskh, uflow;
    static Mfloat xgk[] = {
	0.991455371120812639206854697526e0,
	0.949107912342758524526189684048e0,
	0.864864423359769072789712788641e0,
	0.741531185599394439863864773281e0,
	0.586087235467691130294144838259e0,
	0.405845151377397166906606412077e0,
	0.207784955007898467600689403773e0,
    0.0e0};
    static Mfloat wgk[] = {
	0.022935322010529224963732008059e0,
	0.0630920926299785532907006631892e0,
	0.104790010322250183839876322542e0,
	0.14065325971552591874518959051e0,
	0.169004726639267902826583426599e0,
	0.190350578064785409913256402421e0,
	0.204432940075298892414161999235e0,
    0.209482141084727828012999174892e0};
    static Mfloat wg[] = {
	0.129484966168869693270611432679e0,
	0.279705391489276667901467771424e0,
	0.381830050505118944950369775489e0,
    0.417959183673469387755102040816e0};
    imsl_q4ng (&epmach, &uflow, &oflow);

    centr = F_HALF * (*a + *b);
    hlgth = F_HALF * (*b - *a);
    dhlgth = fabs (hlgth);
    /*
     * COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ABSOLUTE ERROR.
     */
    imsl_e1usr ("ON");
    fc = (*f) (centr);
    imsl_e1usr ("OFF");
    resg = fc * wg[3];
    resk = fc * wgk[7];
    *resabs = fabs (resk);
    for (j = 1; j <= 3; j++) {
	jtw = j * 2;
	absc = hlgth * xgk[jtw - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (fabs (fval1) + fabs (fval2));
    }
    for (j = 1; j <= 4; j++) {
	jtwm1 = j * 2 - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = wgk[7] * fabs (fc - reskh);
    for (j = 1; j <= 7; j++) {
	*resasc += wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
		reskh));
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs ((resk - resg) * hlgth);
    if (*resasc != F_ZERO && *abserr != F_ZERO)
	*abserr = *resasc * imsl_f_min (F_ONE, pow (2.0e02 ** abserr / *resasc, 1.5e00));
    if (*resabs > uflow / (5.0e01 * epmach))
	*abserr = imsl_f_max ((epmach * 5.0e01) ** resabs, *abserr);
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:17:55 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:17:51
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q5AG/DQ5AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q5AG (F, A, B, RESULT, ABSERR, RESABS, RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q5AG
             INTEGRATION RULES
                STANDARD FORTRAN SUBROUTINE
                REAL VERSION

   2.        PURPOSE
                TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR
                               ESTIMATE
                           J = INTEGRAL OF ABS(F) OVER (A,B)

   3.        CALLING SEQUENCE
                CALL Q5AG(F,A,B,RESULT,ABSERR,RESABS,RESASC)

             PARAMETERS
              ON ENTRY
                F      - REAL
                         FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                         FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                         DECLARED E X T E R N A L IN THE CALLING PROGRAM.

                A      - REAL
                         LOWER LIMIT OF INTEGRATION

                B      - REAL
                         UPPER LIMIT OF INTEGRATION

              ON RETURN
                RESULT - REAL
                         APPROXIMATION TO THE INTEGRAL I
                         RESULT IS COMPUTED BY APPLYING THE 31-POINT
                         GAUSS-KRONROD RULE (RESK), OBTAINED BY OPTIMAL
                         ADDITION OF ABSCISSAE TO THE 15-POINT GAUSS
                         RULE (RESG).

                ABSERR - REAL
                         ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                         WHICH SHOULD NOT EXCEED ABS(I-RESULT)

                RESABS - REAL
                         APPROXIMATION TO THE INTEGRAL J

                RESASC - REAL
                         APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
                         OVER (A,B)

   4.        SUBROUTINES OR FUNCTIONS NEEDED
                   - F (USER-PROVIDED FUNCTION)
                   - Q4NG
                   - FORTRAN ABS, AMIN1, AMAX1

       ..................................................................



             THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
             BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
             CORRESPONDING WEIGHTS ARE GIVEN.

             XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
                      XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
                      GAUSS RULE
                      XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
                      ADDED TO THE 15-POINT GAUSS RULE

             WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE

             WG     - WEIGHTS OF THE 15-POINT GAUSS RULE



             LIST OF MAJOR VARIABLES
             -----------------------
             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC   - ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
             RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
             RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
                      I.E. TO I/(B-A)

             MACHINE DEPENDENT CONSTANTS
             ---------------------------
             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q5ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc)
#else
static void /* FUNCTION */ l_q5ag (f, a, b, result, abserr, resabs, resasc)
    Mfloat      (*f) (), *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    Mint        j, jtw, jtwm1;
    Mfloat      _f0, absc, centr, dhlgth, epmach, fc, fsum, fv1[15], fv2[15],
                fval1, fval2, hlgth, oflow, resg, resk, reskh, uflow;
    static Mfloat xgk[] = {
	0.998002298693397060285172840152e0,
	0.987992518020485428489565718587e0,
	0.967739075679139134257347978784e0,
	0.93727339240070590430775894771e0,
	0.897264532344081900882509656455e0,
	0.848206583410427216200648320774e0,
	0.790418501442465932967649294818e0,
	0.724417731360170047416186054614e0,
	0.650996741297416970533735895313e0,
	0.570972172608538847537226737254e0,
	0.485081863640239680693655740232e0,
	0.394151347077563369897207370981e0,
	0.299180007153168812166780024266e0,
	0.201194093997434522300628303395e0,
	0.101142066918717499027074231447e0,
    0.0e0};
    static Mfloat wgk[] = {
	0.00537747987292334898779205143013e0,
	0.0150079473293161225383747630758e0,
	0.0254608473267153201868740010197e0,
	0.0353463607913758462220379484784e0,
	0.0445897513247648766082272993733e0,
	0.0534815246909280872653431472394e0,
	0.0620095678006706402851392309608e0,
	0.0698541213187282587095200770992e0,
	0.0768496807577203788944327774827e0,
	0.0830805028231330210382892472861e0,
	0.0885644430562117706472754436938e0,
	0.0931265981708253212254868727474e0,
	0.0966427269836236785051799076276e0,
	0.0991735987217919593323931734846e0,
	0.100769845523875595044946662618e0,
    0.101330007014791549017374792768e0};
    static Mfloat wg[] = {
	0.0307532419961172683546283935772e0,
	0.0703660474881081247092674164507e0,
	0.107159220467171935011869546686e0,
	0.139570677926154314447804794511e0,
	0.166269205816993933553200860481e0,
	0.186161000015562211026800561866e0,
	0.198431485327111576456118326444e0,
    0.202578241925561272880620199968e0};
    imsl_q4ng (&epmach, &uflow, &oflow);
    centr = F_HALF * (*a + *b);
    hlgth = F_HALF * (*b - *a);
    dhlgth = fabs (hlgth);
    /*
     * COMPUTE THE 31-POINT KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ABSOLUTE ERROR
     */
    imsl_e1usr ("ON");
    fc = (*f) (centr);
    imsl_e1usr ("OFF");
    resg = wg[7] * fc;
    resk = wgk[15] * fc;
    *resabs = fabs (resk);
    for (j = 1; j <= 7; j++) {
	jtw = j * 2;
	absc = hlgth * xgk[jtw - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (fabs (fval1) + fabs (fval2));
    }
    for (j = 1; j <= 8; j++) {
	jtwm1 = j * 2 - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = wgk[15] * fabs (fc - reskh);
    for (j = 1; j <= 15; j++) {
	*resasc += wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
		reskh));
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs ((resk - resg) * hlgth);
    if (*resasc != F_ZERO && *abserr != F_ZERO)
	*abserr = *resasc * imsl_f_min (F_ONE, pow (2.0e02 ** abserr / *resasc, 1.5e00));
    if (*resabs > uflow / (5.0e01 * epmach))
	*abserr = imsl_f_max ((epmach * 5.0e01) ** resabs, *abserr);
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:19:22 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:19:18
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q6AG/DQ6AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q6AG (F, A, B, RESULT, ABSERR, RESABS, RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q6AG
             INTEGRATION RULES
                STANDARD FORTRAN SUBROUTINE
                REAL VERSION

   2.        PURPOSE
                TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR
                               ESTIMATE
                           J = INTEGRAL OF ABS(F) OVER (A,B)

   3.        CALLING SEQUENCE
                CALL Q6AG(F,A,B,RESULT,ABSERR,RESABS,RESASC)

             PARAMETERS
              ON ENTRY
                F      - REAL
                         FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                         FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                         DECLARED E X T E R N A L IN THE CALLING PROGRAM.

                A      - REAL
                         LOWER LIMIT OF INTEGRATION

                B      - REAL
                         UPPER LIMIT OF INTEGRATION

              ON RETURN
                RESULT - REAL
                         APPROXIMATION TO THE INTEGRAL I
                         RESULT IS COMPUTED BY APPLYING THE 41-POINT
                         GAUSS-KRONROD RULE (RESK) OBTAINED BY OPTIMAL
                         ADDITION OF ABSCISSAE TO THE 20-POINT GAUSS
                         RULE (RESG).

                ABSERR - REAL
                         ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                         WHICH SHOULD NOT EXCEED ABS(I-RESULT)

                RESABS - REAL
                         APPROXIMATION TO THE INTEGRAL J

                RESASC - REAL
                         APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
                         OVER (A,B)

   4.        SUBROUTINES OR FUNCTIONS NEEDED
                   - F (USER-PROVIDED FUNCTION)
                   - Q4NG
                   - FORTRAN ABS, AMAX1, AMIN1

       ..................................................................



             THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
             BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
             CORRESPONDING WEIGHTS ARE GIVEN.

             XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
                      XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT
                      GAUSS RULE
                      XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
                      ADDED TO THE 20-POINT GAUSS RULE

             WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE

             WG     - WEIGHTS OF THE 20-POINT GAUSS RULE



             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC   - ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 20-POINT GAUSS FORMULA
             RESK   - RESULT OF THE 41-POINT KRONROD FORMULA
             RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E.
                      TO I/(B-A)

             MACHINE DEPENDENT CONSTANTS
             ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q6ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc)
#else
static void /* FUNCTION */ l_q6ag (f, a, b, result, abserr, resabs, resasc)
    Mfloat      (*f) (), *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    Mint        j, jtw, jtwm1;
    Mfloat      _f0, absc, centr, dhlgth, epmach, fc, fsum, fv1[20], fv2[20],
                fval1, fval2, hlgth, oflow, resg, resk, reskh, uflow;
    static Mfloat xgk[] = {
	0.998859031588277663838315576546e0,
	0.993128599185094924786122388471e0,
	0.98150787745025025919334299472e0,
	0.963971927277913791267666131197e0,
	0.940822633831754753519982722212e0,
	0.912234428251325905867752441203e0,
	0.878276811252281976077442995113e0,
	0.839116971822218823394529061702e0,
	0.795041428837551198350638833273e0,
	0.746331906460150792614305070356e0,
	0.693237656334751384805490711846e0,
	0.636053680726515025452836696226e0,
	0.575140446819710315342946036586e0,
	0.510867001950827098004364050955e0,
	0.443593175238725103199992213493e0,
	0.373706088715419560672548177025e0,
	0.301627868114913004320555356859e0,
	0.227785851141645078080496195369e0,
	0.152605465240922675505220241023e0,
	0.0765265211334973337546404093988e0,
    0.0e0};
    static Mfloat wgk[] = {
	0.00307358371852053150121829324603e0,
	0.0086002698556429421986617879501e0,
	0.0146261692569712529837879603089e0,
	0.0203883734612665235980102314328e0,
	0.0258821336049511588345050670962e0,
	0.0312873067770327989585431193238e0,
	0.0366001697582007980305572407072e0,
	0.0416688733279736862637883059369e0,
	0.0464348218674976747202318809261e0,
	0.0509445739237286919327076700504e0,
	0.0551951053482859947448323724198e0,
	0.0591114008806395723749672206486e0,
	0.0626532375547811680258701221743e0,
	0.0658345971336184221115635569694e0,
	0.0686486729285216193456234118854e0,
	0.0710544235534440683057903617232e0,
	0.0730306903327866674951894176589e0,
	0.0745828754004991889865814183625e0,
	0.0757044976845566746595427753766e0,
	0.0763778676720807367055028350381e0,
    0.0766007119179996564450499015301e0};
    static Mfloat wg[] = {
	0.0176140071391521183118619623519e0,
	0.0406014298003869413310399522749e0,
	0.062672048334109063569506535187e0,
	0.0832767415767047487247581432221e0,
	0.10193011981724043503675013548e0,
	0.118194531961518417312377377711e0,
	0.131688638449176626898494499748e0,
	0.142096109318382051329298325067e0,
	0.149172986472603746787828737002e0,
    0.152753387130725850698084331955e0};
    imsl_q4ng (&epmach, &uflow, &oflow);

    centr = F_HALF * (*a + *b);
    hlgth = F_HALF * (*b - *a);
    dhlgth = fabs (hlgth);
    /*
     * COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ABSOLUTE ERROR
     */
    resg = F_ZERO;
    imsl_e1usr ("ON");
    fc = (*f) (centr);
    imsl_e1usr ("OFF");
    resk = wgk[20] * fc;
    *resabs = fabs (resk);
    for (j = 1; j <= 10; j++) {
	jtw = j * 2;
	absc = hlgth * xgk[jtw - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (fabs (fval1) + fabs (fval2));
    }
    for (j = 1; j <= 10; j++) {
	jtwm1 = j * 2 - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = wgk[20] * fabs (fc - reskh);
    for (j = 1; j <= 20; j++) {
	*resasc += wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
		reskh));
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs ((resk - resg) * hlgth);
    if (*resasc != F_ZERO && *abserr != F_ZERO)
	*abserr = *resasc * imsl_f_min (F_ONE, pow (2.0e02 ** abserr / *resasc, 1.5e00));
    if (*resabs > uflow / (5.0e01 * epmach))
	*abserr = imsl_f_max ((epmach * 5.0e01) ** resabs, *abserr);
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:20:56 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:20:52
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q7AG/DQ7AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q7AG (F, A, B, RESULT, ABSERR, RESABS, RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q7AG
             INTEGRATION RULES
               STANDARD FORTRAN SUBROUTINE
               REAL VERSION

   2.        PURPOSE
                TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR
                               ESTIMATE
                           J = INTEGRAL OF ABS(F) OVER (A,B)

   3.        CALLING SEQUENCE
                CALL Q7AG(F,A,B,RESULT,ABSERR,RESABS,RESASC)

             PARAMETERS
              ON ENTRY
                F      - REAL
                         FUNCTION SUBROUTINE DEFINING THE INTEGRAND
                         FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                         DECLARED E X T E R N A L IN THE CALLING PROGRAM.

                A      - REAL
                         LOWER LIMIT OF INTEGRATION

                B      - REAL
                         UPPER LIMIT OF INTEGRATION

              ON RETURN
                RESULT - REAL
                         APPROXIMATION TO THE INTEGRAL I
                         RESULT IS COMPUTED BY APPLYING THE 51-POINT
                         KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
                         OF ABSCISSAE TO THE 25-POINT GAUSS RULE (RESG).

                ABSERR - REAL
                         ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                         WHICH SHOULD NOT EXCEED ABS(I-RESULT)

                RESABS - REAL
                         APPROXIMATION TO THE INTEGRAL J

                RESASC - REAL
                         APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
                         OVER (A,B)

   4.        SUBROUTINES OR FUNCTIONS NEEDED
                   - F (USER-PROVIDED FUNCTION)
                   - Q4NG
                   - FORTRAN ABS, AMAX1, AMIN1

       ..................................................................



             THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
             BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
             CORRESPONDING WEIGHTS ARE GIVEN.

             XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE
                      XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT
                      GAUSS RULE
                      XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
                      ADDED TO THE 25-POINT GAUSS RULE

             WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE

             WG     - WEIGHTS OF THE 25-POINT GAUSS RULE



             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC   - ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 25-POINT GAUSS FORMULA
             RESK   - RESULT OF THE 51-POINT KRONROD FORMULA
             RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
                      I.E. TO I/(B-A)

             MACHINE DEPENDENT CONSTANTS
             ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q7ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc)
#else
static void /* FUNCTION */ l_q7ag (f, a, b, result, abserr, resabs, resasc)
    Mfloat      (*f) (), *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    Mint        j, jtw, jtwm1;
    Mfloat      _f0, absc, centr, dhlgth, epmach, fc, fsum, fv1[25], fv2[25],
                fval1, fval2, hlgth, oflow, resg, resk, reskh, uflow;
    static Mfloat xgk[] = {
	0.99926210499260983419345748654e0,
	0.995556969790498097908784946894e0,
	0.988035794534077247637331014577e0,
	0.97666392145951751149831538648e0,
	0.96161498642584251241813003366e0,
	0.942974571228974339414011169659e0,
	0.920747115281701561746346084546e0,
	0.894991997878275368851042006783e0,
	0.865847065293275595448996969588e0,
	0.833442628760834001421021108694e0,
	0.797873797998500059410410904994e0,
	0.759259263037357630577282865204e0,
	0.717766406813084388186654079773e0,
	0.673566368473468364485120633248e0,
	0.626810099010317412788122681625e0,
	0.577662930241222967723689841613e0,
	0.526325284334719182599623778158e0,
	0.473002731445714960522182115009e0,
	0.417885382193037748851814394595e0,
	0.361172305809387837735821730128e0,
	0.30308953893110783016747890998e0,
	0.243866883720988432045190362798e0,
	0.18371893942104889201596988876e0,
	0.122864692610710396387359818808e0,
	0.0615444830056850788865463923668e0,
    0.0e0};
    static Mfloat wgk[] = {
	0.00198738389233031592650785188284e0,
	0.00556193213535671375804023690107e0,
	0.00947397338617415160720771052366e0,
	0.013236229195571674813656405847e0,
	0.0168478177091282982315166675363e0,
	0.0204353711458828354565682922359e0,
	0.0240099456069532162200924891649e0,
	0.0274753175878517378029484555178e0,
	0.0307923001673874888911090202152e0,
	0.0340021302743293378367487952296e0,
	0.0371162714834155435603306253676e0,
	0.0400838255040323820748392844671e0,
	0.0428728450201700494768957924395e0,
	0.0455029130499217889098705847527e0,
	0.0479825371388367139063922557569e0,
	0.0502776790807156719633252594334e0,
	0.0523628858064074758643667121379e0,
	0.0542511298885454901445433704599e0,
	0.0559508112204123173082406863828e0,
	0.0574371163615678328535826939395e0,
	0.0586896800223942079619741758568e0,
	0.0597203403241740599790992919326e0,
	0.0605394553760458629453602675176e0,
	0.0611285097170530483058590304163e0,
	0.0614711898714253166615441319653e0,
    0.0615808180678329350787598242401e0};
    static Mfloat wg[] = {
	0.0113937985010262879479029641132e0,
	0.0263549866150321372619018152953e0,
	0.0409391567013063126556234877117e0,
	0.0549046959758351919259368915405e0,
	0.0680383338123569172071871856567e0,
	0.0801407003350010180132349596691e0,
	0.0910282619829636498114972207029e0,
	0.100535949067050644202206890393e0,
	0.10851962447426365311609395705e0,
	0.11485825914571164833932554587e0,
	0.119455763535784772228178126513e0,
	0.122242442990310041688959518946e0,
    0.123176053726715451203902873079e0};
    imsl_q4ng (&epmach, &uflow, &oflow);

    centr = F_HALF * (*a + *b);
    hlgth = F_HALF * (*b - *a);
    dhlgth = fabs (hlgth);
    /*
     * COMPUTE THE 51-POINT KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ABSOLUTE ERROR
     */
    imsl_e1usr ("ON");
    fc = (*f) (centr);
    imsl_e1usr ("OFF");
    resg = wg[12] * fc;
    resk = wgk[25] * fc;
    *resabs = fabs (resk);
    for (j = 1; j <= 12; j++) {
	jtw = j * 2;
	absc = hlgth * xgk[jtw - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (fabs (fval1) + fabs (fval2));
    }
    for (j = 1; j <= 13; j++) {
	jtwm1 = j * 2 - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = wgk[25] * fabs (fc - reskh);
    for (j = 1; j <= 25; j++) {
	*resasc += wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
		reskh));
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs ((resk - resg) * hlgth);
    if (*resasc != F_ZERO && *abserr != F_ZERO)
	*abserr = *resasc * imsl_f_min (F_ONE, pow (2.0e02 ** abserr / *resasc, 1.5e00));
    if (*resabs > uflow / (5.0e01 * epmach))
	*abserr = imsl_f_max ((epmach * 5.0e01) ** resabs, *abserr);
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:22:34 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:22:30
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q8AG/DQ8AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q8AG (F, A, B, RESULT, ABSERR, RESABS, RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


  .......................................................................

   1.     Q8AG
          INTEGRATION RULE
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR
                            ESTIMATE
                        J = INTEGRAL OF ABS(F) OVER (A,B)

   3.     CALLING SEQUENCE
             CALL Q8AG(F,A,B,RESULT,ABSERR,RESABS,RESASC)

          PARAMETERS
           ON ENTRY
             F      - REAL
                      FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                      FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                      DECLARED E X T E R N A L IN THE CALLING PROGRAM.

             A      - REAL
                      LOWER LIMIT OF INTEGRATION

             B      - REAL
                      UPPER LIMIT OF INTEGRATION

           ON RETURN
             RESULT - REAL
                      APPROXIMATION TO THE INTEGRAL I
                      RESULT IS COMPUTED BY APPLYING THE 61-POINT
                      KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
                      OF ABSCISSAE TO THE 30-POINT GAUSS RULE (RESG).

             ABSERR - REAL
                      ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                      WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

             RESABS - REAL
                      APPROXIMATION TO THE INTEGRAL J

             RESASC - REAL
                      APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))


   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - F (USER-PROVIDED FUNCTION)
                - Q4NG
                - FORTRAN ABS, AMAX1, AMIN1

  .....................................................................



             THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
             INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
             ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.

             XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
                     XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
                     GAUSS RULE
                     XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
                     TO THE 30-POINT GAUSS RULE

             WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE

             WG    - WEIGTHS OF THE 30-POINT GAUSS RULE


             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC   - ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 30-POINT GAUSS RULE
             RESK   - RESULT OF THE 61-POINT KRONROD RULE
             RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
                      OVER (A,B), I.E. TO I/(B-A)

             MACHINE DEPENDENT CONSTANTS
             ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q8ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc)
#else
static void /* FUNCTION */ l_q8ag (f, a, b, result, abserr, resabs, resasc)
    Mfloat      (*f) (), *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    Mint        j, jtw, jtwm1;
    Mfloat      _f0, absc, centr, dhlgth, epmach, fc, fsum, fv1[30], fv2[30],
                fval1, fval2, hlgth, oflow, resg, resk, reskh, uflow;
    static Mfloat xgk[] = {
	0.999484410050490637571325895706e0,
	0.996893484074649540271630050919e0,
	0.99163099687040459485862836611e0,
	0.983668123279747209970032581606e0,
	0.973116322501126268374693868424e0,
	0.960021864968307512216871025582e0,
	0.944374444748559979415831324037e0,
	0.926200047429274325879324277081e0,
	0.905573307699907798546522558926e0,
	0.88256053579205268154311646253e0,
	0.857205233546061098958658510659e0,
	0.829565762382768397442898119733e0,
	0.799727835821839083013668942323e0,
	0.767777432104826194917977340975e0,
	0.73379006245322680472617113137e0,
	0.697850494793315796932292388027e0,
	0.660061064126626961370053668149e0,
	0.620526182989242861140477556431e0,
	0.579345235826361691756024932173e0,
	0.536624148142019899264169793311e0,
	0.492480467861778574993693061208e0,
	0.447033769538089176780609900323e0,
	0.400401254830394392535476211543e0,
	0.352704725530878113471037207089e0,
	0.304073202273625077372677107199e0,
	0.254636926167889846439805129818e0,
	0.204525116682309891438957671002e0,
	0.153869913608583546963794672743e0,
	0.102806937966737030147096751318e0,
	0.0514718425553176958330252131667e0,
    0.0e0};
    static Mfloat wgk[] = {
	0.00138901369867700762455159122676e0,
	0.00389046112709988405126720184452e0,
	0.00663070391593129217331982636975e0,
	0.00927327965951776342844114689202e0,
	0.0118230152534963417422328988533e0,
	0.0143697295070458048124514324436e0,
	0.0169208891890532726275722894203e0,
	0.0194141411939423811734089510501e0,
	0.0218280358216091922971674857383e0,
	0.0241911620780806013656863707252e0,
	0.0265099548823331016106017093351e0,
	0.0287540487650412928439787853543e0,
	0.0309072575623877624728842529431e0,
	0.0329814470574837260318141910169e0,
	0.0349793380280600241374996707315e0,
	0.0368823646518212292239110656171e0,
	0.0386789456247275929503486515323e0,
	0.0403745389515359591119952797525e0,
	0.041969810215164246147147541286e0,
	0.0434525397013560693168317281171e0,
	0.0448148001331626631923555516167e0,
	0.0460592382710069881162717355594e0,
	0.0471855465692991539452614781811e0,
	0.0481858617570871291407794922983e0,
	0.0490554345550297788875281653672e0,
	0.0497956834270742063578115693799e0,
	0.0504059214027823468408930856536e0,
	0.0508817958987496064922974730498e0,
	0.0512215478492587721706562826049e0,
	0.0514261285374590259338628792158e0,
    0.0514947294294515675583404336471e0};
    static Mfloat wg[] = {
	0.00796819249616660561546588347467e0,
	0.0184664683110909591423021319121e0,
	0.0287847078833233693497191796113e0,
	0.0387991925696270495968019364464e0,
	0.0484026728305940529029381404228e0,
	0.0574931562176190664817216894021e0,
	0.065974229882180495128128515116e0,
	0.0737559747377052062682438500222e0,
	0.0807558952294202153546949384605e0,
	0.0868997872010829798023875307151e0,
	0.0921225222377861287176327070876e0,
	0.0963687371746442596394686263518e0,
	0.0995934205867952670627802821036e0,
	0.101762389748405504596428952169e0,
    0.102852652893558840341285636705e0};
    imsl_q4ng (&epmach, &uflow, &oflow);

    centr = F_HALF * (*b + *a);
    hlgth = F_HALF * (*b - *a);
    dhlgth = fabs (hlgth);
    /*
     * COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ABSOLUTE ERROR
     */
    resg = F_ZERO;
    imsl_e1usr ("ON");
    fc = (*f) (centr);
    imsl_e1usr ("OFF");
    resk = wgk[30] * fc;
    *resabs = fabs (resk);
    for (j = 1; j <= 15; j++) {
	jtw = j * 2;
	absc = hlgth * xgk[jtw - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (fabs (fval1) + fabs (fval2));
    }
    for (j = 1; j <= 15; j++) {
	jtwm1 = j * 2 - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = wgk[30] * fabs (fc - reskh);
    for (j = 1; j <= 30; j++) {
	*resasc += wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
		reskh));
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs ((resk - resg) * hlgth);
    if (*resasc != F_ZERO && *abserr != F_ZERO)
	*abserr = *resasc * imsl_f_min (F_ONE, pow (2.0e02 ** abserr / *resasc, 1.5e00));
    if (*resabs > uflow / (5.0e01 * epmach))
	*abserr = imsl_f_max ((epmach * 5.0e01) ** resabs, *abserr);
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:24:16 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:24:11
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q9AG/DQ9AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q9AG (F, A, B, RESULT, ABSERR, RESABS, RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q9AG
             INTEGRATION RULES
                STANDARD FORTRAN SUBROUTINE
                REAL VERSION

   2.        PURPOSE
                TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR
                               ESTIMATE
                           J = INTEGRAL OF ABS(F) OVER (A,B)

   3.        CALLING SEQUENCE
                CALL Q9AG(F,A,B,RESULT,ABSERR,RESABS,RESASC)

             PARAMETERS
              ON ENTRY
                F      - REAL
                         FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                         FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                         DECLARED E X T E R N A L IN THE CALLING PROGRAM.

                A      - REAL
                         LOWER LIMIT OF INTEGRATION

                B      - REAL
                         UPPER LIMIT OF INTEGRATION

              ON RETURN
                RESULT - REAL
                         APPROXIMATION TO THE INTEGRAL I
                         RESULT IS COMPUTED BY APPLYING THE 21-POINT
                         KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
                         OF ABSCISSAE TO THE 10-POINT GAUSS RULE (RESG).

                ABSERR - REAL
                         ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                         WHICH SHOULD NOT EXCEED ABS(I-RESULT)

                RESABS - REAL
                         APPROXIMATION TO THE INTEGRAL J

                RESASC - REAL
                         APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
                         OVER (A,B)

   4.        SUBROUTINES OR FUNCTIONS NEEDED
                   - F (USER-PROVIDED FUNCTION)
                   - Q4NG
                   - FORTRAN ABS, AMAX1, AMIN1

       ..................................................................



             THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
             BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
             CORRESPONDING WEIGHTS ARE GIVEN.

             XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
                      XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
                      GAUSS RULE
                      XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
                      ADDED TO THE 10-POINT GAUSS RULE

             WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE

             WG     - WEIGHTS OF THE 10-POINT GAUSS RULE



             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC   - ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
             RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
             RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
                      I.E. TO I/(B-A)


             MACHINE DEPENDENT CONSTANTS
             ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_q9ag (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc)
#else
void        imsl_q9ag (f, a, b, result, abserr, resabs, resasc)
    Mfloat      (*f) (), *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    Mint        j, jtw, jtwm1;
    Mfloat      _f0, absc, centr, dhlgth, epmach, fc, fsum, fv1[10], fv2[10],
                fval1, fval2, hlgth, oflow, resg, resk, reskh, uflow;
    static Mfloat xgk[] = {
	0.995657163025808080735527280689e0,
	0.973906528517171720077964012085e0,
	0.93015749135570822600120718006e0,
	0.865063366688984510732096688424e0,
	0.780817726586416897063717578345e0,
	0.679409568299024406234327365115e0,
	0.562757134668604683339000099273e0,
	0.433395394129247190799265943166e0,
	0.294392862701460198131126603104e0,
	0.14887433898163121088482600113e0,
    0.0e0};
    static Mfloat wgk[] = {
	0.0116946388673718742780643960622e0,
	0.0325581623079647274788189724594e0,
	0.0547558965743519960313813002446e0,
	0.0750396748109199527670431409162e0,
	0.0931254545836976055350654650834e0,
	0.109387158802297641899210590326e0,
	0.123491976262065851077958109831e0,
	0.134709217311473325928054001772e0,
	0.142775938577060080797094273139e0,
	0.147739104901338491374841515972e0,
    0.14944555400291690566493646839e0};
    static Mfloat wg[] = {
	0.0666713443086881375935688098933e0,
	0.149451349150580593145776339658e0,
	0.219086362515982043995534934228e0,
	0.26926671930999635509122692157e0,
    0.295524224714752870173892994651e0};
    imsl_q4ng (&epmach, &uflow, &oflow);

    centr = F_HALF * (*a + *b);
    hlgth = F_HALF * (*b - *a);
    dhlgth = fabs (hlgth);
    /*
     * COMPUTE THE 21-POINT KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ABSOLUTE ERROR
     */
    resg = F_ZERO;
    imsl_e1usr ("ON");
    fc = (*f) (centr);
    imsl_e1usr ("OFF");
    resk = wgk[10] * fc;
    *resabs = fabs (resk);
    for (j = 1; j <= 5; j++) {
	jtw = 2 * j;
	absc = hlgth * xgk[jtw - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (fabs (fval1) + fabs (fval2));
    }
    for (j = 1; j <= 5; j++) {
	jtwm1 = 2 * j - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	imsl_e1usr ("ON");
/*		fval1 = (*f) (ADR(_f0, centr - absc));
		fval2 = (*f) (ADR(_f0, centr + absc)); */
	_f0 = centr - absc;
	fval1 = (*f) (_f0);
	_f0 = centr + absc;
	fval2 = (*f) (_f0);
	imsl_e1usr ("OFF");
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = wgk[10] * fabs (fc - reskh);
    for (j = 1; j <= 10; j++) {
	*resasc += wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
		reskh));
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs ((resk - resg) * hlgth);
    if (*resasc != F_ZERO && *abserr != F_ZERO)
	*abserr = *resasc * imsl_f_min (F_ONE, pow (2.0e02 ** abserr / *resasc, 1.5e00));
    if (*resabs > uflow / (5.0e01 * epmach))
	*abserr = imsl_f_max ((epmach * 5.0e01) ** resabs, *abserr);
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:25:34 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:25:31
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q10G/DQ10G (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q10G (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD,
                           NRMAX)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q10G
             ORDERING ROUTINE
                STANDARD FORTRAN SUBROUTINE
                REAL VERSION

   2.        PURPOSE
                THIS ROUTINE MAINTAINS THE DESCENDING ORDERING
                IN THE LIST OF THE LOCAL ERROR ESTIMATES RESULTING FROM
                THE INTERVAL SUBDIVISION PROCESS. AT EACH CALL TWO ERROR
                ESTIMATES ARE INSERTED USING THE SEQUENTIAL SEARCH
                METHOD, TOP-DOWN FOR THE LARGEST ERROR ESTIMATE
                AND BOTTOM-UP FOR THE SMALLEST ERROR ESTIMATE.

   3.        CALLING SEQUENCE
                CALL Q10G(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)

             PARAMETERS (MEANING AT OUTPUT)
                LIMIT  - INTEGER
                         MAXIMUM NUMBER OF ERROR ESTIMATES THE LIST
                         CAN CONTAIN

                LAST   - INTEGER
                         NUMBER OF ERROR ESTIMATES CURRENTLY IN THE LIST

                MAXERR - INTEGER
                         MAXERR POINTS TO THE NRMAX-TH LARGEST ERROR
                         ESTIMATE CURRENTLY IN THE LIST

                ERMAX  - REAL
                         NRMAX-TH LARGEST ERROR ESTIMATE
                         ERMAX = ELIST(MAXERR)

                ELIST  - REAL
                         VECTOR OF DIMENSION LAST CONTAINING THE ERROR
                         ESTIMATES

                IORD   - INTEGER
                         VECTOR OF DIMENSION LAST, THE FIRST K ELEMENTS
                         OF WHICH CONTAIN POINTERS TO THE ERROR
                         ESTIMATES, SUCH THAT
                         ELIST(IORD(1)),...,  ELIST(IORD(K))
                         FORM A DECREASING SEQUENCE, WITH
                         K = LAST IF LAST.LE.(LIMIT/2+2), AND
                         K = LIMIT+1-LAST OTHERWISE

                NRMAX  - INTEGER
                         MAXERR = IORD(NRMAX)

   4.        NO SUBROUTINES OR FUNCTIONS NEEDED

       ..................................................................


             CHECK WHETHER THE LIST CONTAINS MORE THAN
             TWO ERROR ESTIMATES.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_q10g (Mint *limit, Mint *last, Mint *maxerr, Mfloat *ermax,
                Mfloat elist[], Mint iord[], Mint *nrmax)
#else
void 	    imsl_q10g (limit, last, maxerr, ermax, elist, iord, nrmax)
    Mint       *limit, *last, *maxerr;
    Mfloat     *ermax, elist[];
    Mint        iord[], *nrmax;
#endif
{
    Mint        i, ibeg, ido, isucc, j, jbnd, jupbn, k;
    Mfloat      errmax, errmin;


    if (*last > 2)
	goto L_10;
    /*
     * THIS ROUTINE IS A NUCLEI OF QDAG, BUT IS ALSO USED BY Q3AGI, Q3AGP,
     * Q3AGS, Q3AWC, Q3AWO, AND Q3AWS
     */
    iord[0] = 1;
    iord[1] = 2;
    goto L_90;
    /*
     * THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A DIFFICULT
     * INTEGRAND, SUBDIVISION INCREASED THE ERROR ESTIMATE.  IN THE NORMAL
     * CASE THE INSERT PROCEDURE SHOULD START AFTER THE NRMAX-TH LARGEST
     * ERROR ESTIMATE
     */
L_10:
    errmax = elist[*maxerr - 1];
    if (*nrmax == 1)
	goto L_30;
    ido = *nrmax - 1;
    for (i = 1; i <= ido; i++) {
	isucc = iord[*nrmax - 2];
	/* JUMP OUT OF DO-LOOP */
	if (errmax <= elist[isucc - 1])
	    goto L_30;
	iord[*nrmax - 1] = isucc;
	*nrmax -= 1;
    }
    /*
     * COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED IN
     * DESCENDING ORDER.  THIS NUMBER DEPENDS ON THE NUMBER OF SUBDIVISIONS
     * STILL ALLOWED
     */
L_30:
    jupbn = *last;
    if (*last > (*limit / 2 + 2))
	jupbn = *limit + 3 - *last;
    errmin = elist[*last - 1];
    /*
     * INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN, STARTING COMPARISON
     * FROM THE ELEMENT ELIST(IORD(NRMAX+1))
     */
    jbnd = jupbn - 1;
    ibeg = *nrmax + 1;
    if (ibeg > jbnd)
	goto L_50;
    for (i = ibeg; i <= jbnd; i++) {
	isucc = iord[i - 1];
	/* JUMP OUT OF DO-LOOP */
	if (errmax >= elist[isucc - 1])
	    goto L_60;
	iord[i - 2] = isucc;
    }
L_50:
    iord[jbnd - 1] = *maxerr;
    iord[jupbn - 1] = *last;
    goto L_90;
    /*
     * INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP
     */
L_60:
    iord[i - 2] = *maxerr;
    k = jbnd;
    for (j = i; j <= jbnd; j++) {
	isucc = iord[k - 1];
	/* JUMP OUT OF DO-LOOP */
	if (errmin < elist[isucc - 1])
	    goto L_80;
	iord[k] = isucc;
	k -= 1;
    }
    iord[i - 1] = *last;
    goto L_90;
L_80:
    iord[k] = *last;
    /* SET MAXERR AND ERMAX. */
L_90:
    *maxerr = iord[*nrmax - 1];
    *ermax = elist[*maxerr - 1];
    return;
}				/* end of function */











/*Translated by FOR_C++, v0.1, on 05/01/90 at 14:26:45 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 05/01/90 at 14:26:43
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q4NG/DQ4NG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function using a non-adaptive rule.

    Usage:      CALL Q4NG (EPMACH, UFLOW, OFLOW)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  .......................................................................

    REAL MACHINE CONSTANTS

    EPMACH = THE LARGEST RELATIVE SPACING

    UFLOW  = THE SMALLEST POSITIVE MAGNITUDE

    OFLOW  = THE LARGEST MAGNITUDE

  .......................................................................

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_q4ng (Mfloat *epmach, Mfloat *uflow, Mfloat *oflow)
#else
void        imsl_q4ng (epmach, uflow, oflow)
    Mfloat     *epmach, *uflow, *oflow;
#endif
{
    *epmach = imsl_amach (4);
    *uflow = imsl_amach (1);
    *oflow = imsl_amach (2);
    return;
}				/* end of function */
