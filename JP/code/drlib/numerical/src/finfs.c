#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn_sing (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                va_list argptr);
static void l_q2ags (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mfloat *errabs, Mfloat *errrel, Mfloat *result,
                Mfloat *errest, Mint *maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[]);
static void l_q3ags (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mfloat *epsabs, Mfloat *epsrel, Mint *limit,
                Mfloat *result, Mfloat *abserr, Mint *neval,
                Mint *ier, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint *last);
#else
static VA_LIST_HACK l_int_fcn_sing ();
static void l_q2ags ();
static void l_q3ags ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn_sing (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,...)
#else
Mfloat      imsl_f_int_fcn_sing (fcn, a, b, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, b);

    E1PSH ("imsl_f_int_fcn_sing", "imsl_d_int_fcn_sing");
    lv_value = F_ZERO;
    IMSL_CALL (l_int_fcn_sing (fcn, a, b, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn_sing", "imsl_d_int_fcn_sing");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn_sing (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn_sing (fcn, a, b, argptr)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    va_list     argptr;
#endif
{
    Mfloat      a_float;
    Mfloat      b_float;
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
    Mint        user_err_list = 0;
    Mint        user_err_order = 0;
    Mint        i;
    Mfloat      temp_err_est;
    Mint        temp_n_subinter;
    Mint        temp_n_evals;


    err_abs = sqrt( (double) imsl_amach(4));
    err_rel = sqrt( (double) imsl_amach(4));


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
    
    l_q2ags (fcn, &a_float, &b_float, &err_abs, &err_rel, &lv_value, err_est,
	&max_subinter, n_evals, n_subinter, alist, blist, rlist,
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
    if (imsl_n1rty (0) > 3)
	lv_value = imsl_amach(6);
    return (argptr);
}







/*Translated by FOR_C++, v0.1, on 08/06/90 at 14:21:34 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef INCL_PROTOTYPES
#include "/home/usr2/imsl1/clib/newclib/include/imsl_int.h"
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/06/90 at 14:21:32
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q2AGS/DQ2AGS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function with endpoint singularities.

    Usage:      CALL Q2AGS (F, A, B, ERRABS, ERRREL, RESULT, ERREST,
                            MAXSUB, NEVAL, NSUBIN, ALIST, BLIST, RLIST,
                            ELIST, IORD)

    Arguments:  (See QDAGS)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q2ags (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mfloat *errabs, Mfloat *errrel, Mfloat *result,
                Mfloat *errest, Mint *maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[])
#else
static void l_q2ags (f, a, b, errabs, errrel, result, errest, maxsub,
                neval, nsubin, alist, blist, rlist, elist, iord)
    Mfloat      (*f) (), *a, *b, *errabs, *errrel, *result, *errest;
    Mint        *maxsub, *neval, *nsubin;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint         iord[];
#endif
{
    Mint         ier;


    imsl_e1psh ("Q2AGS ");
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

    l_q3ags (f, a, b, errabs, errrel, maxsub, result, errest, neval,
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
    else if (ier == 4) {
	imsl_e1str (1, *errabs);
	imsl_e1str (2, *errrel);

	/*
	 * (3, 4, "Roundoff error has been detected in the extrapolation
	 * table.  The requested tolerances, err_abs = %(r1) and err_REL =
	 * %(r2) cannot be reached.");
	 */
	imsl_ermes (IMSL_WARNING, IMSL_EXTRAPOLATION_ROUNDOFF);
    }
    else if (ier == 5) {

	/*
	 * (4, 5, "Integral is probably divergent or slowly convergent.");
	 */
	imsl_ermes (IMSL_FATAL, IMSL_DIVERGENT);
    }
L_9000:
    ;
    imsl_e1pop ("Q2AGS ");
    return;
}				/* end of function */












/*Translated by FOR_C++, v0.1, on 08/06/90 at 14:26:03 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef INCL_PROTOTYPES
#include "/home/usr2/imsl1/clib/newclib/include/imsl_int.h"
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/06/90 at 14:25:58
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3AGS/DQ3AGS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q3AGS (F, A, B, EPSABS, EPSREL, LIMIT, RESULT,
                            ABSERR, NEVAL, IER, ALIST, BLIST, RLIST,
                            ELIST, IORD, LAST)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


   ......................................................................

   1.     Q3AGS
          COMPUTATION OF A DEFINITE INTEGRAL
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
             DEFINITE INTEGRAL   I = INTEGRAL OF  F  OVER (A,B),
             HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

   3.     CALLING SEQUENCE
             CALL Q3AGS(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER)

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
                           = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
                               HAS BEEN ACHIEVED. ONE CAN ALLOW MORE SUB-
                               DIVISIONS BY INCREASING THE DATA VALUE OF
                               LIMIT IN Q3AGS (AND TAKING THE ACCORDING
                               DIMENSION ADJUSTMENTS INTO ACCOUNT).
                               HOWEVER, IF THIS YIELDS NO IMPROVEMENT
                               IT IS ADVISED TO ANALYZE THE INTEGRAND
                               IN ORDER TO DETERMINE THE INTEGRATION
                               DIFFICULTIES. IF THE POSITION OF A
                               LOCAL DIFFICULTY CAN BE DETERMINED (E.G.
                               SINGULARITY, DISCONTINUITY WITHIN THE
                               INTERVAL) ONE WILL PROBABLY GAIN FROM
                               SPLITTING UP THE INTERVAL AT THIS POINT
                               AND CALLING THE INTEGRATOR ON THE SUB-
                               RANGES. IF POSSIBLE, AN APPROPRIATE
                               SPECIAL-PURPOSE INTEGRATOR SHOULD BE USED,
                               WHICH IS DESIGNED FOR HANDLING THE TYPE
                               OF DIFFICULTY INVOLVED.
                           = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS DETEC-
                               TED, WHICH PREVENTS THE REQUESTED
                               TOLERANCE FROM BEING ACHIEVED.
                               THE ERROR MAY BE UNDER-ESTIMATED.
                           = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR
                               OCCURS AT SOME  POINTS OF THE INTEGRATION
                               INTERVAL.
                           = 4 THE ALGORITHM DOES NOT CONVERGE. ROUNDOFF
                               ERROR IS DETECTED IN THE EXTRAPOLATION
                               TABLE. IT IS PRESUMED THAT THE REQUESTED
                               TOLERANCE CANNOT BE ACHIEVED, AND THAT THE
                               RETURNED RESULT IS THE BEST WHICH CAN BE
                               OBTAINED.
                           = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR
                               SLOWLY CONVERGENT. IT MUST BE NOTED
                               THAT DIVERGENCE CAN OCCUR WITH ANY OTHER
                               VALUE OF IER.
                           = 6 THE INPUT IS INVALID, BECAUSE
                               EPSABS.LT.0 AND EPSREL.LT.0,
                               RESULT, ABSERR, NEVAL ARE SET TO ZERO.

   4.      SUBROUTINES OR FUNCTIONS NEEDED
                - Q9AG
                - Q10G
                - Q4AWO
                - F (USER-PROVIDED FUNCTION)
                - Q4NG
                - FORTRAN ABS, AMAX1, AMIN1

   ......................................................................



              THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
              LIMEXP IN SUBROUTINE Q4AWO (RLIST2 SHOULD BE OF
              DIMENSION (LIMEXP+2) AT LEAST).


              LIMIT IS THE MAXIMUM NUMBER OF SUBINTERVALS ALLOWED IN THE
              SUBDIVISION PROCESS OF Q3AGS. TAKE CARE THAT LIMIT.GE.1.

              LIST OF MAJOR VARIABLES
              -----------------------

             ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
                         CONSIDERED UP TO NOW
             BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
                         CONSIDERED UP TO NOW
             RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
                         (ALIST(I),BLIST(I))
             RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
                         CONTAINING THE PART OF THE EPSILON TABLE
                         WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS
             ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
             MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
                         ESTIMATE
             ERRMAX    - ELIST(MAXERR)
             ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
                         (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
             AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
             ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
             ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
                         ABS(RESULT))
             *****1    - VARIABLE FOR THE LEFT INTERVAL
             *****2    - VARIABLE FOR THE RIGHT INTERVAL
             LAST      - INDEX FOR SUBDIVISION
             NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
             NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
                         APPROPRIATE APPROXIMATION TO THE COMPOUNDED
                         INTEGRAL HAS BEEN OBTAINED IT IS PUT IN
                         RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
                         BY ONE.
             SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED
                         UP TO NOW, MULTIPLIED BY 1.5
             ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
                         THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
             EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
                         IS ATTEMPTING TO PERFORM EXTRAPOLATION
                         I.E. BEFORE SUBDIVIDING THE SMALLEST INTERVAL
                         WE TRY TO DECREASE THE VALUE OF ERLARG.
             NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
                         IS NO LONGER ALLOWED (TRUE VALUE)

              MACHINE DEPENDENT CONSTANTS
              ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST POSITIVE MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q3ags (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mfloat *epsabs, Mfloat *epsrel, Mint *limit,
                Mfloat *result, Mfloat *abserr, Mint *neval,
                Mint *ier, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint *last)
#else
static void l_q3ags (f, a, b, epsabs, epsrel, limit, result, abserr,
                neval, ier, alist, blist, rlist, elist, iord, last)
    Mfloat      (*f) (), *a, *b, *epsabs, *epsrel;
    Mint        *limit;
    Mfloat     *result, *abserr;
    Mint        *neval, *ier;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint         iord[], *last;
#endif
{
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
    Mlong        extrap, noext;
    Mint         id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin,
                maxerr, nres, nrmax, numrl2;
    Mfloat      a1, a2, abseps, area, area1, area12, area2, b1, b2, correc,
                defab1, defab2, defabs, dres, epmach, erlarg, erlast, errbnd,
                errmax, erro12, error1, error2, errsum, ertest, oflow,
                res3la[3], resabs, reseps, rlist2[52], small, uflow;


    imsl_q4ng (&epmach, &uflow, &oflow);
    /* TEST ON VALIDITY OF PARAMETERS */
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = F_ZERO;
    *abserr = F_ZERO;
    alist[0] = *a;
    blist[0] = *b;
    rlist[0] = F_ZERO;
    elist[0] = F_ZERO;
    if (*epsabs < F_ZERO && *epsrel < F_ZERO)
	*ier = 6;
    if (*ier == 6)
	goto L_180;
    /* FIRST APPROXIMATION TO THE INTEGRAL */
    ierro = 0;
    imsl_q9ag (f, a, b, result, abserr, &defabs, &resabs);

    /* TEST ON ACCURACY. */
    dres = fabs (*result);
    errbnd = imsl_f_max (*epsabs, *epsrel * dres);
    *last = 1;
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 1;
    if (*abserr <= 1.0e02 * epmach * defabs && *abserr > errbnd)
	*ier = 2;
    if (*limit == 1)
	*ier = 1;
    if ((*ier != 0 || (*abserr <= errbnd && *abserr != resabs)) ||
	*abserr == F_ZERO)
	goto L_170;
    /* INITIALIZATION */
    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    numrl2 = 2;
    ktmin = 0;
    extrap = FALSE;
    noext = FALSE;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (F_ONE - 5.0e01 * epmach) * defabs)
	ksgn = 1;

    /* MAIN DO-LOOP */
    for (*last = 2; *last <= *limit; (*last)++) {
	/*
	 * BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR ESTIMATE.
	 */
	a1 = alist[maxerr - 1];
	b1 = F_HALF * (alist[maxerr - 1] + blist[maxerr - 1]);
	a2 = b1;
	b2 = blist[maxerr - 1];
	erlast = errmax;
	imsl_q9ag (f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	imsl_q9ag (f, &a2, &b2, &area2, &error2, &resabs, &defab2);

	/*
	 * IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR AND TEST FOR
	 * ACCURACY.
	 */
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum += erro12 - errmax;
	area += area12 - rlist[maxerr - 1];
	if (defab1 == error1 || defab2 == error2)
	    goto L_20;
	if (fabs (rlist[maxerr - 1] - area12) > 1.0e-05 * fabs (area12) ||
	    erro12 < 9.9e-01 * errmax)
	    goto L_10;
	if (extrap)
	    iroff2 += 1;
	if (!extrap)
	    iroff1 += 1;
L_10:
	if (*last > 10 && erro12 > errmax)
	    iroff3 += 1;
L_20:
	rlist[maxerr - 1] = area1;
	rlist[*last - 1] = area2;
	errbnd = imsl_f_max (*epsabs, *epsrel * fabs (area));
	/*
	 * TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
	 */
	if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
	    *ier = 2;
	if (iroff2 >= 5)
	    ierro = 3;
	/*
	 * SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS EQUALS
	 * LIMIT.
	 */
	if (*last == *limit)
	    *ier = 1;
	/*
	 * SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR AT A POINT
	 * OF THE INTEGRATION RANGE.
	 */
	if (imsl_f_max (fabs (a1), fabs (b2)) <= (F_ONE + 1.0e03 * epmach) *
	    (fabs (a2) + 1.0e03 * uflow))
	    *ier = 4;
	/*
	 * APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
	 */
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
	if (errsum <= errbnd)
	    goto L_140;
	/* JUMP OUT OF DO-LOOP */
	if (*ier != 0)
	    goto L_110;
	if (*last == 2)
	    goto L_90;
	if (noext)
	    goto L_100;
	erlarg -= erlast;
	if (fabs (b1 - a1) > small)
	    erlarg += erro12;
	if (extrap)
	    goto L_50;
	/*
	 * TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE SMALLEST
	 * INTERVAL.
	 */
	if (fabs (blist[maxerr - 1] - alist[maxerr - 1]) > small)
	    goto L_100;
	extrap = TRUE;
	nrmax = 2;
L_50:
	if (ierro == 3 || erlarg <= ertest)
	    goto L_70;

	/*
	 * THE SMALLEST INTERVAL HAS THE LARGEST ERROR. BEFORE BISECTING
	 * DECREASE THE SUM OF THE ERRORS OVER THE LARGER INTERVALS (ERLARG)
	 * AND PERFORM EXTRAPOLATION.
	 */
	id = nrmax;
	jupbnd = *last;
	if (*last > (2 + *limit / 2))
	    jupbnd = *limit + 3 - *last;
	for (k = id; k <= jupbnd; k++) {
	    maxerr = iord[nrmax - 1];
	    errmax = elist[maxerr - 1];
	    /*
	     * JUMP OUT OF DO-LOOP
	     */
	    if (fabs (blist[maxerr - 1] - alist[maxerr - 1]) > small)
		goto L_100;
	    nrmax += 1;
	}
	/* PERFORM EXTRAPOLATION. */
L_70:
	numrl2 += 1;
	rlist2[numrl2 - 1] = area;
	imsl_q4awo (&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	ktmin += 1;
	if (ktmin > 5 && *abserr < 1.0e-03 * errsum)
	    *ier = 5;
	if (abseps >= *abserr)
	    goto L_80;
	ktmin = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
	ertest = imsl_f_max (*epsabs, *epsrel * fabs (reseps));
	/* JUMP OUT OF DO-LOOP */
	if (*abserr <= ertest)
	    goto L_110;
	/*
	 * PREPARE BISECTION OF THE SMALLEST INTERVAL.
	 */
L_80:
	if (numrl2 == 1)
	    noext = TRUE;
	if (*ier == 5)
	    goto L_110;
	maxerr = iord[0];
	errmax = elist[maxerr - 1];
	nrmax = 1;
	extrap = FALSE;
	small *= F_HALF;
	erlarg = errsum;
	goto L_100;
L_90:
	small = fabs (*b - *a) * 3.75e-01;
	erlarg = errsum;
	ertest = errbnd;
	rlist2[1] = area;
L_100:
	;
    }
    /* SET FINAL RESULT AND ERROR ESTIMATE. */
L_110:
    if (*abserr == oflow)
	goto L_140;
    if (*ier + ierro == 0)
	goto L_130;
    if (ierro == 3)
	*abserr += correc;
    if (*ier == 0)
	*ier = 3;
    if (*result != F_ZERO && area != F_ZERO)
	goto L_120;
    if (*abserr > errsum)
	goto L_140;
    if (area == F_ZERO)
	goto L_160;
    goto L_130;
L_120:
    if (*abserr / fabs (*result) > errsum / fabs (area))
	goto L_140;

    /*
     * TEST ON DIVERGENCE.
     */
L_130:
    if (ksgn == (-1) && imsl_f_max (fabs (*result), fabs (area)) <= defabs *
	1.0e-02)
	goto L_160;
    if ((1.0e-02 > (*result / area) || (*result / area) > 1.0e02) || errsum >
	fabs (area))
	*ier = 6;
    goto L_160;
    /* COMPUTE GLOBAL INTEGRAL SUM. */
L_140:
    *result = F_ZERO;
    for (k = 1; k <= *last; k++) {
	*result += rlist[k - 1];
    }
    *abserr = errsum;
L_160:
    if (*ier > 2)
	*ier -= 1;
L_170:
    *neval = 42 ** last - 21;
L_180:
    return;
}				/* end of function */
