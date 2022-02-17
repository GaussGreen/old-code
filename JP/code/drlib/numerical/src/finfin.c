#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn_inf (Mfloat (*fcn) (Mfloat), Mfloat bound,
                Imsl_quad interval, va_list argptr);
static void l_q2agi (Mfloat (*f) (Mfloat), Mfloat *bound, Mint *interv,
                Mfloat *errabs, Mfloat *errrel, Mfloat *result,
                Mfloat *errest, Mint *maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[]);
static void l_q4agi (Mfloat (*f) (Mfloat), Mfloat *boun, Mint *inf,
                Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc);
#else
static VA_LIST_HACK l_int_fcn_inf ();
static void l_q2agi ();
static void l_q4agi ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn_inf (Mfloat (*fcn) (Mfloat), Mfloat bound, Imsl_quad interval,
    ...)
#else
Mfloat      imsl_f_int_fcn_inf (fcn, bound, interval, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      bound;
    Imsl_quad        interval;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, interval);
    E1PSH ("imsl_f_int_fcn_inf", "imsl_d_int_fcn_inf");
    lv_value = F_ZERO;
    IMSL_CALL (l_int_fcn_inf (fcn, bound, interval, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn_inf", "imsl_d_int_fcn_inf");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn_inf (Mfloat (*fcn) (Mfloat), Mfloat bound,
                Imsl_quad interval, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn_inf (fcn, bound, interval, argptr)
    Mfloat      (*fcn) ();
    Mfloat      bound;
    Imsl_quad        interval;
    va_list     argptr;
#endif
{
    Mfloat      bound_float;
    Mint        code;
    Mint        arg_number = 3;
    Mint 	l_interval = 0;
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

    bound_float = bound;
    if (interval == IMSL_INF_BOUND)
	l_interval = -1;
    if (interval == IMSL_BOUND_INF)
	l_interval = 1;
    if (interval == IMSL_INF_INF)
	l_interval = 2;
    l_q2agi (fcn, &bound_float, &l_interval, &err_abs, &err_rel, &lv_value, err_est,
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
    if (imsl_n1rty(0) > 3)
	lv_value = imsl_amach(6);
    return (argptr);
}












/*Translated by FOR_C++, v0.1, on 08/07/90 at 09:24:18 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/07/90 at 09:24:15
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q2AGI/DQ2AGI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 29, 1988

    Purpose:    Integrate a function over an infinite or semi-infinite
                interval.

    Usage:      CALL Q2AGI (F, BOUND, INTERV, ERRABS, ERRREL, RESULT,
                            ERREST, MAXSUB, NEVAL, NSUBIN, ALIST, BLIST,
                            RLIST, ELIST, IORD)

    Arguments:  (See QDAGI)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q2agi (Mfloat (*f) (Mfloat), Mfloat *bound, Mint *interv,
                Mfloat *errabs, Mfloat *errrel, Mfloat *result,
                Mfloat *errest, Mint *maxsub, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[])
#else
static void l_q2agi (f, bound, interv, errabs, errrel, result,
                errest, maxsub, neval, nsubin, alist, blist, rlist, elist, iord)
    Mfloat      (*f) (), *bound;
    Mint        *interv;
    Mfloat     *errabs, *errrel, *result, *errest;
    Mint        *maxsub, *neval, *nsubin;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint         iord[];
#endif
{
    Mint         ier;
    Mfloat      aa, bb, cc, dd, tt, uu;


    imsl_e1psh ("Q2AGI ");
    /*
     * CHECK INTERV
     */
    if ((*interv != -1 && *interv != 1) && *interv != 2) {
	imsl_e1sti (1, *interv);

	/*
	 * (5, 2, "The chosen interval, interval, is equal to '%(i1)'.  Valid
	 * choices for the variable interval include IMSL_INF_BOUND,
	 * IMSL_BOUND_INF, and IMSL_INF_INF.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_INTERVAL_BOUNDS);
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

    imsl_q3agi (f, bound, interv, errabs, errrel, maxsub, result, errest,
	neval, &ier, alist, blist, rlist, elist, iord, nsubin);

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
	/*
	 * Report original interval, not transformed interval.
	 */
	tt = alist[iord[0] - 1];
	uu = blist[iord[0] - 1];
	if (*interv == 1) {
	    cc = *bound + (F_ONE - tt) / tt;
	    dd = *bound + (F_ONE - uu) / uu;
	    aa = imsl_f_min (cc, dd);
	    bb = imsl_f_max (cc, dd);
	}
	else if (*interv == -1) {
	    cc = *bound + (tt - F_ONE) / tt;
	    dd = *bound + (uu - F_ONE) / uu;
	    aa = imsl_f_min (cc, dd);
	    bb = imsl_f_max (cc, dd);
	}
	else {
	    cc = fabs ((F_ONE - tt) / tt);
	    dd = fabs ((F_ONE - uu) / uu);
	    aa = imsl_f_min (cc, dd);
	    bb = imsl_f_max (cc, dd);
	    cc = -bb;
	    dd = -aa;
	}
	imsl_e1str (1, aa);
	imsl_e1str (2, bb);
	if (abs (*interv) == 1) {
	    imsl_ermes (IMSL_WARNING, IMSL_PRECISION_DEGRADATION);
	}
	else {
	    imsl_e1str (3, cc);
	    imsl_e1str (4, dd);
	    imsl_ermes (IMSL_WARNING, IMSL_PRECISION_DEGRADATION);
	}
    }
    else if (ier == 4) {
	imsl_e1str (1, *errabs);
	imsl_e1str (2, *errrel);

	/*
	 * (3, 4, "Roundoff error has been detected in the extrapolation
	 * table.  The requested tolerances, err_abs = %(r1) and err_rel =
	 * %(r2) cannot be reached.");
	 */
	imsl_ermes (IMSL_WARNING, IMSL_EXTRAPOLATION_ROUNDOFF);
    }
    else if (ier == 5) {

	/*
	 * (4, 5, "The integral is probably divergent or slowly
	 * convergent.");
	 */
	imsl_ermes (IMSL_FATAL, IMSL_DIVERGENT);
    }
L_9000:
    ;
    imsl_e1pop ("Q2AGI ");
    return;
}				/* end of function */


















/*Translated by FOR_C++, v0.1, on 08/07/90 at 09:25:20 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/07/90 at 09:25:16
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3AGI/DQ3AGI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function over a semi-infinite interval.

    Usage:      CALL Q3AGI (F, BOUND, INF, EPSABS, EPSREL, LIMIT, RESULT,
                            ABSERR, NEVAL, IER, ALIST, BLIST, RLIST,
                            ELIST, IORD, LAST)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


   ......................................................................

   1.     Q3AGI
          INTEGRATION OVER INFINITE INTERVALS
             STANDARD FORTRAN SUBROUTINE

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
             INTEGRAL    I = INTEGRAL OF  F  OVER (BOUND,+INFINITY)
                      OR I = INTEGRAL OF  F  OVER (-INFINITY,BOUND)
                      OR I = INTEGRAL OF  F  OVER (-INFINITY,+INFINITY),
             HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

   3.     CALLING SEQUENCE
             CALL Q3AGI(F,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,
                        NEVAL,IER)

          PARAMETERS
           ON ENTRY
              F      - REAL
                       FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                       FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                       DECLARED E X T E R N A L IN THE DRIVER PROGRAM.

              BOUND  - REAL
                       FINITE BOUND OF INTEGRATION RANGE
                       (HAS NO MEANING IF INTERVAL IS DOUBLY-INFINITE)

              INF    - REAL
                       INDICATING THE KIND OF INTEGRATION RANGE INVOLVED
                       INF = 1 CORRESPONDS TO  (BOUND,+INFINITY),
                       INF = -1            TO  (-INFINITY,BOUND),
                       INF = 2             TO (-INFINITY,+INFINITY).

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
                       IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE. THE
                               ESTIMATES FOR RESULT AND ERROR ARE LESS
                               RELIABLE. IT IS ASSUMED THAT THE REQUESTED
                               ACCURACY HAS NOT BEEN ACHIEVED.
                       IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
                               HAS BEEN ACHIEVED. ONE CAN ALLOW MORE
                               SUBDIVISIONS BY INCREASING THE DATA VALUE
                               OF LIMIT IN Q3AGI(AND TAKING THE ACCORDING
                               DIMENSION ADJUSTMENTS INTO ACCOUNT).
                               HOWEVER, IF THIS YIELDS NO IMPROVEMENT
                               IT IS ADVISED TO ANALYZE THE INTEGRAND
                               IN ORDER TO DETERMINE THE INTEGRATION
                               DIFFICULTIES. IF THE POSITION OF A LOCAL
                               DIFFICULTY CAN BE DETERMINED (E.G.
                               SINGULARITY, DISCONTINUITY WITHIN THE
                               INTERVAL) ONE WILL PROBABLY GAIN FROM
                               SPLITTING UP THE INTERVAL AT THIS POINT
                               AND CALLING THE INTEGRATOR ON THE
                               SUBRANGES. IF POSSIBLE, AN APPROPRIATE
                               SPECIAL-PURPOSE INTEGRATOR SHOULD BE USED,
                               WHICH IS DESIGNED FOR HANDLING THE TYPE
                               OF DIFFICULTY INVOLVED.
                           = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS
                               DETECTED, WHICH PREVENTS THE REQUESTED
                               TOLERANCE FROM BEING ACHIEVED.
                               THE ERROR MAY BE UNDER-ESTIMATED.
                           = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS
                               AT SOME POINTS OF THE INTEGRATION
                               INTERVAL.
                           = 4 THE ALGORITHM DOES NOT CONVERGE. ROUNDOFF
                               ERROR IS DETECTED IN THE EXTRAPOLATION
                               TABLE. IT IS ASSUMED THAT THE REQUESTED
                               TOLERANCE CANNOT BE ACHIEVED, AND THAT THE
                               RETURNED RESULT IS THE BEST WHICH CAN BE
                               OBTAINED.
                           = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR
                               SLOWLY CONVERGENT. IT MUST BE NOTED THAT
                               DIVERGENCE CAN OCCUR WITH ANY OTHER VALUE
                               OF IER.
                           = 6 THE INPUT IS INVALID, BECAUSE
                               INF.NE.1.AND.INF.NE.-1.AND.INF.NE.2,
                               EPSABS.LT.0 AND EPSREL.LT.0,
                               RESULT, ABSERR, NEVAL ARE SET TO ZERO.

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q4AGI
                - Q10G
                - Q4AWO
                - F (USER-PROVIDED FUNCTION)
                - Q4NG
                - FORTRAN ABS, AMAX1

   ......................................................................



              THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
              LIMEXP IN SUBROUTINE Q4AWO.


             LIMIT IS THE MAXIMUM NUMBER OF SUBINTERVALS ALLOWED IN THE
             SUBDIVISION PROCESS OF Q3AGI. TAKE CARE THAT LIMIT.GE.1.

              LIST OF MAJOR VARIABLES
              -----------------------

             ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
                         CONSIDERED UP TO NOW
             BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
                         CONSIDERED UP TO NOW
             RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
                         (ALIST(I),BLIST(I))
             RLIST2    - ARRAY OF DIMENSION AT LEAST (LIMEXP+2),
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
             *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
             *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
             LAST      - INDEX FOR SUBDIVISION
             NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
             NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
                         APPROPRIATE APPROXIMATION TO THE COMPOUNDED
                         INTEGRAL HAS BEEN OBTAINED, IT IS PUT IN
                         RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
                         BY ONE.
             SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
                         TO NOW, MULTIPLIED BY 1.5
             ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
                         THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
             EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
                         IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
                         BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
                         TRY TO DECREASE THE VALUE OF ERLARG.
             NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
                         IS NO LONGER ALLOWED (TRUE-VALUE)

              MACHINE DEPENDENT CONSTANTS
              ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST POSITIVE MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_q3agi (Mfloat (*f) (Mfloat), Mfloat *bound, Mint *inf,
                Mfloat *epsabs, Mfloat *epsrel, Mint *limit,
                Mfloat *result, Mfloat *abserr, Mint *neval,
                Mint *ier, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint *last)
#else
void        imsl_q3agi (f, bound, inf, epsabs, epsrel, limit, result,
                abserr, neval, ier, alist, blist, rlist, elist, iord, last)
    Mfloat      (*f) (), *bound;
    Mint        *inf;
    Mfloat     *epsabs, *epsrel;
    Mint        *limit;
    Mfloat     *result, *abserr;
    Mint        *neval, *ier;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint         iord[], *last;
#endif
{
    Mlong        extrap, noext;
    Mint         id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin,
                maxerr, nres, nrmax, numrl2;
    Mint         IMSLTRUE = 1, IMSLFALSE = 0;
    Mfloat      _f0, _f1, a1, a2, abseps, area, area1, area12, area2, b1,
                b2, boun, correc, defab1, defab2, defabs, dres, epmach,
                erlarg, erlast, errbnd, errmax, erro12, error1, error2,
                errsum, ertest, oflow, res3la[3], resabs, reseps, rlist2[52],
                small, uflow;


    imsl_q4ng (&epmach, &uflow, &oflow);
    /* TEST ON VALIDITY OF PARAMETERS */
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = F_ZERO;
    *abserr = F_ZERO;
    alist[0] = F_ZERO;
    blist[0] = F_ONE;
    rlist[0] = F_ZERO;
    elist[0] = F_ZERO;
    iord[0] = 0;
    if (*epsabs < F_ZERO && *epsrel < F_ZERO)
	*ier = 6;
    if ((*inf != 1 && *inf != (-1)) && *inf != 2)
	*ier = 6;
    if (*ier == 6)
	goto L_170;
    /*
     * FIRST APPROXIMATION TO THE INTEGRAL - DETERMINE THE INTERVAL TO BE
     * MAPPED ONTO (0,1). IF INF = 2 THE INTEGRAL IS COMPUTED AS I = I1+I2,
     * WHERE I1 = INTEGRAL OF F OVER (-INFINITY,0), I2 = INTEGRAL OF F OVER
     * (0,+INFINITY)
     */
    boun = *bound;
    if (*inf == 2)
	boun = F_ZERO;
    _f0 = F_ZERO;
    _f1 = F_ONE;
    l_q4agi (f, &boun, inf, &_f0, &_f1, result,
	abserr, &defabs, &resabs);
    /* TEST ON ACCURACY */
    *last = 1;
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 1;
    dres = fabs (*result);
    errbnd = imsl_f_max (*epsabs, *epsrel * dres);
    if (*abserr <= 1.0e02 * epmach * defabs && *abserr > errbnd)
	*ier = 2;
    if (*limit == 1)
	*ier = 1;
    if ((*ier != 0 || (*abserr <= errbnd && *abserr != resabs)) ||
	*abserr == F_ZERO)
	goto L_160;
    /* INITIALIZATION */
    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    ktmin = 0;
    numrl2 = 2;
    extrap = IMSLFALSE;
    noext = IMSLFALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (F_ONE - 5.0e01 * epmach) * defabs)
	ksgn = 1;

    /* MAIN DO-LOOP */
    for (*last = 2; *last <= *limit; (*last)++) {
	/*
	 * BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE
	 */
	a1 = alist[maxerr - 1];
	b1 = F_HALF * (alist[maxerr - 1] + blist[maxerr - 1]);
	a2 = b1;
	b2 = blist[maxerr - 1];
	erlast = errmax;
	/*
	 * SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR AT SOME
	 * POINTS OF THE INTEGRATION RANGE
	 */
	if (*last != 2) {
	    if (imsl_f_max (fabs (a1), fabs (b2)) <= (F_ONE + 1.0e03 *
		    epmach) * (fabs (a2) + 1.0e03 * uflow)) {
		*ier = 4;
		goto L_140;
	    }
	}
	l_q4agi (f, &boun, inf, &a1, &b1, &area1, &error1, &resabs,
	    &defab1);
	l_q4agi (f, &boun, inf, &a2, &b2, &area2, &error2, &resabs,
	    &defab2);

	/*
	 * IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR AND TEST FOR
	 * ACCURACY
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
	 * TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG
	 */
	if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
	    *ier = 2;
	if (iroff2 >= 5)
	    ierro = 3;
	/*
	 * SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS EQUALS
	 * LIMIT
	 */
	if (*last == *limit)
	    *ier = 1;
	/*
	 * APPEND THE NEWLY-CREATED INTERVALS TO THE LIST
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
	 * LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT)
	 */
L_40:
	imsl_q10g (limit, last, &maxerr, &errmax, elist, iord, &nrmax);
	if (errsum <= errbnd)
	    goto L_140;
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
	 * INTERVAL
	 */
	if (fabs (blist[maxerr - 1] - alist[maxerr - 1]) > small)
	    goto L_100;
	extrap = IMSLTRUE;
	nrmax = 2;
L_50:
	if (ierro == 3 || erlarg <= ertest)
	    goto L_70;
	/*
	 * THE SMALLEST INTERVAL HAS THE LARGEST ERROR. BEFORE BISECTING
	 * DECREASE THE SUM OF THE ERRORS OVER THE LARGER INTERVALS (ERLARG)
	 * AND PERFORM EXTRAPOLATION
	 */
	id = nrmax;
	jupbnd = *last;
	if (*last > (2 + *limit / 2))
	    jupbnd = *limit + 3 - *last;
	for (k = id; k <= jupbnd; k++) {
	    maxerr = iord[nrmax - 1];
	    errmax = elist[maxerr - 1];
	    if (fabs (blist[maxerr - 1] - alist[maxerr - 1]) > small)
		goto L_100;
	    nrmax += 1;
	}
	/* PERFORM EXTRAPOLATION */
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
	if (*abserr <= ertest)
	    goto L_110;
	/*
	 * PREPARE BISECTION OF THE SMALLEST INTERVAL
	 */
L_80:
	if (numrl2 == 1)
	    noext = IMSLTRUE;
	if (*ier == 5)
	    goto L_110;
	maxerr = iord[0];
	errmax = elist[maxerr - 1];
	nrmax = 1;
	extrap = IMSLFALSE;
	small *= F_HALF;
	erlarg = errsum;
	goto L_100;
L_90:
	small = 3.75e-01;
	erlarg = errsum;
	ertest = errbnd;
	rlist2[1] = area;
L_100:
	;
    }
    /* SET FINAL RESULT AND ERROR ESTIMATE */
L_110:
    if (*abserr == oflow)
	goto L_140;
    if ((*ier + ierro) == 0)
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
     * TEST ON DIVERGENCE
     */
L_130:
    if (ksgn == (-1) && imsl_f_max (fabs (*result), fabs (area)) <= defabs *
	1.0e-02)
	goto L_160;
    if ((1.0e-02 > (*result / area) || (*result / area) > 1.0e02) || errsum >
	fabs (area))
	*ier = 6;
    goto L_160;
    /* COMPUTE GLOBAL INTEGRAL SUM */
L_140:
    *result = F_ZERO;
    if (*ier == 4)
	*last -= 1;
    for (k = 1; k <= *last; k++) {
	*result += rlist[k - 1];
    }
    *abserr = errsum;
L_160:
    *neval = 30 ** last - 15;
    if (*inf == 2)
	*neval *= 2;
    if (*ier > 2)
	*ier -= 1;
L_170:
    return;
}				/* end of function */

















/*Translated by FOR_C++, v0.1, on 08/07/90 at 09:26:15 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/07/90 at 09:26:12
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q4AGI/DQ4AGI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function over a (semi)-infinite interval.

    Usage:      CALL Q4AGI (F, BOUN, INF, A, B, RESULT, ABSERR, RESABS,
                            RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q4AGI
             INTEGRATION RULE
                STANDARD FORTRAN SUBROUTINE
                REAL VERSION

   2.        PURPOSE
                THE ORIGINAL (INFINITE) INTEGRATION RANGE IS MAPPED
                ONTO THE INTERVAL (0,1) AND (A,B) IS A PART OF (0,1).
                IT IS THE PURPOSE TO COMPUTE
                I = INTEGRAL OF TRANSFORMED INTEGRAND OVER (A,B),
                J = INTEGRAL OF ABS(TRANSFORMED INTEGRAND) OVER (A,B).

   3.        CALLING SEQUENCE
               CALL Q4AGI(F,BOUN,INF,A,B,RESULT,ABSERR,RESABS,RESASC)

             PARAMETERS
              ON ENTRY
                F      - REAL
                         FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                         FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                         DECLARED E X T E R N A L IN THE CALLING PROGRAM.

                BOUN   - REAL
                         FINITE BOUND OF ORIGINAL INTEGRATION RANGE
                         (SET TO ZERO IF INF = +2)

                INF    - INTEGER
                         IF INF = -1, THE ORIGINAL INTERVAL IS
                                     (-INFINITY,BOUND),
                         IF INF = +1, THE ORIGINAL INTERVAL IS
                                     (BOUND,+INFINITY),
                         IF INF = +2, THE ORIGINAL INTERVAL IS
                                     (-INFINITY,+INFINITY) AND
                         THE INTEGRAL IS COMPUTED AS THE SUM OF TWO
                         INTEGRALS, ONE OVER (-INFINITY,0)
                         AND ONE OVER (0,+INFINITY).

                A      - REAL
                         LOWER LIMIT FOR INTEGRATION OVER SUBRANGE
                         OF (0,1)

                B      - REAL
                         UPPER LIMIT FOR INTEGRATION OVER SUBRANGE
                         OF (0,1)

              ON RETURN
                RESULT - REAL
                         APPROXIMATION TO THE INTEGRAL I
                         RESULT IS COMPUTED BY APPLYING THE 15-POINT
                         KRONROD RULE(RESK) OBTAINED BY OPTIMAL ADDITION
                         OF ABSCISSAE TO THE 7-POINT GAUSS RULE (RESG).

                ABSERR - REAL
                         ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                         WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

                RESABS - REAL
                         APPROXIMATION TO THE INTEGRAL J

                RESASC - REAL
                         APPROXIMATION TO THE INTEGRAL OF
                         ABS((TRANSFORMED INTEGRAND)-I/(B-A)) OVER (A,B)

   4.        SUBROUTINES OR FUNCTIONS NEEDED
                   - F (USER-PROVIDED FUNCTION)
                   - Q4NG
                   - FORTRAN ABS, AMAX1, AMIN1, MIN0

       ..................................................................



             THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
             (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
             THEIR CORRESPONDING WEIGHTS ARE GIVEN.

             XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
                      XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
                      GAUSS RULE
                      XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
                      ADDED TO THE 7-POINT GAUSS RULE

             WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE

             WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING
                      TO THE ABSCISSAE XGK(2), XGK(4), ...
                      WG(1), WG(3), ... ARE SET TO ZERO.



             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC*  - ABSCISSA
             TABSC* - TRANSFORMED ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
             RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
             RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
                      INTEGRAND OVER (A,B), I.E. TO I/(B-A)

             MACHINE DEPENDENT CONSTANTS
             ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */

static Mfloat lv_xgk[] = {
    0.991455371120812639206854697526e0,
    0.949107912342758524526189684048e0,
    0.864864423359769072789712788641e0,
    0.741531185599394439863864773281e0,
    0.586087235467691130294144838259e0,
    0.405845151377397166906606412077e0,
    0.207784955007898467600689403773e0,
    0.0e0
};

static Mfloat lv_wgk[] = {
    0.022935322010529224963732008059e0,
    0.0630920926299785532907006631892e0,
    0.104790010322250183839876322542e0,
    0.14065325971552591874518959051e0,
    0.169004726639267902826583426599e0,
    0.190350578064785409913256402421e0,
    0.204432940075298892414161999235e0,
    0.209482141084727828012999174892e0
};

static Mfloat lv_wg[] = {
    0.0e0,
    0.129484966168869693270611432679e0,
    0.0e0,
    0.279705391489276667901467771424e0,
    0.0e0,
    0.381830050505118944950369775489e0,
    0.0e0,
    0.417959183673469387755102040816e0
};
#ifdef ANSI
static void l_q4agi (Mfloat (*f) (Mfloat), Mfloat *boun, Mint *inf,
                Mfloat *a, Mfloat *b, Mfloat *result,
                Mfloat *abserr, Mfloat *resabs, Mfloat *resasc)
#else
static void l_q4agi (f, boun, inf, a, b, result, abserr, resabs,
                resasc)
    Mfloat      (*f) (), *boun;
    int        *inf;
    Mfloat     *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    int         j;
    Mfloat      absc, absc1, absc2, centr, dinf, epmach, fc, fsum,
                fv1[7], fv2[7], fval1, fval2, hlgth, oflow, resg, resk,
                reskh, tabsc1, tabsc2, uflow;
    imsl_q4ng (&epmach, &uflow, &oflow);
    dinf = imsl_i_min (1, *inf);

    centr = F_HALF * (*a + *b);
    hlgth = F_HALF * (*b - *a);
    tabsc1 = *boun + dinf * (F_ONE - centr) / centr;
    imsl_e1usr ("ON");
    fval1 = (*f) (tabsc1);
    imsl_e1usr ("OFF");
    imsl_e1usr ("ON");
    if (*inf == 2)
	fval1 += (*f) (-tabsc1);
    imsl_e1usr ("OFF");
    fc = (fval1 / centr) / centr;
    /*
     * COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ERROR.
     */
    resg = lv_wg[7] * fc;
    resk = lv_wgk[7] * fc;
    *resabs = fabs (resk);
    for (j = 1; j <= 7; j++) {
	absc = hlgth * lv_xgk[j - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	tabsc1 = *boun + dinf * (F_ONE - absc1) / absc1;
	tabsc2 = *boun + dinf * (F_ONE - absc2) / absc2;
	imsl_e1usr ("ON");
	fval1 = (*f) (tabsc1);
	imsl_e1usr ("OFF");
	imsl_e1usr ("ON");
	fval2 = (*f) (tabsc2);
	imsl_e1usr ("OFF");
	imsl_e1usr ("ON");
	if (*inf == 2)
	    fval1 += (*f) (-tabsc1);
	imsl_e1usr ("OFF");
	imsl_e1usr ("ON");
	if (*inf == 2)
	    fval2 += (*f) (-tabsc2);
	imsl_e1usr ("OFF");
	fval1 = (fval1 / absc1) / absc1;
	fval2 = (fval2 / absc2) / absc2;
	fv1[j - 1] = fval1;
	fv2[j - 1] = fval2;
	fsum = fval1 + fval2;
	resg += lv_wg[j - 1] * fsum;
	resk += lv_wgk[j - 1] * fsum;
	*resabs += lv_wgk[j - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = lv_wgk[7] * fabs (fc - reskh);
    for (j = 1; j <= 7; j++) {
	*resasc += lv_wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
		reskh));
    }
    *result = resk * hlgth;
    *resasc *= hlgth;
    *resabs *= hlgth;
    *abserr = fabs ((resk - resg) * hlgth);
    if (*resasc != F_ZERO && *abserr != F_ZERO)
	*abserr = *resasc * imsl_f_min (F_ONE, pow (2.0e02 ** abserr / *resasc, 1.5e00));
    if (*resabs > uflow / (5.0e01 * epmach))
	*abserr = imsl_f_max ((epmach * 5.0e01) ** resabs, *abserr);
    return;
}				/* end of function */
