#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI

static VA_LIST_HACK l_int_fcn_sing_pts (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,

                Mint npoints, Mfloat *points, va_list argptr);

static void l_q2agp (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mint *npts,

                Mfloat *points, Mfloat *errabs, Mfloat *errrel,

                Mfloat *result, Mfloat *errest, Mint *maxsub,

                Mint *neval, Mint *nsubin, Mfloat *alist,

                Mfloat *blist, Mfloat *rlist, Mfloat *elist,

                Mint *iord, Mint *level, Mfloat *wk, Mint *iwk);

static void l_q3agp (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mint *npts2,

                Mfloat points[], Mfloat *epsabs, Mfloat *epsrel,

                Mint *limit, Mfloat *result, Mfloat *abserr,

                Mint *neval, Mint *ier, Mfloat alist[], Mfloat blist[],

                Mfloat rlist[], Mfloat elist[], Mint iord[],

                Mint level[], Mint *last, Mfloat pts[], Mint ndin[]);

#else

static VA_LIST_HACK l_int_fcn_sing_pts ();

static void l_q2agp ();

static void l_q3agp ();

#endif



static Mfloat lv_value;

#ifdef ANSI

Mfloat      imsl_f_int_fcn_sing_pts (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,

                Mint npoints, Mfloat *points,...)

#else

Mfloat      imsl_f_int_fcn_sing_pts (fcn, a, b, npoints, points, va_alist)

    Mfloat      (*fcn) ();

    Mfloat      a;

    Mfloat      b;

    Mint        npoints;

    Mfloat     *points;

va_dcl

#endif

{

    va_list     argptr;

    VA_START (argptr, points);



    E1PSH ("imsl_f_int_fcn_sing_pts", "imsl_d_int_fcn_sing_pts");

    lv_value = 0.0;

    IMSL_CALL (l_int_fcn_sing_pts (fcn, a, b, npoints, points, argptr));

    va_end (argptr);

    E1POP ("imsl_f_int_fcn_sing_pts", "imsl_d_int_fcn_sing_pts");

    return lv_value;

}







#ifdef ANSI

static VA_LIST_HACK l_int_fcn_sing_pts (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,

                Mint npoints, Mfloat *points, va_list argptr)

#else

static VA_LIST_HACK l_int_fcn_sing_pts (fcn, a, b, npoints, points, argptr)

    Mfloat      (*fcn) ();

    Mfloat      a;

    Mfloat      b;

    Mint        npoints;

    Mfloat     *points;

    va_list     argptr;

#endif



{

    Mint        code;

    Mint        arg_number = 5;

    Mfloat      temp_a;

    Mfloat      temp_b;

    Mfloat      err_abs;

    Mfloat      err_rel;

    Mfloat     *err_est = NULL;

    Mfloat     *alist = NULL;

    Mfloat     *blist = NULL;

    Mfloat     *rlist = NULL;

    Mfloat     *elist = NULL;

    Mfloat     *work = NULL;

    Mint       *iwork = NULL;

    Mint       *level = NULL;

    Mint        max_subinter = 500;

    Mint       *n_subinter = NULL;

    Mint       *n_evals = NULL;

    Mint       *iord = NULL;

    Mint        user_err_est = 0;

    Mint        user_n_subinter = 0;

    Mint        user_n_evals = 0;

    Mint        user_err_list = 0;

    Mint        user_err_order = 0;

    Mint        i;

    err_abs = sqrt( (double) imsl_amach (4));

    err_rel = sqrt( (double) imsl_amach (4));



    code = 1;

    while (code > 0) {

	code = va_arg (argptr, int);

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

	    err_abs = (Mfloat) va_arg (argptr, double);

	    arg_number++;

	    break;

	case IMSL_ERR_REL:

	    err_rel = (Mfloat) va_arg (argptr, double);

	    arg_number++;

	    break;

	case IMSL_ERR_EST:

	    user_err_est = 1;

	    err_est = va_arg (argptr, Mfloat *);

	    arg_number++;

	    break;

	case IMSL_MAX_SUBINTER:

	    max_subinter = va_arg (argptr, Mint);

	    arg_number++;

	    break;

	case IMSL_N_SUBINTER:

	    user_n_subinter = 1;

	    n_subinter = va_arg (argptr, Mint *);

	    arg_number++;

	    break;

	case IMSL_N_EVALS:

	    user_n_evals = 1;

	    n_evals = va_arg (argptr, Mint *);

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



    if (max_subinter < 1) {

	imsl_e1sti (1, max_subinter);

	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_SUBINTER_SMALL);

    }



    if (npoints < 0) {

	imsl_e1sti (1, npoints);

	imsl_e1stl (1, "npoints");

	/*

	 * (5, 1, "The number of break points %(l1) = %(i1).  It must be at

	 * least zero.");

	 */

	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_BREAK_POINTS);

    }



    if (imsl_n1rty (0))

	goto RETURN;



    if (fcn == NULL) {

	imsl_e1stl (1, "fcn");

	/* (5, 1, "The required argument %(L1) is NULL."); */

	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);

    }



    if (points == NULL) {

	imsl_e1stl (1, "points");

	/* (5, 1, "The required argument %(L1) is NULL."); */

	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);

    }

    if (imsl_n1rty (0))

	goto RETURN;





    if (!user_n_subinter)

	n_subinter = (Mint *) imsl_malloc (sizeof (*n_subinter));



    if (!user_err_est)

	err_est = (Mfloat *) imsl_malloc (sizeof (*err_est));



    if (!user_n_evals)

	n_evals = (Mint *) imsl_malloc (sizeof (*n_evals));



    if (!user_err_list)

	elist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*elist));



    if (!user_err_order)

	iord = (Mint *) imsl_malloc (max_subinter * sizeof (*iord));



    alist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*alist));

    blist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*blist));

    rlist = (Mfloat *) imsl_malloc (max_subinter * sizeof (*rlist));

    level = (Mint *) imsl_malloc (max_subinter * sizeof (*level));

    if (npoints >= 0) {

	work = (Mfloat *) imsl_malloc ((npoints + 2) * sizeof (*work));

	iwork = (Mint *) imsl_malloc ((npoints + 2) * sizeof (*iwork));

    }



    if (n_subinter == NULL || err_est == NULL || n_evals == NULL || elist == NULL || iord == NULL || alist == NULL || blist == NULL || rlist == NULL) {

	imsl_e1sti (1, max_subinter);

	/* (5, 1, "Not enough memory.  max_subinter = %(i1)"); */

        imsl_e1stl (1, "max_subinter");

	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

	goto FREE_SPACE;

    }



/*    imsl_q2ag (fcn, a, b, err_abs, err_rel, rule, &lv_value, err_est, max_subinter,

	       n_evals, n_subinter, alist, blist, rlist, elist, iord);

*/



    temp_a = a;

    temp_b = b;

    l_q2agp (fcn, &temp_a, &temp_b, &npoints, points, &err_abs, &err_rel, &lv_value, err_est,

	&max_subinter, n_evals, n_subinter, alist, blist, rlist, elist,

	iord, level, work, iwork);

FREE_SPACE:

    if (!user_n_subinter && n_subinter != NULL)

	imsl_free (n_subinter);

    if (!user_err_est && err_est != NULL)

	imsl_free (err_est);

    if (!user_n_evals && n_evals != NULL)

	imsl_free (n_evals);

    if (!user_err_list && elist != NULL)

	imsl_free (elist);

    if (!user_err_order && iord != NULL)

	imsl_free (iord);

    if (level != NULL)

	imsl_free (level);

    if (work != NULL)

	imsl_free (work);

    if (iwork != NULL)

	imsl_free (iwork);

    if (alist != NULL)

	imsl_free (alist);

    if (blist != NULL)

	imsl_free (blist);

    if (rlist != NULL)

	imsl_free (rlist);

RETURN:

    if (imsl_n1rty (0) > 3)

	lv_value = imsl_amach(6);

    return (argptr);

}





















/* Structured by FOR_STRUCT, v0.2, on 07/17/90 at 10:27:11

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Q2AGP/DQ2AGP (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 29, 1985



    Purpose:    Integrate a function with singularity points given.



    Usage:      CALL Q2AGP (F, A, B, NPTS, POINTS, ERRABS, ERRREL,

                            RESULT, ERREST, MAXSUB, NEVAL, NSUBIN, ALIST,

                            BLIST, RLIST, ELIST, IORD, LEVEL, WK, IWK)



    Arguments:  (See QDAGP)



    Chapter:    MATH/LIBRARY Integration and Differentiation



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_q2agp (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mint *npts,

                Mfloat *points, Mfloat *errabs, Mfloat *errrel,

                Mfloat *result, Mfloat *errest, Mint *maxsub,

                Mint *neval, Mint *nsubin, Mfloat *alist,

                Mfloat *blist, Mfloat *rlist, Mfloat *elist,

                Mint *iord, Mint *level, Mfloat *wk, Mint *iwk)

#else

static void l_q2agp (f, a, b, npts, points, errabs, errrel, result,

                errest, maxsub, neval, nsubin, alist, blist, rlist, elist, iord,

                level, wk, iwk)

    Mfloat      (*f) (), *a, *b;

    Mint        *npts;

    Mfloat      points[], *errabs, *errrel, *result, *errest;

    Mint        *maxsub, *neval, *nsubin;

    Mfloat      alist[], blist[], rlist[], elist[];

    Mint         iord[], level[];

    Mfloat      wk[];

    Mint         iwk[];

#endif

{

    Mint         ier = 0, npts2 = 0;

/*

    Mint         ier, npts2;

*/

    imsl_e1psh ("l_q2agp");

    /* CHECK NPTS */

    if (*npts < 0) {

	imsl_e1sti (1, *npts);

	imsl_e1stl (1, "npts");

	/*

	 * (5, 1, "The number of break points %(l1) = %(i1).  It must be at

	 * least zero.");

	 */

	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_BREAK_POINTS);

    }

    /* CHECK ERRABS */

    if (*errabs < 0.0e0) {

	imsl_e1str (1, *errabs);

	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_ABS_SMALL);

    }

    /* CHECK ERRREL */

    if (*errrel < 0.0e0) {

	imsl_e1str (1, *errrel);

	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_SMALL);

    }

    /*

     * CHECK ERRABS AND ERRREL

     */

    if (*errabs == 0.0e0 && *errrel == 0.0e0) {

	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_TOL_ZERO);

    }

    /* CHECK ERRREL .GE. 1 */

    if (*errrel >= 1.0) {

	imsl_e1str (1, *errrel);

	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_BIG);

    }

    if (imsl_n1rty (0) != 0)

	goto L_9000;



    npts2 = *npts + 2;



    l_q3agp (f, a, b, &npts2, points, errabs, errrel, maxsub, result,

	errest, neval, &ier, alist, blist, rlist, elist, iord, level,

	nsubin, wk, iwk);



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

	 * table.  The requested tolerances, ERRABS = %(r1) and ERRREL =

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

    imsl_e1pop ("l_q2agp");

    return;

}				/* end of function */

















/*Translated by FOR_C++, v0.1, on 07/17/90 at 10:28:48 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 07/17/90 at 10:28:30

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Q3AGP/DQ3AGP (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    September 6, 1988



    Purpose:    Integrate a function with points of singularity given.



    Usage:      CALL Q3AGP (F, A, B, NPTS2, POINTS, EPSABS, EPSREL,

                            LIMIT, RESULT, ABSERR, NEVAL, IER, ALIST,

                            BLIST, RLIST, ELIST, IORD, LEVEL, LAST,

                            PTS, NDIN)



    Arguments:  (See comment block below)



    Chapter:    MATH/LIBRARY Integration and Differentiation



    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.





   ......................................................................



   1.     Q3AGP

          COMPUTATION OF A DEFINITE INTEGRAL

            STANDARD FORTRAN SUBROUTINE

            REAL VERSION



   2.     PURPOSE

             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN

             DEFINITE INTEGRAL   I = INTEGRAL OF  F  OVER (A,B),

             HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY

             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

             INTERIOR BREAK POINTS  OF THE INTEGRATION INTERVAL, WHERE

             LOCAL DIFFICULTIES OF THE INTEGRAND MAY OCCUR (E.G.

             SINGULARITIES, DISCONTINUITIES), ARE PROVIDED BY THE USER.



   3.     CALLING SEQUENCE

             CALL Q3AGP(F,A,B,NPTS2,POINTS,EPSABS,EPSREL,RESULT,

                        ABSERR,NEVAL,IER)



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



              NPTS2  - INTEGER

                       NUMBER EQUAL TO TWO MORE THAN THE NUMBER

                       OF USER-SUPPLIED BREAK POINTS WITHIN THE

                       INTEGRATION RANGE, NPTS2.GE.2.

                       IF NPTS2.LT.2, THE ROUTINE WILL END WITH IER = 6.



              POINTS - REAL

                       VECTOR OF DIMENSION NPTS2, THE FIRST (NPTS2-2)

                       ELEMENTS OF WHICH ARE THE USER PROVIDED INTERIOR

                       BREAK POINTS. IF THESE POINTS DO NOT

                       CONSTITUTE AN ASCENDING SEQUENCE THERE WILL BE AN

                       AUTOMATIC SORTING.



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

                       IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE.

                               THE ESTIMATES FOR INTEGRAL AND ERROR ARE

                               LESS RELIABLE. IT IS ASSUMED THAT THE

                               REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.

                       IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED

                               HAS BEEN ACHIEVED. ONE CAN ALLOW MORE

                               SUBDIVISIONS BY INCREASING THE DATA VALUE

                               OF LIMIT IN Q3AGP(AND TAKING THE ACCORDING

                               DIMENSION ADJUSTMENTS INTO ACCOUNT).

                               HOWEVER, IF THIS YIELDS NO IMPROVEMENT

                               IT IS ADVISED TO ANALYZE THE INTEGRAND

                               IN ORDER TO DETERMINE THE INTEGRATION

                               DIFFICULTIES. IF THE POSITION OF A LOCAL

                               DIFFICULTY CAN BE DETERMINED (I.E.

                               SINGULARITY, DISCONTINUITY WITHIN THE

                               INTERVAL), IT SHOULD BE SUPPLIED TO THE

                               ROUTINE AS AN ELEMENT OF THE VECTOR

                               POINTS. IF NECESSARY, AN APPROPRIATE

                               SPECIAL-PURPOSE INTEGRATOR MUST BE USED,

                               WHICH IS DESIGNED FOR HANDLING THE TYPE

                               OF DIFFICULTY INVOLVED.

                           = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS

                               DETECTED, WHICH PREVENTS THE REQUESTED

                               TOLERANCE FROM BEING ACHIEVED.

                               THE ERROR MAY BE UNDER-ESTIMATED.

                           = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS

                               AT SOME POINTS OF THE INTEGRATION

                               INTERVAL.

                           = 4 THE ALGORITHM DOES NOT CONVERGE.  ROUNDOFF

                               ERROR IS DETECTED IN THE EXTRAPOLATION

                               TABLE. IT IS PRESUMED THAT THE REQUESTED

                               TOLERANCE CANNOT BE ACHIEVED, AND THAT

                               THE RETURNED RESULT IS THE BEST WHICH

                               CAN BE OBTAINED.

                           = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR

                               SLOWLY CONVERGENT. IT MUST BE NOTED THAT

                               DIVERGENCE CAN OCCUR WITH ANY OTHER VALUE

                               OF IER.GT.0.

                           = 6 THE INPUT IS INVALID BECAUSE

                               NPTS2.LT.2 OR

                               BREAK POINTS ARE SPECIFIED OUTSIDE

                               THE INTEGRATION RANGE OR

                               EPSABS.LT.0 AND EPSREL.LT.0,

                               OR LIMIT.LT.NPTS2.

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

              SUBDIVISION PROCESS OF Q3AGP. TAKE CARE THAT LIMIT.GE.2.



              LIST OF MAJOR VARIABLES

              -----------------------



             ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS

                         CONSIDERED UP TO NOW

             BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS

                         CONSIDERED UP TO NOW

             RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER

                         (ALIST(I),BLIST(I))

             RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2

                         CONTAINING THE PART OF THE EPSILON TABLE WHICH

                         IS STILL NEEDED FOR FURTHER COMPUTATIONS

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

             NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN

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

             NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS

                         NO LONGER ALLOWED (TRUE-VALUE)



              MACHINE DEPENDENT CONSTANTS

              ---------------------------



             EPMACH IS THE LARGEST RELATIVE SPACING.

             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.

             OFLOW IS THE LARGEST POSITIVE MAGNITUDE.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_q3agp (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, Mint *npts2,

                Mfloat points[], Mfloat *epsabs, Mfloat *epsrel,

                Mint *limit, Mfloat *result, Mfloat *abserr,

                Mint *neval, Mint *ier, Mfloat alist[], Mfloat blist[],

                Mfloat rlist[], Mfloat elist[], Mint iord[],

                Mint level[], Mint *last, Mfloat pts[], Mint ndin[])

#else

static void l_q3agp (f, a, b, npts2, points, epsabs, epsrel, limit,

                result, abserr, neval, ier, alist, blist, rlist, elist, iord,

                level, last, pts, ndin)

    Mfloat      (*f) (), *a, *b;

    Mint        *npts2;

    Mfloat      points[], *epsabs, *epsrel;

    Mint        *limit;

    Mfloat     *result, *abserr;

    Mint        *neval, *ier;

    Mfloat      alist[], blist[], rlist[], elist[];

    Mint         iord[], level[], *last;

    Mfloat      pts[];

    Mint         ndin[];

#endif

{

/*	LOGICAL32       extrap, noext; */

    Mlong        extrap, noext;

    Mint 	i, id, ierro, ind1, ind2, ip1, iroff1, iroff2, iroff3, j,

                jlow, jupbnd, k=0, ksgn, ktmin, levcur, levmax, maxerr, nintp1,

                nintt, npts, nres, nrmax, numrl2;

    Mfloat      a1, a2, abseps, area, area1, area12, area2, b1, b2, correc,

                defab1, defab2, defabs, dres, epmach, erlarg, erlast, errbnd,

                errmax, erro12, error1, error2, errsum, ertest, oflow,

                res3la[3], resa, resabs, reseps, rlist2[52], siign, temp,

                uflow;

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



    imsl_q4ng (&epmach, &uflow, &oflow);

    /* TEST ON VALIDITY OF PARAMETERS */

    *ier = 0;

    *neval = 0;

    *last = 0;

    *result = 0.0e00;

    *abserr = 0.0e00;

    alist[0] = *a;

    blist[0] = *b;

    rlist[0] = 0.0e00;

    elist[0] = 0.0e00;

    iord[0] = 0;

    level[0] = 0;

    npts = *npts2 - 2;

    if ((*npts2 < 2 || *limit < npts) || (*epsabs < 0.0e00 && *epsrel <

	    0.0e00))

	*ier = 6;

    if (*ier == 6)

	goto L_260;

    /*

     * IF ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN ASCENDING

     * SEQUENCE.

     */

    siign = 1.0e00;

    if (*a > *b)

	siign = -1.0e00;

    pts[0] = imsl_f_min (*a, *b);

    if (npts == 0)

	goto L_20;

    for (i = 1; i <= npts; i++) {

	pts[i] = points[i - 1];

    }

L_20:

    pts[npts + 1] = imsl_f_max (*a, *b);

    nintt = npts + 1;

    a1 = pts[0];

    if (npts == 0)

	goto L_40;

    nintp1 = nintt + 1;

    for (i = 1; i <= nintt; i++) {

	ip1 = i + 1;

	for (j = ip1; j <= nintp1; j++) {

	    if (pts[i - 1] <= pts[j - 1])

		goto L_30;

	    temp = pts[i - 1];

	    pts[i - 1] = pts[j - 1];

	    pts[j - 1] = temp;

    L_30:

	    ;

	}

    }

    if (pts[0] != imsl_f_min (*a, *b) || pts[nintp1 - 1] != imsl_f_max (*a, *b))

	*ier = 6;

    if (*ier == 6)

	goto L_260;

    /*

     * COMPUTE FIRST INTEGRAL AND ERROR APPROXIMATIONS.

     */

L_40:

    resabs = 0.0e00;

    for (i = 1; i <= nintt; i++) {

	b1 = pts[i];

	imsl_q9ag (f, &a1, &b1, &area1, &error1, &defabs, &resa);

	*abserr += error1;

	*result += area1;

	ndin[i - 1] = 0;

	if (error1 == resa && error1 != 0.0e00)

	    ndin[i - 1] = 1;

	resabs += defabs;

	level[i - 1] = 0;

	elist[i - 1] = error1;

	alist[i - 1] = a1;

	blist[i - 1] = b1;

	rlist[i - 1] = area1;

	iord[i - 1] = i;

	a1 = b1;

    }

    errsum = 0.0e00;

    for (i = 1; i <= nintt; i++) {

	if (ndin[i - 1] == 1)

	    elist[i - 1] = *abserr;

	errsum += elist[i - 1];

    }

    /* TEST ON ACCURACY. */

    *last = nintt;

    *neval = 21 * nintt;

    dres = fabs (*result);

    errbnd = imsl_f_max (*epsabs, *epsrel * dres);

    if (*abserr <= 1.0e02 * epmach * resabs && *abserr > errbnd)

	*ier = 2;

    if (nintt == 1)

	goto L_90;

    for (i = 1; i <= npts; i++) {

	jlow = i + 1;

	ind1 = iord[i - 1];

	for (j = jlow; j <= nintt; j++) {

	    ind2 = iord[j - 1];

	    if (elist[ind1 - 1] > elist[ind2 - 1])

		goto L_70;

	    ind1 = ind2;

	    k = j;

    L_70:

	    ;

	}

	if (ind1 == iord[i - 1])

	    goto L_80;

	iord[k - 1] = iord[i - 1];

	iord[i - 1] = ind1;

L_80:

	;

    }

    if (*limit < *npts2)

	*ier = 1;

L_90:

    if (*ier != 0 || *abserr <= errbnd)

	goto L_260;

    /* INITIALIZATION */

    rlist2[0] = *result;

    maxerr = iord[0];

    errmax = elist[maxerr - 1];

    area = *result;

    nrmax = 1;

    nres = 0;

    numrl2 = 1;

    ktmin = 0;

    extrap = FALSE;

    noext = FALSE;

    erlarg = errsum;

    ertest = errbnd;

    levmax = 1;

    iroff1 = 0;

    iroff2 = 0;

    iroff3 = 0;

    ierro = 0;

    *abserr = oflow;

    ksgn = -1;

    if (dres >= (1.0e00 - 5.0e01 * epmach) * resabs)

	ksgn = 1;



    /* MAIN DO-LOOP */

    for (*last = *npts2; *last <= *limit; (*last)++) {

	/*

	 * BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR ESTIMATE.

	 */

	levcur = level[maxerr - 1] + 1;

	a1 = alist[maxerr - 1];

	b1 = 5.0e-01 * (alist[maxerr - 1] + blist[maxerr - 1]);

	a2 = b1;

	b2 = blist[maxerr - 1];

	erlast = errmax;

	imsl_q9ag (f, &a1, &b1, &area1, &error1, &resa, &defab1);

	imsl_q9ag (f, &a2, &b2, &area2, &error2, &resa, &defab2);



	/*

	 * IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR AND TEST FOR

	 * ACCURACY.

	 */

	*neval += 42;

	area12 = area1 + area2;

	erro12 = error1 + error2;

	errsum += erro12 - errmax;

	area += area12 - rlist[maxerr - 1];

	if (defab1 == error1 || defab2 == error2)

	    goto L_110;

	if (fabs (rlist[maxerr - 1] - area12) > 1.0e-05 * fabs (area12) ||

	    erro12 < 9.9e-01 * errmax)

	    goto L_100;

	if (extrap)

	    iroff2 += 1;

	if (!extrap)

	    iroff1 += 1;

L_100:

	if (*last > 10 && erro12 > errmax)

	    iroff3 += 1;

L_110:

	level[maxerr - 1] = levcur;

	level[*last - 1] = levcur;

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

	 * OF THE INTEGRATION RANGE

	 */

	if (imsl_f_max (fabs (a1), fabs (b2)) <= (1.0e00 + 1.0e03 * epmach) *

	    (fabs (a2) + 1.0e03 * uflow))

	    *ier = 4;

	/*

	 * APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.

	 */

	if (error2 > error1)

	    goto L_120;

	alist[*last - 1] = a2;

	blist[maxerr - 1] = b1;

	blist[*last - 1] = b2;

	elist[maxerr - 1] = error1;

	elist[*last - 1] = error2;

	goto L_130;

L_120:

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

L_130:

	imsl_q10g (limit, last, &maxerr, &errmax, elist, iord, &nrmax);



	/* JUMP OUT OF DO-LOOP */

	if (errsum <= errbnd)

	    goto L_230;

	/* JUMP OUT OF DO-LOOP */

	if (*ier != 0)

	    goto L_200;

	if (noext)

	    goto L_190;

	erlarg -= erlast;

	if (levcur + 1 <= levmax)

	    erlarg += erro12;

	if (extrap)

	    goto L_140;

	/*

	 * TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE SMALLEST

	 * INTERVAL.

	 */

	if (level[maxerr - 1] + 1 <= levmax)

	    goto L_190;

	extrap = TRUE;

	nrmax = 2;

L_140:

	if (ierro == 3 || erlarg <= ertest)

	    goto L_160;



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

	    /* JUMP OUT OF DO-LOOP */

	    if (level[maxerr - 1] + 1 <= levmax)

		goto L_190;

	    nrmax += 1;

	}

	/* PERFORM EXTRAPOLATION. */

L_160:

	numrl2 += 1;

	rlist2[numrl2 - 1] = area;

	if (numrl2 <= 2)

	    goto L_180;

	imsl_q4awo (&numrl2, rlist2, &reseps, &abseps, res3la, &nres);

	ktmin += 1;

	if (ktmin > 5 && *abserr < 1.0e-03 * errsum)

	    *ier = 5;

	if (abseps >= *abserr)

	    goto L_170;

	ktmin = 0;

	*abserr = abseps;

	*result = reseps;

	correc = erlarg;

	ertest = imsl_f_max (*epsabs, *epsrel * fabs (reseps));

	/* JUMP OUT OF DO-LOOP */

	if (*abserr < ertest)

	    goto L_200;

	/*

	 * PREPARE BISECTION OF THE SMALLEST INTERVAL.

	 */

L_170:

	if (numrl2 == 1)

	    noext = TRUE;

	if (*ier >= 5)

	    goto L_200;

L_180:

	maxerr = iord[0];

	errmax = elist[maxerr - 1];

	nrmax = 1;

	extrap = FALSE;

	levmax += 1;

	erlarg = errsum;

L_190:

	;

    }

    /* SET THE FINAL RESULT. */

L_200:

    if (*abserr == oflow)

	goto L_230;

    if ((*ier + ierro) == 0)

	goto L_220;

    if (ierro == 3)

	*abserr += correc;

    if (*ier == 0)

	*ier = 3;

    if (*result != 0.0e00 && area != 0.0e00)

	goto L_210;

    if (*abserr > errsum)

	goto L_230;

    if (area == 0.0e00)

	goto L_250;

    goto L_220;

L_210:

    if (*abserr / fabs (*result) > errsum / fabs (area))

	goto L_230;



    /*

     * TEST ON DIVERGENCE.

     */

L_220:

    if (ksgn == (-1) && imsl_f_max (fabs (*result), fabs (area)) <= resabs *

	1.0e-02)

	goto L_250;

    if ((1.0e-02 > (*result / area) || (*result / area) > 1.0e02) || errsum >

	fabs (area))

	*ier = 6;

    goto L_250;

    /* COMPUTE GLOBAL INTEGRAL SUM. */

L_230:

    *result = 0.0e00;

    for (k = 1; k <= *last; k++) {

	*result += rlist[k - 1];

    }

    *abserr = errsum;

L_250:

    if (*ier > 2)

	*ier -= 1;

    *result *= siign;

L_260:

    return;

}				/* end of function */













/*Translated by FOR_C++, v0.1, on 07/17/90 at 14:19:04 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 07/17/90 at 14:18:59

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  Q4AWO/DQ4AWO (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 29, 1985



    Purpose:    Integrate a function.



    Usage:      CALL Q4AWO (N, EPSTAB, RESULT, ABSERR, RES3LA, NRES)



    Arguments:  (See comment block below)



    Chapter:    MATH/LIBRARY Integration and Differentiation



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.





       ................................................................



   1.        Q4AWO

             EPSILON ALGORITHM

                STANDARD FORTRAN SUBROUTINE

                REAL VERSION



   2.        PURPOSE

                THE ROUTINE DETERMINES THE LIMIT OF A GIVEN SEQUENCE OF

                APPROXIMATIONS, BY MEANS OF THE EPSILON ALGORITHM

                OF P. WYNN.

                AN ESTIMATE OF THE ABSOLUTE ERROR IS ALSO GIVEN.

                THE CONDENSED EPSILON TABLE IS COMPUTED. ONLY THOSE

                ELEMENTS NEEDED FOR THE COMPUTATION OF THE NEXT DIAGONAL

                ARE PRESERVED.



   3.        CALLING SEQUENCE

                CALL Q4AWO(N,EPSTAB,RESULT,ABSERR,RES3LA,NRES)



             PARAMETERS

                N      - INTEGER

                         EPSTAB(N) CONTAINS THE NEW ELEMENT IN THE

                         FIRST COLUMN OF THE EPSILON TABLE.



                EPSTAB - REAL

                         VECTOR OF DIMENSION 52 CONTAINING THE ELEMENTS

                         OF THE TWO LOWER DIAGONALS OF THE

                         TRIANGULAR EPSILON TABLE

                         THE ELEMENTS ARE NUMBERED STARTING AT THE

                         RIGHT-HAND CORNER OF THE TRIANGLE.



                RESULT - REAL

                         RESULTING APPROXIMATION TO THE INTEGRAL



                ABSERR - REAL

                         ESTIMATE OF THE ABSOLUTE ERROR COMPUTED FROM

                         RESULT AND THE 3 PREVIOUS RESULTS



                RES3LA - REAL

                         VECTOR OF DIMENSION 3 CONTAINING THE LAST 3

                         RESULTS



                NRES   - INTEGER

                         NUMBER OF CALLS TO THE ROUTINE

                         (SHOULD BE ZERO AT FIRST CALL)



   4.        SUBROUTINES OR FUNCTIONS NEEDED

                       - Q4NG

                       - FORTRAN ABS, AMAX1



       ..................................................................





             LIST OF MAJOR VARIABLES

             -----------------------



             E0     - THE 4 ELEMENTS ON WHICH THE

             E1       COMPUTATION OF A NEW ELEMENT IN

             E2       THE EPSILON TABLE IS BASED

             E3                 E0

                          E3    E1    NEW

                                E2

             NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW

                      DIAGONAL

             ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)

             RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE

                      OF ERROR



             MACHINE DEPENDENT CONSTANTS

             ---------------------------



             EPMACH IS THE LARGEST RELATIVE SPACING.

             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.

             OFLOW IS THE LARGEST POSITIVE MAGNITUDE.

             LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON

             TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER

             DIAGONAL OF THE EPSILON TABLE IS DELETED.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

void 	    imsl_q4awo (Mint *n, Mfloat epstab[], Mfloat *result,

			Mfloat *abserr, Mfloat res3la[], Mint *nres)

#else

void        imsl_q4awo (n, epstab, result, abserr, res3la, nres)

    Mint        *n;

    Mfloat      epstab[], *result, *abserr, res3la[];

    Mint        *nres;

#endif

{

    Mint         i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num;

    Mfloat      delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epmach,

                epsinf, err1, err2, err3, error, oflow, res, ss, tol1,

                tol2, tol3, uflow;





    imsl_q4ng (&epmach, &uflow, &oflow);

    *nres += 1;

    *abserr = oflow;

    *result = epstab[*n - 1];

    if (*n < 3)

	goto L_100;

    limexp = 50;

    epstab[*n + 1] = epstab[*n - 1];

    newelm = (*n - 1) / 2;

    epstab[*n - 1] = oflow;

    num = *n;

    k1 = *n;

    for (i = 1; i <= newelm; i++) {

	k2 = k1 - 1;

	k3 = k1 - 2;

	res = epstab[k1 + 1];

	e0 = epstab[k3 - 1];

	e1 = epstab[k2 - 1];

	e2 = res;

	e1abs = fabs (e1);

	delta2 = e2 - e1;

	err2 = fabs (delta2);

	tol2 = imsl_f_max (fabs (e2), e1abs) * epmach;

	delta3 = e1 - e0;

	err3 = fabs (delta3);

	tol3 = imsl_f_max (e1abs, fabs (e0)) * epmach;

	if (err2 > tol2 || err3 > tol3)

	    goto L_10;

	/*

	 * IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE ACCURACY, CONVERGENCE

	 * IS ASSUMED. RESULT = E2 ABSERR = ABS(E1-E0)+ABS(E2-E1)

	 */

	*result = res;

	*abserr = err2 + err3;

	/* JUMP OUT OF DO-LOOP */

	goto L_100;

L_10:

	e3 = epstab[k1 - 1];

	epstab[k1 - 1] = e1;

	delta1 = e1 - e3;

	err1 = fabs (delta1);

	tol1 = imsl_f_max (e1abs, fabs (e3)) * epmach;

	/*

	 * IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT A PART OF THE

	 * TABLE BY ADJUSTING THE VALUE OF N

	 */

	if ((err1 <= tol1 || err2 <= tol2) || err3 <= tol3)

	    goto L_20;

	ss = 1.0e00 / delta1 + 1.0e00 / delta2 - 1.0e00 / delta3;

	epsinf = fabs (ss * e1);

	/*

	 * TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND EVENTUALLY

	 * OMIT A PART OF THE TABLE ADJUSTING THE VALUE OF N.

	 */

	if (epsinf > 1.0e-04)

	    goto L_30;

L_20:

	*n = i + i - 1;

	/* JUMP OUT OF DO-LOOP */

	goto L_50;

	/*

	 * COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST THE VALUE OF RESULT.

	 */

L_30:

	res = e1 + 1.0e00 / ss;

	epstab[k1 - 1] = res;

	k1 -= 2;

	error = err2 + fabs (res - e2) + err3;

	if (error > *abserr)

	    goto L_40;

	*abserr = error;

	*result = res;

L_40:

	;

    }

    /* SHIFT THE TABLE. */

L_50:

    if (*n == limexp)

	*n = 2 * (limexp / 2) - 1;

    ib = 1;

    if ((num / 2) * 2 == num)

	ib = 2;

    ie = newelm + 1;

    for (i = 1; i <= ie; i++) {

	ib2 = ib + 2;

	epstab[ib - 1] = epstab[ib2 - 1];

	ib = ib2;

    }

    if (num == *n)

	goto L_80;

    indx = num - *n + 1;

    for (i = 1; i <= *n; i++) {

	epstab[i - 1] = epstab[indx - 1];

	indx += 1;

    }

L_80:

    if (*nres >= 4)

	goto L_90;

    res3la[*nres - 1] = *result;

    *abserr = oflow;

    goto L_100;

    /*

     * COMPUTE ERROR ESTIMATE

     */

L_90:

    *abserr = fabs (*result - res3la[2]) + fabs (*result - res3la[1]) +

	fabs (*result - res3la[0]);

    res3la[0] = res3la[1];

    res3la[1] = res3la[2];

    res3la[2] = *result;

L_100:

    *abserr = imsl_f_max (*abserr, 5.0e00 * epmach * fabs (*result));

    return;

}				/* end of function */

