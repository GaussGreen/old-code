#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn_fourier (Mfloat (*fcn) (Mfloat), Mfloat a,
                Imsl_quad weight, Mfloat omega, va_list argptr);
static void l_q2awf (Mfloat (*f) (Mfloat), Mfloat *a, Mint *iweigh,
                Mfloat *omega, Mfloat *errabs, Mfloat *result,
                Mfloat *errest, Mint *maxcyl, Mint *maxsub,
                Mint *maxcby, Mint *neval, Mint *ncycle,
                Mfloat rslist[], Mfloat erlist[], Mint ierlst[],
                Mint *nsubin, Mfloat wk[], Mint *iwk);
static void l_q3awf (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *omega,
                Mint *integr, Mfloat *epsabs, Mint *limlst,
                Mint *limit, Mint *maxp1, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mint *ier,
                Mfloat rslst[], Mfloat erlst[], Mint ierlst[],
                Mint *lst, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint nnlog[], Mfloat *chebmo);
#else
static VA_LIST_HACK l_int_fcn_fourier ();
static void l_q2awf ();
static void l_q3awf ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn_fourier (Mfloat (*fcn) (Mfloat), Mfloat a,
                Imsl_quad weight, Mfloat omega,...)
#else
Mfloat      imsl_f_int_fcn_fourier (fcn, a, weight, omega, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Imsl_quad        weight;
    Mfloat      omega;
va_dcl
#endif
{
    va_list argptr;
    VA_START (argptr, omega);
    E1PSH ("imsl_f_int_fcn_fourier", "imsl_d_int_fcn_fourier");
    lv_value = F_ZERO;
    IMSL_CALL (l_int_fcn_fourier (fcn, a, weight, omega, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn_fourier", "imsl_d_int_fcn_fourier");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn_fourier (Mfloat (*fcn) (Mfloat), Mfloat a,
                Imsl_quad weight, Mfloat omega, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn_fourier (fcn, a, weight, omega, argptr)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Imsl_quad        weight;
    Mfloat      omega;
    va_list     argptr;
#endif
{
    Mfloat      a_float;
    Mfloat      omega_float;
    Mint        code;
    Mint        arg_number = 4;
    Mint 	l_weight = 0;
    Mfloat      err_abs;
    Mfloat     *rslist = NULL;
    Mfloat     *erlist = NULL;
    Mint       *ierlst = NULL;
    Mfloat     *err_est = NULL;
    Mint        max_subinter = 500;
    Mint        maxcyl = 50;
    Mint        maxcby = 21;
    Mint       *n_subinter = NULL;
    Mint       *n_evals = NULL;
    Mfloat     *wk = NULL;
    Mint       *iwk = NULL;
    Mint        user_err_list = 0;
    Mfloat      temp_err_est;
    Mint        temp_n_subinter;
    Mint        temp_n_evals;
    Mint        temp_ncycle;
    Mint       *n_cycle = NULL;
    err_abs = sqrt((double) imsl_amach (4));

    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
        case IMSL_ERR_ABS_ADR:
            arg_number++;
            err_abs = *(va_arg (argptr, Mfloat *));
            break;
	case IMSL_ERR_ABS:
	    arg_number++;
	    err_abs = (Mfloat) va_arg (argptr, Mdouble);
	    break;
	case IMSL_ERR_EST:
	    arg_number++;
	    err_est = va_arg (argptr, Mfloat *);
	    break;
	case IMSL_MAX_SUBINTER:
	    arg_number++;
	    max_subinter = va_arg (argptr, Mint);
	    break;
	case IMSL_N_EVALS:
	    arg_number++;
	    n_evals = va_arg (argptr, Mint *);
	    break;
	case IMSL_MAX_MOMENTS:
	    arg_number++;
	    maxcby = va_arg (argptr, Mint);
	    break;
	case IMSL_MAX_CYCLES:
	    arg_number++;
	    maxcyl = va_arg (argptr, Mint);
	    break;
	case IMSL_N_CYCLES:
	    arg_number++;
	    n_cycle = va_arg (argptr, Mint *);
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
	imsl_e1sti (1, max_subinter);
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_SUBINTER_SMALL);
    }
    if (imsl_n1rty (0))
	goto RETURN;
    if (maxcby < 1) {
	imsl_e1sti (1, maxcby);
	/*
	 * (5, 3, "The maximum number of moments maxcby = %(i1).  It must be
	 * at least 1.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_MOMENTS_REACHED);
    }
    if (imsl_n1rty (0))
	goto RETURN;
    if (maxcyl < 3) {
	imsl_e1sti (1, maxcyl);
	/*
	 * (5, 2, "The maximum number of cycles maxcyl = %(i1).  It must be
	 * at least 3.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_CYCLES_REACHED);
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
	erlist = (Mfloat *) imsl_malloc (maxcyl * sizeof (*erlist));

    ierlst = (Mint *) imsl_malloc (maxcyl * sizeof (*ierlst));

    rslist = (Mfloat *) imsl_malloc (maxcyl * sizeof (*rslist));
    wk = (Mfloat *) imsl_malloc ((4 * max_subinter + 25 * maxcby) * sizeof (*wk));
    iwk = (Mint *) imsl_malloc (2 * max_subinter * sizeof (*iwk));

    if (erlist == NULL || ierlst == NULL || rslist == NULL
	|| wk == NULL || iwk == NULL) {
	/* Not enough memory, with %(L1) = %(I1). */
	imsl_e1sti (1, max_subinter);
	imsl_e1sti (2, maxcyl);
	imsl_e1stl (1, "max_subinter");
	imsl_e1stl (2, "max_cycles");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	goto FREE_SPACE;
    }

    if (err_est == NULL)
	err_est = &temp_err_est;
    if (n_subinter == NULL)
	n_subinter = &temp_n_subinter;
    if (n_evals == NULL)
	n_evals = &temp_n_evals;
    if (n_cycle == NULL)
	n_cycle = &temp_ncycle;

    a_float = a;
    omega_float = omega;

    if (weight == IMSL_COS)
	l_weight = 1;
    if (weight == IMSL_SIN)
	l_weight = 2;
    l_q2awf (fcn, &a_float, &l_weight, &omega_float, &err_abs, &lv_value,
	err_est, &maxcyl, &max_subinter, &maxcby, n_evals, n_cycle,
	rslist, erlist, ierlst, n_subinter, wk, iwk);
FREE_SPACE:
    ;
    if (!user_err_list && erlist != NULL)
	imsl_free (erlist);
    if (ierlst != NULL)
	imsl_free (ierlst);
    if (rslist != NULL)
	imsl_free (rslist);
    if (wk != NULL)
	imsl_free (wk);
    if (iwk != NULL)
	imsl_free (iwk);
RETURN:
    ;
    if (imsl_n1rty(0) > 3)
	lv_value = imsl_amach(6);
    return (argptr);
}








/*Translated by FOR_C++, v0.1, on 08/14/90 at 11:24:12 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/14/90 at 11:24:10
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q2AWF/DQ2AWF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Compute a Fourier integral.

    Usage:      CALL Q2AWF (F, A, IWEIGH, OMEGA, ERRABS, RESULT, ERREST,
                            MAXCYL, MAXSUB, MAXCBY, NEVAL, NCYCLE,
                            RSLIST, ERLIST, IERLST, NSUBIN, WK, IWK)

    Arguments:  (See QDAWF)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* NSUBIN is not used here, but leave the calling sequence intact. */
#ifdef ANSI
static void l_q2awf (Mfloat (*f) (Mfloat), Mfloat *a, Mint *iweigh,
                Mfloat *omega, Mfloat *errabs, Mfloat *result,
                Mfloat *errest, Mint *maxcyl, Mint *maxsub,
                Mint *maxcby, Mint *neval, Mint *ncycle,
                Mfloat rslist[], Mfloat erlist[], Mint ierlst[],
                Mint *nsubin, Mfloat wk[], Mint *iwk)
#else
static void l_q2awf (f, a, iweigh, omega, errabs, result, errest,
                maxcyl, maxsub, maxcby, neval, ncycle, rslist, erlist, ierlst,
                nsubin, wk, iwk)
    Mfloat      (*f) (), *a;
    Mint       *iweigh;
    Mfloat     *omega, *errabs, *result, *errest;
    Mint       *maxcyl, *maxsub, *maxcby, *neval, *ncycle;
    Mfloat      rslist[], erlist[];
    Mint        ierlst[], *nsubin;
    Mfloat      wk[];
    Mint       *iwk;
#endif
{
    Mint        ier;


    imsl_e1psh ("Q2AWF ");
    if (*maxsub < 1) {
	imsl_e1sti (1, *maxsub);
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_SUBINTER_SMALL);
    }
    if (*maxcby < 1) {
	imsl_e1sti (1, *maxcby);
	/*
	 * (5, 3, "The maximum number of moments maxcby = %(i1).  It must be
	 * at least 1.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_MOMENTS_REACHED);
    }
    if (*maxcyl < 3) {
	imsl_e1sti (1, *maxcyl);
	/*
	 * (5, 2, "The maximum number of cycles maxcyl = %(i1).  It must be
	 * at least 3.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_CYCLES_REACHED);
    }
    /* CHECK IWEIGH */
    if (*iweigh < 1 || *iweigh > 2) {
	imsl_e1sti (1, *iweigh);

	/*
	 * (5, 4, "The weight choice, weight, is equal to ' %(i1)'.  Valid
	 * choices for the variable weight include IMSL_COS and IMSL_SIN.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_WEIGHT_CHOICE);
    }
    /* CHECK ERRABS */
    if (*errabs <= F_ZERO) {
	imsl_e1str (1, *errabs);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_ABS_SMALL);
    }
    if (imsl_n1rty (0) == 5)
	goto L_9000;

    l_q3awf (f, a, omega, iweigh, errabs, maxcyl, maxsub, maxcby, result,
	errest, neval, &ier, rslist, erlist, ierlst, ncycle, &wk[0],
	&wk[*maxsub], &wk[*maxsub * 2], &wk[*maxsub * 3], iwk,
	iwk + (*maxsub), &wk[*maxsub * 4]);

    if (ier == 7) {

	/*
	 * (3, 1, "Bad integrand behavior occurs within one or more of the
	 * cycles.");
	 */
	imsl_ermes (IMSL_WARNING, IMSL_BAD_INTEGRAND_BEHAVIOR);
    }
    else if (ier == 8) {
	imsl_e1sti (1, *maxcyl);

	/*
	 * (4, 2, "Maximum number of cycles allowed, maxcyl = %(i1), has been
	 * achieved.");
	 */
	imsl_ermes (IMSL_FATAL, IMSL_MAX_CYCLES_ACHIEVED);
    }
    else if (ier == 9) {
	imsl_e1str (1, *errabs);
	/*
	 * (3, 3, "Extrapolation table constructed for convergence
	 * acceleration of the series formed by the integral contributions of
	 * the cycles, does not converge to the requested accuracy, err_abs =
	 * %(r1).");
	 */
	imsl_ermes (IMSL_WARNING, IMSL_EXTRAPOLATION_PROBLEMS);
    }
L_9000:
    ;
    imsl_e1pop ("Q2AWF ");
    return;
}				/* end of function */




























/*Translated by FOR_C++, v0.1, on 08/14/90 at 14:27:01 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/14/90 at 14:26:59
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3AWF/DQ3AWF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q3AWF (F, A, OMEGA, INTEGR, EPSABS, LIMLST, LIMIT,
                            MAXP1, RESULT, ABSERR, NEVAL, IER, RSLST,
                            ERLST, IERLST, LST, ALIST, BLIST, RLIST,
                            ELIST, IORD, NNLOG, CHEBMO)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


   ......................................................................

   1.     Q3AWF
          COMPUTATION OF FOURIER INTEGRALS
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A
             GIVEN FOURIER INTEGRAL
                   I = INTEGRAL OF F(X)*W(X) OVER (A,INFINITY)
                       WHERE W(X) = COS(OMEGA*X)
                          OR W(X) = SIN(OMEGA*X),
             HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.EPSABS.

   3.     CALLING SEQUENCE
             CALL Q3AWF(F,A,OMEGA,INTEGR,EPSABS,LIMLST,LIMIT,MAXP1,
                  RESULT,ABSERR,NEVAL,IER,RSLST,ERLST,IERLST,LST,ALIST,
                  BLIST,RLIST,ELIST,IORD,NNLOG,CHEBMO)

          PARAMETERS
           ON ENTRY
              F      - REAL
                       FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                       FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO
                       BE DECLARED E X T E R N A L IN THE DRIVER
                       PROGRAM.

              A      - REAL
                       LOWER LIMIT OF INTEGRATION

              OMEGA  - REAL
                       PARAMETER IN THE WEIGHT FUNCTION

              INTEGR - INTEGER
                       INDICATES WHICH WEIGHT FUNCTION IS USED
                       INTEGR = 1      W(X) = COS(OMEGA*X)
                       INTEGR = 2      W(X) = SIN(OMEGA*X)
                       IF INTEGR.NE.1.AND.INTEGR.NE.2, THE ROUTINE WILL
                       END WITH IER = 6.

              EPSABS - REAL
                       ABSOLUTE ACCURACY REQUESTED, EPSABS.GT.0
                       IF EPSABS.LE.0, THE ROUTINE WILL END WITH IER = 6.

              LIMLST - INTEGER
                       LIMLST GIVES AN UPPER BOUND ON THE NUMBER OF
                       CYCLES, LIMLST.GE.1.  IF LIMLST.LT.3,
                       THE ROUTINE WILL END WITH IER = 6.

              LIMIT  - INTEGER
                       GIVES AN UPPER BOUND ON THE NUMBER OF SUBINTERVALS
                       ALLOWED IN THE PARTITION OF EACH CYCLE,
                       LIMIT.GE.1.

              MAXP1  - INTEGER
                       GIVES AN UPPER BOUND ON THE NUMBER OF
                       CHEBYSHEV MOMENTS WHICH CAN BE STORED, I.E.
                       FOR THE INTERVALS OF LENGTHS ABS(B-A)*2**(-L),
                       L=0,1, ..., MAXP1-2, MAXP1.GE.1

           ON RETURN
              RESULT - REAL
                       APPROXIMATION TO THE INTEGRAL

              ABSERR - REAL
                       ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                       WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

              NEVAL  - INTEGER
                       NUMBER OF INTEGRAND EVALUATIONS

              IER    - IER = 0 NORMAL AND RELIABLE TERMINATION OF
                               THE ROUTINE. IT IS ASSUMED THAT THE
                               REQUESTED ACCURACY HAS BEEN ACHIEVED.
                       IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE
                               THE ESTIMATES FOR INTEGRAL AND ERROR
                               ARE LESS RELIABLE. IT IS ASSUMED THAT
                               THE REQUESTED ACCURACY HAS NOT BEEN
                               ACHIEVED.
                      IF OMEGA.NE.0
                       IER = 6 THE INPUT IS INVALID BECAUSE
                               (INTEGR.NE.1 AND INTEGR.NE.2) OR
                                EPSABS.LE.0 OR LIMLST.LT.3.
                                RESULT, ABSERR, NEVAL, LST ARE SET
                                TO ZERO.
                           = 7 BAD INTEGRAND BEHAVIOUR OCCURS WITHIN
                               ONE OR MORE OF THE CYCLES. LOCATION
                               AND TYPE OF THE DIFFICULTY INVOLVED
                               CAN BE DETERMINED FROM THE VECTOR IERLST.
                               HERE LST IS THE NUMBER OF CYCLES ACTUALLY
                               NEEDED (SEE BELOW).
                           = 8 MAXIMUM NUMBER OF CYCLES ALLOWED
                               HAS BEEN ACHIEVED., I.E. OF SUBINTERVALS
                               (A+(K-1)C,A+KC) WHERE
                               C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
                               FOR K = 1, 2, ..., LST.
                               ONE CAN ALLOW MORE CYCLES BY INCREASING
                               THE VALUE OF LIMLST (AND TAKING THE
                               ACCORDING DIMENSION ADJUSTMENTS INTO
                               ACCOUNT).
                               EXAMINE THE ARRAY IERLST WHICH CONTAINS
                               THE ERROR FLAGS ON THE CYCLES, IN ORDER
                               TO EVENTUAL LOOK FOR LOCAL INTEGRATION
                               DIFFICULTIES.
                               IF THE POSITION OF A LOCAL DIFFICULTY CAN
                               BE DETERMINED (E.G. SINGULARITY,
                               DISCONTINUITY WITHIN THE INTERVAL)
                               ONE WILL PROBABLY GAIN FROM SPLITTING UP
                               THE INTERVAL AT THIS POINT AND
                               CALLING APPROPRIATE INTEGRATORS ON
                               THE SUBRANGES.
                           = 9 THE EXTRAPOLATION TABLE CONSTRUCTED FOR
                               CONVERGENCE ACCELERATION OF THE SERIES
                               FORMED BY THE INTEGRAL CONTRIBUTIONS
                               OVER THE CYCLES, DOES NOT CONVERGE TO
                               WITHIN THE REQUIRED ACCURACY.
                               AS IN THE CASE OF IER = 8, IT IS ADVISED
                               TO EXAMINE THE ARRAY IERLST WHICH
                               CONTAINS THE ERROR FLAGS ON THE CYCLES.
                               IERLST(K) = 1 THE MAXIMUM NUMBER OF
                                             SUBDIVISIONS (= LIMIT)
                                             HAS BEEN ACHIEVED ON THE
                                             K TH CYCLE.
                                         = 2 OCCURRENCE OF ROUNDOFF
                                             ERROR IS DETECTED AND
                                             PREVENTS THE TOLERANCE
                                             IMPOSED ON THE KTH CYCLE,
                                             FROM BEING ACHIEVED.
                                         = 3 EXTREMELY BAD INTEGRAND
                                             BEHAVIOUR OCCURS AT SOME
                                             POINTS OF THE K TH CYCLE.
                                         = 4 THE INTEGRATION PROCEDURE
                                             OVER THE K TH CYCLE DOES
                                             NOT CONVERGE (TO WITHIN THE
                                             REQUIRED ACCURACY) DUE TO
                                             ROUNDOFF IN THE
                                             EXTRAPOLATION PROCEDURE
                                             INVOKED ON THIS CYCLE. IT
                                             IS ASSUMED THAT THE RESULT
                                             ON THIS INTERVAL IS THE
                                             BEST WHICH CAN BE OBTAINED.
                                         = 5 THE INTEGRAL OVER THE K TH
                                             CYCLE IS PROBABLY DIVERGENT
                                             OR SLOWLY CONVERGENT. IT
                                             MUST BE NOTED THAT
                                             DIVERGENCE CAN OCCUR WITH
                                             ANY OTHER VALUE OF
                                             IERLST(K).
                      IF OMEGA = 0 AND INTEGR = 1,
                      THE INTEGRAL IS CALCULATED BY MEANS OF Q3AGI
                      AND IER = IERLST(1) (WITH MEANING AS DESCRIBED
                      FOR IERLST(K), K = 1).

              RSLST  - REAL
                       VECTOR OF DIMENSION AT LEAST LIMLST
                       RSLST(K) CONTAINS THE INTEGRAL CONTRIBUTION
                       OVER THE INTERVAL (A+(K-1)C,A+KC) WHERE
                       C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
                       K = 1, 2, ..., LST.
                       NOTE THAT, IF OMEGA = 0, RSLST(1) CONTAINS
                       THE VALUE OF THE INTEGRAL OVER (A,INFINITY).

              ERLST  - REAL
                       VECTOR OF DIMENSION AT LEAST LIMLST
                       ERLST(K) CONTAINS THE ERROR ESTIMATE
                       CORRESPONDING WITH RSLST(K).

              IERLST - INTEGER
                       VECTOR OF DIMENSION AT LEAST LIMLST
                       IERLST(K) CONTAINS THE ERROR FLAG CORRESPONDING
                       WITH RSLST(K). FOR THE MEANING OF THE LOCAL ERROR
                       FLAGS SEE DESCRIPTION OF OUTPUT PARAMETER IER.

              LST    - INTEGER
                       NUMBER OF SUBINTERVALS NEEDED FOR THE INTEGRATION
                       IF OMEGA = 0 THEN LST IS SET TO 1.

              ALIST, BLIST, RLIST, ELIST - REAL
                       VECTOR OF DIMENSION AT LEAST LIMIT,

              IORD, NNLOG - INTEGER
                       VECTOR OF DIMENSION AT LEAST LIMIT, PROVIDING
                       SPACE FOR THE QUANTITIES NEEDED IN THE SUBDIVISION
                       PROCESS OF EACH CYCLE

              CHEBMO - REAL
                       ARRAY OF DIMENSION AT LEAST (MAXP1,25), PROVIDING
                       SPACE FOR THE CHEBYSHEV MOMENTS NEEDED WITHIN THE
                       CYCLES

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q3AWO
                - Q3AGI
                - Q4AWO
                - Q4AGI
                - Q0AG
                - Q5AWO
                - Q8AWO
                - Q7AWO
                - F (USER-PROVIDED FUNCTION)
                - Q6AWO
                - Q4NG
                - FORTRAN ABS, AMAX1, AMIN1

   ......................................................................





              THE DIMENSION OF  PSUM  IS DETERMINED BY THE VALUE OF
              LIMEXP IN SUBROUTINE Q4AWO (PSUM MUST BE
              OF DIMENSION (LIMEXP+2) AT LEAST).

             LIST OF MAJOR VARIABLES
             -----------------------

             C1, C2    - END POINTS OF SUBINTERVAL (OF LENGTH CYCLE)
             CYCLE     - (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA)
             PSUM      - VECTOR OF DIMENSION AT LEAST (LIMEXP+2)
                         (SEE ROUTINE Q4AWO)
                         PSUM CONTAINS THE PART OF THE EPSILON TABLE
                         WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS.
                         EACH ELEMENT OF PSUM IS A PARTIAL SUM OF
                         THE SERIES WHICH SHOULD SUM TO THE VALUE OF
                         THE INTEGRAL.
             ERRSUM    - SUM OF ERROR ESTIMATES OVER THE
                         SUBINTERVALS, CALCULATED CUMULATIVELY
             EPSA      - ABSOLUTE TOLERANCE REQUESTED OVER CURRENT
                         SUBINTERVAL
             CHEBMO    - ARRAY CONTAINING THE MODIFIED CHEBYSHEV
                         MOMENTS (SEE ALSO ROUTINE Q5AWO)


             MACHINE DEPENDENT CONSTANTS
             --------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q3awf (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *omega,
                Mint *integr, Mfloat *epsabs, Mint *limlst,
                Mint *limit, Mint *maxp1, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mint *ier,
                Mfloat rslst[], Mfloat erlst[], Mint ierlst[],
                Mint *lst, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint nnlog[], Mfloat *chebmo)
#else
static void l_q3awf (f, a, omega, integr, epsabs, limlst, limit,
                maxp1, result, abserr, neval, ier, rslst, erlst, ierlst, lst,
                alist, blist, rlist, elist, iord, nnlog, chebmo)
    Mfloat      (*f) (), *a, *omega;
    Mint       *integr;
    Mfloat     *epsabs;
    Mint       *limlst, *limit, *maxp1;
    Mfloat     *result, *abserr;
    Mint       *neval, *ier;
    Mfloat      rslst[], erlst[];
    Mint        ierlst[], *lst;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[], nnlog[];
    Mfloat     *chebmo;
#endif
{
#define CHEBMO(I_,J_)	(chebmo+(I_)*(*maxp1)+(J_))
    Mint        _l0, ktmin, l, last, ll, momcom, nev, nres, numrl2;
    Mfloat      _f0, _f1, abseps, c1, c2, correc, dl, drl, ep, epmach,
                eps, epsa, errsum, fact, oflow, p1, psum[52], res3la[3],
                reseps, tcycle, uflow;
    static Mfloat p = 9.0e-01;
    static Mfloat pi = 3.14159265358979323846264338328e00;



    imsl_q4ng (&epmach, &uflow, &oflow);
    /* TEST ON VALIDITY OF PARAMETERS */
    *result = F_ZERO;
    *abserr = F_ZERO;
    *neval = 0;
    *lst = 0;
    *ier = 0;
    if (((*integr != 1 && *integr != 2) || *epsabs <= F_ZERO) || *limlst <
	3)
	*ier = 6;
    if (*ier == 6)
	goto L_90;
    if (*omega != F_ZERO)
	goto L_10;
    /*
     * INTEGRATION BY Q3AGI IF OMEGA IS ZERO
     */
    if (*integr == 1)
	_f0 = F_ZERO;
    _l0 = 1;
    _f1 = F_ZERO;
    imsl_q3agi (f, &_f0, &_l0, epsabs, &_f1,
	limit, result, abserr, neval, ier, alist, blist, rlist, elist,
	iord, &last);
    rslst[0] = *result;
    erlst[0] = *abserr;
    ierlst[0] = *ier;
    *lst = 1;
    goto L_90;
    /* INITIALIZATIONS */
L_10:
    l = fabs (*omega);
    dl = 2 * l + 1;
    tcycle = dl * pi;
    tcycle /= fabs (*omega);
    *ier = 0;
    ktmin = 0;
    *neval = 0;
    numrl2 = 0;
    nres = 0;
    c1 = *a;
    c2 = tcycle + *a;
    p1 = F_ONE - p;
    eps = *epsabs;
    if (*epsabs > uflow / p1)
	eps = *epsabs * p1;
    ep = eps;
    fact = F_ONE;
    correc = F_ZERO;
    *abserr = F_ZERO;
    errsum = F_ZERO;
    /* MAIN DO-LOOP */
    for (*lst = 1; *lst <= *limlst; (*lst)++) {
	/* INTEGRATE OVER CURRENT SUBINTERVAL. */
	epsa = eps * fact;
	_f0 = F_ZERO;
	imsl_q3awo (f, &c1, &c2, omega, integr, &epsa, &_f0,
	    limit, lst, maxp1, &rslst[*lst - 1], &erlst[*lst - 1], &nev,
	    &ierlst[*lst - 1], alist, blist, rlist, elist, iord, nnlog,
	    &last, &momcom, chebmo);
	*neval += nev;
	fact *= p;
	errsum += erlst[*lst - 1];
	drl = 5.0e01 * fabs (rslst[*lst - 1]);
	/*
	 * TEST ON ACCURACY WITH PARTIAL SUM
	 */
	if ((errsum + drl) <= *epsabs && *lst >= 6)
	    goto L_80;
	correc = imsl_f_max (correc, erlst[*lst - 1]);
	if (ierlst[*lst - 1] != 0)
	    eps = imsl_f_max (ep, correc * p1);
	if (ierlst[*lst - 1] != 0)
	    *ier = 7;
	if ((*ier == 7 && (errsum + drl) <= correc * F_TEN) && *lst >
	    5)
	    goto L_80;
	numrl2 += 1;
	if (*lst > 1)
	    goto L_20;
	psum[0] = rslst[0];
	goto L_40;
L_20:
	psum[numrl2 - 1] = psum[ll - 1] + rslst[*lst - 1];
	if (*lst == 2)
	    goto L_40;
	/*
	 * TEST ON MAXIMUM NUMBER OF SUBINTERVALS
	 */
	if (*lst == *limlst)
	    *ier = 8;
	/*
	 * PERFORM NEW EXTRAPOLATION
	 */
	imsl_q4awo (&numrl2, psum, &reseps, &abseps, res3la, &nres);

	/*
	 * TEST WHETHER EXTRAPOLATED RESULT IS INFLUENCED BY ROUNDOFF
	 */
	ktmin += 1;
	if (ktmin >= 15 && *abserr <= 1.0e-03 * (errsum + drl))
	    *ier = 9;
	if (abseps > *abserr && *lst != 3)
	    goto L_30;
	*abserr = abseps;
	*result = reseps;
	ktmin = 0;
	/*
	 * IF IER IS NOT 0, CHECK WHETHER DIRECT RESULT (PARTIAL SUM) OR
	 * EXTRAPOLATED RESULT YIELDS THE BEST INTEGRAL APPROXIMATION
	 */
	if ((*abserr + F_TEN * correc) <= *epsabs || (*abserr <= *epsabs &&
		F_TEN * correc >= *epsabs))
	    goto L_60;
L_30:
	if (*ier != 0 && *ier != 7)
	    goto L_60;
L_40:
	ll = numrl2;
	c1 = c2;
	c2 += tcycle;
    }
    /* SET FINAL RESULT AND ERROR ESTIMATE */
L_60:
    *abserr += F_TEN * correc;
    if (*ier == 0)
	goto L_90;
    if (*result != F_ZERO && psum[numrl2 - 1] != F_ZERO)
	goto L_70;
    if (*abserr > errsum)
	goto L_80;
    if (psum[numrl2 - 1] == F_ZERO)
	goto L_90;
L_70:
    if (*abserr / fabs (*result) > (errsum + drl) / fabs (psum[numrl2 - 1]))
	goto L_80;
    if (*ier >= 1 && *ier != 7)
	*abserr += drl;
    goto L_90;
L_80:
    *result = psum[numrl2 - 1];
    *abserr = errsum + drl;
L_90:
    return;
}				/* end of function */

#undef CHEBMO
