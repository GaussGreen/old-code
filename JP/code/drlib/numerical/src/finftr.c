#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn_trig (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Imsl_quad weight, Mfloat omega, va_list argptr);
static void l_q2awo (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mint *iweigh, Mfloat *omega, Mfloat *errabs,
                Mfloat *errrel, Mfloat *result, Mfloat *errest,
                Mint *maxsub, Mint *maxcby, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint nnlog[], Mfloat wk[]);
static void l_q5awo (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mfloat *omega, Mint *integr, Mint *nrmom,
                Mint *maxp1, Mint *ksave, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mfloat *resabs,
                Mfloat *resasc, Mint *momcom, Mfloat *chebmo);
static Mfloat l_q6awo (Mfloat *x, Mfloat *omega, Mfloat *p2, Mfloat *p3,
                Mfloat *p4, Mint *integr);
#else
static VA_LIST_HACK l_int_fcn_trig ();
static void l_q2awo ();
static void l_q5awo ();
static Mfloat l_q6awo ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn_trig (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Imsl_quad weight, Mfloat omega,...)
#else
Mfloat      imsl_f_int_fcn_trig (fcn, a, b, weight, omega, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    Imsl_quad        weight;
    Mfloat      omega;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, omega);
    E1PSH ("imsl_f_int_fcn_trig", "imsl_d_int_fcn_trig");
    lv_value = F_ZERO;
    IMSL_CALL (l_int_fcn_trig (fcn, a, b, weight, omega, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn_trig", "imsl_d_int_fcn_trig");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn_trig (Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                Imsl_quad weight, Mfloat omega, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn_trig (fcn, a, b, weight, omega, argptr)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
    Imsl_quad        weight;
    Mfloat      omega;
    va_list     argptr;
#endif
{
    Mfloat      a_float;
    Mfloat      b_float;
    Mfloat      omega_float;
    Mint        code;
    Mint 	l_weight=0;
    Mint        arg_number = 5;
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
    Mfloat     *wk = NULL;
    Mint       *nnlog = NULL;
    Mint        maxcby = 21;
    err_abs = sqrt( (double) imsl_amach (4));
    err_rel = sqrt( (double) imsl_amach (4));

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
	case IMSL_MAX_MOMENTS:
	    arg_number++;
	    maxcby = va_arg (argptr, Mint);
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

    if (maxcby < 1) {
	imsl_e1sti (1, maxcby);
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_MOMENTS_REACHED);
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
    nnlog = (Mint *) imsl_malloc (max_subinter * sizeof (*nnlog));
    wk = (Mfloat *) imsl_malloc (25 * maxcby * sizeof (*wk));

    if (elist == NULL || iord == NULL || alist == NULL || blist == NULL ||
	rlist == NULL || nnlog == NULL || wk == NULL) {
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
    omega_float = omega;

    if (weight == IMSL_COS)
	l_weight = 1;
    if (weight == IMSL_SIN)
	l_weight = 2;

    l_q2awo (fcn, &a_float, &b_float, &l_weight, &omega_float, &err_abs,
	&err_rel, &lv_value, err_est, &max_subinter, &maxcby,
	n_evals, n_subinter, alist, blist, rlist, elist, iord,
	nnlog, wk);
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
    if (nnlog != NULL)
	imsl_free (nnlog);
    if (wk != NULL)
	imsl_free (wk);
RETURN:
    ;
    if (imsl_n1rty (0) > 3)
	lv_value = imsl_amach(6);
    return (argptr);
}





/*Translated by FOR_C++, v0.1, on 08/10/90 at 09:54:32 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 09:54:28
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q2AWO/DQ2AWO (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate an oscillatory function.

    Usage:      CALL Q2AWO (F, A, B, IWEIGH, OMEGA, ERRABS, ERRREL,
                            RESULT, ERREST, MAXSUB, MAXCBY, NEVAL,
                            NSUBIN, ALIST, BLIST, RLIST, ELIST, IORD,
                            NNLOG, WK)

    Arguments:  (See QDAWO)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_q2awo (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mint *iweigh, Mfloat *omega, Mfloat *errabs,
                Mfloat *errrel, Mfloat *result, Mfloat *errest,
                Mint *maxsub, Mint *maxcby, Mint *neval,
                Mint *nsubin, Mfloat alist[], Mfloat blist[],
                Mfloat rlist[], Mfloat elist[], Mint iord[],
                Mint nnlog[], Mfloat wk[])
#else
static void l_q2awo (f, a, b, iweigh, omega, errabs, errrel, result,
                errest, maxsub, maxcby, neval, nsubin, alist, blist, rlist, elist,
                iord, nnlog, wk)
    Mfloat      (*f) (), *a, *b;
    Mint       *iweigh;
    Mfloat     *omega, *errabs, *errrel, *result, *errest;
    Mint       *maxsub, *maxcby, *neval, *nsubin;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[], nnlog[];
    Mfloat      wk[];
#endif
{
    Mint        _l0, ier, momcom;


    imsl_e1psh ("l_q2awo");

    /* CHECK MAXCBY */
    if (*maxcby < 1) {
	imsl_e1sti (1, *maxcby);
	/*
	 * (5, 2, "The maximum number of moments, maxcby, is equal to %(i1).
	 * It must be at least 1.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_MAX_MOMENTS_REACHED);
    }
    /* CHECK IWEIGH */
    if (*iweigh < 1 || *iweigh > 2) {
	imsl_e1sti (1, *iweigh);
	/*
	 * (5, 3, "Weight choice is '%(i1)'.  Valid choices for the variable
	 * weight include IMSL_COS and IMSL_SIN.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_WEIGHT_CHOICE);
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
    /* CHECK ERRABS AND ERRREL */
    if (*errabs == F_ZERO && *errrel == F_ZERO) {
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_TOL_ZERO);
    }
    /* CHECK ERRREL .GE. 1 */
    if (*errrel >= F_ONE) {
	imsl_e1str (1, *errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_BIG);
    }
    if (imsl_n1rty (0) == 5)
	goto L_9000;

    _l0 = 1;
    imsl_q3awo (f, a, b, omega, iweigh, errabs, errrel, maxsub, &_l0,
	maxcby, result, errest, neval, &ier, alist, blist, rlist, elist,
	iord, nnlog, nsubin, &momcom, wk);

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
	 * (4, 5, "The integral is probably divergent or slowly
	 * convergent.");
	 */
	imsl_ermes (IMSL_FATAL, IMSL_DIVERGENT);
    }
L_9000:
    ;
    imsl_e1pop ("l_q2awo");
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 09:55:39 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 09:55:34
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3AWO/DQ3AWO (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q3AWO (F, A, B, OMEGA, INTEGR, EPSABS, EPSREL,
                            LIMIT, ICALL, MAXP1, RESULT, ABSERR, NEVAL,
                            IER, ALIST, BLIST, RLIST, ELIST, IORD, NNLOG,
                            LAST, MOMCOM, CHEBMO)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


   ......................................................................

   1.     Q3AWO
          COMPUTATION OF OSCILLATORY INTEGRALS
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
             DEFINITE INTEGRAL
                I = INTEGRAL OF F(X)*W(X) OVER (A,B)
                    WHERE W(X) = COS(OMEGA*X)
                       OR W(X) = SIN(OMEGA*X),
             HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
             Q3AWO IS CALLED BY DQAWO AND DQAWF.
             HOWEVER, IT CAN ALSO BE CALLED DIRECTLY IN A USER-WRITTEN
             PROGRAM. IN THE LATTER CASE IT IS POSSIBLE FOR THE USER
             TO DETERMINE THE FIRST DIMENSION OF ARRAY CHEBMO(MAXP1,25).
             SEE ALSO PARAMETER DESCRIPTION OF MAXP1. ADDITIONALLY SEE
             PARAMETER DESCRIPTION OF ICALL FOR EVENTUALLY RE-USING
             CHEBYSHEV MOMENTS COMPUTED DURING FORMER CALL ON SUBINTERVAL
             OF EQUAL LENGTH ABS(B-A).

   3.     CALLING SEQUENCE
             CALL  Q3AWO(F,A,B,OMEGA,INTEGR,EPSABS,EPSREL,LIMIT,ICALL,
                         MAXP1,RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,
                         ELIST,IORD,NNLOG,MOMCOM,CHEBMO)

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

              OMEGA  - REAL
                       PARAMETER IN THE INTEGRAND WEIGHT FUNCTION

              INTEGR - INTEGER
                       INDICATES WHICH OF THE WEIGHT FUNCTIONS IS TO BE
                       USED
                       INTEGR = 1      W(X) = COS(OMEGA*X)
                       INTEGR = 2      W(X) = SIN(OMEGA*X)
                       IF INTEGR.NE.1 AND INTEGR.NE.2, THE ROUTINE
                       WILL END WITH IER = 6.

              EPSABS - REAL
                       ABSOLUTE ACCURACY REQUESTED
              EPSREL - REAL
                       RELATIVE ACCURACY REQUESTED
                       IF  EPSABS.LT.0 AND EPSREL.LT.0,
                       THE ROUTINE WILL END WITH IER = 6.

              LIMIT  - INTEGER
                       GIVES AN UPPER BOUND ON THE NUMBER OF SUBDIVISIONS
                       IN THE PARTITION OF (A,B), LIMIT.GE.1.

              ICALL  - INTEGER
                       IF Q3AWO IS TO BE USED ONLY ONCE, ICALL MUST
                       BE SET TO 1.  ASSUME THAT DURING THIS CALL, THE
                       CHEBYSHEV MOMENTS (FOR CLENSHAW-CURTIS INTEGRATION
                       OF DEGREE 24) HAVE BEEN COMPUTED FOR INTERVALS OF
                       LENGHTS (ABS(B-A))*2**(-L), L=0,1,2,...MOMCOM-1.
                       THE CHEBYSHEV MOMENTS ALREADY COMPUTED CAN BE
                       RE-USED IN SUBSEQUENT CALLS, IF Q3AWO MUST BE
                       CALLED TWICE OR MORE TIMES ON INTERVALS OF THE
                       SAME LENGTH ABS(B-A). FROM THE SECOND CALL ON,
                       ONE HAS TO PUT THEN ICALL.GT.1.
                       IF ICALL.LT.1, THE ROUTINE WILL END WITH IER = 6.

              MAXP1  - INTEGER
                       GIVES AN UPPER BOUND ON THE NUMBER OF
                       CHEBYSHEV MOMENTS WHICH CAN BE STORED, I.E.
                       FOR THE INTERVALS OF LENGHTS ABS(B-A)*2**(-L),
                       L=0,1, ..., MAXP1-2, MAXP1.GE.1.
                       IF MAXP1.LT.1, THE ROUTINE WILL END WITH IER = 6.
                       INCREASING (DECREASING) THE VALUE OF MAXP1
                       DECREASES (INCREASES) THE COMPUTATIONAL TIME BUT
                       INCREASES (DECREASES) THE REQUIRED MEMORY SPACE.

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
                               ROUTINE. IT IS ASSUMED THAT THE
                               REQUESTED ACCURACY HAS BEEN ACHIEVED.
                     - IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE.
                               THE ESTIMATES FOR INTEGRAL AND ERROR ARE
                               LESS RELIABLE. IT IS ASSUMED THAT THE
                               REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
                       IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
                               HAS BEEN ACHIEVED. ONE CAN ALLOW MORE
                               SUBDIVISIONS BY INCREASING THE VALUE
                               OF LIMIT (AND TAKING ACCORDING
                               DIMENSION ADJUSTMENTS INTO ACCOUNT).
                               HOWEVER, IF THIS YIELDS NO IMPROVEMENT IT
                               IS ADVISED TO ANALYZE THE INTEGRAND,
                               IN ORDER TO DETERMINE THE INTEGRATION
                               DIFFICULTIES. IF THE POSITION OF A LOCAL
                               DIFFICULTY CAN BE DETERMINED (E.G.
                               SINGULARITY, DISCONTINUITY WITHIN THE
                               INTERVAL) ONE WILL PROBABLY GAIN FROM
                               SPLITTING UP THE INTERVAL AT THIS POINT
                               AND CALLING THE INTEGRATOR ON THE
                               SUBRANGES. IF POSSIBLE, AN APPROPRIATE
                               SPECIAL-PURPOSE INTEGRATOR SHOULD BE
                               USED WHICH IS DESIGNED FOR HANDLING THE
                               TYPE OF DIFFICULTY INVOLVED.
                           = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS
                               DETECTED, WHICH PREVENTS THE REQUESTED
                               TOLERANCE FROM BEING ACHIEVED.
                               THE ERROR MAY BE UNDER-ESTIMATED.
                           = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS
                               AT SOME POINTS OF THE INTEGRATION
                               INTERVAL.
                           = 4 THE ALGORITHM DOES NOT CONVERGE. ROUNDOFF
                               ERROR IS DETECTED IN THE EXTRAPOLATION
                               TABLE. IT IS PRESUMED THAT THE REQUESTED
                               TOLERANCE CANNOT BE ACHIEVED DUE TO
                               ROUNDOFF IN THE EXTRAPOLATION TABLE,
                               AND THAT THE RETURNED RESULT IS THE
                               BEST WHICH CAN BE OBTAINED.
                           = 5 THE INTEGRAL IS PROBABLY DIVERGENT, OR
                               SLOWLY CONVERGENT. IT MUST BE NOTED THAT
                               DIVERGENCE CAN OCCUR WITH ANY OTHER VALUE
                               OF IER.GT.0.
                           = 6 THE INPUT IS INVALID, BECAUSE
                               EPSABS.LT.0 AND EPSREL.LT.0,
                               OR (INTEGR.NE.1 AND INTEGR.NE.2) OR
                               ICALL.LT.1 OR MAXP1.LT.1.
                               RESULT, ABSERR, NEVAL, LAST, RLIST(1),
                               ELIST(1), IORD(1) AND NNLOG(1) ARE SET
                               TO ZERO. ALIST(1) AND BLIST(1) ARE SET
                               TO A AND B RESPECTIVELY.

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
                       VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                        LAST  ELEMENTS OF WHICH ARE THE INTEGRAL
                       APPROXIMATIONS ON THE SUBINTERVALS

              ELIST  - REAL
                       VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
                        LAST  ELEMENTS OF WHICH ARE THE MODULI OF THE
                       ABSOLUTE ERROR ESTIMATES ON THE SUBINTERVALS

              IORD   - INTEGER
                       VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST K
                       ELEMENTS OF WHICH ARE POINTERS TO THE ERROR
                       ESTIMATES OVER THE SUBINTERVALS, SUCH THAT
                       ELIST(IORD(1)), ..., ELIST(IORD(K)) FORM A
                       DECREASING SEQUENCE, WITH
                       K = LAST IF LAST.LE.(LIMIT/2+2), AND
                       K = LIMIT+1-LAST OTHERWISE.

              NNLOG  - INTEGER
                       VECTOR OF DIMENSION AT LEAST LIMIT, CONTAINING THE
                       SUBDIVISION LEVELS OF THE SUBINTERVALS, I.E.
                       IWORK(I) = L MEANS THAT THE SUBINTERVAL
                       NUMBERED I IS OF LENGTH ABS(B-A)*2**(1-L)

           ON ENTRY AND RETURN
              MOMCOM - INTEGER
                       INDICATING THAT THE CHEBYSHEV MOMENTS
                       HAVE BEEN COMPUTED FOR INTERVALS OF LENGTHS
                       (ABS(B-A))*2**(-L), L=0,1,2, ..., MOMCOM-1,
                       MOMCOM.LT.MAXP1

              CHEBMO - REAL
                       ARRAY OF DIMENSION (MAXP1,25) CONTAINING THE
                       CHEBYSHEV MOMENTS

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q5AWO
                - Q10G
                - Q4AWO
                - Q8AWO
                - Q4NG
                - F(USER PROVIDED FUNCTION)
                - FORTRAN ABS, AMAX1, AMIN1
                - Q6AWO
                - Q7AWO

   ......................................................................




              THE DIMENSION OF RLIST2 IS DETERMINED BY  THE VALUE OF
              LIMEXP IN SUBROUTINE  Q4AWO (RLIST2 SHOULD BE OF
              DIMENSION (LIMEXP+2) AT LEAST).

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
             MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
                         ERROR ESTIMATE
             ERRMAX    - ELIST(MAXERR)
             ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
             AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
             ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
             ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
                         ABS(RESULT))
             *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
             *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
             LAST      - INDEX FOR SUBDIVISION
             NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
             NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN APPROPRIATE
                         APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS
                         BEEN OBTAINED IT IS PUT IN RLIST2(NUMRL2) AFTER
                         NUMRL2 HAS BEEN INCREASED BY ONE
             SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED
                         UP TO NOW, MULTIPLIED BY 1.5
             ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
                         THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
             EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS
                         ATTEMPTING TO PERFORM EXTRAPOLATION, I.E. BEFORE
                         SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO
                         DECREASE THE VALUE OF ERLARG
             NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
                         IS NO LONGER ALLOWED (TRUE  VALUE)

              MACHINE DEPENDENT CONSTANTS
              ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST POSITIVE MAGNITUDE.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_q3awo (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mfloat *omega, Mint *integr, Mfloat *epsabs,
                Mfloat *epsrel, Mint *limit, Mint *icall,
                Mint *maxp1, Mfloat *result, Mfloat *abserr,
                Mint *neval, Mint *ier, Mfloat alist[],
                Mfloat blist[], Mfloat rlist[], Mfloat elist[],
                Mint iord[], Mint nnlog[], Mint *last, Mint *momcom,
                Mfloat *chebmo)
#else
void        imsl_q3awo (f, a, b, omega, integr, epsabs, epsrel, limit,
                icall, maxp1, result, abserr, neval, ier, alist, blist, rlist,
                elist, iord, nnlog, last, momcom, chebmo)
    Mfloat      (*f) (), *a, *b, *omega;
    Mint       *integr;
    Mfloat     *epsabs, *epsrel;
    Mint       *limit, *icall, *maxp1;
    Mfloat     *result, *abserr;
    Mint       *neval, *ier;
    Mfloat      alist[], blist[], rlist[], elist[];
    Mint        iord[], nnlog[], *last, *momcom;
    Mfloat     *chebmo;
#endif
{
#define CHEBMO(I_,J_)	(chebmo+(I_)*(*maxp1)+(J_))
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
#ifdef COMPUTER_DECOSF
    int         extall, extrap, noext;
#else
    long        extall, extrap, noext;
#endif
    Mint        _l0, id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn,
                ktmin, maxerr, nev, nres, nrmax, nrmom, numrl2;
    Mfloat      a1, a2, abseps, area, area1, area12, area2, b1, b2, correc,
                defab1, defab2, defabs, domega, dres, epmach, erlarg, erlast,
                errbnd, errmax, erro12, error1, error2, errsum, ertest,
                oflow, resabs, reseps, rlist2[52], small, uflow, width;
    Mfloat      res3la[3] = {0.0, 0.0, 0.0}; /* init. (Jim M.) */

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
    iord[0] = 0;
    nnlog[0] = 0;
    if ((((*integr != 1 && *integr != 2) || (*epsabs < F_ZERO && *epsrel <
		    F_ZERO)) || *icall < 1) || *maxp1 < 1)
	*ier = 6;
    if (*ier == 6)
	goto L_240;
    /* FIRST APPROXIMATION TO THE INTEGRAL */
    domega = fabs (*omega);
    nrmom = 0;
    if (*icall > 1)
	goto L_10;
    *momcom = 0;
L_10:
    _l0 = 0;
    l_q5awo (f, a, b, &domega, integr, &nrmom, maxp1, &_l0, result,
	abserr, neval, &defabs, &resabs, momcom, chebmo);
    /* TEST ON ACCURACY. */
    dres = fabs (*result);
    errbnd = imsl_f_max (*epsabs, *epsrel * dres);
    rlist[0] = *result;
    elist[0] = *abserr;
    iord[0] = 1;
    if (*abserr <= 1.0e02 * epmach * defabs && *abserr > errbnd)
	*ier = 2;
    if (*limit == 1)
	*ier = 1;
    if (*ier != 0 || *abserr <= errbnd)
	goto L_230;
    /* INITIALIZATIONS */
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktmin = 0;
    small = fabs (*b - *a) * 7.5e-01;
    nres = 0;
    numrl2 = 0;
    extall = FALSE;
    if (F_HALF * fabs (*b - *a) * domega > F_TWO)
	goto L_20;
    numrl2 = 1;
    extall = TRUE;
    rlist2[0] = *result;
L_20:
    if (2.5e-01 * fabs (*b - *a) * domega <= F_TWO)
	extall = TRUE;
    ksgn = -1;
    if (dres >= (F_ONE - 5.0e01 * epmach) * defabs)
	ksgn = 1;

    /* MAIN DO-LOOP */
    for (*last = 2; *last <= *limit; (*last)++) {
	/*
	 * BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR ESTIMATE.
	 */
	nrmom = nnlog[maxerr - 1] + 1;
	a1 = alist[maxerr - 1];
	b1 = F_HALF * (alist[maxerr - 1] + blist[maxerr - 1]);
	a2 = b1;
	b2 = blist[maxerr - 1];
	erlast = errmax;
	_l0 = 0;
	l_q5awo (f, &a1, &b1, &domega, integr, &nrmom, maxp1, &_l0,
	    &area1, &error1, &nev, &resabs, &defab1, momcom, chebmo);
	*neval += nev;
	_l0 = 1;
	l_q5awo (f, &a2, &b2, &domega, integr, &nrmom, maxp1, &_l0,
	    &area2, &error2, &nev, &resabs, &defab2, momcom, chebmo);
	*neval += nev;
	/*
	 * IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR AND TEST FOR
	 * ACCURACY.
	 */
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum += erro12 - errmax;
	area += area12 - rlist[maxerr - 1];
	if (defab1 == error1 || defab2 == error2)
	    goto L_40;
	if (fabs (rlist[maxerr - 1] - area12) > 1.0e-05 * fabs (area12) ||
	    erro12 < 9.9e-01 * errmax)
	    goto L_30;
	if (extrap)
	    iroff2 += 1;
	if (!extrap)
	    iroff1 += 1;
L_30:
	if (*last > 10 && erro12 > errmax)
	    iroff3 += 1;
L_40:
	rlist[maxerr - 1] = area1;
	rlist[*last - 1] = area2;
	nnlog[maxerr - 1] = nrmom;
	nnlog[*last - 1] = nrmom;
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
	 * LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL WITH NRMAX-TH
	 * LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
	 */
L_60:
	imsl_q10g (limit, last, &maxerr, &errmax, elist, iord, &nrmax);

	/* JUMP OUT OF DO-LOOP */
	if (errsum <= errbnd)
	    goto L_200;
	if (*ier != 0)
	    goto L_170;
	if (*last == 2 && extall)
	    goto L_140;
	if (noext)
	    goto L_160;
	if (!extall)
	    goto L_70;
	erlarg -= erlast;
	if (fabs (b1 - a1) > small)
	    erlarg += erro12;
	if (extrap)
	    goto L_90;
	/*
	 * TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE SMALLEST
	 * INTERVAL.
	 */
L_70:
	width = fabs (blist[maxerr - 1] - alist[maxerr - 1]);
	if (width > small)
	    goto L_160;
	if (extall)
	    goto L_80;
	/*
	 * TEST WHETHER WE CAN START WITH THE EXTRAPOLATION PROCEDURE (WE DO
	 * THIS IF WE INTEGRATE OVER THE NEXT INTERVAL WITH USE OF A
	 * GAUSS-KRONROD RULE - SEE SUBROUTINE Q5AWO).
	 */
	small *= F_HALF;
	if (2.5e-01 * width * domega > F_TWO)
	    goto L_160;
	extall = TRUE;
	goto L_150;
L_80:
	extrap = TRUE;
	nrmax = 2;
L_90:
	if (ierro == 3 || erlarg <= ertest)
	    goto L_110;

	/*
	 * THE SMALLEST INTERVAL HAS THE LARGEST ERROR. BEFORE BISECTING
	 * DECREASE THE SUM OF THE ERRORS OVER THE LARGER INTERVALS (ERLARG)
	 * AND PERFORM EXTRAPOLATION.
	 */
	jupbnd = *last;
	if (*last > (*limit / 2 + 2))
	    jupbnd = *limit + 3 - *last;
	id = nrmax;
	for (k = id; k <= jupbnd; k++) {
	    maxerr = iord[nrmax - 1];
	    errmax = elist[maxerr - 1];
	    if (fabs (blist[maxerr - 1] - alist[maxerr - 1]) > small)
		goto L_160;
	    nrmax += 1;
	}
	/* PERFORM EXTRAPOLATION. */
L_110:
	numrl2 += 1;
	rlist2[numrl2 - 1] = area;
	if (numrl2 < 3)
	    goto L_130;
	imsl_q4awo (&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	ktmin += 1;
	if (ktmin > 5 && *abserr < 1.0e-03 * errsum)
	    *ier = 5;
	if (abseps >= *abserr)
	    goto L_120;
	ktmin = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
	ertest = imsl_f_max (*epsabs, *epsrel * fabs (reseps));
	/* JUMP OUT OF DO-LOOP */
	if (*abserr <= ertest)
	    goto L_170;
	/*
	 * PREPARE BISECTION OF THE SMALLEST INTERVAL.
	 */
L_120:
	if (numrl2 == 1)
	    noext = TRUE;
	if (*ier == 5)
	    goto L_170;
L_130:
	maxerr = iord[0];
	errmax = elist[maxerr - 1];
	nrmax = 1;
	extrap = FALSE;
	small *= F_HALF;
	erlarg = errsum;
	goto L_160;
L_140:
	small *= F_HALF;
	numrl2 += 1;
	rlist2[numrl2 - 1] = area;
L_150:
	ertest = errbnd;
	erlarg = errsum;
L_160:
	;
    }
    /* SET THE FINAL RESULT. */
L_170:
    if (*abserr == oflow || nres == 0)
	goto L_200;
    if (*ier + ierro == 0)
	goto L_190;
    if (ierro == 3)
	*abserr += correc;
    if (*ier == 0)
	*ier = 3;
    if (*result != F_ZERO && area != F_ZERO)
	goto L_180;
    if (*abserr > errsum)
	goto L_200;
    if (area == F_ZERO)
	goto L_220;
    goto L_190;
L_180:
    if (*abserr / fabs (*result) > errsum / fabs (area))
	goto L_200;

    /*
     * TEST ON DIVERGENCE.
     */
L_190:
    if (ksgn == (-1) && imsl_f_max (fabs (*result), fabs (area)) <= defabs *
	1.0e-02)
	goto L_220;
    if ((1.0e-02 > (*result / area) || (*result / area) > 1.0e02) || errsum >=
	fabs (area))
	*ier = 6;
    goto L_220;
    /* COMPUTE GLOBAL INTEGRAL SUM. */
L_200:
    *result = F_ZERO;
    for (k = 1; k <= *last; k++) {
	*result += rlist[k - 1];
    }
    *abserr = errsum;
L_220:
    if (*ier > 2)
	*ier -= 1;
L_230:
    if (*integr == 2 && *omega < F_ZERO)
	*result = -*result;
L_240:
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 09:56:46 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 09:56:41
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q5AWO/DQ5AWO (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q5AWO (F, A, B, OMEGA, INTEGR, NRMOM, MAXP1, KSAVE,
                            RESULT, ABSERR, NEVAL, RESABS, RESASC,
                            MOMCOM, CHEBMO)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


   ......................................................................

   1.     Q5AWO
          INTEGRATION RULES FOR FUNCTIONS WITH COS OR SIN FACTOR
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             TO COMPUTE  THE  INTEGRAL
                I = INTEGRAL OF F(X)*W(X) OVER (A,B)
                    WHERE W(X) = COS(OMEGA*X)
                       OR W(X) = SIN(OMEGA*X),
             AND TO COMPUTE J = INTEGRAL OF ABS(F) OVER (A,B).
             FOR SMALL VALUES OF OMEGA OR SMALL INTERVALS (A,B) THE
             15-POINT GAUSS-KRONROD RULE IS USED. IN ALL OTHER CASES A
             GENERALIZED CLENSHAW-CURTIS METHOD IS USED, I.E. A
             TRUNCATED CHEBYSHEV EXPANSION OF THE FUNCTION F IS COMPUTED
             ON (A,B), SO THAT THE INTEGRAND CAN BE WRITTEN AS A SUM OF
             TERMS OF THE FORM W(X)T(K,X), WHERE T(K,X) IS THE CHEBYSHEV
             POLYNOMIAL OF DEGREE K. THE CHEBYSHEV MOMENTS ARE COMPUTED
             WITH USE OF A LINEAR RECURRENCE RELATION.

   3.     CALLING SEQUENCE
             CALL Q5AWO(F,A,B,OMEGA,INTEGR,NRMOM,MAXP1,KSAVE,RESULT,
                         ABSERR,NEVAL,RESABS,RESASC,MOMCOM,CHEBMO)

          PARAMETERS
           ON ENTRY
             F      - REAL
                      FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                      FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO
                      BE DECLARED E X T E R N A L IN THE CALLING PROGRAM.

             A      - REAL
                      LOWER LIMIT OF INTEGRATION

             B      - REAL
                      UPPER LIMIT OF INTEGRATION

             OMEGA  - REAL
                      PARAMETER IN THE WEIGHT FUNCTION

             INTEGR - INTEGER
                      INDICATES WHICH WEIGHT FUNCTION IS TO BE USED
                         INTEGR = 1   W(X) = COS(OMEGA*X)
                         INTEGR = 2   W(X) = SIN(OMEGA*X)

             NRMOM  - INTEGER
                      THE LENGTH OF INTERVAL (A,B) IS EQUAL TO THE LENGTH
                      OF THE ORIGINAL INTEGRATION INTERVAL DIVIDED BY
                      2**NRMOM (WE SUPPOSE THAT THE ROUTINE IS USED IN AN
                      ADAPTIVE INTEGRATION PROCESS, OTHERWISE SET
                      NRMOM = 0). NRMOM MUST BE ZERO AT THE FIRST CALL.

             MAXP1  - INTEGER
                      GIVES AN UPPER BOUND ON THE NUMBER OF CHEBYSHEV
                      MOMENTS WHICH CAN BE STORED, I.E. FOR THE INTERVALS
                      OF LENGTHS ABS(BB-AA)*2**(-L), L = 0,1,2, ...,
                      MAXP1-2.

             KSAVE  - INTEGER
                      KEY WHICH IS ONE WHEN THE MOMENTS FOR THE
                      CURRENT INTERVAL HAVE BEEN COMPUTED

           ON RETURN
             RESULT - REAL
                      APPROXIMATION TO THE INTEGRAL I

             ABSERR - REAL
                      ESTIMATE OF THE MODULUS OF THE ABSOLUTE
                      ERROR, WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

             NEVAL  - INTEGER
                      NUMBER OF INTEGRAND EVALUATIONS

             RESABS - REAL
                      APPROXIMATION TO THE INTEGRAL J

             RESASC - REAL
                      APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))

           ON ENTRY AND RETURN
             MOMCOM - INTEGER
                      FOR EACH INTERVAL LENGTH WE NEED TO COMPUTE
                      THE CHEBYSHEV MOMENTS. MOMCOM COUNTS THE NUMBER
                      OF INTERVALS FOR WHICH THESE MOMENTS HAVE ALREADY
                      BEEN COMPUTED. IF NRMOM.LT.MOMCOM OR KSAVE = 1,
                      THE CHEBYSHEV MOMENTS FOR THE INTERVAL (A,B)
                      HAVE ALREADY BEEN COMPUTED AND STORED, OTHERWISE
                      WE COMPUTE THEM AND WE INCREASE MOMCOM.

             CHEBMO - REAL
                      ARRAY OF DIMENSION AT LEAST (MAXP1,25) CONTAINING
                      THE MODIFIED CHEBYSHEV MOMENTS FOR THE FIRST MOMCOM
                      INTERVAL LENGTHS

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                - Q8AWO
                - Q7AWO
                - Q4NG
                - F(USER-PROVIDED FUNCTION)
                - Q6AWO
                - FORTRAN ABS, COS, AMAX1, AMIN1, SIN

   ......................................................................




             THE DATA VALUE OF MAXP1 GIVES AN UPPER BOUND
             ON THE NUMBER OF CHEBYSHEV MOMENTS WHICH CAN BE
             COMPUTED, I.E. FOR THE INTERVAL (BB-AA), ...,
             (BB-AA)/2**(MAXP1-2).
             SHOULD THIS NUMBER BE ALTERED, THE FIRST DIMENSION OF
             CHEBMO NEEDS TO BE ADAPTED.


             THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
             K = 1, ...,11, TO BE USED FOR THE CHEBYSHEV EXPANSION OF F



             LIST OF MAJOR VARIABLES
             ----------------------
             CENTR  - MID POINT OF THE INTEGRATION INTERVAL
             HLGTH  - HALF LENGTH OF THE INTEGRATION INTERVAL
             FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
                      (B-A)*0.5*COS(K*PI/12) + (B+A)*0.5
                      K = 0, ...,24
             CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
                      OF DEGREE 12, FOR THE FUNCTION F, IN THE
                      INTERVAL (A,B)
             CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
                      OF DEGREE 24, FOR THE FUNCTION F, IN THE
                      INTERVAL (A,B)
             RESC12 - APPROXIMATION TO THE INTEGRAL OF
                      COS(0.5*(B-A)*OMEGA*X)*F(0.5*(B-A)*X+0.5*(B+A))
                      OVER (-1,+1), USING THE CHEBYSHEV SERIES
                      EXPANSION OF DEGREE 12
             RESC24 - APPROXIMATION TO THE SAME INTEGRAL, USING THE
                      CHEBYSHEV SERIES EXPANSION OF DEGREE 24
             RESS12 - THE ANALOGUE OF RESC12 FOR THE SINE
             RESS24 - THE ANALOGUE OF RESC24 FOR THE SINE


             MACHINE DEPENDENT CONSTANT
             --------------------------

             OFLOW IS THE LARGEST POSITIVE MAGNITUDE.


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
static void l_q5awo (Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b,
                Mfloat *omega, Mint *integr, Mint *nrmom,
                Mint *maxp1, Mint *ksave, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mfloat *resabs,
                Mfloat *resasc, Mint *momcom, Mfloat *chebmo)
#else
static void l_q5awo (f, a, b, omega, integr, nrmom, maxp1, ksave,
                result, abserr, neval, resabs, resasc, momcom, chebmo)
    Mfloat      (*f) (), *a, *b, *omega;
    Mint       *integr, *nrmom, *maxp1, *ksave;
    Mfloat     *result, *abserr;
    Mint       *neval;
    Mfloat     *resabs, *resasc;
    Mint       *momcom;
    Mfloat     *chebmo;
#endif
{
#define CHEBMO(I_,J_)	(chebmo+(I_)*(*maxp1)+(J_))
    Mint        i, isym, j, k, m, noeq1, noequ;
    Mfloat      ac, an, an2, as, asap, ass, centr, cheb12[13], cheb24[25],
                conc, cons, cospar, d[28], d1[28], d2[28], d3[28], epmach,
                estc, ests, fval[25], hlgth, oflow, p2, p3, p4, par2, par22,
                parint, resc12, resc24, ress12, ress24, sinpar, uflow, v[28];
    Mfloat      temp_float;
    Mfloat      temp2_float;
    static Mint nmac = 28;
    imsl_q4ng (&epmach, &uflow, &oflow);

    centr = F_HALF * (*b + *a);
    hlgth = F_HALF * (*b - *a);
    parint = *omega * hlgth;
    /*
     * COMPUTE THE INTEGRAL USING THE 15-POINT GAUSS-KRONROD FORMULA IF THE
     * VALUE OF THE PARAMETER IN THE INTEGRAND IS SMALL OR IF THE LENGTH OF
     * THE INTEGRATION INTERVAL IS LESS THAN (BB-AA)/2**(MAXP1-2), WHERE
     * (AA,BB) IS THE ORIGINAL INTEGRATION INTERVAL
     */
    if (fabs (parint) > F_TWO)
	goto L_10;
    imsl_q8awo (f, l_q6awo, omega, &p2, &p3, &p4, integr, a, b, result, abserr,
	resabs, resasc);
    *neval = 15;
    goto L_190;
    /*
     * COMPUTE THE INTEGRAL USING THE GENERALIZED CLENSHAW- CURTIS METHOD
     */
L_10:
    conc = hlgth * cos (centr ** omega);
    cons = hlgth * sin (centr ** omega);
    *resasc = oflow;
    *neval = 25;
    /*
     * CHECK WHETHER THE CHEBYSHEV MOMENTS FOR THIS INTERVAL HAVE ALREADY
     * BEEN COMPUTED
     */
    if (*nrmom < *momcom || *ksave == 1)
	goto L_140;
    /*
     * COMPUTE A NEW SET OF CHEBYSHEV MOMENTS
     */
    m = *momcom + 1;
    par2 = parint * parint;
    par22 = par2 + F_TWO;
    sinpar = sin (parint);
    cospar = cos (parint);
    /*
     * COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO COSINE
     */
    v[0] = F_TWO * sinpar / parint;
    v[1] = (F_EIGHT * cospar + (par2 + par2 - F_EIGHT) * sinpar / parint) /
	par2;
    v[2] = (3.2e01 * (par2 - 1.2e01) * cospar + (F_TWO * ((par2 - 8.0e01) *
		par2 + 1.92e02) * sinpar) / parint) / (par2 * par2);
    ac = F_EIGHT * cospar;
    as = 2.4e01 * parint * sinpar;
    if (fabs (parint) > 2.4e01)
	goto L_70;
    /*
     * COMPUTE THE CHEBYSHEV MOMENTS AS THE SOLUTIONS OF A BOUNDARY VALUE
     * PROBLEM WITH 1 INITIAL VALUE (V(3)) AND 1 END VALUE (COMPUTED USING AN
     * ASYMPTOTIC FORMULA)
     */
    noequ = nmac - 3;
    noeq1 = noequ - 1;
    an = F_SIX;
    for (k = 1; k <= noeq1; k++) {
	an2 = an * an;
	d[k - 1] = -F_TWO * (an2 - F_FOUR) * (par22 - an2 - an2);
	d2[k - 1] = (an - F_ONE) * (an - F_TWO) * par2;
	d1[k - 1] = (an + F_THREE) * (an + F_FOUR) * par2;
	v[k + 2] = as - (an2 - F_FOUR) * ac;
	an += F_TWO;
    }
    an2 = an * an;
    d[noequ - 1] = -F_TWO * (an2 - F_FOUR) * (par22 - an2 - an2);
    v[noequ + 2] = as - (an2 - F_FOUR) * ac;
    v[3] += -5.6e01 * par2 * v[2];
    ass = parint * sinpar;
    asap = (((((2.1e02 * par2 - F_ONE) * cospar - (1.05e02 * par2 - 6.3e01) *
		    ass) / an2 - (F_ONE - 1.5e01 * par2) * cospar + 1.5e01 * ass) / an2 -
	    cospar + F_THREE * ass) / an2 - cospar) / an2;
    v[noequ + 2] += -F_TWO * asap * par2 * (an - F_ONE) * (an - F_TWO);
    /*
     * SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN ELIMINATION WITH
     * PARTIAL PIVOTING
     */
    for (i = 1; i <= noequ; i++) {
	d3[i - 1] = F_ZERO;
    }
    d2[noequ - 1] = F_ZERO;
    for (i = 1; i <= noeq1; i++) {
	if (fabs (d1[i - 1]) <= fabs (d[i - 1]))
	    goto L_40;
	an = d1[i - 1];
	d1[i - 1] = d[i - 1];
	d[i - 1] = an;
	an = d2[i - 1];
	d2[i - 1] = d[i];
	d[i] = an;
	d3[i - 1] = d2[i];
	d2[i] = F_ZERO;
	an = v[i + 3];
	v[i + 3] = v[i + 2];
	v[i + 2] = an;
L_40:
	d[i] += -d2[i - 1] * d1[i - 1] / d[i - 1];
	d2[i] += -d3[i - 1] * d1[i - 1] / d[i - 1];
	v[i + 3] += -v[i + 2] * d1[i - 1] / d[i - 1];
    }
    v[noequ + 2] /= d[noequ - 1];
    v[noequ + 1] = (v[noequ + 1] - d2[noeq1 - 1] * v[noequ + 2]) / d[noeq1 - 1];
    for (i = 2; i <= noeq1; i++) {
	k = noequ - i;
	v[k + 2] = (v[k + 2] - d3[k - 1] * v[k + 4] - d2[k - 1] * v[k + 3]) /
	    d[k - 1];
    }
    goto L_90;
    /*
     * COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD RECURSION
     */
L_70:
    an = F_FOUR;
    for (i = 4; i <= 13; i++) {
	an2 = an * an;
	v[i - 1] = ((an2 - F_FOUR) * (F_TWO * (par22 - an2 - an2) * v[i - 2] -
		ac) + as - par2 * (an + F_ONE) * (an + F_TWO) * v[i - 3]) / (par2 *
	    (an - F_ONE) * (an - F_TWO));
	an += F_TWO;
    }
L_90:
    for (j = 1; j <= 13; j++) {
	*CHEBMO (j * 2 - 2, m - 1) = v[j - 1];
    }
    /*
     * COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO SINE
     */
    v[0] = F_TWO * (sinpar - parint * cospar) / par2;
    v[1] = (1.8e01 - 4.8e01 / par2) * sinpar / par2 + (-F_TWO + 4.8e01 /
	par2) * cospar / parint;
    ac = -2.4e01 * parint * cospar;
    as = -F_EIGHT * sinpar;
    *CHEBMO (1, m - 1) = v[0];
    *CHEBMO (3, m - 1) = v[1];
    if (fabs (parint) > 2.4e01)
	goto L_120;
    for (k = 3; k <= 12; k++) {
	an = k;
	*CHEBMO (k * 2 - 1, m - 1) = -sinpar / (an * (F_TWO * an - F_TWO)) -
	    2.5e-01 * parint * (v[k] / an - v[k - 1] / (an - F_ONE));
    }
    goto L_140;
    /*
     * COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD RECURSION
     */
L_120:
    an = F_THREE;
    for (i = 3; i <= 12; i++) {
	an2 = an * an;
	v[i - 1] = ((an2 - F_FOUR) * (F_TWO * (par22 - an2 - an2) * v[i - 2] +
		as) + ac - par2 * (an + F_ONE) * (an + F_TWO) * v[i - 3]) / (par2 *
	    (an - F_ONE) * (an - F_TWO));
	an += F_TWO;
	*CHEBMO (i * 2 - 1, m - 1) = v[i - 1];
    }
L_140:
    if (*nrmom < *momcom)
	m = *nrmom + 1;
    if (*momcom < (*maxp1 - 1) && *nrmom >= *momcom)
	*momcom += 1;

    /*
     * COMPUTE THE COEFFICIENTS OF THE CHEBYSHEV EXPANSIONS OF DEGREES 12 AND
     * 24 OF THE FUNCTION F
     */
    imsl_e1usr ("ON");
    temp_float = centr + hlgth;
    fval[0] = F_HALF * (*f) (temp_float);
/*	fval[0] = 5.0e-01 * (*f) (ADR(_f0, centr + hlgth)); */
    fval[12] = (*f) (centr);
/*	fval[12] = (*f) (&centr); */
    temp_float = centr - hlgth;
    fval[24] = F_HALF * (*f) (temp_float);
/*	fval[24] = 5.0e-01 * (*f) (ADR(_f0, centr - hlgth)); */
    for (i = 2; i <= 12; i++) {
	isym = 26 - i;
	temp_float = hlgth * lv_x[i - 2] + centr;
	fval[i - 1] = (*f) (temp_float);
/*		fval[i - 1] = (*f) (ADR(_f0, hlgth * lv_x[i - 2] + centr)); */
	temp2_float = centr - hlgth * lv_x[i - 2];
	fval[isym - 1] = (*f) (temp2_float);
/*		fval[isym - 1] = (*f) (ADR(_f0, centr - hlgth * lv_x[i - 2])); */
    }
    imsl_e1usr ("OFF");
    imsl_q7awo (lv_x, fval, cheb12, cheb24);
    /*
     * COMPUTE THE INTEGRAL AND ERROR ESTIMATES
     */
    resc12 = cheb12[12] ** CHEBMO (12, m - 1);
    ress12 = F_ZERO;
    estc = fabs (cheb24[24] ** CHEBMO (24, m - 1)) + fabs ((cheb12[12] -
	    cheb24[12]) ** CHEBMO (12, m - 1));
    ests = F_ZERO;
    k = 11;
    for (j = 1; j <= 6; j++) {
	resc12 += cheb12[k - 1] ** CHEBMO (k - 1, m - 1);
	ress12 += cheb12[k] ** CHEBMO (k, m - 1);
	estc += fabs ((cheb12[k - 1] - cheb24[k - 1]) ** CHEBMO (k - 1, m - 1));
	ests += fabs ((cheb12[k] - cheb24[k]) ** CHEBMO (k, m - 1));
	k -= 2;
    }
    resc24 = cheb24[24] ** CHEBMO (24, m - 1);
    ress24 = F_ZERO;
    *resabs = fabs (cheb24[24]);
    k = 23;
    for (j = 1; j <= 12; j++) {
	resc24 += cheb24[k - 1] ** CHEBMO (k - 1, m - 1);
	ress24 += cheb24[k] ** CHEBMO (k, m - 1);
	*resabs += fabs (cheb24[k - 1]) + fabs (cheb24[k]);
	if (j <= 5)
	    estc += fabs (cheb24[k - 1] ** CHEBMO (k - 1, m - 1));
	if (j <= 5)
	    ests += fabs (cheb24[k] ** CHEBMO (k, m - 1));
	k -= 2;
    }
    *resabs *= fabs (hlgth);
    if (*integr == 2)
	goto L_180;
    *result = conc * resc24 - cons * ress24;
    *abserr = fabs (conc * estc) + fabs (cons * ests);
    goto L_190;
L_180:
    *result = conc * ress24 + cons * resc24;
    *abserr = fabs (conc * ests) + fabs (cons * estc);
L_190:
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 09:59:52 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 09:59:50
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q6AWO/DQ6AWO (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      Q6AWO(X, OMEGA, P2, P3, P4, INTEGR)

    Arguments:
       X      - Scalar.  (Input)
       OMEGA  - Scalar.  (Input)
       P2     - Scalar.  (Input)
       P3     - Scalar.  (Input)
       P4     - Scalar.  (Input)
       INTEGR - Scalar.  (Input)
       Q6AWO  - Function that takes either the SIN or COS of OMEGA*X.

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* P2, P3, and P4 are not used , but leave the calling sequence intact. */
#ifdef ANSI
static Mfloat l_q6awo (Mfloat *x, Mfloat *omega, Mfloat *p2, Mfloat *p3,
                Mfloat *p4, Mint *integr)
#else
static Mfloat l_q6awo (x, omega, p2, p3, p4, integr)
    Mfloat     *x, *omega, *p2, *p3, *p4;
    Mint       *integr;
#endif
{
    Mfloat      omx, q6awo_v;


    omx = *omega ** x;
    if (*integr == 2)
	goto L_20;
    q6awo_v = cos (omx);
    goto L_30;
L_20:
    q6awo_v = sin (omx);
L_30:
    return (q6awo_v);
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 09:57:51 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 09:57:48
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q7AWO/DQ7AWO (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q7AWO (X, FVAL, CHEB12, CHEB24)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


  .......................................................................

   1.     Q7AWO
          CHEBYSHEV SERIES EXPANSION
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THIS ROUTINE COMPUTES THE CHEBYSHEV SERIES EXPANSION
             OF DEGREES 12 AND 24 OF A FUNCTION USING A
             FAST FOURIER TRANSFORM METHOD
             F(X) = SUM(K=1, ...,13) (CHEB12(K)*T(K-1,X)),
             F(X) = SUM(K=1, ...,25) (CHEB24(K)*T(K-1,X)),
             WHERE T(K,X) IS THE CHEBYSHEV POLYNOMIAL OF DEGREE K.

   3.     CALLING SEQUENCE
             CALL Q7AWO(X,FVAL,CHEB12,CHEB24)

          PARAMETERS
            ON ENTRY
             X      - REAL
                      VECTOR OF DIMENSION 11 CONTAINING THE VALUES
                      COS(K*PI/24), K = 1, ..., 11

             FVAL   - REAL
                      VECTOR OF DIMENSION 25 CONTAINING THE
                      FUNCTION VALUES AT THE POINTS
                      (B+A+(B-A)*COS(K*PI/24))/2, K = 0, ...,24,
                      WHERE (A,B) IS THE APPROXIMATION INTERVAL.
                      FVAL(1) AND FVAL(25) ARE DIVIDED BY TWO
                      (THESE VALUES ARE DESTROYED AT OUTPUT).

            ON RETURN
             CHEB12 - REAL
                      VECTOR OF DIMENSION 13 CONTAINING THE
                      CHEBYSHEV COEFFICIENTS FOR DEGREE 12

             CHEB24 - REAL
                      VECTOR OF DIMENSION 25 CONTAINING THE
                      CHEBYSHEV COEFFICIENTS FOR DEGREE 24

   4.     NO SUBROUTINES OR FUNCTIONS NEEDED

  .......................................................................



  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_q7awo (Mfloat x[], Mfloat fval[], Mfloat cheb12[],
                Mfloat cheb24[])
#else
void        imsl_q7awo (x, fval, cheb12, cheb24)
    Mfloat      x[], fval[], cheb12[], cheb24[];
#endif
{
    Mint        i, j;
    Mfloat      alam, alam1, alam2, part1, part2, part3, v[12];


    for (i = 1; i <= 12; i++) {
	j = 26 - i;
	v[i - 1] = fval[i - 1] - fval[j - 1];
	fval[i - 1] += fval[j - 1];
    }
    alam1 = v[0] - v[8];
    alam2 = x[5] * (v[2] - v[6] - v[10]);
    cheb12[3] = alam1 + alam2;
    cheb12[9] = alam1 - alam2;
    alam1 = v[1] - v[7] - v[9];
    alam2 = v[3] - v[5] - v[11];
    alam = x[2] * alam1 + x[8] * alam2;
    cheb24[3] = cheb12[3] + alam;
    cheb24[21] = cheb12[3] - alam;
    alam = x[8] * alam1 - x[2] * alam2;
    cheb24[9] = cheb12[9] + alam;
    cheb24[15] = cheb12[9] - alam;
    part1 = x[3] * v[4];
    part2 = x[7] * v[8];
    part3 = x[5] * v[6];
    alam1 = v[0] + part1 + part2;
    alam2 = x[1] * v[2] + part3 + x[9] * v[10];
    cheb12[1] = alam1 + alam2;
    cheb12[11] = alam1 - alam2;
    alam = x[0] * v[1] + x[2] * v[3] + x[4] * v[5] + x[6] * v[7] + x[8] * v[9] +
	x[10] * v[11];
    cheb24[1] = cheb12[1] + alam;
    cheb24[23] = cheb12[1] - alam;
    alam = x[10] * v[1] - x[8] * v[3] + x[6] * v[5] - x[4] * v[7] + x[2] * v[9] -
	x[0] * v[11];
    cheb24[11] = cheb12[11] + alam;
    cheb24[13] = cheb12[11] - alam;
    alam1 = v[0] - part1 + part2;
    alam2 = x[9] * v[2] - part3 + x[1] * v[10];
    cheb12[5] = alam1 + alam2;
    cheb12[7] = alam1 - alam2;
    alam = x[4] * v[1] - x[8] * v[3] - x[0] * v[5] - x[10] * v[7] + x[2] * v[9] +
	x[6] * v[11];
    cheb24[5] = cheb12[5] + alam;
    cheb24[19] = cheb12[5] - alam;
    alam = x[6] * v[1] - x[2] * v[3] - x[10] * v[5] + x[0] * v[7] - x[8] * v[9] -
	x[4] * v[11];
    cheb24[7] = cheb12[7] + alam;
    cheb24[17] = cheb12[7] - alam;
    for (i = 1; i <= 6; i++) {
	j = 14 - i;
	v[i - 1] = fval[i - 1] - fval[j - 1];
	fval[i - 1] += fval[j - 1];
    }
    alam1 = v[0] + x[7] * v[4];
    alam2 = x[3] * v[2];
    cheb12[2] = alam1 + alam2;
    cheb12[10] = alam1 - alam2;
    cheb12[6] = v[0] - v[4];
    alam = x[1] * v[1] + x[5] * v[3] + x[9] * v[5];
    cheb24[2] = cheb12[2] + alam;
    cheb24[22] = cheb12[2] - alam;
    alam = x[5] * (v[1] - v[3] - v[5]);
    cheb24[6] = cheb12[6] + alam;
    cheb24[18] = cheb12[6] - alam;
    alam = x[9] * v[1] - x[5] * v[3] + x[1] * v[5];
    cheb24[10] = cheb12[10] + alam;
    cheb24[14] = cheb12[10] - alam;
    for (i = 1; i <= 3; i++) {
	j = 8 - i;
	v[i - 1] = fval[i - 1] - fval[j - 1];
	fval[i - 1] += fval[j - 1];
    }
    cheb12[4] = v[0] + x[7] * v[2];
    cheb12[8] = fval[0] - x[7] * fval[2];
    alam = x[3] * v[1];
    cheb24[4] = cheb12[4] + alam;
    cheb24[20] = cheb12[4] - alam;
    alam = x[7] * fval[1] - fval[3];
    cheb24[8] = cheb12[8] + alam;
    cheb24[16] = cheb12[8] - alam;
    cheb12[0] = fval[0] + fval[2];
    alam = fval[1] + fval[3];
    cheb24[0] = cheb12[0] + alam;
    cheb24[24] = cheb12[0] - alam;
    cheb12[12] = v[0] - v[2];
    cheb24[12] = cheb12[12];
    alam = F_ONE / F_SIX;
    for (i = 2; i <= 12; i++) {
	cheb12[i - 1] *= alam;
    }
    alam *= F_HALF;
    cheb12[0] *= alam;
    cheb12[12] *= alam;
    for (i = 2; i <= 24; i++) {
	cheb24[i - 1] *= alam;
    }
    cheb24[0] *= F_HALF * alam;
    cheb24[24] *= F_HALF * alam;
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/10/90 at 09:58:43 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/10/90 at 09:58:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q8AWO/DQ8AWO (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function.

    Usage:      CALL Q8AWO (F, W, P1, P2, P3, P4, KP, A, B, RESULT,
                            ABSERR, RESABS, RESASC)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


       ..................................................................

   1.        Q8AWO
             INTEGRATION RULES
                STANDARD FORTRAN SUBROUTINE
                REAL VERSION

   2.        PURPOSE
                TO COMPUTE I = INTEGRAL OF F*W OVER (A,B), WITH ERROR
                               ESTIMATE
                           J = INTEGRAL OF ABS(F*W) OVER (A,B)

   3.        CALLING SEQUENCE
                CALL Q8AWO(F,W,P1,P2,P3,P4,KP,A,B,RESULT,ABSERR,
                            RESABS,RESASC)

             PARAMETERS
               ON ENTRY
                F      - REAL
                         FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                         FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
                         DECLARED E X T E R N A L IN THE CALLING PROGRAM.

                W      - REAL
                         FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
                         WEIGHT FUNCTION W(X). THE ACTUAL NAME FOR W
                         NEEDS TO BE DECLARED E X T E R N A L IN THE
                         CALLING PROGRAM.

                P1, P2, P3, P4 - REAL
                         PARAMETERS IN THE WEIGHT FUNCTION

                KP     - INTEGER
                         KEY FOR INDICATING THE TYPE OF WEIGHT FUNCTION

                A      - REAL
                         LOWER LIMIT OF INTEGRATION

                B      - REAL
                         UPPER LIMIT OF INTEGRATION

              ON RETURN
                RESULT - REAL
                         APPROXIMATION TO THE INTEGRAL I
                         RESULT IS COMPUTED BY APPLYING THE 15-POINT
                         KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
                         OF ABSCISSAE TO THE 7-POINT GAUSS RULE (RESG).

                ABSERR - REAL
                         ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                         WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

                RESABS - REAL
                         APPROXIMATION TO THE INTEGRAL OF ABS(F)

                RESASC - REAL
                         APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))

   4.        SUBROUTINES OR FUNCTIONS NEEDED
                    - F (USER-PROVIDED FUNCTION)
                    - W (Q5AWC, Q6AWO OR Q6AWS, DEPENDING ON THE
                         CALLING ROUTINE)
                    - Q4NG
                    - FORTRAN ABS, AMAX1, AMIN1


       ..................................................................



             THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
             BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
             CORRESPONDING WEIGHTS ARE GIVEN.

             XGK    - ABSCISSAE OF THE 15-POINT GAUSS-KRONROD RULE
                      XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
                      GAUSS RULE
                      XGK(1), XGK(3), ... ABSCISSAE WHICH ARE OPTIMALLY
                      ADDED TO THE 7-POINT GAUSS RULE

             WGK    - WEIGHTS OF THE 15-POINT GAUSS-KRONROD RULE

             WG     - WEIGHTS OF THE 7-POINT GAUSS RULE



             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTERVAL
             ABSC*  - ABSCISSA
             FVAL*  - FUNCTION VALUE
             RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
             RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
             RESKH  - APPROXIMATION TO THE MEAN VALUE OF F*W OVER (A,B),
                      I.E. TO I/(B-A)

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
    0.129484966168869693270611432679e0,
    0.279705391489276667901467771424e0,
    0.381830050505118944950369775489e0,
    0.417959183673469387755102040816e0
};


#ifdef ANSI
void        imsl_q8awo (Mfloat (*f) (Mfloat), Mfloat (*w) (Mfloat*, 
		Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mint*), Mfloat *p1,
                Mfloat *p2, Mfloat *p3, Mfloat *p4, Mint *kp,
                Mfloat *a, Mfloat *b, Mfloat *result, Mfloat *abserr,
                Mfloat *resabs, Mfloat *resasc)
#else
void        imsl_q8awo (f, w, p1, p2, p3, p4, kp, a, b, result, abserr,
                resabs, resasc)
    Mfloat      (*f) (), (*w) (), *p1, *p2, *p3, *p4;
    Mint       *kp;
    Mfloat     *a, *b, *result, *abserr, *resabs, *resasc;
#endif
{
    Mint        j, jtw, jtwm1;
    Mfloat      absc, absc1, absc2, centr, dhlgth, epmach, fc, fsum, fv1[7],
                fv2[7], fval1, fval2, hlgth, oflow, resg, resk, reskh,
                uflow;
    imsl_q4ng (&epmach, &uflow, &oflow);

    centr = F_HALF * (*a + *b);
    hlgth = F_HALF * (*b - *a);
    dhlgth = fabs (hlgth);
    /*
     * COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE INTEGRAL, AND
     * ESTIMATE THE ERROR
     */
    imsl_e1usr ("ON");
    fc = (*f) (centr) * (*w) (&centr, p1, p2, p3, p4, kp);
    imsl_e1usr ("OFF");
    resg = lv_wg[3] * fc;
    resk = lv_wgk[7] * fc;
    *resabs = fabs (resk);
    for (j = 1; j <= 3; j++) {
	jtw = j * 2;
	absc = hlgth * lv_xgk[jtw - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	imsl_e1usr ("ON");
	fval1 = (*f) (absc1) * (*w) (&absc1, p1, p2, p3, p4, kp);
	fval2 = (*f) (absc2) * (*w) (&absc2, p1, p2, p3, p4, kp);
	imsl_e1usr ("OFF");
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += lv_wg[j - 1] * fsum;
	resk += lv_wgk[jtw - 1] * fsum;
	*resabs += lv_wgk[jtw - 1] * (fabs (fval1) + fabs (fval2));
    }
    for (j = 1; j <= 4; j++) {
	jtwm1 = j * 2 - 1;
	absc = hlgth * lv_xgk[jtwm1 - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	imsl_e1usr ("ON");
	fval1 = (*f) (absc1) * (*w) (&absc1, p1, p2, p3, p4, kp);
	fval2 = (*f) (absc2) * (*w) (&absc2, p1, p2, p3, p4, kp);
	imsl_e1usr ("OFF");
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += lv_wgk[jtwm1 - 1] * fsum;
	*resabs += lv_wgk[jtwm1 - 1] * (fabs (fval1) + fabs (fval2));
    }
    reskh = resk * F_HALF;
    *resasc = lv_wgk[7] * fabs (fc - reskh);
    for (j = 1; j <= 7; j++) {
	*resasc += lv_wgk[j - 1] * (fabs (fv1[j - 1] - reskh) + fabs (fv2[j - 1] -
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
