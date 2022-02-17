#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK l_int_fcn_smooth (Mfloat (*fcn) (Mfloat), Mfloat a,
                Mfloat b, va_list argptr);
static void l_qdng (Mfloat (*f) (Mfloat), Mfloat a, Mfloat b,
                Mfloat errabs, Mfloat errrel, Mfloat *result,
                Mfloat *errest);
static void l_q3ng (Mfloat (*f) (Mfloat), Mfloat a, Mfloat b,
                Mfloat epsabs, Mfloat epsrel, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mint *ier);
#else
static VA_LIST_HACK l_int_fcn_smooth ();
static void l_qdng ();
static void l_q3ng ();
#endif

static Mfloat lv_value;
#ifdef ANSI
Mfloat      imsl_f_int_fcn_smooth (Mfloat (*fcn) (Mfloat), Mfloat a,
                Mfloat b,...)
#else
Mfloat      imsl_f_int_fcn_smooth (fcn, a, b, va_alist)
    Mfloat      (*fcn) ();
    Mfloat      a;
    Mfloat      b;
va_dcl
#endif
{
    va_list argptr;
    VA_START (argptr, b);
    E1PSH ("imsl_f_int_fcn_smooth", "imsl_d_int_fcn_smooth");
    lv_value = 0.0;
    IMSL_CALL (l_int_fcn_smooth (fcn, a, b, argptr));
    va_end (argptr);
    E1POP ("imsl_f_int_fcn_smooth", "imsl_d_int_fcn_smooth");
    return lv_value;
}



#ifdef ANSI
static VA_LIST_HACK l_int_fcn_smooth (Mfloat (*fcn) (Mfloat), Mfloat a,
                Mfloat b, va_list argptr)
#else
static VA_LIST_HACK l_int_fcn_smooth (fcn, a, b, argptr)
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
    Mfloat      temp_err_est;


    err_abs = sqrt ( (double)imsl_amach (4));
    err_rel = sqrt ( (double)imsl_amach (4));


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
	case IMSL_ERR_REL:
	    arg_number++;
	    err_rel = (Mfloat) va_arg (argptr, Mdouble);
	    break;
	case IMSL_ERR_ABS:
	    arg_number++;
	    err_abs = (Mfloat) va_arg (argptr, Mdouble);
	    break;
	case IMSL_ERR_EST:
	    arg_number++;
	    err_est = va_arg (argptr, Mfloat *);
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

    if (imsl_n1rty (0))
	goto RETURN;
    if (*fcn == NULL) {
	imsl_e1stl (1, "fcn");
	/* (5, 1, "The required argument %(L1) is NULL."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
    }
    if (imsl_n1rty (0))
	goto RETURN;


    if (err_est == NULL)
	err_est = &temp_err_est;

    a_float = a;
    b_float = b;

    l_qdng (fcn, a_float, b_float, err_abs, err_rel, &lv_value,
	err_est);
RETURN:
    ;
    if (imsl_n1rty (0) > 3)
	lv_value = imsl_amach(6);
    return (argptr);
}











/*Translated by FOR_C++, v0.1, on 08/16/90 at 13:20:13 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 13:20:11
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  QDNG/DQDNG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a smooth function using a nonadaptive rule.

    Usage:      CALL QDNG (F, A, B, ERRABS, ERRREL, RESULT, ERREST)

    Arguments:
       F      - User-supplied FUNCTION to be integrated.  The form is
                F(X), where
                X      - Independent variable.  (Input)
                F      - The function value.  (Output)
                F must be declared EXTERNAL in the calling program.
       A      - Lower limit of integration.  (Input)
       B      - Upper limit of integration.  (Input)
       ERRABS - Absolute accuracy desired.  (Input)
       ERRREL - Relative accuracy desired.  (Input)
       RESULT - Estimate of the integral from A to B of F.  (Output)
       ERREST - Estimate of the absolute value of the error.  (Output)

    Remarks:
    1. Informational error
       Type Code
         4   1  The maximum number of steps allowed have been taken.
                The integral is too difficult for QDNG.

    2. If EXACT is the exact value, QDNG attempts to find RESULT
       such that ABS(EXACT-RESULT) .LE. MAX(ERRABS,ERRREL*ABS(EXACT)).
       To specify only a relative error, set ERRABS to zero.  Similarly,
       to specify only an absolute error, set ERRREL to zero.

    3. This routine is designed for efficiency, not robustness.  If
       the above error is encountered, try QDAGS.

    Keywords:   Univariate quadrature; Numerical integration

    GAMS:       H2a1a1

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_qdng (Mfloat (*f) (Mfloat), Mfloat a, Mfloat b,
                Mfloat errabs, Mfloat errrel, Mfloat *result,
                Mfloat *errest)
#else
static void l_qdng (f, a, b, errabs, errrel, result, errest)
    Mfloat      (*f) (), a, b, errabs, errrel, *result, *errest;
#endif
{
    Mint        ier, neval;


    imsl_e1psh ("QDNG  ");
    /* CHECK ERRABS */
    if (errabs < 0.0e0) {
	imsl_e1str (1, errabs);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_ABS_SMALL);
    }
    /* CHECK ERRREL */
    if (errrel < 0.0e0) {
	imsl_e1str (1, errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_SMALL);
    }
    /*
     * CHECK ERRABS AND ERRREL
     */
    if (errabs == 0.0e0 && errrel == 0.0e0) {
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_TOL_ZERO);
    }
    /* CHECK ERRREL .GE. 1 */
    if (errrel >= 1.0) {
	imsl_e1str (1, errrel);
	imsl_ermes (IMSL_TERMINAL, IMSL_ERR_REL_BIG);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;

    l_q3ng (f, a, b, errabs, errrel, result, errest, &neval, &ier);

    if (ier == 1) {

	/*
	 * (4, 1, "The maximum number of steps allowed has been executed.
	 * The integral is too difficult for imsl_f_int_fcn_smooth; try
	 * imsl_f_int_fcn_sing.");
	 */
	imsl_ermes (IMSL_FATAL, IMSL_MAX_STEPS_ALLOWED);
    }
L_9000:
    ;
    imsl_e1pop ("QDNG  ");
    return;
}				/* end of function */



























/*Translated by FOR_C++, v0.1, on 08/16/90 at 13:21:34 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 13:21:30
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  Q3NG/DQ3NG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Integrate a function using a non-adaptive rule.

    Usage:      CALL Q3NG (F, A, B, EPSABS, EPSREL, RESULT, ABSERR,
                           NEVAL, IER)

    Arguments:  (See comment block below)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.


  .......................................................................

   1.     Q3NG
          NON-ADAPTIVE INTEGRATION
             STANDARD FORTRAN SUBROUTINE
             REAL VERSION

   2.     PURPOSE
             THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO
             A GIVEN DEFINITE INTEGRAL  I = INTEGRAL OF  F  OVER (A,B),
             HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).

   3.     CALLING SEQUENCE
             CALL Q3NG(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER)

          PARAMETERS
           ON ENTRY
             F      - REAL
                      FUNCTION SUBPROGRAM DEFINING THE INTEGRAND FUNCTION
                      F(X). THE ACTUAL NAME FOR F NEEDS TO BE DECLARED
                      E X T E R N A L IN THE DRIVER PROGRAM.

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
                      APPROXIMATION TO THE INTEGRAL I
                      RESULT IS OBTAINED BY APPLYING THE 21-POINT
                      GAUSS-KRONROD RULE (RES21) OBTAINED BY OPTIMAL
                      ADDITION OF ABSCISSAE TO THE 10-POINT GAUSS RULE
                      (RES10), OR BY APPLYING THE 43-POINT RULE (RES43)
                      OBTAINED BY OPTIMAL ADDITION OF ABSCISSAE TO THE
                      21-POINT GAUSS-KRONROD RULE, OR BY APPLYING THE
                      87-POINT RULE (RES87) OBTAINED BY OPTIMAL ADDITION
                      OF ABSCISSAE TO THE 43-POINT RULE.

             ABSERR - REAL
                      ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
                      WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)

             NEVAL  - INTEGER
                      NUMBER OF INTEGRAND EVALUATIONS

             IER    - IER = 0 NORMAL AND RELIABLE TERMINATION OF THE
                              ROUTINE. IT IS ASSUMED THAT THE REQUESTED
                              ACCURACY HAS BEEN ACHIEVED.
                      IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE. IT IS
                              ASSUMED THAT THE REQUESTED ACCURACY HAS
                              NOT BEEN ACHIEVED.
                      IER = 1 THE MAXIMUM NUMBER OF STEPS HAS BEEN
                              EXECUTED. THE INTEGRAL IS PROBABLY TOO
                              DIFFICULT TO BE CALCULATED BY Q3NG.
                          = 6 THE INPUT IS INVALID, BECAUSE
                              EPSABS.LT.0 AND EPSREL.LT.0,
                              RESULT, ABSERR AND NEVAL ARE SET TO ZERO.

   4.     SUBROUTINES OR FUNCTIONS NEEDED
                 - F (USER-PROVIDED FUNCTION)
                 - Q4NG
                 - FORTRAN ABS, AMAX1, AMIN1

  .......................................................................



             THE FOLLOWING DATA STATEMENTS CONTAIN THE
             ABSCISSAE AND WEIGHTS OF THE INTEGRATION RULES USED.

             X1      ABSCISSAE COMMON TO THE 10-, 21-, 43-
                     AND 87-POINT RULE
             X2      ABSCISSAE COMMON TO THE 21-, 43- AND 87-POINT RULE
             X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT RULE
             X4      ABSCISSAE OF THE 87-POINT RULE
             W10     WEIGHTS OF THE 10-POINT FORMULA
             W21A    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X1
             W21B    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X2
             W43A    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X1, X3
             W43B    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X3
             W87A    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X1,
                     X2, X3
             W87B    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X4


             LIST OF MAJOR VARIABLES
             -----------------------

             CENTR  - MID POINT OF THE INTEGRATION INTERVAL
             HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
             FCENTR - FUNCTION VALUE AT MID POINT
             ABSC   - ABSCISSA
             FVAL   - FUNCTION VALUE
             SAVFUN - ARRAY OF FUNCTION VALUES WHICH
                      HAVE ALREADY BEEN COMPUTED
             RES10  - 10-POINT GAUSS RESULT
             RES21  - 21-POINT KRONROD RESULT
             RES43  - 43-POINT RESULT
             RES87  - 87-POINT RESULT
             RESABS - APPROXIMATION TO THE INTEGRAL OF ABS(F)
             RESASC - APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))

             MACHINE DEPENDENT CONSTANTS
             ---------------------------

             EPMACH IS THE LARGEST RELATIVE SPACING.
             UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
             OFLOW IS THE LARGEST MAGNITUDE.

  -----------------------------------------------------------------------
 */

static Mfloat lv_x1[] = {
    0.973906528517171720077964012085e0,
    0.865063366688984510732096688424e0,
    0.679409568299024406234327365115e0,
    0.433395394129247190799265943166e0,
    0.14887433898163121088482600113e0
};

static Mfloat lv_w10[] = {
    0.0666713443086881375935688098933e0,
    0.149451349150580593145776339658e0,
    0.219086362515982043995534934228e0,
    0.26926671930999635509122692157e0,
    0.295524224714752870173892994651e0
};

static Mfloat lv_x2[] = {
    0.995657163025808080735527280689e0,
    0.93015749135570822600120718006e0,
    0.780817726586416897063717578345e0,
    0.562757134668604683339000099273e0,
    0.294392862701460198131126603104e0
};

static Mfloat lv_w21a[] = {
    0.0325581623079647274788189724594e0,
    0.0750396748109199527670431409162e0,
    0.109387158802297641899210590326e0,
    0.134709217311473325928054001772e0,
    0.147739104901338491374841515972e0
};

static Mfloat lv_w21b[] = {
    0.0116946388673718742780643960622e0,
    0.0547558965743519960313813002446e0,
    0.0931254545836976055350654650834e0,
    0.123491976262065851077958109831e0,
    0.142775938577060080797094273139e0,
    0.14944555400291690566493646839e0
};

static Mfloat lv_x3[] = {
    0.99933336090193208139409932392e0,
    0.987433402908088869795961478381e0,
    0.954807934814266299257919200291e0,
    0.900148695748328293625099494069e0,
    0.825198314983114150847066732589e0,
    0.732148388989304982612354848756e0,
    0.622847970537725238641159120344e0,
    0.4994795740710564999522148855e0,
    0.364901661346580768043989548503e0,
    0.222254919776601296498260928066e0,
    0.0746506174613833220439144357965e0
};

static Mfloat lv_w43a[] = {
    0.0162967342896665649242819746177e0,
    0.0375228761208695014616137958981e0,
    0.054694902058255442147212685465e0,
    0.0673554146094780860755531663022e0,
    0.0738701996323939534321406952514e0,
    0.00576855605976979618418432790866e0,
    0.0273718905932488420812760692892e0,
    0.0465608269104288307433391544338e0,
    0.0617449952014425644962403360309e0,
    0.0713872672686933977685591144255e0
};

static Mfloat lv_w43b[] = {
    0.00184447764021241410038910655297e0,
    0.0107986895858916517404654067413e0,
    0.0218953638677954281025231230752e0,
    0.0325974639753456894438822225261e0,
    0.042163137935191811847627924328e0,
    0.0507419396001845777801890200921e0,
    0.0583793955426192483754753693302e0,
    0.0647464049514458855446892595175e0,
    0.0695661979123564845286333150384e0,
    0.0728244414718332081509395351928e0,
    0.0745077510141751182735718138429e0,
    0.0747221475174030055944251682804e0
};

static Mfloat lv_x4[] = {
    0.999902977262729234490529830592e0,
    0.997989895986678745427496322366e0,
    0.992175497860687222808523352251e0,
    0.981358163572712773571916941624e0,
    0.965057623858384619128284110608e0,
    0.943167613133670596816416634507e0,
    0.91580641468550720959182643072e0,
    0.883221657771316501372117548744e0,
    0.845710748462415666605902011505e0,
    0.803557658035230982788739474981e0,
    0.757005730685495558328942793432e0,
    0.706273209787321819824094274741e0,
    0.651589466501177922534422205017e0,
    0.593223374057961088875273770349e0,
    0.531493605970831932285268948563e0,
    0.466763623042022844871966781659e0,
    0.399424847859218804732101665818e0,
    0.329874877106188288265053371825e0,
    0.258503559202161551802280975429e0,
    0.185695396568346652015917141168e0,
    0.111842213179907468172398359241e0,
    0.0373521233946198708149981654377e0
};

static Mfloat lv_w87a[] = {
    0.00814837738414917290000287844819e0,
    0.0187614382015628222439350590038e0,
    0.0273474510500522861615828297413e0,
    0.0336777073116379300465810569576e0,
    0.0369350998204279076145895867425e0,
    0.0028848724302115305013341562487e0,
    0.0136859460227127018889500352731e0,
    0.0232804135028883111234092910304e0,
    0.0308724976117133586754663941264e0,
    0.035693633639418770719351355457e0,
    0.000915283345202241360843392549948e0,
    0.00539928021930047136773874339105e0,
    0.0109476796011189311343278268568e0,
    0.0162987316967873352626657032233e0,
    0.0210815688892038351124330601882e0,
    0.0253709697692538272434679998317e0,
    0.0291896977564757525014461540849e0,
    0.0323732024672027896857881948896e0,
    0.0347830989503651427507819979496e0,
    0.0364122207313517875628011636876e0,
    0.0372538755030477085395920011912e0
};

static Mfloat lv_w87b[] = {
    0.000274145563762072350016527092881e0,
    0.00180712415505794294834131175325e0,
    0.00409686928275916486445807068348e0,
    0.00675829005184737869981657789742e0,
    0.00954995767220164653605358132538e0,
    0.0123294476522448536946266399638e0,
    0.0150104473463889523766972860419e0,
    0.0175489679862431910996653529259e0,
    0.0199380377864408882022781927307e0,
    0.0221949359610122867963321029595e0,
    0.0243391471260008054703606470415e0,
    0.0263745054148392072415037865526e0,
    0.028286910788771200659968002988e0,
    0.0300525811280926953225211103473e0,
    0.0316467513714399294045860510789e0,
    0.0330504134199785032907859448627e0,
    0.0342550997042260617870828210468e0,
    0.0352624126601566810337827179984e0,
    0.0360769896228887011855003180039e0,
    0.0366986044984560944980180474411e0,
    0.0371205492698325761141199584136e0,
    0.0373342287519350403212354490947e0,
    0.0373610737626790234103212417666e0
};


#ifdef ANSI
static void l_q3ng (Mfloat (*f) (Mfloat), Mfloat a, Mfloat b,
                Mfloat epsabs, Mfloat epsrel, Mfloat *result,
                Mfloat *abserr, Mint *neval, Mint *ier)
#else
static void l_q3ng (f, a, b, epsabs, epsrel, result, abserr, neval,
                ier)
    Mfloat      (*f) (), a, b, epsabs, epsrel, *result, *abserr;
    Mint       *neval, *ier;
#endif
{
    Mint        ipx, k, l;
    Mfloat      absc, centr, dhlgth, epmach, fcentr, fv1[5],
                fv2[5], fv3[5], fv4[5], fval, fval1, fval2, hlgth, oflow,
                res10, res21, res43, res87, resabs, resasc, reskh, savfun[21],
                uflow;
    imsl_q4ng (&epmach, &uflow, &oflow);
    /* TEST ON VALIDITY OF PARAMETERS */
    *result = 0.0e00;
    *abserr = 0.0e00;
    *neval = 0;
    *ier = 6;
    if (epsabs < 0.0e00 && epsrel < 0.0e00)
	goto L_130;
    hlgth = 5.0e-01 * (b - a);
    dhlgth = fabs (hlgth);
    centr = 5.0e-01 * (b + a);
    imsl_e1usr ("ON");
    fcentr = (*f) (centr);
    imsl_e1usr ("OFF");
    *neval = 21;
    *ier = 1;
    /*
     * COMPUTE THE INTEGRAL USING THE 10- 10- AND 21-POINT FORMULA
     */
    for (l = 1; l <= 3; l++) {
	if (l == 2)
	    goto L_50;
	if (l == 3)
	    goto L_80;

	res10 = 0.0e00;
	res21 = lv_w21b[5] * fcentr;
	resabs = lv_w21b[5] * fabs (fcentr);
	for (k = 1; k <= 5; k++) {
	    absc = hlgth * lv_x1[k - 1];
	    imsl_e1usr ("ON");
	    fval1 = (*f) (centr + absc);
	    fval2 = (*f) (centr - absc);
	    imsl_e1usr ("OFF");
	    fval = fval1 + fval2;
	    res10 += lv_w10[k - 1] * fval;
	    res21 += lv_w21a[k - 1] * fval;
	    resabs += lv_w21a[k - 1] * (fabs (fval1) + fabs (fval2));
	    savfun[k - 1] = fval;
	    fv1[k - 1] = fval1;
	    fv2[k - 1] = fval2;
	}
	ipx = 5;
	for (k = 1; k <= 5; k++) {
	    ipx += 1;
	    absc = hlgth * lv_x2[k - 1];
	    imsl_e1usr ("ON");
	    fval1 = (*f) (centr + absc);
	    fval2 = (*f) (centr - absc);
	    imsl_e1usr ("OFF");
	    fval = fval1 + fval2;
	    res21 += lv_w21b[k - 1] * fval;
	    resabs += lv_w21b[k - 1] * (fabs (fval1) + fabs (fval2));
	    savfun[ipx - 1] = fval;
	    fv3[k - 1] = fval1;
	    fv4[k - 1] = fval2;
	}
	/* TEST FOR CONVERGENCE. */
	*result = res21 * hlgth;
	resabs *= dhlgth;
	reskh = 5.0e-01 * res21;
	resasc = lv_w21b[5] * fabs (fcentr - reskh);
	for (k = 1; k <= 5; k++) {
	    resasc += lv_w21a[k - 1] * (fabs (fv1[k - 1] - reskh) + fabs (fv2[k - 1] -
		    reskh)) + lv_w21b[k - 1] * (fabs (fv3[k - 1] - reskh) +
		fabs (fv4[k - 1] - reskh));
	}
	*abserr = fabs ((res21 - res10) * hlgth);
	resasc *= dhlgth;
	goto L_110;
	/*
	 * COMPUTE THE INTEGRAL USING THE 43-POINT FORMULA
	 */
L_50:
	res43 = lv_w43b[11] * fcentr;
	*neval = 43;
	for (k = 1; k <= 10; k++) {
	    res43 += savfun[k - 1] * lv_w43a[k - 1];
	}
	for (k = 1; k <= 11; k++) {
	    ipx += 1;
	    absc = hlgth * lv_x3[k - 1];
	    imsl_e1usr ("ON");
	    fval = (*f) (absc + centr) + (*f) (centr - absc);
	    imsl_e1usr ("OFF");
	    res43 += fval * lv_w43b[k - 1];
	    savfun[ipx - 1] = fval;
	}
	/* TEST FOR CONVERGENCE. */
	*result = res43 * hlgth;
	*abserr = fabs ((res43 - res21) * hlgth);
	goto L_110;
	/*
	 * COMPUTE THE INTEGRAL USING THE 87-POINT FORMULA
	 */
L_80:
	res87 = lv_w87b[22] * fcentr;
	*neval = 87;
	for (k = 1; k <= 21; k++) {
	    res87 += savfun[k - 1] * lv_w87a[k - 1];
	}
	for (k = 1; k <= 22; k++) {
	    absc = hlgth * lv_x4[k - 1];
	    imsl_e1usr ("ON");
	    res87 += lv_w87b[k - 1] * ((*f) (absc + centr) +
		(*f) (centr - absc));
	    imsl_e1usr ("OFF");
	}
	*result = res87 * hlgth;
	*abserr = fabs ((res87 - res43) * hlgth);
L_110:
	if (resasc != 0.0e00 && *abserr != 0.0e00)
	    *abserr = resasc * imsl_f_min (1.0e00, pow (2.0e02 ** abserr / resasc, 1.5e00));
	if (resabs > uflow / (5.0e01 * epmach))
	    *abserr = imsl_f_max ((epmach * 5.0e01) * resabs, *abserr);
	if (*abserr <= imsl_f_max (epsabs, epsrel * fabs (*result)))
	    *ier = 0;
	/* JUMP OUT OF DO-LOOP */
	if (*ier == 0)
	    goto L_130;
    }
L_130:
    return;
}				/* end of function */
