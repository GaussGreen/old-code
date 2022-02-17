/*
 * $Log: armdef.h,v $
 * Revision 1.64  2004/01/21 13:14:10  mab
 * #define ARM_NB_TERMS          400
 *
 * Revision 1.63  2004/01/08 17:29:25  emezzine
 * added ARM_MAX_ITER
 *
 * Revision 1.62  2003/12/11 14:56:21  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.61  2003/12/09 09:47:06  emezzine
 * Rename a constant
 *
 * Revision 1.60  2003/12/09 09:36:36  emezzine
 * added two define to parameter the extrapolation method
 *
 * Revision 1.59  2003/11/26 18:10:16  jpriaudel
 * added flag in refvalue for credit
 *
 * Revision 1.58  2003/11/14 11:10:32  mab
 * Added: K_FX_VOL_SP_INTERP K_FX_VOL_TRIANG_INT
 *
 * Revision 1.57  2003/11/07 11:39:57  ebenhamou
 * added no base
 *
 * Revision 1.56  2003/10/23 16:56:58  ykhlif
 * ajout des modes de generation de portfolio mode SUMMIT et mode MANUAL
 *
 * Revision 1.55  2003/10/16 07:39:31  rguillemot
 * Double Barrier Pricing with Mixture
 *
 * Revision 1.54  2003/09/22 07:57:10  mab
 * Added: #define K_FX_VOL            7
 *
 * Revision 1.53  2003/09/09 11:40:34  ebenhamou
 * added flags TimeToStart and TimeToExpiry
 *
 * Revision 1.52  2003/09/04 17:32:48  ykhlif
 *  ajout de ModelTypeTable
 *
 * Revision 1.51  2003/09/03 10:24:38  ykhlif
 * ajout de SecurityTypeTable
 *
 * Revision 1.50  2003/08/26 11:58:54  ebenhamou
 * added flag for tenor and expiry
 *
 * Revision 1.49  2003/08/21 06:08:57  ebenhamou
 * added flag K_NX_ASINFEND
 *
 * Revision 1.48  2003/08/19 12:29:07  jpriaudel
 * Add flag for SABR
 *
 * Revision 1.47  2003/08/13 08:45:02  mab
 * Correction in PRCS Flags definition : K_OPT_NO K_OPT_FWD K_OPT
 *
 * Revision 1.46  2003/08/13 05:52:02  ebenhamou
 * change flag for inflation notional
 *
 * Revision 1.45  2003/08/04 17:58:04  ekoehler
 * definition of string constants S_ for armdef.cpp
 *
 * Revision 1.44  2003/07/31 12:40:30  ebenhamou
 * added flag for inflation leg
 *
 * Revision 1.43  2003/07/30 11:58:03  ekoehler
 * added some comment as well as a tag for daycount
 *
 * $Log: armdef.h,v $
 * Revision 1.64  2004/01/21 13:14:10  mab
 * #define ARM_NB_TERMS          400
 *
 * Revision 1.63  2004/01/08 17:29:25  emezzine
 * added ARM_MAX_ITER
 *
 * Revision 1.62  2003/12/11 14:56:21  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.61  2003/12/09 09:47:06  emezzine
 * Rename a constant
 *
 * Revision 1.60  2003/12/09 09:36:36  emezzine
 * added two define to parameter the extrapolation method
 *
 * Revision 1.59  2003/11/26 18:10:16  jpriaudel
 * added flag in refvalue for credit
 *
 * Revision 1.58  2003/11/14 11:10:32  mab
 * Added: K_FX_VOL_SP_INTERP K_FX_VOL_TRIANG_INT
 *
 * Revision 1.57  2003/11/07 11:39:57  ebenhamou
 * added no base
 *
 * Revision 1.56  2003/10/23 16:56:58  ykhlif
 * ajout des modes de generation de portfolio mode SUMMIT et mode MANUAL
 *
 * Revision 1.55  2003/10/16 07:39:31  rguillemot
 * Double Barrier Pricing with Mixture
 *
 * Revision 1.54  2003/09/22 07:57:10  mab
 * Added: #define K_FX_VOL            7
 *
 * Revision 1.53  2003/09/09 11:40:34  ebenhamou
 * added flags TimeToStart and TimeToExpiry
 *
 * Revision 1.52  2003/09/04 17:32:48  ykhlif
 *  ajout de ModelTypeTable
 *
 * Revision 1.51  2003/09/03 10:24:38  ykhlif
 * ajout de SecurityTypeTable
 *
 * Revision 1.50  2003/08/26 11:58:54  ebenhamou
 * added flag for tenor and expiry
 *
 * Revision 1.49  2003/08/21 06:08:57  ebenhamou
 * added flag K_NX_ASINFEND
 *
 * Revision 1.48  2003/08/19 12:29:07  jpriaudel
 * Add flag for SABR
 *
 * Revision 1.47  2003/08/13 08:45:02  mab
 * Correction in PRCS Flags definition : K_OPT_NO K_OPT_FWD K_OPT
 *
 * Revision 1.46  2003/08/13 05:52:02  ebenhamou
 * change flag for inflation notional
 *
 * Revision 1.45  2003/08/04 17:58:04  ekoehler
 * definition of string constants S_ for armdef.cpp
 *
 * Revision 1.44  2003/07/31 12:40:30  ebenhamou
 * added flag for inflation leg
 *
 * Revision 1.41  2003/07/16 06:57:25  ebenhamou
 * added flag for inflation
 *
 * Revision 1.40  2003/06/30 16:17:26  ebenhamou
 * added interpolation type for inflation
 *
 * Revision 1.39  2003/06/11 09:12:52  emezzine
 * Cut comment
 *
 * Revision 1.38  2003/06/11 08:28:42  emezzine
 * Added #define K_FRM_TOL       1.0E-12  /* double precision tolerance 
 *
 * Revision 1.37  2003/05/26 12:13:52  emezzine
 * Added FRMVOL_LAG_THRESHOLD
 *
 * Revision 1.36  2003/05/26 09:16:41  mab
 * Added : #define ARM_CALIB_GLOB      0
 * and #define ARM_CALIB_BOOTSTRAP 1
 *
 * Revision 1.35  2003/05/07 15:54:57  mab
 * Formatting
 *
 * Revision 1.34  2003/04/28 15:24:25  emezzine
 * Ajout de deux cst K_BOX_MULLER et INV_NORM_CUM
 *
 * Revision 1.33  2003/04/14 16:09:25  jmprie
 * ajout de K_MARKOVTREE_MARKOVIANDRIFT
 *
 * Revision 1.32  2003/04/07 08:36:26  jmprie
 * ajout de K_MARKOVTREE_SPOTDRIFTCALIB et K_MARKOVTREE_MARKOVIANDRIFTCALIB
 *
 * Revision 1.31  2003/03/18 17:05:12  mab
 * Added : Flags for PRCS K_OPT_NO, K_OPT_YES, K_OPT_FWD
 *
 * Revision 1.30  2003/03/07 17:22:42  mab
 * defines added for FRMMarkovTree : Drift Mode
 *
 * Revision 1.29  2003/03/07 12:56:52  mab
 * adedd : #define ZERO_PREC 1e-16
 *
 * Revision 1.28  2003/02/25 12:28:41  sgasquet
 * Modif K_STD par K_SPREAD
 *
 * Revision 1.27  2003/02/25 12:20:46  sgasquet
 * Ajout K_STD et K_DIGITALE pour corridor dans LDAna
 *
 * Revision 1.26  2003/02/10 16:53:47  mab
 * added : ARM_DEF_MAX_ITER, UPPER_INFINITE_BOUND, LOW_INFINITE_BOUND,
 * ONE_MONTH, MAX_PARAMS_FIT, ARM_SIGMA, ARM_RHO, ARM_NU
 * K_STEPUP_LEFT, K_GLOBAL_CALIB, K_BOOTSTRAP_CALIB
 *
 * Revision 1.25  2002/11/15 10:08:27  mab
 * SwapLeg Shedules flags added : #define K_START_DATES     0
 * .... #define K_NOTIONAL_VALUES 15
 *
 * Revision 1.24  2002/10/09 09:29:26  mab
 * Added : #define ARM_NB_MAX_CHAR_TERMS 30
 * ARM_NB_TERMS changed to 300
 *
 * Revision 1.23  2002/08/08 09:53:20  mab
 * Added : #define K_YES     1 And #define K_NO      0
 *
 * Revision 1.22  2002/05/23 14:16:14  sgasquet
 * Ajout typage modele pour Spread option
 *
 * Revision 1.21  2002/03/28 13:12:29  mab
 * Rajout de : #define K_AMORT_MORTGAGE                6
 *
 * Revision 1.20  2002/03/01 14:35:15  mab
 * Rajout de : #define K_ACCRUE
 *
 * Revision 1.19  2002/02/25 13:40:43  mab
 * Rajout de : K_AMORT_FIXED K_AMORT_FIXEDEND K_AMORT_PERCENTAGE K_AMORT_ANNUITY
 * K_AMORT_ANNUITY_REDUCED K_NX_NONE K_NX_START K_NX_END K_NX_BOTH
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : armdef.h                                                     */
/*                                                                            */
/* DESCRIPTION : global utilities                                             */
/*                                                                            */
/* DATE        : Tue Apr 22 1997                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*
 * do not add any #ifdef or #ifndef as they are not supported by applix
 * therefore do not add either any variable definition.
 */

/*
 * procedure to add new define variables
 *		1) add the appropriate flag name with its corresponding value
 *			something like #define FLAGNAME    1
 *
 *		2) check whether the category already exists
 *			if not, add the corresponding category type
 *			something like #define S_CATEGORY "CategoryName"
 *
 *		3) add the corresponding Table defined in armdef.cpp
 *			if not defined, 
 *				- add something a new line to ParamMappingTable
 *				with something like S_CATEGORY ,	 OurCategoryTable
 *				- define the corresponding table
 *				something like 
 *				tagArrayRow OurCategoryTable[]=
 *				{
 *					/// methodFlag			methodName
 *					{	KFLAGNAME,		"What it corresponds to" 	}, 
 *					/// do not forget to put this end of line
 *					{	PARAMVIEW_ENDOFLINE, PARAMVIEW_ENDOFLINE_CHAR }
 *				};
 *	
 *				many examples can be found in armdef.cpp
 *			else
 *				- just add one line to the corresponding category table		
 * 
 *		THAT IS IT!
 */

#ifdef WIN32
#error "ce fichier ne doit plus etre inclus"
#endif

#define ARM_DEFAULT_DATE "01/01/1981"
#define ARM_CONV_ERR  -9999
 
#define RET_OK    0
#define RET_KO    -1

#define K_YES     1
#define K_NO      0 


#define ARM_NB_TERMS          400
#define ARM_NB_MAX_CHAR_TERMS 30

 
#define ZERO_PREC 1e-16



/*---- For Calibration ----*/

#define ARM_DEF_MAX_ITER    1000
#define UPPER_INFINITE_BOUND 1e+20
#define ARM_MAX_ITER         100
#define LOW_INFINITE_BOUND   -UPPER_INFINITE_BOUND
#define ONE_MONTH            (1.0/12.0)

#define MAX_PARAMS_FIT       150


#define ARM_CALIB_GLOB      0
#define ARM_CALIB_BOOTSTRAP 1

/*-------FRMVOL-----------*/
 
#define FRMVOL_LAG_THRESHOLD 7 


#define K_FRM_TOL       1.0E-12  

/*---- BS Smiled parameters names ----*/
#define ARM_SIGMA    1
#define ARM_RHO      2
#define ARM_NU       3


/*---- For PRCS ----*/

#define K_OPT_NO      0
#define K_OPT_FWD     1
#define K_OPT         2

/*---- For SABR ----*/

#define K_LD          0
#define K_SABR_GEO    1
#define K_SABR_ARITH  2

/*---- Booleans ----*/
#define K_FALSE			0
#define K_TRUE			1

#define K_ROW			0
#define K_DIAG			1



/*---- Swap leg schedules ----*/
#define S_SWAPLEG_SCHEDULE  "SwapLeg Schedule"
#define K_START_DATES     0
#define K_END_DATES       1
#define K_RESET_DATES     2
#define K_PAY_DATES       3
#define K_FWD_START_DATES 4
#define K_FWD_END_DATES   5
#define K_AMORT_DATES     6


#define K_FWD_RATES       10
#define K_INT_DAYS        11
#define K_INT_TERMS       12
#define K_FLOW_VALUES     13
#define K_AMORT_VALUES    14
#define K_NOTIONAL_VALUES 15



/*---- RandomGenerators -------*/
/*---- Gaussian generator ----*/
#define S_GAUSSIAN_GENERATOR "Random Generator"
#define K_BOX_MULLER		0
#define K_INV_NORM_CUM	1
/*---- MC generator ----*/
#define S_MC_GENERATOR  "Monte Carlo Generator"
#define K_MC_SCRAMBLE		0
#define K_MC_FAURE		1
#define K_MC_SIMPLE		2


 
#define ARM_NULL_OBJECT -11111
 
#define K_YEAR_LEN         365        /* double precision tolerance */
#define K_NEW_DOUBLE_TOL   1.0E-08
#define K_DOUBLE_TOL       2.776E-17  /* double precision tolerance */
#define K_HUGE_DOUBLE      9.0E37     /* huge doubles */
#define K_LOG_HUGE_DOUBLE  87.51 /* log (base e) of limit for floating overflow */
#define K_TINY                  1.0e-20
 
#define YC_INSTANTANEOUS_MAT    1.0E-6
 
#define S_OPTION_TYPE "Option Type"
#define K_CALL                  1
#define K_PUT                   -1
 
#define S_CAPLIKE_TYPE          "Cap Like Type"
#define K_CAP                   1
#define K_FLOOR                 -1
#define K_CAPFLOOR				0
#define K_NO_CAPFLOOR           2
#define K_ACCRUE				3

#define K_IN_FINE		0
#define K_TRI			1


/*---- swap leg type ----*/
#define S_LEGTYPE               "Swap Leg Type"
#define K_FLOATING_LEG         -1
#define K_MARKET_FIXED_LEG      0
#define K_FIXED_LEG             1
#define K_ZEROCOUPON_LEG        2
#define K_YEARTOYEAR_LEG		3
#define K_OATTYPE_LEG			4
#define K_GENERICINF_LEG		5


#define K_ZEROCURVE        10

#define K_UP                  1
#define K_DOWN               -1
 
#define K_IN				1
/* #define K_OUT               0 */
#define K_OUT				-1


/* Exercise styles */
#define S_EXERCISE_STYLE         "Exercise Style"
#define K_EUROPEAN              0
#define K_AMERICAN              1
#define K_BERMUDAN              2
#define K_AUTOMATIC             3
#define K_CUSTOMIZED           -1
 
#define S_COMPUT_TYPE_RESULT    "Computing Type Result"
#define K_PRICE                 1
#define K_YIELD                 0
#define K_INDEX_RATE            2

#define K_PRIX                  0
#define K_TAUX                  1
#define K_PRICE_TO_RATE         1
#define K_RATE_TO_PRICE         2
#define K_KNOCK_IN              1
/* #define K_KNOCK_OUT             0 */
#define S_BARRIER               "Barrier"
#define K_KNOCK_OUT             -1
#define K_UP                    1
#define K_DOWN                 -1
#define K_DOUBLE				2
 
 
#define ARM_MIN_INT             INT_MIN
 
#define ARM_MAX_DOUBLE          1.0E+306
 
 
/* Forward Rules */
#define S_FORWARD_RULES      "Forward Rules"
#define K_PREVIOUS         -1
#define K_MOD_PREVIOUS     -2
#define K_FOLLOWING         1
#define K_MOD_FOLLOWING     2
 
/* Interest Rules */
#define S_INTEREST_RULES "Interest Rules"
#define K_ADJUSTED          1
#define K_UNADJUSTED        0
 
/* Stub Rules */
#define S_STUB_RULES         "Stub Rules"
#define K_SHORTSTART        1
#define K_LONGSTART         2
#define K_SHORTEND          3
#define K_LONGEND           4
 
/* Amortization Type */
#define S_AMORTIZATION_TYPE              "Amortization Type"
#define K_AMORT_FIXED                   1
#define K_AMORT_FIXEDEND                2
#define K_AMORT_PERCENTAGE              3
#define K_AMORT_ANNUITY                 4
#define K_AMORT_ANNUITY_REDUCED         5
#define K_AMORT_MORTGAGE                6


/* Nominal Exchange type */
#define S_NOMINAL_EXCHANGE_TYPE         "Nominal Exchange type"
#define K_NX_NONE                       0
#define K_NX_START                      1
#define K_NX_END                        2
#define K_NX_BOTH                       3
#define K_NX_INFEND                     4
#define K_NX_ASINFEND					5


/* Frequency */
#define S_FREQUENCY     "Frequency"
#define K_ANNUAL       1
#define K_SEMIANNUAL   2
#define K_QUARTERLY    4
#define K_BIMONTHLY    6
#define K_MONTHLY     12
#define K_WEEKLY      52
#define K_DAILY      365
#define K_ZEROCOUPON   0
#define K_DEF_FREQ    -1 /* default value used differently in various contexts */
 
/* Date types */
#define S_DATE_TYPES         "Date Types"
#define K_FLOW_PAYMENT       1
#define K_INDEX_RESET        2
#define K_START_INTEREST     3
#define K_END_INTEREST       4

/* Calculation Method for ref Value */
#define S_CALCULATION_METHOD     "Calculation Method for ref Value"
#define K_CONSTANT_REF           0
#define K_LININTERPOL_REF        1
#define K_DISCRETE_REF           2
#define K_ZEROCOUPON_REF         3
#define K_STEPUP_RIGHT           4
#define K_STEPUP_LEFT            5 
#define K_PERFECT_DISCRETE_REF   6


/* for spread option*/
#define K_2LOG		            0
#define K_NOR	                1
#define K_LOG				    2

#define K_INPUTED               0
#define K_COMPUTED              1
/* end spread option*/



/* Interpolating Method */
/* for Curves and for Corridors with MC */
/* and inflation curve type */
#define S_INTERPOL_METHOD   "Interpolation Method"
#define K_CONTINUOUS		0
#define K_LINEAR			1
#define K_SQRT				2
#define K_SLOPELIN			3
#define K_SLOPESQRT			4

#define K_CPILINEAR			5 /* LINEAR sur les CPI */
#define K_CPISTEPWISE		6 /* StepWise sur les CPI */
#define K_ZCLINEAR			7 /* Linear sur les ZC */
#define K_ZCCTFWD			8 /* Utilise Ct Fwd ZC */

#define K_CPISTEPWISESTART	9 /* StepWiseStart sur les CPI */
#define K_CPISTEPWISEEND	10 /* StepWiseEnd sur les CPI */
#define K_CPISTEPWISEMIDDLE	11 /* StepWiseMiddle sur les CPI */

/* Extrapolation method*/
#define S_EXTRAPOLATION_METHOD    "Extrapolation Method"
#define K_LASTTWO		21
#define K_MIDDLE		22
#define K_FIRSTTWO		23

#define K_LINEAR_EXTRAPOL		0 /* Linear out the boundaries*/
#define K_CONSTANT_EXTRAPOL		1 /* constant out the boundaries*/

/* for pricing corridor with LDAna */
#define K_SPREAD	0
#define K_DIGITALE	1


/* Missing Maturity Filling Method */
#define S_MISSING_MATURITY_FILLING   "Curve Mod"
#define K_PAR			0
#define K_RAW			1

/* Market Segment */ 
#define S_MARKET_SEGMENT    "Curve Priorities"
#define K_FUT			0
#define K_MM			1
#define K_SWAP			2
#define K_BOND			3

/* Rate types */
#define K_FLOAT_RATE        -1000000.0
#define K_FIXED_RATE        0.0
#define K_MARKET_RATE       -1.0
#define K_ATMF_CAPLETS      -2.0
 
/* Index Types */
#define S_INDEX_TYPES   "Index Types"
#define    K_FIXED     0
#define    K_LIBOR1M   1
#define    K_LIBOR2M   2
#define    K_LIBOR3M   3
#define    K_LIBOR6M   4
#define    K_LIBOR1Y   5
#define    K_PIBOR1M   6
#define    K_PIBOR2M   7
#define    K_PIBOR3M   8
#define    K_PIBOR6M   9
#define    K_PIBOR1Y   10
#define    K_EURIBOR1M   11
#define    K_EURIBOR2M   12
#define    K_EURIBOR3M   13
#define    K_EURIBOR6M   14
#define    K_EURIBOR1Y   15
#define    K_CMT1      16
#define    K_CMT2      17
#define    K_CMT5      18
#define    K_CMT10     19
#define    K_CMT15     20
#define    K_CMT20     21
#define    K_CMT30     22
#define    K_CMS1      23
#define    K_CMS2      24
#define    K_CMS3      25	
#define    K_CMS4      26
#define    K_CMS5      27
#define    K_CMS6      28
#define    K_CMS7      29
#define    K_CMS8      30
#define    K_CMS9      31
#define    K_CMS10     32
#define    K_CMS11     33
#define    K_CMS12     34
#define    K_CMS13     35
#define    K_CMS14     36
#define    K_CMS15     37
#define    K_CMS16     38 
#define    K_CMS17     39
#define    K_CMS18     40
#define    K_CMS19     41
#define    K_CMS20     42
#define    K_CMS21     43
#define    K_CMS22     44
#define    K_CMS23     45
#define    K_CMS24     46
#define    K_CMS25     47
#define    K_CMS26     48
#define    K_CMS27     49
#define    K_CMS28     50
#define    K_CMS29     51
#define    K_CMS30     52
#define    K_TEC5      53
#define    K_TEC10     54
#define    K_T4M       55
#define    K_T4M_FIXED 56
#define    K_TAM       57
#define    K_TAG       58
#define    K_TMP       59

#define K_STD           0
#define K_CMSLIBOR      1
#define K_CMSPIBOR      2
#define K_CMT           3
#define K_TEC           4
 
 
 
/* Receive or Pay */
#define S_RECEIVE_PAY  "Receive or Pay"
#define K_RCV               1
#define K_PAY               -1
 
#define K_DOMESTIC          1
#define K_FOREIGN           0
 
/* Rate Compounding methods */
#define S_RATE_COMPOUNDING  "Rate Compounding methods"
#define K_COMP_CONT         0
#define K_COMP_PROP         -1
#define K_COMP_ANNUAL       1
#define K_COMP_SEMIANNUAL   2
#define K_COMP_QUARTERLY    4
#define K_COMP_MONTHLY      12
#define K_COMP_BIMONTHLY    6
#define K_COMP_DAILY_360    360
#define K_COMP_DAILY_365    365
/* fin de categorie*/ 

/* Rate Type : Money Market, Actuarial (or Bond) and Continuous Compounding*/
#define K_MM_RATE         -1
#define K_ACT_RATE       1
#define K_CONT_RATE       0
 
 
 
/* reset and payment rules */
#define S_TIMING_MOD         "Timing Mod"
#define K_ADVANCE           1
#define K_ARREARS           -1
 
 
 
/* days of the week */
#define S_WEEKDAY           "Week Day"
#define K_SUNDAY            0
#define K_MONDAY            1
#define K_TUESDAY           2
#define K_WEDNESDAY         3
#define K_THURSDAY          4
#define K_FRIDAY            5
#define K_SATURDAY          6
 
 
 
 
/* Basis rules for day count */
#define S_DAYCOUNT				"DayCount" 
#define KACTUAL_ACTUAL          1
#define KACTUAL_365             2
#define KACTUAL_360             3
#define K30_360                 4
#define KACTUAL_REAL            5
#define KACTUAL_FEB29           6
#define KACTUAL_ISMA            7
#define KNOBASE			        8
 
 
 
/*--- strike type ---*/
 
#define K_STK_TYPE_YIELD        0
#define K_STK_TYPE_PRICE        1
#define K_STK_TYPE_DELTA        2
#define K_STK_TYPE_MATURITY     3
#define K_STK_TYPE_MATU_SMILE   4
 


/*--- Vol type ----*/
 
#define K_ATMF_VOL			0
#define K_SMILE_VOL			1
#define K_CUBE_VOL			2
#define K_TENOR				3
#define K_EXPIRY			4
#define K_TIMETOSTART		5
#define K_TIMETOEXPIRY		6
#define K_FX_VOL            7
#define K_FX_VOL_SP_INTERP  8
#define K_FX_VOL_TRIANG_INT 9

 
/* Notional calculation (zero coupon swap) */
 
#define K_NTL_FIXED          0
#define K_NTL_IMPLIED        1
 
#define K_SVITY_SPOT            1
#define K_SVITY_SPOT2           2
#define K_SVITY_DRATE           3
#define K_SVITY_TIME            4
#define K_SVITY_SIGMA           5
#define K_SVITY_SPOTSIGMA       6
 
#define ARM_VIEW_FILE			"VF"
#define ARM_CREDIT_FILE			"SC"

#define K_H                     0.00001



/*---- For TAM/T4M Curve ----*/

#define FIRST_DF_FOR_NEWTON     0.75
#define PRECISION               1e-10
#define H_FOR_DERIVATIVE        1.0e-8
#define TAM_PERIOD_IN_MONTH     12
#define MID_MONTH_DAY           15
#define ISPRESENT               0.0001




/* Mode de depart d'un swap T4M/TAM */
 
#define K_BWD                   0
#define K_FWD                   1



/* Mode de depart d'un swap T4M/TAM */
 
#define K_T4M_INDEX             0
#define K_TAM_INDEX             1



/* Longstaff Monte Carlo Mode       */
#define S_MC_MOD                "Monte Carlo Mod"
#define K_SCRMC                 0
#define K_FAURE                 1
#define K_SIMMC                 2


/* Tree fine step modes             */
#define K_L_NONE                  0
#define K_L_DEFAULT               1
#define K_L_QUICK                 2    
#define K_L_2                     3
#define K_L_4                     4
#define K_L_8                     5
#define K_L_12                    6
#define K_L_16                    7
#define K_L_24                    8
#define K_L_7DAY                  9
#define K_L_5DAY                  10
#define K_L_3DAY                  11
#define K_L_2DAY                  12
#define K_L_1DAY                  13
#define K_L_A                     14
#define K_L_B                     15
#define K_L_C                     16

/*---- FRM autocalibration modes ----*/
#define S_FRM_AUTOCALIBRATION_MOD   "FRM autocalibration modes"
#define K_SWOPT                     0
#define K_IRG                       1
#define K_SWOPT_IRG                 2
#define K_IRG_SWOPT                 3


/*---- Calibration Type ----*/
#define S_CALIBRATION_TYPE          "Calibration Type"
#define K_DIAG_TYPE                10


/*---- Calibration Methode ----*/
#define S_CALIBRATION_METHOD    "Calibration Method"
#define K_GLOBAL_CALIB       0
#define K_BOOTSTRAP_CALIB    1


/* FRMMarkovTree Drift Mode */
#define S_FRMMARKOVTREE_DRIFTMOD           "FRMMarkovTree Drift Mode"
#define K_MARKOVTREE_DRIFTAPPROX            0
#define K_MARKOVTREE_DRIFTCALIB             1
#define K_MARKOVTREE_SPOTDRIFTCALIB         2
#define K_MARKOVTREE_MARKOVIANDRIFT         3

#define K_MARKOVTREE_MARKOVIANDRIFTCALIB    4

#define S_SPACING                        "Tree Spacing"
#define K_CONSTSTANTESPACING             0

/* Security Types */
#define S_SECURITY_TYPES			"Security Types"
#define S_MODEL_TYPES				"Model Types"

#define K_SUMMIT				0
#define K_MANUAL				1

/* Replic Mode */
#define S_REPLIC_MODE_TYPES			"Replic Mode"
#define K_CONST_STEP			0
#define K_CONST_PREC			1
#define K_ALTN_PREC				2

/* Replic Mode */
#define S_REPLIC_STOP_MODE_TYPES	"Replic Stop Mode"
#define K_STOP_PRICE				0
#define K_STOP_PRICE_WEIGHT			1

/*----------------------------------------------------------------------------*/
/*---- End of File ----*/
