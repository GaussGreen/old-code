/*
 * Copyright CDC IXIS CM Paris July 2003
 * $Log: armdef.cpp,v $
 * Revision 1.18  2004/06/03 16:45:21  mab
 * added:  {       K_LIVRET_A,                     "LIVRETA"}
 *
 * Revision 1.17  2004/01/21 12:06:13  emezzine
 * remove ARM_FRMANALYTIC
 *
 * Revision 1.16  2003/12/11 14:57:06  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.15  2003/11/07 11:39:41  ebenhamou
 * added no base
 *
 * Revision 1.14  2003/10/16 07:44:31  rguillemot
 * Double Barrier Pricing with Mixture
 *
 * Revision 1.13  2003/09/04 17:31:40  ykhlif
 * ajout de ModelTypeTable
 *
 * Revision 1.12  2003/09/03 10:23:32  ykhlif
 * ajout de SecurityTypeTable
 *
 * Revision 1.11  2003/08/21 06:09:15  ebenhamou
 * added flag K_NX_ASINFEND
 *
 * Revision 1.10  2003/08/15 07:05:11  ebenhamou
 * formatting for better reading
 *
 * Revision 1.9  2003/08/13 11:18:05  ekoehler
 * GaussianGenaratorTable instead of RandomGeneratorTable
 *
 * Revision 1.8  2003/08/13 06:43:57  ebenhamou
 * add the tag for inflation notional
 *
 * Revision 1.7  2003/08/04 17:58:41  ekoehler
 * Definition of ARM_TagArrayRow variables and extension of ARM_ParamMappingRow
 *
 * Revision 1.4  2003/07/30 12:20:36  ekoehler
 * added end of line
 *
 * Revision 1.3  2003/07/30 12:19:32  ekoehler
 * added endofline and a lne at the end
 *
 * Revision 1.2  2003/07/30 12:16:23  ekoehler
 * added log and copyright
 *
 *
 */
#include "armdef.h"
#include "paramview.h"





ARM_TagArrayRow SwapLegScheduleTable[]=
	{
		/* methodFlag				methodName */
		{	K_START_DATES,			"START_DATES"	}, 
		{	K_END_DATES,			"END_DATES"	}, 
		{	K_RESET_DATES,			"RESET_DATES"	}, 
		{	K_PAY_DATES,			"PAY_DATES"	}, 
		{	K_FWD_START_DATES,		"FWD_START_DATES"	}, 
		{	K_FWD_END_DATES,		"FWD_END_DATES"	}, 
		{	K_AMORT_DATES,			"AMORT_DATES"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};



ARM_TagArrayRow GaussianGeneratorTable[]=
	{
		/* methodFlag				methodName */
		{	K_BOX_MULLER,			"BOX_MULLER	"	}, 
		{	K_INV_NORM_CUM,			"INV_NORM_CUM"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow MonteCarloGeneratorTable[]=
	{
		/* methodFlag				methodName */
		{	K_MC_SCRAMBLE,			"MC_SCRAMBLE	"	}, 
		{	K_MC_FAURE,				"MC_FAURE	"	}, 
		{	K_MC_SIMPLE,			"MC_SIMPLE	"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow OptionTypeTable[]=
	{
		/* methodFlag				methodName */
		{	K_CALL,					"CALL"	}, 
		{	K_PUT,					"PUT"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow CapLikeTypeTable[]=
	{
		/* methodFlag				methodName */
		{	K_CAP,					"CAP"	}, 
		{	K_FLOOR,				"FLOOR"	}, 
		{	K_CAPFLOOR,				"CAPFLOOR"	}, 
		{	K_NO_CAPFLOOR,			"NO_CAPFLOOR"	}, 
		{	K_ACCRUE,				"ACCRUE"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow SwapLegTypeTable[]=
	{
		/* methodFlag				methodName */
		{	K_FLOATING_LEG,			"FLOATING_LEG "	}, 
		{	K_MARKET_FIXED_LEG,		"MARKET_FIXED_LEG"	}, 
		{	K_FIXED_LEG,			"FIXED_LEG"	}, 
		{	K_ZEROCOUPON_LEG,		"ZEROCOUPON_LEG"	}, 
		{	K_YEARTOYEAR_LEG,		"YEARTOYEAR_LEG"	}, 
		{	K_OATTYPE_LEG,			"OATTYPE_LEG"	}, 
		{	K_GENERICINF_LEG,		"GENERICINF_LEG"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow ExerciseStyleTable[]=
	{
		/* methodFlag				methodName */
		{	K_EUROPEAN,				"EUROPEAN"	}, 
		{	K_AMERICAN,				"AMERICAN"	}, 
		{	K_BERMUDAN,				"BERMUDAN"	}, 
		{	K_AUTOMATIC,			"AUTOMATIC"	}, 
		{	K_CUSTOMIZED,			"CUSTOMIZED"	}, 
		/* do not forget to put this following end of line
        used as an equivalent of '\0' to stop a loop*/
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow ComputingTypeResultTable[]=
	{
		/* methodFlag				methodName */
		{	K_PRICE,				"PRICE"	}, 
		{	K_YIELD,				"YIELD"	}, 
		{	K_INDEX_RATE,			"INDEX_RATE"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow BarrierTable[]=
	{
		/* methodFlag				methodName */
		{	K_KNOCK_IN,				"KNOCK_IN"	},
		{	K_KNOCK_OUT,			"KNOCK_OUT"	},
		{	K_UP,					"UP"	},
		{	K_DOWN,					"DOWN"	},
		{	K_DOUBLE,				"DOUBLE" },
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow FwdRulesTable[]=
	{
		/* methodFlag				methodName */
		{	K_PREVIOUS,				"PREVIOUS"	}, 
		{	K_MOD_PREVIOUS,			"MOD_PREVIOUS"	}, 
		{	K_FOLLOWING,			"FOLLOWING"	}, 
		{	K_MOD_FOLLOWING,		"MOD_FOLLOWING"	}, 
		/* do not forget to put this following end of line
        used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow InterestRulesTable[]=
	{
		/* methodFlag				methodName */
		{	K_ADJUSTED,				"ADJUSTED"	}, 
		{	K_UNADJUSTED,			"UNADJUSTED"	}, 
		{	K_MATUNADJUSTED,		"MATUNADJUSTED"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow StubRulesTable[]=
	{
		/* methodFlag				methodName */
		{	K_SHORTSTART,			"SHORTSTART"	}, 
		{	K_LONGSTART,			"LONGSTART"	}, 
		{	K_SHORTEND,				"SHORTEND"	}, 
		{	K_LONGEND,				"LONGEND"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow AmortizationTypeTable[]=
	{
		/* methodFlag				methodName */
		{	K_AMORT_FIXED,			"AMORT_FIXED"	}, 
		{	K_AMORT_FIXEDEND,		"AMORT_FIXEDEND"	}, 
		{	K_AMORT_PERCENTAGE,		"AMORT_PERCENTAGE"	}, 
		{	K_AMORT_ANNUITY,		"AMORT_ANNUITY"	}, 
		{	K_AMORT_ANNUITY_REDUCED,"AMORT_ANNUITY_REDUCED "	}, 
		{	K_AMORT_MORTGAGE,		"AMORT_MORTGAGE"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow NominalExchangeTypeTable[]=
	{
		/* methodFlag		methodName */
		{	K_NX_NONE,		"NX_NONE"		}, 
		{	K_NX_START,		"NX_START"		}, 
		{	K_NX_END,		"NX_END"		}, 
		{	K_NX_BOTH,		"NX_BOTH"		}, 
		{	K_NX_INFEND,	"NX_INF_END"	}, 
		{	K_NX_ASINFEND,	"NX_AS_END_INF"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow FrequencyTable[]=
	{
		/* methodFlag		methodName */
		{	K_ANNUAL,		"ANNUAL"	}, 
		{	K_SEMIANNUAL,	"SEMIANNUAL"	}, 
		{	K_QUARTERLY,	"QUARTERLY"	}, 
		{	K_BIMONTHLY,	"BIMONTHLY"	}, 
		{	K_MONTHLY,		"MONTHLY"	}, 
		{	K_WEEKLY,		"WEEKLY"	}, 
		{	K_DAILY,		"DAILY"	}, 
		{	K_ZEROCOUPON,	"ZEROCOUPON"	}, 
		{	K_DEF_FREQ,		"DEF_FREQ, default value used differently in various contexts "	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow DateTypesTable[]=
	{
		/* methodFlag		methodName */
		{	K_FLOW_PAYMENT,		"FLOW_PAYMENT"	}, 
		{	K_INDEX_RESET,		"INDEX_RESET"	}, 
		{	K_START_INTEREST,	"START_INTEREST"	}, 
		{	K_END_INTEREST,		"END_INTEREST"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow CalculationMethodTable[]=
	{
		/* methodFlag		methodName */
		{	K_CONSTANT_REF,		"CONSTANT_REF"	}, 
		{	K_LININTERPOL_REF,	"LININTERPOL_REF"	}, 
		{	K_DISCRETE_REF,		"DISCRETE_REF"	}, 
		{	K_ZEROCOUPON_REF,	"ZEROCOUPON_REF"	}, 
		{	K_STEPUP_RIGHT,		"STEPUP_RIGHT"	}, 
		{	K_STEPUP_LEFT,		"STEPUP_LEFT"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};



ARM_TagArrayRow InterpolationMethodTable[]=
	{
		/* methodFlag				methodName */
		{	K_CONTINUOUS,			"CONTINUOUS"	}, 
		{	K_LINEAR,				"LINEAR"	}, 
		{	K_SQRT,					"SQRT"	}, 
		{	K_SLOPELIN,				"SLOPELIN"	}, 
		{	K_SLOPESQRT,			"SLOPESQRT"	}, 
		{	K_CPILINEAR,			"LINEAR on the CPI"	}, 
		{	K_CPISTEPWISE,			"StepWise on the CPI "	}, 
		{	K_ZCLINEAR,				"Linear on the Zero Coupon"	}, 
		{	K_ZCCTFWD,				"Uses Constant t Forward on the Zero Coupon"	}, 
		{	K_CPISTEPWISESTART,		"StepWiseStart on the CPI"	}, 
		{	K_CPISTEPWISEEND,		"StepWiseEnd on the CPI"	}, 
		{	K_CPISTEPWISEMIDDLE,	"StepWiseMiddle on the CPI"	}, 
        {   K_SPLINE,               "Spline Interpolation" },
        /* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow ExtrapolationMethodTable[]=
	{
		/* methodFlag		methodName */
		{	K_LASTTWO		,		"LASTTWO		"	}, 
		{	K_MIDDLE		,		"MIDDLE		"	}, 
		{	K_FIRSTTWO		,		"FIRSTTWO		"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow CurveModTable[]=
	{
		/* methodFlag			methodName */
		{	K_PAR		,		"PAR"	}, 
		{	K_RAW		,		"RAW"	}, 
		{	K_FORWARD	,		"FORWARD"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow CurvePrioritiesTable[]=
	{
		/* methodFlag			methodName */
		{	K_FUT		,		"FUT		"	}, 
		{	K_MM		,		"MM		"	}, 
		{	K_SWAP		,		"SWAP		"	}, 
		{	K_BOND		,		"BOND		"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow IndexTypesTable[]=
	{
		/* methodFlag			methodName */
		{	K_FIXED,			"FIXED"	}, 
		{	K_LIBOR1M,			"LIBOR1M"	}, 
		{	K_LIBOR2M,			"LIBOR2M"	}, 
		{	K_LIBOR3M,			"LIBOR3M"	}, 
		{	K_LIBOR6M,			"LIBOR6M"	}, 
		{	K_LIBOR1Y,			"LIBOR1Y"	}, 
		{	K_PIBOR1M,			"PIBOR1M"	}, 
		{	K_PIBOR2M,			"PIBOR2M"	}, 
		{	K_PIBOR3M,			"PIBOR3M"	}, 
		{	K_PIBOR6M,			"PIBOR6M"	}, 
		{	K_PIBOR1Y,			"PIBOR1Y"	}, 
		{	K_EURIBOR1M,		"EURIBOR1M"	}, 
		{	K_EURIBOR2M,		"EURIBOR2M"	}, 
		{	K_EURIBOR3M,		"EURIBOR3M"	}, 
		{	K_EURIBOR6M,		"EURIBOR6M"	}, 
		{	K_EURIBOR1Y,		"EURIBOR1Y"	}, 
		{	K_CMT1,				"CMT1"	}, 
		{	K_CMT2,				"CMT2"	}, 
		{	K_CMT5,				"CMT5"	}, 
		{	K_CMT10,			"CMT10"	}, 
		{	K_CMT15,			"CMT15"	}, 
		{	K_CMT20,			"CMT20"	}, 
		{	K_CMT30,			"CMT30"	}, 
		{	K_CMS1,				"CMS1"	}, 
		{	K_CMS2,				"CMS2"	}, 
		{	K_CMS3,				"CMS3 "	}, 
		{	K_CMS4,				"CMS4"	}, 
		{	K_CMS5,				"CMS5"	}, 
		{	K_CMS6,				"CMS6"	}, 
		{	K_CMS7,				"CMS7"	}, 
		{	K_CMS8,				"CMS8"	}, 
		{	K_CMS9,				"CMS9"	}, 
		{	K_CMS10,			"CMS10"	}, 
		{	K_CMS11,			"CMS11"	}, 
		{	K_CMS12,			"CMS12"	}, 
		{	K_CMS13,			"CMS13"	}, 
		{	K_CMS14,			"CMS14"	}, 
		{	K_CMS15,			"CMS15"	}, 
		{	K_CMS16,			"CMS16"	}, 
		{	K_CMS17,			"CMS17"	}, 
		{	K_CMS18,			"CMS18"	}, 
		{	K_CMS19,			"CMS19"	}, 
		{	K_CMS20,			"CMS20"	}, 
		{	K_CMS21,			"CMS21"	}, 
		{	K_CMS22,			"CMS22"	}, 
		{	K_CMS23,			"CMS23"	}, 
		{	K_CMS24,			"CMS24"	}, 
		{	K_CMS25,			"CMS25"	}, 
		{	K_CMS26,			"CMS26"	}, 
		{	K_CMS27,			"CMS27"	}, 
		{	K_CMS28,			"CMS28"	}, 
		{	K_CMS29,			"CMS29"	}, 
		{	K_CMS30,			"CMS30"	}, 
		{	K_TEC5,				"TEC5"	}, 
		{	K_TEC10,			"TEC10"	}, 
		{	K_T4M,				"T4M"	}, 
		{	K_T4M_FIXED,		"T4M_FIXED"	}, 
		{	K_TAM,				"TAM"	}, 
		{	K_TAG,				"TAG"	}, 
		{	K_EONIA,			"EONIA"	}, 
		{	K_STD,				"STD"	}, 
		{	K_CMSLIBOR,			"CMSLIBOR"	}, 
		{	K_CMSPIBOR,			"CMSPIBOR"	}, 
		{	K_CMT,				"CMT"	}, 
		{	K_TEC,				"TEC"	}, 
		{	K_LIVRET_A,			"LIVRETA"},
		{	K_EUR3M,			"EUR3M"},
		{	K_EUR12,			"EUR12"},
		{	K_EUR1M,			"EUR1M"},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow ReceiveOrPayTable[]=
	{
		/* methodFlag		methodName */
		{	K_RCV,		"RCV"	}, 
		{	K_PAY,		"PAY"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow RateCompoundingMethodsTable[]=
	{
		/* methodFlag			methodName */
		{	K_COMP_CONT,		"COMPOUNDING_CONT"	}, 
		{	K_COMP_PROP,		"COMPOUNDING_PROP"	}, 
		{	K_COMP_ANNUAL,		"COMPOUNDING_ANNUAL"	}, 
		{	K_COMP_SEMIANNUAL,	"COMPOUNDING_SEMIANNUAL"	}, 
		{	K_COMP_QUARTERLY,	"COMPOUNDING_QUARTERLY"	}, 
		{	K_COMP_MONTHLY,		"COMPOUNDING_MONTHLY"	}, 
		{	K_COMP_BIMONTHLY,	"COMPOUNDING_BIMONTHLY"	}, 
		{	K_COMP_DAILY_360,	"COMPOUNDING_DAILY_360"	}, 
		{	K_COMP_DAILY_365,	"COMPOUNDING_DAILY_365"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow TimingModTable[]=
	{
		/* methodFlag		methodName */
		{	K_ADVANCE,		"ADVANCE"	}, 
		{	K_ARREARS,		"ARREARS"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow WeekDayTable[]=
	{
		/* methodFlag		methodName */
		{	K_SUNDAY,		"SUNDAY"	}, 
		{	K_MONDAY,		"MONDAY"	}, 
		{	K_TUESDAY,		"TUESDAY"	}, 
		{	K_WEDNESDAY,	"WEDNESDAY"	}, 
		{	K_THURSDAY,		"THURSDAY"	}, 
		{	K_FRIDAY,		"FRIDAY"	}, 
		{	K_SATURDAY,		"SATURDAY"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
};


ARM_TagArrayRow DayCountTable[]=
	{
		/* methodFlag			methodName      */
		{	KACTUAL_ACTUAL,		"ACTUAL ACTUAL"	}, 
		{	KACTUAL_365,		"ACTUAL 365"	},  
		{	KACTUAL_360,		"ACTUAL 360"	},  
		{	K30_360,			"30 360"		},  
		{	KACTUAL_REAL,		"ACTUAL REAL"	},  
		{	KACTUAL_FEB29,		"ACTUAL FEB29"	},  
		{	KACTUAL_ISMA,		"ACTUAL ISMA"	},
		{	KNOBASE,			"NO BASE"		},
		{	K30_360E,			"30 360E"		},
		/* do not forget to put this following end of line
        used as an equivalent of '\0' to stop a loop*/
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow MonteCarloModTable[]=
	{
		/* methodFlag		methodName */
		{	K_SCRMC,		"SCRMC"	}, 
		{	K_FAURE,		"FAURE"	}, 
		{	K_SIMMC,		"SIMMC"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow FRMAutoCalibrationModesTable[]=
	{
		/* methodFlag			methodName */
		{	K_SWOPT,			"SWOPT"	}, 
		{	K_IRG,				"IRG"	}, 
		{	K_SWOPT_IRG,		"SWOPT_IRG"	}, 
		{	K_IRG_SWOPT,		"IRG_SWOPT"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow CalibrationTypeTable[]=
	{
		/* methodFlag				methodName */
		{	K_DIAG_TYPE,			"DIAG_TYPE"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow CalibrationMethodTable[]=
	{
		/* methodFlag				methodName */
		{	K_GLOBAL_CALIB,			"GLOBAL_CALIB"	}, 
		{	K_BOOTSTRAP_CALIB,		"BOOTSTRAP_CALIB"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow FRMMarkovTreeDriftModeTable[]=
{
	/* methodFlag							methodName */
	{	K_MARKOVTREE_DRIFTAPPROX,			"MARKOVTREE_DRIFTAPPROX"	}, 
	{	K_MARKOVTREE_DRIFTCALIB,			"MARKOVTREE_DRIFTCALIB"	}, 
	{	K_MARKOVTREE_SPOTDRIFTCALIB,		"MARKOVTREE_SPOTDRIFTCALIB"	}, 
	{	K_MARKOVTREE_MARKOVIANDRIFT,		"MARKOVTREE_MARKOVIANDRIFT"	}, 
	/* do not forget to put this following end of line
	used as an equivalent of '\0' to stop a loop */
	{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow TreeSpacingTable[]=
	{
		/* meth2odFlag						methodName */
		{	K_CONSTSTANTESPACING,			"CONSTANTE_SPACING"	}, 
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};
ARM_TagArrayRow SecurityTypeTable[]=
	{
		/* Security Types Flag             Securiy Type Name */
		{	ARM_SWAPLEG,				"SWAPLEG"				},
		{	ARM_FIXEDLEG,				"FIXEDLEG"				},
		{	ARM_CMSLEG,					"CMSLEG"				},
		{	ARM_CMTLEG,					"CMTLEG"				},
		{	ARM_TMLEG,					"TMLEG"					},
		{	ARM_SPREADLEG,				"SPREADLEG"				},
		{	ARM_SPREADOPTION,			"SPREADOPTION"			},	
		{	ARM_CORRIDORLEG,			"CORRIDORLEG"			},
		{	ARM_REVERSESTICKYLEG,		"REVERSESTICKYLEG"		},
		{   ARM_SWAP,					"SWAP"					},
		{   ARM_SWAPTION,				"SWAPTION"				},
		{   ARM_SWAPTION_CAPFLOOR,		"SWAPTION_CAPFLOOR"		},
		{   ARM_BARRIER,				"BARRIER"				},
		{   ARM_CAPFLOOR,				"CAPFLOOR"				},
		{   ARM_RATCHET,				"RATCHET"				},
		{   ARM_STICKY,					"STICKY"				},
		{   ARM_CALLSTRUCT,				"CALLSTRUCT"			},
		{   ARM_DIGITAL,				"DIGITAL"				},
		{   ARM_FLEXIBLECAPFLOOR,		"FLEXIBLECAPFLOOR"		},
		{   ARM_RESTRIKABLECAPFLOOR,	"RESTRIKABLECAPFLOOR"	},
		{   ARM_RANGENOTE,				"RANGENOTE"				},
		{   ARM_OPTION,					"OPTION"				},
		{   ARM_OPTIONPORTFOLIO,		"OPTIONPORTFOLIO"		},
		{   ARM_IDXAMORTSEC,			"IDXAMORTSEC"			},
		{   ARM_REVERSE,				"REVERSE"				},
		{   ARM_POWERREVERSE,			"POWERREVERSE"			},
		{   ARM_FOREX,			        "FOREX"			        },
		{   ARM_SUMOPTION,			    "SUMOPTION"			    },
		{   ARM_SMILEDSWAPTION,			"SMILEDSWAPTION"		},
		{   ARM_PORTFOLIO,				"PORTFOLIO"				},
		{   ARM_STRIPOPTION,			"STRIPOPTION"			},
		{   ARM_STRIPDIGITALOPTION,		"STRIPDIGITALOPTION"	},
		{   ARM_FXSPREADSTRIPOPTION,	"FXSPREADSTRIPOPTION"	},
		{   ARM_CORRIDORDBLCONDITION,	"RA2CONDITION"			},

		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};

ARM_TagArrayRow ModelTypeTable[]=
	{
		/* Model Types Flag             Model Type Name */
		{	ARM_MODEL,							"MODEL"						},					
		{	ARM_BSMODEL,						"BSMODEL"					},
		{	ARM_BSSMILEDMODEL,					"BSSMILEDMODEL"				},	
		{	ARM_INFBSSMILEDMODEL,				"INFBSSMILEDMODEL"			},		
		{	ARM_BSCORRMODEL,					"BSCORRMODEL"				},
		{	ARM_DFBSMODEL,						"DFBSMODEL"					},
		{	ARM_YCMODEL,						"YCMODEL"					},
		{	ARM_Y2CMODEL,						"Y2CMODEL"					},
		{	ARM_GYCMODEL,						"GYCMODEL"					},
		{	ARM_DFGYCMODEL,						"DFGYCMODEL"				},
		{	ARM_GYCLSMODEL,						"GYCLSMODEL"				},	
		{	ARM_HWTWOFACTORMODEL,				"HWTWOFACTORMODEL"			},
		{	ARM_TWOFACTORIRTREE,				"TWOFACTORIRTREE"			},
		{	ARM_HWTWOFACTORIRTREE,				"HWTWOFACTORIRTREE"			},
		{	ARM_HWSIGVARTWOFACTOR,				"HWSIGVARTWOFACTOR"			},
		{	ARM_HWSIGVARTWOFACTORMODEL,			"HWSIGVARTWOFACTORMODEL"	},
		{	ARM_PDESOLVER,						"PDESOLVER"					},
		{	ARM_TREESOLVER,						"TREESOLVER"				},	
		{	ARM_CRRTREE,						"CRRTREE"					},
		{	ARM_FRMMARKOVTREE,					"FRMMARKOVTREE"				},				
		{	ARM_RSIRTREE,						"RSIRTREE"					},	
		{	ARM_IRTREE,							"IRTREE"					},
		{	ARM_IRTREEHW,						"IRTREEHW"					},	
		{	ARM_BKIRTREE,						"BKIRTREE"					},
		{	ARM_HWIRTREE,						"HWIRTREE"					},
		{	ARM_DFIRTREE,						"DFIRTREE"					},	
		{	ARM_DFIRTREEHW,						"DFIRTREEHW"				},
		{	ARM_HWSIGCSTTREE,					"HWSIGCSTTREE"				},	
		{	ARM_HWSIGVARTREE,					"HWSIGVARTREE"				},
		{	ARM_DFHWSIGVARTREE,					"DFHWSIGVARTREE"			},
		{	ARM_HULLWHITE1SIGVAR,				"HULLWHITE1SIGVAR"			},	
		{	ARM_DFHULLWHITESIGVAR,				"DFHULLWHITESIGVAR"			},
		{	ARM_DFGYCSIGVARMODEL,				"DFGYCSIGVARMODEL"			},
		{	ARM_GYCSIGVARMODEL,					"GYCSIGVARMODEL"			},
		{	ARM_MONTECARLO,						"MONTECARLO"				},
		{	ARM_HWMONTECARLO,					"HWMONTECARLO"				},
		{	ARM_BS1MONTECARLO,					"BS1MONTECARLO"				},
		{	ARM_BSNMONTECARLO,					"BSNMONTECARLO"				},
		{	ARM_G2YCMODEL,						"G2YCMODEL"					},
		{	ARM_LG2FDIFFUSION,					"LG2FDIFFUSION"				},	
		{	ARM_TREE3F,							"TREE3F"					},
		{	ARM_MCMODEL,						"MCMODEL"					},	
		{	ARM_MCFNHW1FMODEL,					"MCFNHW1FMODEL"				},
		{	ARM_MCFNHW1FSIGVAR,					"MCFNHW1FSIGVAR"			},	
		{	ARM_MCFNHW2FSIGVAR,					"MCFNHW2FSIGVAR"			},
		{	ARM_MCRNLOGDECAL,					"MCRNLOGDECAL"				},
		{	ARM_SMILEDMCRNLDC,					"SMILEDMCRNLDC"				},	
		{	ARM_MCRNQUANTOLDC,					"MCRNQUANTOLDC"				},
		{	ARM_HW1F,							"HW1F"						},
		{	ARM_LOGDECAL,						"LOGDECAL"					},	
		{	ARM_LOGDECALANA,					"LOGDECALANA"				},
		{	ARM_SMILEDLDCANA,					"SMILEDLDCANA"				},
		{	ARM_FRM,							"FRM"						},
		{	ARM_FRM_ANA,						"RM_ANA"					},
		{	ARM_SMCFRM,							"SMCFRM"					},
		{	ARM_BMCFRM,							"BMCFRM"					},
		{	ARM_FRM_TREE,						"FRM_TREE"					},
		{	ARM_BSFRM_TREE,						"BSFRM_TREE"				},
		{	ARM_FRMMODEL,						"FRMMODEL"					},
		{	ARM_FRMMCMODEL,						"FRMMCMODEL"				},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};


ARM_TagArrayRow ReplicModeTable[]=
	{
		/* Model Types Flag             Model Type Name */
		{	K_CONST_STEP,							"CONST_STEP"						},
		{	K_CONST_PREC,							"CONST_PREC"						},
		{	K_ALTN_PREC,							"ALTN_PREC"							},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};

ARM_TagArrayRow ReplicStopModeTable[]=
	{
		/* Model Types Flag             Model Type Name */
		{	K_STOP_PRICE,						"PRICE"								},
		{	K_STOP_PRICE_WEIGHT,				"PRICE_WEIGHT"						},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};

ARM_TagArrayRow LongShortTable[]=
	{
		/* Model Types Flag             Model Type Name */
		{	K_RCV,								"LONG"								},
		{	K_PAY,								"SHORT"								},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};

ARM_TagArrayRow MaturityCapModeTable[]=
	{
		/* Model Types Flag             Model Type Name */
		{	K_INFINE,							"INFINE"							},
		{	K_TRI,								"TRI"								},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};

ARM_TagArrayRow MaturityCapCalibModeTable[]=
	{
		/* Model Types Flag             Model Type Name */
		{	K_MCAP_ATM,						"ATM"									},
		{	K_MCAP_FLAT,					"FLAT"									},
		{	K_MCAP_EX_BOUNDARY,				"EX_BOUNDARY"							},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};

ARM_TagArrayRow CompoundingSpreadTable[]=
	{
		/* Model Types Flag             Model Type Name */
		{	K_COMP_NONE,					"NONE"									},
		{	K_SPREAD_INC,					"SPREAD INCLUDED"						},
		/* do not forget to put this following end of line
		used as an equivalent of '\0' to stop a loop */
		{	ARM_PARAMVIEW_ENDOFLINE, ARM_PARAMVIEW_ENDOFLINE_CHAR }
	};
			
/* table to get the name that corresponds to a constant */
ARM_ParamMappingRow ARM_ParamView::itsParamMappingTable[ARM_PARAMTABLESIZE] =
{
		S_SWAPLEG_SCHEDULE,			SwapLegScheduleTable,
		S_GAUSSIAN_GENERATOR,		GaussianGeneratorTable,
		S_MC_GENERATOR,				MonteCarloGeneratorTable,
		S_OPTION_TYPE,				OptionTypeTable,
		S_CAPLIKE_TYPE,				CapLikeTypeTable,
		S_LEGTYPE,					SwapLegTypeTable,
		S_EXERCISE_STYLE,			ExerciseStyleTable,
		S_COMPUT_TYPE_RESULT,		ComputingTypeResultTable,
		S_BARRIER,					BarrierTable,
		S_FORWARD_RULES,			FwdRulesTable,
		S_INTEREST_RULES,			InterestRulesTable,
		S_STUB_RULES,				StubRulesTable,
		S_AMORTIZATION_TYPE,		AmortizationTypeTable,
		S_NOMINAL_EXCHANGE_TYPE,	NominalExchangeTypeTable,
		S_FREQUENCY,				FrequencyTable,
		S_DATE_TYPES,				DateTypesTable,
		S_CALCULATION_METHOD,		CalculationMethodTable,
		S_INTERPOL_METHOD,			InterpolationMethodTable,
		S_EXTRAPOLATION_METHOD,		ExtrapolationMethodTable,
		S_MISSING_MATURITY_FILLING,	CurveModTable,
		S_MARKET_SEGMENT,			CurvePrioritiesTable,
		S_INDEX_TYPES,				IndexTypesTable,
		S_RECEIVE_PAY,				ReceiveOrPayTable,
		S_RATE_COMPOUNDING,			RateCompoundingMethodsTable,
		S_TIMING_MOD,				TimingModTable,
		S_WEEKDAY,					WeekDayTable,
		S_DAYCOUNT,					DayCountTable,
		S_MC_MOD,					MonteCarloModTable,
		S_FRM_AUTOCALIBRATION_MOD,	FRMAutoCalibrationModesTable,
		S_CALIBRATION_TYPE,			CalibrationTypeTable,
		S_CALIBRATION_METHOD,		CalibrationMethodTable,
		S_FRMMARKOVTREE_DRIFTMOD,	FRMMarkovTreeDriftModeTable,
		S_SPACING,					TreeSpacingTable,
		S_SECURITY_TYPES,			SecurityTypeTable,
		S_MODEL_TYPES,				ModelTypeTable,
		S_REPLIC_MODE_TYPES,		ReplicModeTable,
		S_REPLIC_STOP_MODE_TYPES,	ReplicStopModeTable,
		S_LONG_SHORT_TYPES,			LongShortTable,
		S_MCAP_MODE_TYPES,			MaturityCapModeTable,
		S_MCAP_CALIB_MODE_TYPES,	MaturityCapCalibModeTable,
		S_COMPOUNDING_SPREAD_TYPES,	CompoundingSpreadTable,

};




/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


