 /*****************************************************************************
  *    SCCS Keyword Information
  *    ------------------------
  *    Module name    :   enumeration.h
  *    Version        :   2.175
  *    Extracted      :  11/20/97 @ 10:30:50
  *    Last Updated   :   11/20/97 @ 10:30:48
  ****************************************************************************/
/*M***************************************************************************

    CLASS(ES):     

    NOTE(S):       

*****************************************************************************/
#ifndef ENUMERATION
#define ENUMERATION

#ifndef __boolean_h
#include "boolean.h"
#endif

typedef enum{
	CAP,
	FLOOR
} cap_or_floor;

typedef enum{
        USE_DEFAULT_VOLS,
	USE_BASE_VOLS,
	USE_SWAPTION_VOLS,
	USE_DIAG_VOLS
} sticky_vols_to_use;


typedef enum{
	EARLIER,
	LATER,
	ON_OR_EARLIER,
	ON_OR_LATER
} relative_date_type;


typedef enum{
	NO_EXCESS_SERVICING,
	PREVIOUS_EXCESS_SERVICING,
	MAXIMUM_EXCESS_SERVICING
} excess_servicing_type;


typedef enum {
        NO_SORT,
	SORT_ASCENDING,
	SORT_DESCENDING
} sort_order_type;


typedef enum {
	NULL_TRANSACTION,
	BUY,
	SELL
} transaction_type;

typedef enum {
	FI_INST,
	PORT_INST,
	TSY_INST,
	DSCNT_INST,
	MTG_INST,
	ARM_MTG_INST,
	ETP_INST,
	RULE_INST
} instrument_type;


typedef enum {
	RAMP,
	STEP,
	IMPULSE,
	QUADRATIC,
	FLAT
} interpolation_type;

/* SC4.0 chokes on the word GLOBAL hence the name change. */
#ifdef SC4_0
typedef enum {
	GLOBAL_SCOPE,
	LOCAL_SCOPE
} scope;
#else
typedef enum {
	GLOBAL,
	LOCAL
} scope;
#endif


typedef enum {
	DAILY	   = -1,
	WEEKLY     = -2,
	BIWEEKLY   = -3,
	MONTHLY    = -4,
	QUARTERLY  = -5,
	SEMIANNUAL = -6,
	ANNUAL     = -7,
	EVERY2YRS  = -8,
	EVERY3YRS  = -9,
	EVERY5YRS  = -10,
	NEVER	   =  -19,
	IMM	   = -20,
	DISCOUNT   = -21,
	SIMPLE	   = -22,
	IRREGULAR  = -23,
	ACTUAL	   = -24
} frequency_type;


typedef enum {
	BACKWARD,	/* go backward if not business day */
	FORWARD, 	/* go forward if not business day */
	MODFOL, 	/* Modified Following: go forward unless you go */
			/* over the end of month, then go backwards */
			/* e.g., starting from 30 (not bday), then 31 (not */
			/* bday), then use 29 ... */
	NOADJUST	/* Don't do any adjustments */
} date_adjustment;


typedef enum {
	MAXIMUM,
	MINIMUM
} maxmin_type;

typedef enum {
	NON_ECON,
	ECONO,
	MEMBER
} tr_grp_typ;

/* copy type is being phased out, see migration documentation 10/28/91 mtl */

#ifdef __cplusplus
typedef char* DUMMY_TYPE;
#define DUMMY (char*)0
#else
typedef enum {
	SHARED  = -1,
	PRIVATE	= -2
} copy_type;
#endif

typedef enum {
	CACHE_TYPE_FIFO,
	CACHE_TYPE_LRU,
	CACHE_TYPE_LFU
} cache_type;


typedef enum {
    ORIGINAL_STATE,
    ACTUAL_STATE,
    CURRENT_STATE
}   state_type;

#ifdef ETP_READY
typedef enum {
    PRICE,
    YIELD,
    SPREAD,
    OAS,
    OAY,
    NULL_PRICING,
    PV,
    DBMARGIN,
    DISCMGN,
    ROLLSPRD,
    BEYSPRD,
    DISCNT,
    MMKTYLD,
    MATRIX,
    YCMARGIN,
    YTC,
    YTP,
    YTM,
    RTM,
    USWPSPRD,
    AVGLIFE,
    MODDUR,
    SWPSPRD,
    PROCEEDS,
    TSYROLL,
    UNCAPPDM,
    FUTURES,
    TBASPRD,
    TSYOAS,
    LIBOROAS,
    DB_PRICE,
    SWAPTION,
    PTTABLE,
    ADJOAS,
    BOND_OPTION
}   price_yield_type;
#else
typedef enum {
    PRICE,
    YIELD,
    SPREAD,
    OAS,
    OAY,
    NULL_PRICING,
    PV,
    DBMARGIN,
    DISCMGN,
    ROLLSPRD,
    BEYSPRD,
    DISCNT,
    MMKTYLD,
    MATRIX,
    YCMARGIN,
    YTC,
    YTP,
    YTM,
    RTM,
    USWPSPRD,
    AVGLIFE,
    MODDUR,
    SWPSPRD,
    PROCEEDS,
    TSYROLL,
    UNCAPPDM,
    TBASPRD,
    TSYOAS,
    LIBOROAS,
    PTTABLE,
    BOND_OPTION
}   price_yield_type;
#endif

typedef  enum {
    NDX,
    YC,
    NO_SELECTOR
} rate_selector;


typedef enum {
    NO_SOURCE,
    GIVEN,
    CALCULATED
} source_type;


typedef enum {
    ECONOMIC,
    PRICING
} environment_type;


typedef enum {
    INCLUSIVE,
    EXCLUSIVE
}   boundry_type;


typedef enum {
    BEY,
    MTG
}   yield_basis_type;


typedef enum{
	STANDARD_REPORT,
	LOTUS_REPORT,
	BALANCE_REPORT
} report_format;

typedef enum
{
	PRICE_INPUT,
	YIELD_INPUT,
	OAS_INPUT,
	COUPON_INPUT,
	SPREAD_INPUT,
	EFF_YIELD_INPUT,
	S_PRICE_INPUT,
	S_YIELD_INPUT,
	S_OAS_INPUT,
	S_SPREAD_INPUT,
	NULL_INPUT,	/* No analytics_input specified */
	PV_INPUT
} analytics_input;

typedef enum
{
	NEW_DEAL,
	ALL_POOLS,
	WGT_AVG,
	DET_POOLS,
	CNVTL_AVG,
	CLUSTERED
} cmo_coll_mode;

typedef enum{
	GT,
	GE,
	EQ,
	LE,
	LT
} rel_op_type;

typedef enum 
{
	PERIODS,
	DATES
} sched_mode;


# define NO_TRACE	0
# define NOTRACE	0
# define TRACE		1

# define SUCCESS	0
# define WARNING	1
# define PRIMUS_FAILURE	2
# define GTO_SUCCESS	0
# define GTO_FAILURE	-1

typedef enum {
    NO_ROUNDING,
    ROUND_UP,
    ROUND_DOWN,
    ROUND_TO_NEAREST
}  rounding_type;

typedef enum {
	NEED_INTEREST,
	NEED_TOTAL_PRINCIPAL,
	NEED_TOTAL_CASHFLOW
} need_cashflow_type;

typedef enum {
	ACTUAL_DATE,
	FUDGED_DATE
} cmo_state_date_mode;

typedef enum {
	POSTSCRIPT,
	HP,
	COLORPS,
	FAXSVR
} printer_type;

#define MAX_TRANCHES 100

typedef enum {
	ORIGINAL_FACE,
	CURRENT_FACE,
	SHARE_FACE,
	PERCENTAGE_FACE,
	COST_FACE,
	SETTLEMENT_FACE,
	REFORIG_FACE,
	REFCURR_FACE,
	NCONTRACTS_FACE
} FaceAmountType;

typedef enum {
	INVALID_BOND_TYPE,
	T_BOND,
	AGENCY,
	COUP_STRIP,
	PRIN_STRIP,
	REFCORP,
	RFCRP_STRP,
	RFCRP_PRIN,
	T_BILL,
	DSC_AGENCY,
	CORPORATE
} BondIdType;

typedef enum {
	NO_ANALYTICS,
	OPENING_ANALYTICS,
	CLOSING_ANALYTICS,
	OPENING_AND_CLOSING_ANALYTICS
} AnalyticsType;
	
typedef enum {
                LEVEL_CLASS,
                GPM_CLASS,
		TPM_CLASS,
                ARM_CLASS,
		BALLOON_CLASS,
		INVALID
             }  mortgage_type;

typedef enum {
	MOVE_RIGHT,
	MOVE_LEFT,
	MOVE_NEXT,
	MOVE_PREVIOUS
} MoveDirection;

typedef enum {
	NULL_DURATION,
	IMPLIED_DURATION,
	CACHED_DURATION,
	OAS_DURATION,
	MODIFIED_DURATION
} DurationType;

typedef enum {
	AMERICAN,
	EUROPEAN,
	MULTI_EURO
} AmerEuro;

typedef enum {
	NULL_PC,
	PUT,
	CALL,
	P_WORST,
	P_MATURITY
} PutCall;

typedef enum {
	AVLF_TSY,
	BASE_TSY,
	SPFD_TSY
} spread_tsy_type;

#define NULL_SPREAD_TSY_TYPE AVLF_TSY

typedef enum	{
	CREDIT,
	DEBIT
}   account_transaction_type;

typedef enum {
		BARCUNDEF,
		GNMAISF30,
		FNMASF30,
		FHLMCSF30
} barc_type;

typedef enum{
    SIMPLE_CONTEXT,
    PORTFOLIO_CONTEXT,
    DERIVATIVE_CONTEXT
} fi_instrument_context;

typedef enum {
    /*	valid for derivative fi_instrument_context */
    GOTO_COMMAND_INTERPRETER,
    /*  valid for portfolio fi_instrument_context */
    DELETE_OBJECT,
    APPEND_OBJECT,
    GOTO_FIRST,
    GOTO_LAST,
    GOTO_NEXT,
    GOTO_PREVIOUS,
    /*  valid for simple fi_instrument_context */
    GOTO_EXIT
} fi_instrument_edit_command;

typedef enum {
    POLYMORPHIC,
    STANDARD,
    BVALUE,
    ECON,
    QUARTILE
}   fi_instrument_form;

typedef enum   {
   
    LIBOR1MT,   /* ACT/360 SIMPLE */
    LIBOR3MT,	/* ACT/360 SIMPLE */
    LIBOR6MT,	/* ACT/360 SIMPLE */
    LIBOR1YR,   /* ACT/360 SIMPLE */
    CMT1MT,     /* ACTACT  SEMIANNUAL */
    CMT3MT,     /* ACTACT  SEMIANNUAL */
    CMT6MT,     /* ACTACT  SEMIANNUAL */
    CMT1YR,     /* ACTACT  SEMIANNUAL */
    CMT2YR,     /* ACTACT  SEMIANNUAL */
    CMT3YR,	/* ACTACT  SEMIANNUAL */
    CMT5YR,     /* ACTACT  SEMIANNUAL */
    CMT7YR,     /* ACTACT  SEMIANNUAL */
    CMT10YR,    /* ACTACT  SEMIANNUAL */
    CMT30YR,    /* ACTACT  SEMIANNUAL */
    COFI,  	/* FHLB 11th District Cost of Funds (monthly average)	    */
	   	/* 30/360  MONTHLY Source: FNMA numeric subtype:  1.*/ 
    TSY1MT,	
    TSY6MT,     
    TSY10YR,	
    TBILL6MT,
    CMS10YR,
    CMS2YR,
    CMS3YR,
    CMS5YR,
    CMS7YR,
    CMS30YR,
    
    N_INDEX_TYPES  /* This should be the last line in this enumeration */
} index_type;

typedef enum {
	NOSTATUS,
	CURRENT,
	DELIQUENT,
	FORECLOSED
}  mortgage_status;

typedef enum {
    SMM,
    CPR,
    PSA,
    VECTOR_PREPAY, 
    NOPP,	/* For backward compatibility */
    SCPR
} pp_model_type;

typedef enum {
    MATRIX_CPR,
    MATRIX_PSA,
    MGRP,
    JPM,
    ARM_PREPAY,
    USER_PREPAY,
    REGRESSN,
    NULL_PPSRC,
    ADCO,
    FQA,
    DLR_MX_PREPAY,
    DLR_MX_PREPAY_SEAS,
    BASEPP
} PrepaySource;

		/*
		pp_usage determines which prepay model is used by
		a prepaying instrument (e.g., a mortgage) in
		calculating cashflows.  pp_usage just established the
		order in which prepay models are scanned.
		*/
typedef enum {
	LOCAL_YIELD_PREPAY,	
	LOCAL_OAS_PREPAY,
	GLOBAL_YIELD_PREPAY,
	GLOBAL_OAS_PREPAY
} prepay_model_usage;

typedef enum {
	CONST_PYMT,	
	CONST_FCTR,
	CUST_FCTRS
} amortization_mode;

typedef enum {
	CONST_PAYMENT,	
	CONST_FACTOR,
	CUSTOM_FACTORS,
	REFERENCE,
	MAP,
	CMO_SWAP,
	INDEX_AMORT,
	PREPAY_LINKED
} principal_paydown_mode;

typedef enum {
	STANDARD_AMORT,
	STRAIGHT_LINE_AMORT,
	PREPAY_ONLY_AMORT,
	REVERSE_AMORT
} mortgage_amort_mode;

typedef enum {
	FACTORS,
	AVERAGE_LIVES,
        NULL_FACTOR_MAP_TYPE
} factor_map_type;

typedef enum {
	LOGNORMAL,
	TRIMEANRV,
	TRMIMPLPP,
	TWOFACTOR,
	HJM
} oas_model_type;

typedef enum {
        OAS_FLOATER,
        OAS_FIXED,
	OAS_DIAG_VOLS,
	OAS_DIAG_BASE_VOLS,
	OAS_FLAT,
	OAS_FLAT_CMT
} flt_oas_mode;

typedef enum{
	IRS,    /* Interest Rate Swap(vainilla swap)		*/
        PPAS,	/* Prepayment amortizing Swap			*/
        IAS,	/* Index amortizing Swap			*/
	ACAP,   /* Amortizing CAP				*/
	AFLOOR  /* Amortizing FLOOR				*/
} swap_type;

typedef enum 
{
	GNMAI,
	FNMA,
	GOLDPC,
	MIDGET,
	DWARF,
	PC15,
	GP5BAL,
	GP7BAL,
	FN7BAL,
	PC7BAL,
	NPREPAY_AGENCIES /* Should be the last entry in this enum */
} PrepayAgencyType;

typedef enum
{
	NEW_REM,
	SEASONED_REM,
	ALL_REM,
	MOD_REM,
	NREMS /* Should be the last entry in this enum */
} RemRangeType;


typedef enum {
    NULL_ADAPTOR,  /* Should be the first in this enum */
    EVAL_ADAPTOR,
    TXN_ADAPTOR,
    ENV_ADAPTOR,
    PARMS_ADAPTOR,
    PORT_ADAPTOR,
    INST_ADAPTOR
} AdaptorNum;


typedef enum {
    PRICE_FMT,
    DATE_FMT,
    NUMBER_FMT,
    STRING_FMT,
    DOLLAR_FMT
} OutputType;


typedef enum
{
	SEVERE_ICT,
	PREPAY_ICT
} InstChangeType;

typedef enum
{
	DEAL_PPCM,
	POOL_PPCM
} PrepayCollMode;

typedef enum
{
	STEPUP,
	STICKY,
	FWDSTEP,
	STICKYM,
	STICKYA
} PeriodicCapMethod;

typedef enum
{
        JPMDUR,
	USERDUR,
        OASDUR,
	MODELDUR,
	MODFDDUR
} JPMDurType;

typedef enum
{
	NO_LOSS,
	SDA,
	MDR,
	CDR,
	ORIGMDR,
	ORIGCDR,
	PHM
} lossModel;

typedef enum
{
    SAME_AS_OPENING,
    ROLLED_FORWARD
} ClosingYCType;

typedef enum
{
	PREPAY_MAP,	/* Is really CPR_MAP, left for backward compatibility */
	FACTOR_MAP,
	RATE_MAP,
	CPR_MAP,
	PSA_MAP
} IndexMapType;

typedef enum{
	MATRIX_PSA_FILE,
	MATRIX_CPR_FILE,
	REGRESS_PARMS_FILE,
	REGRESS_INITSMM_FILE,
	ARM_PREPAY_FILE,
	DLR_MX_PREPAY_FILE,
	DLR_MX_SEAS_FILE,
	PREPAY_FILES_SIZE	/* Last line in the enumeration */
} PrepayFiles;


typedef enum{
        OAS_PATH_NET_INTEREST,
        OAS_PATH_PV
} oas_path_value;

typedef enum
{
    NO_FACTOR,
    ONE_FACTOR,
    TWO_FACTOR
} FactorModelType;

typedef enum
{
    NO_ZCURVE,
    DISC_ZCURVE,
    INDEX_ZCURVE
} ZeroCurveType;

typedef enum
{
    VOL_LINEAR,
    VOL_CONST_SPOT
} VolInterpType;

typedef enum
{
    FIXED_IO,
    FIXED_PO,
    FLOAT_ONLY,
    SFS_SWAP,
    FIXED_ONLY
} SFSDealPartType;

typedef enum{
	GET_PRIN,
	GET_TOTAL_CF
} prin_or_total_cf;

typedef enum{
	YLDPPEQ,
	AVGPPEQ,
	BALPPEQ
} PPEquivType;

typedef enum{
	SETTLED,
	PENDING,
	UNALLOCATED,
	TO_BE_SETTLED
} SettlementStatus;

typedef enum{
        MARKET,
	MODEL,
	RVG_MARKET,
	RVG_MODEL
} DataSource;

typedef enum{
        FLAT_VOL,
	FLAT_CMT_VOL,    
	BASE_VOL,
	DIAGONAL_BASE_VOL,
	DIAGONAL_CMT_VOL,
        SWAP_VOL
} VolType;

typedef	enum { 
	NYSE = 2, 
	NYBANK = 4, 
	UK = 8,
	JAPAN = 16,
	NO_CALENDAR = 32,
	NYUK = 12	/* NYBANK + UK, here for convenience */
} calendar_type;

typedef enum {
	ED_FUTURE,
	BOND_FUTURE,
	FF_FUTURE,
	NULL_FUTURE,
	LIBOR_FUTURE
} FutureType;

typedef enum {
        ONE_YEAR_HUMP,
	TEN_YEAR_HUMP
} HumpType;

typedef enum {
        OAS_BEHAVE_DUMPRATES = 1,
	OAS_BEHAVE_DUMPVOLS = 2,
	OAS_BEHAVE_DEBUGON = 4
} OASBehaveFlag;

typedef enum {
	SIMPLE_AVG,
	WEIGHTED_AVG,
	COMPOUNDING,
	FLAT_COMPOUNDING
} InterestAvgMethod;

typedef enum {
        TSY_MTG_SPREAD,
	SWAP_MTG_SPREAD
} MtgSpreadType;

#define NUM_YC_PTS 13
#define NUM_PARTIAL_PTS (NUM_YC_PTS + 2)  // 12 and 15 year points as well
#define NSWAP_LIBORS  10  // Overnight, 1M, 2, 3, 4, 5, 6, 9, 12mts, 18mts
#define NSWAP_HUMPS 14
#define NLIBOR_HUMPS NSWAP_LIBORS

#define MAX_NOAS_PATHS_UNDISTRIBUTED   1000
#define MAX_NUMDBL_PER_STRING          350

#endif /* ENUMERATION */
