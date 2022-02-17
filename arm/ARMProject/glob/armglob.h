
#ifndef _ARMGLOB_H
#define _ARMGLOB_H

#include "firsttoinc.h"
#include <stdio.h>
#include <limits.h>
#include "armdef.h"

#ifdef __cplusplus

#include <string>

using std::string;

#endif


#define PI 3.1415926535

#define GETDEFAULTVALUE    -1111
#define GETDEFAULTVALUESTR "DEFAULT"


#define ARM_OPTIMIZ_BUF_DIAG_SIZE 255

#define CV_HISTORY_NAME    "MO"

#define NEAREST_INT(x) ( (x-floor(x)) < (ceil(x)-x) ? floor(x) : ceil(x) )


#define MEMCPY(t, src, sz)    memcpy((char *) t, (char *) src, sz)
#define MEMSET(t, c, sz)    memset(t, c, sz)



/* to avoid warning of redefinition
 * we have take the same definiton 
 * as the one in ournag.h 
 * therefore, with the max superior or equal
 */

#ifndef _SHOW_MAXMIN_WARNING
	#ifndef MIN
		#define MIN(x,y)   (x < y ? x : y)
	#endif

	/* this is superior or equal */
	#ifndef MAX
		#define MAX(x,y)    (x >= y ? x : y) 
	#endif

	#ifndef ROUND
		#define ROUND(x)    (int)(floor((x)+0.5)) // returns nearest integer
	#endif

/* the old stuff
 * with the max min non completely
 * compliant to the NAG def
 */


#else
	#define MIN(x,y)    ( (x) < (y) ? (x) : (y))
	#define MAX(x,y)    ( (x) < (y) ? (y) : (x))
	#define ROUND(x)    ((int) floor((x) + 0.5))  // returns nearest integer

#endif

#define POS(a)       MAX(a, 0.0)              // returns positive part of a
#define SQR(x)        ((x)*(x))



typedef char ARM_CRV_TERMS[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];





typedef enum
{
    ARM_NOR_MOD,
    ARM_MOD_LOGNOR,
    ARM_MOD_SQR
}
ARM_CALC_MOD_TYPE;
 
 
typedef enum
{
    ARM_YIELD,
    ARM_VOL
}
ARM_YIELD_VOL;




typedef enum _COUNTRY_SYMBOL
{
    ARM_ZAR,
    ARM_HUF,
    ARM_PLN,
    ARM_ATS,
    ARM_BEF,
    ARM_CAD,
    ARM_CHF,
    ARM_DEM,
    ARM_DKK,
    ARM_ESP,
    ARM_EUR,
    ARM_FRF,
    ARM_GBP,
    ARM_GRD,
    ARM_ITL,
    ARM_JPY,
    ARM_NLG,
    ARM_NOK,
    ARM_PTE,
    ARM_SEK,
    ARM_USD,
    ARM_XEU,
    ARM_IEP,
    ARM_FIM,
    ARM_AUD,
    ARM_CZK,
    ARM_HKD,
    ARM_SGD,
    ARM_PHP,
    ARM_NZD,
    ARM_SKK,
    ARM_CNY,
    ARM_TWD,
    ARM_KRW,
    ARM_INR
}
ARM_COUNTRY_SYMBOL;


typedef struct _COUNTRY_SYMBOL_ISO
{
    ARM_COUNTRY_SYMBOL symb;

    char isoName[4];
}
ARM_COUNTRY_SYMBOL_ISO;


typedef enum _ARM_PRICER_TYPE
{
    PT_NONE,
    PT_FNMONTECARLO,
    PT_RNMONTECARLO,
    PT_TREE,
    PT_HWTREEQUANTO,
    PT_EXPLICITE,
    PT_TWOPLUSTREE,
    PT_SCRMONTECARLO,
    PT_MTree,
    PT_MARKOVTREE,
    PT_MCCREDIT,
    PT_MCPREPAY,
    PT_MCANTICREDIT
}
ARM_PRICER_TYPE;


typedef enum _PRODUCT_TYPE
{
    PT_YIELD,
    PT_IRG,
    PT_SWOPT,
    PT_SWOPT_IRG,
    PT_IRG_SWOPT,
    PT_SWOPT_SWOPT
}
ARM_PRODUCT_TYPE;
 
 
typedef enum _GEN_MATU
{
        GM_1M,
        GM_2M,
        GM_3M,
        GM_4M,
        GM_5M,
        GM_6M,
        GM_9M,
        GM_12M,
        GM_18M,
        GM_1Y,
        GM_2Y,
        GM_3Y,
        GM_4Y,
        GM_5Y,
        GM_6Y,
        GM_7Y,
        GM_8Y,
        GM_9Y,
        GM_10Y,
        GM_12Y,
        GM_15Y,
        GM_20Y,
        GM_25Y,
        GM_30Y,
        GM_40Y
}
ARM_GEN_MATU;
 
 
typedef int ARM_Booleen;


/*---- Global Tracing variables ----*/

extern int gTrace;


typedef char ARM_Country[10];

extern ARM_Country ARM_DEFAULT_COUNTRY;


extern void InitCurveName(char* cvName);
extern char* ARM_GetDefaultCountry(void);
extern char* ARM_SetDefaultCountry(ARM_Country newCountry);

extern ARM_COUNTRY_SYMBOL ARM_GetCountrySymbol(char* isoName);
extern char* ARM_GetIsoFromCountrySymbol(ARM_COUNTRY_SYMBOL isoSymbol);

extern void CUSTOMISE_CURRENCY(void);

extern void CUSTOMISE_COUNTRY(void);


extern void CUSTOMISE_ENV(void);

#define ARM_DBG_BEGIN(funcName, curModuleTrace) \
    if ((curModuleTrace) || (gTrace)) \
       fprintf(stderr, "\n -----> IN : %s : %s \n", __FILE__, funcName)



#define ARM_BEGIN(curModuleTrace) ARM_DBG_BEGIN(__ARM_FUNC__, curModuleTrace) 


    
#define ARM_DBG_END(funcName, curModuleTrace) \
    if ((curModuleTrace) || (gTrace)) \
       fprintf(stderr, "\n <----- OK : %s : %s \n\n\n", __FILE__, funcName)

#define ARM_END(curModuleTrace) ARM_DBG_END(__ARM_FUNC__, curModuleTrace)

 

extern void ARM_ERROR(char* funcName, int traceModule,...); 

extern void ARM_PRINTF(char* funcName, int traceModule,...);

extern char* ARM_UPPERCASE(char* str);

extern int ARM_GetCalendarFile(char* fileName);

extern void ARM_GetViewFile(char* fileName, char* id, char* fOutName);

extern FILE* ARM_GetOutputHomeFile(char* fileName, char* mode);

extern void ARM_GetTmpAbsFile(char* fileName, char* fOutName);



/* Index types which involve special calculations */


typedef enum _ARM_INDEX_TYPE
{
    IDXFIXED,
    LIBOR1M,
    LIBOR2M,
    LIBOR3M,
    LIBOR6M,
    LIBOR1Y,
    PIBOR1M,
    PIBOR2M,
    PIBOR3M,
    PIBOR6M,
    PIBOR1Y,
    EURIBOR1M,
    EURIBOR2M,
    EURIBOR3M,
    EURIBOR6M,
    EURIBOR1Y,
    CMT1,
    CMT2,
    CMT5,
    CMT10,
    CMT15,
    CMT20,
    CMT30,
    CMS1,
    CMS2,
    CMS3,
    CMS4,
    CMS5,
    CMS6,
    CMS7,
    CMS8,
    CMS9,
    CMS10,
    CMS11,
    CMS12,
    CMS13,
    CMS14,
    CMS15,
    CMS16,
    CMS17,
    CMS18,
    CMS19,
    CMS20,
    CMS21,
    CMS22,
    CMS23,
    CMS24,
    CMS25,
    CMS26,
    CMS27,
    CMS28,
    CMS29,
    CMS30,
    TEC5,
    TEC10,
    T4M,
    T4M_FIXED,
    TAM,
    TAG,
    EONIA, // or TMP
	LIVRETA,
	EUR3M,
	EUR12,
	EUR1M
}
ARM_INDEX_TYPE;

#define TMP EONIA // Changed when needed for UNIX/ARM

typedef enum _ARM_INDEX_STYLE
{
    VANILLA,
    CMS,
    CMT,
    TEC,
    IN_ADVANCE,
    IN_ARREARS,
    AVERAGE
}
ARM_INDEX_STYLE;


/* for inflation 
 * ne pas changer l'ordre car important
 * en effet une fonction avec un numero plus eleve
 * inclut les fonctionalites precedentes
 */

typedef enum _ARM_INF_PRICING_INFO
{ 
	PRICE_NOTHING = 0,
	PRICE_FWD_CPI, 
	PRICE_FWDNOPTION 
} ARM_INF_PRICING_INFO;


/* Market data types for model calibration */

typedef enum _ARM_MARKET_DATA_TYPE
{
    PRICE,
    RATE,
    VOLATILITY,
    CORRELATION
}
ARM_MARKET_DATA_TYPE;



/* Market data types for model calibration */

typedef enum _ARM_FITTING_METHOD
{
    QUADMIN,
    CALIBZERO
}
ARM_FITTING_METHOD;

/* types of Calibration Portfolio */

typedef enum {
 
  GLOBATM =  0,

  DIAGATM =  1,

  DIASMILE = 2


} ARM_Calibration_PORTTYPE;


typedef enum {

  LAUNCHCALAGE = 0,

  INPUTOFTHEMODEL = 1,

  COPYFROMTHESECURITY = 2 

} ARM_DFHWSigVarTree_CalibrationAction;



typedef enum {
  MIN_STRIKE_ATM_CAL = 0,
  STRIKE_CAL = 1, 
  FWD_CAL = 2,
  MAX_STRIKE_ATM_CAL = 3
} ARM_Calibration_TYPE;



typedef enum {

    NOJUMP = 0,
    MAXMINJUMP = 1,
    PREVIOUSTONEXTJUMP = 2,
    BOTHJUMP =3

} ARM_Calibration_CALAGEQUALITY;



/*---- Bond Future type ----*/

typedef enum _BDF_Type
{
    BDF_NOTIONNAL,
    BDF_BUND,
    BDF_GILT
}
ARM_BDFTYPE;



/*---- IR future ----*/

typedef enum _IRF_TYPE
{
    IRF_EUROMARK,
    IRF_EUROLIRE,
    IRF_EUROSWISS,
    IRF_STERLING
}
ARM_IRF_TYPE;



/*---- Market type ----*/

typedef enum _ARM_MARKET
{
    MATIF, /* France  */
    LIFFE, /* GB      */
    DTB    /* Germany */
}
ARM_MARKET;


#ifdef __cplusplus

class ARM_Vector;

#ifndef ARM_PARAM
#define ARM_PARAM ARM_Vector
#endif

#endif


#ifndef NULL
#define NULL 0
#endif


/*!
 * List of all the derived object from ARM_OBJECT
 * in short all the ARM_OBJECT
 */

typedef enum
{
    ARM_OBJECT,
    ARM_ABSTRACTMARKETCLASS,
	ARM_SCALAR_DATA,
    ARM_VAL,
    ARM_DATE,
    ARM_GEN_MATRIX,
    ARM_VECTOR,
    ARM_MATRIX,
    ARM_TDIAG,
	ARM_BASIS_CURVE,
    ARM_ZERO_CURVE,
    ARM_ZERO_INTERPOLATION,
    ARM_ZERO_LIN_INTERPOL,
    ARM_ZERO_CUBDIFF,
    ARM_ZERO_SPLINES,
    ARM_ZERO_SPLICUB,
    ARM_ZERO_VASICEK,
    ARM_ZERO_FLAT,
	ARM_ZERO_SPLSUM,
    ARM_LIVRETACURVE,
    ARM_SMILE_CURVE,
    ARM_VOL_CURVE,
    ARM_VOL_FLAT,
    ARM_VOL_LIN_INTERPOL,
    ARM_VOL_SPLINE_INTERPOL,
    ARM_VOL_CUBE,
	ARM_VOL_SABR,
	ARM_FX_VOLAT,
    ARM_SPARSE_VOL_CUBE,
    ARM_SECURITY,
    ARM_CURRENCY,
    ARM_BOND,
    ARM_RISKYBOND,
    ARM_RISKYBONDWITHCF,
	ARM_BONDTEC,
    ARM_FOREX,
    ARM_BONDFUTURE,
    ARM_IRFUTURE,
    ARM_IRINDEX,
    ARM_FIXEDINDEX,
    ARM_CMSINDEX,
    ARM_CMTINDEX,
    ARM_SWAPLEG,
    ARM_FIXEDLEG,
    ARM_CMSLEG,
    ARM_CMTLEG,
    ARM_TMLEG,
    ARM_SPREADLEG,
    ARM_SPREADOPTION,
    ARM_CORRIDORLEG,
    ARM_REVERSESTICKYLEG,
    ARM_SWAP,
    ARM_BARRIER_SWAP,
    ARM_SWAPTION,
    ARM_SWAPTION_CAPFLOOR,
    ARM_BARRIER_SWAPTION,
    ARM_BARRIER,
    ARM_CAPFLOOR,
	ARM_CORRIDORDBLCONDITION,
    ARM_RATCHET,
    ARM_STICKY,
    ARM_CALLSTRUCT,
    ARM_DIGITAL,
    ARM_FLEXIBLECAPFLOOR,
    ARM_RESTRIKABLECAPFLOOR,
    ARM_RANGENOTE,
    ARM_OPTION,
    ARM_FX_OPTION,
    ARM_OPTIONPORTFOLIO,
    ARM_XOPTION,
    ARM_IDXAMORTSEC,
    ARM_REVERSECOUPON,
    ARM_REVERSE,
    ARM_POWERREVERSE,
    ARM_MODEL,
    ARM_BSMODEL,
    ARM_BSPRICINGMODEL,
    ARM_BSSMILEDMODEL,
    ARM_INFBSSMILEDMODEL,
	ARM_BSCORRMODEL, 
    ARM_Q_MODEL,
    ARM_DFBSMODEL,
    ARM_YCMODEL,
    ARM_Y2CMODEL,
    ARM_GYCMODEL,
    ARM_DFGYCMODEL,
    ARM_GYCLSMODEL,
    ARM_HWTWOFACTORMODEL,
    ARM_TWOFACTORIRTREE,
    ARM_HWTWOFACTORIRTREE,
    ARM_HWSIGVARTWOFACTOR,
    ARM_HWSIGVARTWOFACTORMODEL,
    ARM_PDESOLVER,
    ARM_TREESOLVER,
    ARM_CRRTREE,
    ARM_FRMMARKOVTREE,
    ARM_FRMNODE,
    ARM_FRMSLICE,
    ARM_RSIRTREE,
    ARM_IRTREE,
    ARM_IRTREEHW,
    ARM_BKIRTREE,
    ARM_HWIRTREE,
    ARM_DFIRTREE,
    ARM_DFIRTREEHW,
    ARM_HWSIGCSTTREE,
    ARM_HWSIGVARTREE,
    ARM_DFHWSIGVARTREE,
    ARM_HULLWHITE1SIGVAR,
    ARM_DFHULLWHITESIGVAR,
    ARM_DFGYCSIGVARMODEL,
    ARM_GYCSIGVARMODEL,
    ARM_MONTECARLO,
    ARM_HWMONTECARLO,
    ARM_BS1MONTECARLO,
    ARM_BSNMONTECARLO,
    ARM_G2YCMODEL,
    ARM_LG2FDIFFUSION,
    ARM_TREE3F,
    ARM_LISTEVAL,
    ARM_EXERCISE_STYLE,
    ARM_REFERENCE_VALUE,
    ARM_IAREFVAL,
    ARM_PAYOFF,
    ARM_PORTFOLIO,
    ARM_CONTAINER,
    ARM_BUCKET_SHIFT, 
    ARM_STRUCTURE,
    ARM_GENTREE,
    ARM_1DLATTICE,
    ARM_2DLATTICE,
    ARM_CONSTTREE,
    ARM_2DCONSTLATTICE,
    ARM_2DTDEPLATTICE,
    ARM_GENITER,
    ARM_1DITER,
    ARM_2DITER,
    ARM_HWSIGCSTITER,
    ARM_HWSIGVARITER,
    ARM_DFHWSIGVARITER,
    ARM_MCMODEL,
    ARM_MCFNHW1FMODEL,
    ARM_MCFNHW1FSIGVAR,
    ARM_MCFNHW2FSIGVAR,
    ARM_MCRNLOGDECAL,
    ARM_SMILEDMCRNLDC,
    ARM_MCRNQUANTOLDC,
    ARM_HW1F,
    ARM_LOGDECAL,
    ARM_LOGDECALANA,
    ARM_SMILEDLDCANA,
    ARM_CALIBRATION,
    ARM_MINIMIZERONPF,
    ARM_HWSVCALIBRATION,
    ARM_FRMCALIBRATION,
    ARM_LOGDECCALIBRATION,
    ARM_FRM,
    ARM_FRM_ANA,
    ARM_SMCFRM,
    ARM_BMCFRM,
    ARM_FRM_TREE,
    ARM_BSFRM_TREE,
    ARM_FRM_LATTICE,
    ARM_FRMCALIBRATOR,
    ARM_FRM_FLAT_CAL,
    ARM_FRM_INTERP_CAL,
    ARM_SHAPEVECTOR,
    ARM_FORWARDS,
    ARM_NOISEMATRIX,
    ARM_FAUREGENERATOR,
    ARM_FRM_TREE_NODE,
    ARM_MATCAPFLOOR,
    ARM_MULTIINDEX,
    ARM_MINDEXESLEG,
    ARM_CAPFLOORST,
    ARM_EXOTICSTRUCT,
    ARM_REVERSEFLOATERSWAP,
    ARM_EPOCHS,
    ARM_OPTIONALACCRUALZC,
    ARM_FLEXACCRETSWAPTION,
    ARM_MARKET_DATA,
    ARM_YCS,
    ARM_CYCS,
    ARM_IRMA,
    ARM_GRID_SPEC,
    ARM_CURVE_SHIFT,
    ARM_TK_LEG,
    ARM_TK_OPT,
    ARM_QUANTO_SIDE,
    ARM_QUANTO_COMMON,
    ARM_TKCMS_SPREAD_SWAP,
    ARM_TKCMS_SPREAD_OPTION,
    ARM_FRMMARKOVVOL,
    ARM_FRMHWVOL,
    ARM_FRMMODEL,
    ARM_FRMMCMODEL,
    ARM_FRMSTRUCTURE,
    ARM_FRMVOL,
    ARM_QUASINEWTON,
    ARM_DATESTRIP,
    ARM_FIXZC,
	ARM_CORRELMATRIX,
	ARM_CORRELMANAGER,
	ARM_CORRELATORMANAGER,
	ARM_DUALCAP,
	ARM_CROSSMODEL,
	ARM_FRMMODELMIXTURE,
	ARM_FRMMCMODELMIXTURE,
	ARM_CONVADJ,
    ARM_CALIBRATORSFRM,
    ARM_TRIBSMODEL,
    ARM_TRIBSDUALMODEL,
	ARM_SUMOPTION,
	ARM_SMILEDSWAPTION,
	ARM_HYPER_CUBE,
	ARM_INDEX_INDEX_CORREL_CUBE,
	ARM_SECURITY_FLOWS,
	ARM_STRIPOPTION,
	ARM_STRIPDIGITALOPTION,
	ARM_FXSPREADSTRIPOPTION,
	ARM_GLOBALCAP,
	ARM_TRIXBSMODEL,
	ARM_FIXING_SCHED,

	/* MERCURE SECTION	*/
	ARM_MERCURE_RESULT,
	ARM_MERCURE_MARKETDATAMANAGER,
	ARM_MERCURE_HELP,

    /* INFLATION SECTION */
    ARM_INFCURV,
    ARM_INFIDX,
    ARM_INFLEG,
	ARM_INFCURVMODEL,
	ARM_RESETMANAGER,
	ARM_INFCAPFLOOR,
	ARM_INFCAPFLOORRIELYIELD,
	ARM_INFHYBRIDDIGITAL,
	ARM_INFBSMODEL,
	ARM_SEASONMANAGER,
	ARM_INFSWAPTION,
	ARM_INFMULTIBSMODEL,
	ARM_INFSWWOPTVOL_PRODUCER,
	ARM_INFSWWOPTVOL_STD_PRODUCER,
	ARM_INFSWWOPTVOL_EQUALWEIGHT_PRODUCER,
	ARM_INFSWWOPTVOL_SIMPLEWEIGHT_PRODUCER,
	
	/* GENERIC PRICER SECTION */
	ARM_DEALDES,
	ARM_GENPRICER,
	ARM_GENSECURITY,
	ARM_MODELNAMEMAP,
	ARM_MODELPARAM,
    ARM_CALIBMETHOD,
    ARM_CALIBMETHODS,
	ARM_VANILLADENSITY,
	ARM_DENSITYFUNCTOR,
	ARM_MODELFITTERDES,
	ARM_PRICINGMODEL,
	ARM_IRFWDMOD,
	ARM_NUMMETHOD,
	ARM_BINUMMETHOD,
	ARM_FINUMMETHOD,
	ARM_NUMERAIRE,
	ARM_HW1FMODEL,
	ARM_HW2FMODEL,
	ARM_TREEMETHOD,
	ARM_GRAMHELPER,
	ARM_SFRMMODEL,
	ARM_MCMETHOD,
	ARM_AMCANDERSEN,
	ARM_AMCLONGSTAFFSCHWARTZ,
	ARM_MCSCHEME,
	ARM_GRAMFCTORARGDICT,
	ARM_PRICERINFO,
	ARM_RANDOMGEN,
    ARM_MARGINFWDMODEL,
    ARM_BASISFWDIRMODEL,
	ARM_GENERIC_CURVE,
	ARM_EQSHIFTEDLOG_MODEL,
	ARM_EQHYBRID_MODEL,
	ARM_QGM1F_MODEL,
	ARM_QGM2F_MODEL,
	ARM_GP_VECTOR,
	ARM_GP_MATRIX,
	ARM_GP_TENSOR,
	ARM_TRIGO_CORRELMAT,
	ARM_Q1F_MODEL,
	ARM_GP_VANILLAARG,

	/* GENERIC PRICER CALCULATOR SECTION */
	ARM_GENCALCULATOR,
	ARM_GENCALIBRATOR,
	ARM_MKTDATAMANAGER,
	ARM_CALLREVFLOATER,
	ARM_TARN,
	ARM_CAPTION,
	ARM_DATESTRIP_COMBINER,
	ARM_ERRVIEWER,
	ARM_EVENT_VIEWER,
	ARM_MATURITYCAP,
	ARM_TARN_SNOWBALL,
	ARM_CALLABLE_SNOWBALL,
	ARM_CALLABLE_SPREADOPTION,
	ARM_CALLABLE_CORRIDOR_SPREADOPTION,
	ARM_STRIPPER,
	ARM_GLOBALCAP_CALCULATOR,
	ARM_VOLBOND_CALCULATOR,
	ARM_FXVANILLACALCULATOR,

	/* REPLICATION SECTION */
	ARM_BSCONVADJUST,
	ARM_REPLICCONVADJUST,
	ARM_BSCONVADJUSTREP,
	ARM_MAPCONVADJUST,
	ARM_PAYOFFFORM,
	ARM_FORMOPLET,
	ARM_FORMOPLETSQR,
	ARM_FORMSWAPLET,
	ARM_FORMSWAPLETSQR,
	ARM_FORMSWAPTION,
	ARM_FORMINVANNUITY,
	ARM_REPLICPORT,
	ARM_REPLICMOD,
    
	/* NEW CREDIT SECTION */
    ICM_SECURITY,
    ICM_GENCF,
    ICM_LEG,
    ICM_CDS,
    ICM_BOND,
    ICM_FRN,
    ICM_CLN,
    ICM_FTD,
    ICM_NTD,
    ICM_MEZ,
    ICM_ZCCLN,
    ICM_PF,
    ICM_STOCK,
    ICM_STOCKCALLOPTION,
    ICM_CBOPTION,
    ICM_CONVERTIBLE,
    // ICM_BASKET,
    ICM_CBASW,
    ICM_DIFFSEC,
    ICM_DIFFNULL,
    ICM_DIFFSTOCK,
    ICM_DIFFLEG,
    ICM_DIFFNDLEG,
    ICM_DIFFOPTION,
    ICM_DIFFFORCED,
    ICM_DIFFCOMPLEX,
    ICM_CDO2,
	ICM_CMCDO,
	ICM_CREDIT_INDEX,
	ICM_OPTION,
	ICM_COLLATERAL,
	ICM_CAPFLOORCMCDS,
	ICM_LEGS_BASKET,
	ICM_CPPI,
	ICM_CORRIDORLEG,
	ICM_SPREADOPTION,
	ICM_CUSTOMIZED_CDO,
	// ICM_OPTION_INDEX_GEN,
	// ICM_CORRIDOR_INDEX_GEN,
	ICM_LSS_GAP_OPTION,
	ICM_CPDO,
	ICMOPTION_RESTRIKABLE,

	//credit models
    ICM_YCMODELDEFAULT,
    ICM_YCMODELCREDIT,
    ICM_MODELMULTICURVES,
    ICM_PWCSHORTCURVE,
    ICM_PWCDEFAULT,
    ICM_DEFAULTTREEMODEL,
	ICM_MODELMULTICURVES_SMILE,
	ICM_META_MODEL,
	ICM_MKTDATAMNG,
	ICM_CREDIT_MANAGER,
	// ICM_TREE_MODEL,
	// ICM_BINOMIAL_TREE_MODEL,
	// ICM_TRINOMIAL_TREE_MODEL,
    ICM_CUSTOMIZED_CREDIT_MULTICURVES,
	ICM_IMPLIED_LOSS_TREE,

	//credit pricers
    ICM_PRICER,
    ICM_PRICERDEFAULTLEG,
    ICM_PRICERDEFAULTTWOLEGS,
    // ICM_PRICERTREE,
    ICM_MCPRICER,
    ICM_MC_PRICER_CDO2,
    ICM_DIFFUSIONPRICER,
    ICM_CLDPRICERCDO2,
    ICM_PRICERSECURITY,
    ICM_PRICERHOMOGENEOUS,
    ICM_PRICER_CDS,
	ICM_PRICER_OPTION,
	ICM_PRICER_CDS_INDEX,
    ICM_PRICER_ANALYTIC_CDO2,	
	ICM_PRICER_ANALYTIC_CDO2_STRIKE,
	ICM_PRICER_HOMOGENEOUS_SMILE,
	ICM_PRICER_CDSOPTION,
	ICM_PRICER_INDEXOPTION,
	ICM_PRICER_MONTE_CARLO,
	ICM_PRICER_CAPFLOORCMCDS,
	// ICM_TREE_PRICER,
	// ICM_BINOMIAL_TREE_PRICER,
	// ICM_TRINOMIAL_TREE_PRICER,
	ICM_BINOMIAL_TREE_PRICER_CDS,
	// ICM_TRINOMIAL_TREE_PRICER_CDS,
	ICM_PRICER_CPPI,
	ICM_PRICER_CORRIDOR,
	ICM_BINOMIAL_TREE_PRICER_SPREADOPTION,
	//ICM_PRICER_INDEXOPTION_SMILE,
	//ICM_PRICER_INDEXOPTION_SABR,
	// ICM_PRICER_INDEX_CAPGEN_VANILLA,
	//ICM_PRICER_INDEX_CAPGEN_BLACKSCHOLES,
	//ICM_PRICER_INDEX_CAPGEN_BLACKSCHOLES_STRIKES,
	//ICM_PRICER_INDEX_CAPGEN_SMILE,
	ICM_PRICER_GENERIC_MC_CDO,
	ICM_PRICER_GENERIC_MC_CDO2,
	ICM_PRICER_IR_INFCURVE,
	ICM_PRICER_IR_YCMOD,
	ICM_PRICER_HOMOGENEOUS_SMILE_RF,
	ICM_PRICER_VANILLA_FWD,
	ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD,
	ICM_PRICER_CLN,
	ICM_PRICER_TRANCHE_CLN,
	ICM_IMPLIED_LOSS_TREE_PRICER,
	ICM_PRICER_LSS_GAP_OPTION,
	ICM_PRICER_CPDO,
	ICM_BASEPRICER,
	ICM_PRICER_MC_CDO,
	ICM_PRICER_TREE_HW_CDO,
	ICM_PRICER_TRANCHE_RESTRIKABLE,

	ICM_CORRELATION,
	ICM_FLAT_CORRELATION,
    ICM_CORRMATRIX,
    ICM_BETA_CORRMATRIX,
    ICM_SMILE_CORRMATRIX,
    ARM_CDS_ZERO_CURVE,
    ARM_ZERO_PARYIELD,
	ICM_NODE,
	ICM_BINOMIAL_NODE,
	ICM_TRINOMIAL_NODE,
	ICM_SLICE,
	ICM_BINOMIAL_SLICE,
	ICM_TRINOMIAL_SLICE,
	ICM_INDEX_CORRELATION,
	ICM_CORRELATION_SECTOR,

	//credit curves
    ICM_DEFAULTCURVE,
    ICM_CST_PIECEWISE,
    ICM_LINEAR_PIECEWISE,
    ICM_DEFAULT_CURVE_MODEL,
	// ICM_INTERPOL_DEFCRV,
	// _DEFCURVE_INDEX,
	ICM_DEFCURVE_QUICK,
	ICM_FIXING_CURVE,

	//credit utils
	ICM_MATRIX,
	ICM_QMATRIX,
	ICM_QCUBIX,
	ICM_PARAMETERS,
	ICM_LEVERAGE,
	ICM_VOL_CURVE,
	ICM_SCHEDULE_INFO,
	ICM_INDEX_FIXING_DATE,

	ICM_GENERIC_NAME1,
	ICM_GENERIC_NAME2,
	ICM_GENERIC_NAME3,
	ICM_GENERIC_NAME4,
	ICM_GENERIC_NAME5,
	ICM_GENERIC_NAME6,
	ICM_GENERIC_NAME7,
	ICM_GENERIC_NAME8,
	ICM_GENERIC_NAME9
}
ARM_CLASS_NAME;



#define IS_OBJECT_CLASS_OK(obj, oClassName, objName) \
        { \
          char buf[200]; \
          if ( obj == NULL ) \
          { \
             sprintf(buf, "Unknown or NULL Object : %s expected", \
                            objName); \
             return(GET_ERR(ERR_OBJECT_NULL, buf)); \
          } \
          if ( obj->GetName() != oClassName ) \
          { \
             sprintf(buf, \
               "Object is not instance of expected class : %s ",  \
                  objName); \
             return(GET_ERR(ERR_BAD_OBJECT_CLASS, buf)); \
          } \
        }

#define IS_OBJECT_ROOT_CLASS_OK(obj, oClassName, objName) \
        { \
          char buf[200]; \
          if ( obj == NULL ) \
          { \
             sprintf(buf, "Unknown or NULL Object : %s expected", \
                            objName); \
             return(GET_ERR(ERR_OBJECT_NULL, buf)); \
          } \
          if ( obj->GetRootName() != oClassName ) \
          { \
             sprintf(buf, \
               "Object is not instance of expected class or subclass: %s ",  \
                  objName); \
             return(GET_ERR(ERR_BAD_OBJECT_CLASS, buf)); \
          } \
        }

#define IS_CURVE_CLASS_OK(obj, objName) \
        { \
          char buf[200]; \
          if ( obj == NULL ) \
          { \
             sprintf(buf, "Unknown or NULL Object : %s expected", \
                            objName); \
             return(GET_ERR(ERR_OBJECT_NULL, buf)); \
          } \
          if ( obj->GetRootName() != ARM_ZERO_CURVE && \
               obj->GetRootName() != ARM_VOL_CURVE ) \
          { \
             sprintf(buf, \
               "Object is not instance of expected class or subclass: %s ",  \
                  objName); \
             return(GET_ERR(ERR_BAD_OBJECT_CLASS, buf)); \
          } \
        }


    

#ifdef __cplusplus

namespace ARM{
class ARM_Currency;
};


extern ARM::ARM_Currency* ARM_DEFAULT_CURRENCY;

extern ARM::ARM_Currency* ARM_FRF_CURRENCY;


template<class T>
class ARM_Pointer
{
    private :

        T*  itsPointer;

    public :

    ARM_Pointer(void)
    {
        itsPointer = NULL;
    }

    ARM_Pointer(T* in)
    {
        itsPointer = in;
    }


    inline T* operator ->(void)
    {
        return itsPointer;
    }

    inline ARM_Pointer& operator = (const ARM_Pointer& in)
    {
        itsPointer = in.itsPointer     ;
    }

    inline ARM_Pointer& operator = ( T* in)
    {
        itsPointer = in;
    }

    inline T* GetPointer(void)
    {
        return itsPointer;
    }

    inline void SetPointer(T* in)
    {
        itsPointer = in;
    }


    inline T operator * (void)
    {
        return *itsPointer
    }

    inline void free(void)
    {
        cleanit(itsPointer);
    }

};


/*!
 * ARM_Object is the base of any object in ARM
 * 
 * methods to be supported are
 *		- BitwiseCopy copy the attributes of its argument
 *		- Copy that create an exactly similar object to its argument src,
 *			hence creates a new object and called BitwiseCopy
 *		- Clone that copies this as a new object, 
 *			hence implemented in terms of Copy
 *		- View that is used to print in the ARM menu in XL
 *		- Print method used for console application
 *
 * Note that the name given to a derived class of an ARM_Object
 * has to be registered in ARM_CLASS_NAME
 */

class ARM_Object
{
    private:
	
	    ARM_CLASS_NAME name;
	
    public:

	/*
	 * Because the default constructor
	 * do not take a name
	 * you are forced to use SetName
	 * to change the name of ARM_OBJECT
	 */
	ARM_Object(void)
	{
		name = ARM_OBJECT;
	}
	
	ARM_Object(const ARM_Object& obj)
	{
		name = obj.name;
	}

	/*!
	 * because we force this to be virtual
	 * all the derived class have a virtual
	 * destructor even if it is not explicitly
	 * written
	 */
    virtual ~ARM_Object(void)
	{
	}

virtual int IsMarketData(void)
    {
        return(0);
    }

    void BitwiseCopy(const ARM_Object* src)
    {
		name = src->name;
    }
	
	virtual void Copy(const ARM_Object* src)
    {
		ARM_Object* obj = (ARM_Object *) src;
		
		BitwiseCopy(obj);
	}
	
	/*!
	 * Clone is used to implement 
	 * a virtual copy constructor
	 */
    virtual ARM_Object* Clone(void)
	{
		ARM_Object* theClone = new ARM_Object();
		
		theClone->Copy(this);
		
		return(theClone);
	}
	virtual ARM_Object* Clone(void) const
	{
		return (const_cast<ARM_Object*>(this))->Clone() ;
	}
    virtual ARM_Object& operator = (const ARM_Object& obj)
	{
		this->BitwiseCopy(&obj); /// old code style
		
		return(*this);
	};
	
virtual double CalcNumericalObjectSignature(void);

virtual void View(char* id = NULL, FILE* ficOut = NULL)
    {
		FILE* fOut;
		char fOutName[200];
		//char strDate[20];
		
		
		if ( ficOut == NULL )
		{
			ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
			
			fOut = fopen(fOutName, "w");
		}
		else
		{
			fOut = ficOut;
		}
		
		fprintf(fOut, "\n\n >>>>>>>??? No View() method for this Object \n");
		
                if ( ficOut == NULL )
                {
		    fclose(fOut);
                }
    }
	
	virtual ARM_CLASS_NAME GetRootName(void)
	{
		return(name);
	}
	
	inline ARM_CLASS_NAME GetName(void) const 
	{
		return(name);
	}

	/*!
	 * Useful because of the lack of constructor that
	 * thakes no name as argument
	 */
	inline void SetName(ARM_CLASS_NAME nm)
	{
		name = nm;
	}

    virtual void Print(void)
    {
    }
};



/*******************************************************************************/










/* located here because this is really C++ */

/*!
 * general constants for all curves
 * that are defined on a 100 basis
 */




#ifdef __ARM_NO_NAMESPACE

#define ARM_Constants 

extern double rateBase;
extern double volBase;
extern double correlBase;


#else

namespace ARM_Constants
{
	static const double rateBase = 100.0;
	static const double volBase	 = 100.0;
	static const double correlBase = 100.0;
}

#endif

#endif 


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
