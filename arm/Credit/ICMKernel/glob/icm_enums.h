#ifndef _ICM_ENUMS_H_
#define _ICM_ENUMS_H_

#include "ICMKernel/util/icm_macro.h"

#include <math.h>
#include <map>


// ***********************************************************************************************
// *							DEFINITION DES ENUMS											 *	
// ***********************************************************************************************


#ifndef MIN
# define MIN(a,b) (((a)<(b))?(a):(b)) 
#endif
#ifndef MAX
# define MAX(a,b) (((a)>(b))?(a):(b)) 
#endif

#define MINIMAL_ERROR 1e-9
#define MINIMAL_INCREMENT 5e-4
#define DB_TOL	1.E-8
#define NUMBER12 12
#define PARVALUE 1
#define NBMAXITER 100
#define K_JUNIOR 0.2
#define K_SENIOR 0.4 
#define ISSUER_UNDEFINE "undefine" 
#define ISSUER_IN_DEFAULT "_DF" 

#ifndef ABS
	#define ABS(x) (x>=0 ? x : -(x)) /* Absolute value */
#endif

#ifndef NEG_SQRT
	#define NEG_SQRT(x)   ((x<0) ? (-sqrt(fabs(x))) : (sqrt(x)))
#endif

#ifndef CHECK_EQUAL
	#define CHECK_EQUAL(x,y)   ((fabs(x-y)<DB_TOL) ? 1 : 0)
#endif

#ifndef CHECK_NULL
	#define CHECK_NULL(x)   ((fabs(x)<DB_TOL) ? 1 : 0)
#endif

#ifndef NAG_PARAM_QUAD
  #define NAG_PARAM_QUAD 32
#endif 

#ifndef _size_zclabel_
#define _size_zclabel_ 60
#endif


#define CREDIT_DEFAULT_VALUE -999
#define NBDECIMAL 10
#define LOSS_UNIT_MIN 50000
#define MAX_SIZE_LOSS_UNIT 10000
#define INTEGRATION_STEP 40
#define DEFAULT_LEVEL_IN_BP 10 //1e5 Bp eq to default

#define TRX_EUR_NBNAMES 125
#define TRX_EUR_RECOVERY 0.4
#define TRX_USD_NBNAMES 125
#define TRX_USD_RECOVERY 0.4

#define DEFAULT_CREDIT_LAG 0
#define DEFAULT_CREDIT_LAG_INDX 30

// #ifndef ARM_NB_TERMS
// #define ARM_NB_TERMS 400
// #endif
// #ifndef ARM_NB_MAX_CHAR_TERMS
// #define ARM_NB_MAX_CHAR_TERMS 30
// #endif

// typedef char typechar[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];

#define BC_MIN_STRIKE		0.03
#define BC_MIN_STRIKE_ZERO  0.05

#define INCLUDE_MATURITY true
#define EXCLUDE_MATURITY false
#define DEFAULT_FRQ_DEFLEG -1

//special interpolation type
#define K_ICM_STEPUP_RIGHT_MATURITY 1001
#define K_ICM_STEPUP_LEFT_MATURITY  1002
#define K_ICM_STEPUP_LEFT_LINEAR	1003

//Types used for Diffusion  
typedef enum ENUM_DIFFUSION_TYPE{
		qDIFFUSION_NONE,		//no diffusion process
		qDIFFUSION_SPREADS_LN,	//Spreads log-normal
        qDIFFUSION_OU,			//Orstein- Uhlenbeck intensity
		qDIFFUSION_SSRJ			//SSRJD intensity
} qDIFFUSION_TYPE; 

//Types used for Cds Schedule generation
typedef enum ENUM_CDS_ADJ{
		/*qCredit_Special_None_Date,	 // NoAdj for date
		qCredit_Special_None_YF,	 // NoAdj for year fractions */
        qCredit_Default,			 // Schedule Cds type Summit
		qCredit_Adjust20,			 //	Quaterly Roll 20 mar/jun/sep/dec
		qCredit_CDSDTRX,			 //	Départ au 20 mar 2005	
		qCredit_CDSDIND,			 //	Départ au 20 jun 2005 
		qCredit_CDSINDZ,			 // Départ au 20 dec 2005 
		qCredit_Adjust20N,			 //	Idem qCredit_Adjust20 + 3M
		qCredit_Adjust20SA,			 //	Semi Annual Roll 20 mar/sep
		qCredit_STDINDEX,			 // Semi Annal Index 20 jun/dec changing mar/sep
		qCredit_CDSINDH5,			 // Idem qCredit_CDSDTRX
		qCredit_CDSINDM5,			 // Idem qCredit_CDSDIND
		qCredit_CDSINDZ5,			 // Idem qCredit_CDSINDZ
		qCredit_CDSINDM6,			 // Départ 20 jun 2006 
		qCredit_CDSINDZ6,			 // Départ 20 dec 2006
		qCredit_CDSINDM4,			 // Départ 20 jun 2004
		qCredit_CDSINDU4,			 // Départ 20 sep 2004
		qCredit_CDSINDZ4,			 // Départ 20 dec 2004
		qCredit_CDSINDM7,			 // Départ 20 jun 2007
		qCredit_CDSINDZ7,			 // Départ 20 dec 2007
		qCredit_CDSINDM8,			 // Départ 20 jun 2008
		qCredit_CDSINDZ8			 // Départ 20 dec 2008
} qCDS_ADJ; 


//Type used for different Sensitivty Computation
typedef enum ENUM_SENSITIVITY
{
		 ICMSPREAD_TYPE,					// Default Curve Type  
		 ICMIRCURVE_TYPE,					// Interest Rate Curve Type  
		 ICMRECOVERY_TYPE,					// Recovery Type
		 ICMCORRELATION_TYPE,				// Correlation Type
		 ICMBETA_TYPE,						// Beta Type
		 ICMRECOVERY_BAR_TYPE,				// UNUSED
		 ICM_RECOVERY_DEFAULT_TYPE,			// UNUSED
		 ICM_ISSUER_DEFAULT,				// Default to Recovery
		 xICM_FAST_SPREAD_TYPE,				// Fast Sensi for spread
		 ICM_SAMECORRELATION_TYPE,			// Unique Correlation
		 ICM_SAMEBETA_TYPE,					// Unique Correlation
		 ICMSPRELSHIFT_TYPE,				// Relative shift to spreads
		 ICMIRCURVE_WITHOUT_DEFCURVE_TYPE,	// Interest Rate Curve Type with no impact on decurve calibration  
		 ICMBETA_WITH_SPREAD_SHIFT,			// Sensitivity au beta avec un shift // au spread
		 ICM_THETA_Type,					// Sensi au passage du temps
		 xICM_FAST_RECOVERY_TYPE,			// Fast Sensi for Recovery
		 ICMCORREL_STRIKE_DOWN_TYPE,		// Correlation Strike DOWN Type
		 ICMCORREL_STRIKE_UP_TYPE,			// Correlation Strike UP Type
		 ICM_DTR_TYPE,						// Default to Recovery (Issuer Default)	
 		 ICM_GREEK_DELTA_TYPE,				// Delta  
		 ICM_GREEK_GAMMA_TYPE,				// Gamma 
		 ICM_GREEK_VEGA_TYPE,				// Vega
		 ICM_GREEK_THETA_TYPE,				// Theta
		 ICM_GREEK_RHO_TYPE,				// Rho
		 ICM_SLN_SIGMA,						// SLN Sigma
		 ICM_SLN_M,							// SLN M
		 ICM_SLN_DELTA,						// SLN delta
		 ICM_SABR_ALPHA,					// SABR Alpha
		 ICM_SABR_NU,						// SABR Nu
		 ICM_SABR_RHO,						// SABR Rho
		 ICM_CORREL_BET_CDO_UP_TYPE,		// Correlation inter Cdo UP
		 ICM_CORREL_BET_CDO_DW_TYPE,		// Correlation inter Cdo DW
		 ICM_INFLATION_CPN_CURVE,			// Inflation coupon curve
		 ICM_INTEREST_CPN_CURVE,			// IR coupon curve
		 ICM_INDX_SPREAD_RESCALING,			// Sensi Spread sur indice for rescaling
		 ICM_CM_SPREAD_PARALLEL,			// Credit Manager: Parallel Shift
		 ICM_CM_DEFAULT,					// Credit Manager: Default
		 ICM_CM_RECOVERY_SENS,				// Credit Manager: Recovery Sensitivity
		 ICM_CM_RECOVERY_LOSS,				// Credit Manager: Recovery Loss
		 ICM_CM_RECOVERY_GLOBAL,			// Credit Manager: Recovery Global
		 ICM_CM_CORRELATION,				// Credit Manager: Correlation
		 ICM_GREEK_VEGA_ATM_TYPE,			// Vega by perturbating the ATM credit vol 
		 ICM_GREEK_IRVEGA_ATM_TYPE,			// Vega by perturbating the ATM ir vol 
		 ICM_IRCURVE_WITH_CPN,				// IR Curve and Cpn Curve if same Index
		 ICM_SUBORDINATION	,				// perturbating Lower and Upper Bounds for a Mez ( parallel shift)
		 //	----------------------
		 //	Add here  new measures

		 //	----------------------
		 //	don't change this item (should be the last)
		 //
		 ICM_SENSILAST		
}qSENSITIVITY_TYPE;

//
//		Enumerates the PRICING MEASURES available in the pricers. 
//		PRICING MEASURES are double values
//
//		Pricer classes handles the result = Price(measure name) 
//		that will compute & cache the results. 
//		
typedef enum ENUM_CMPMETH
{	
		qCMPPRICE ,				//Net Present Value
		qCMPSPREAD,				//Break even Spread
		qCMPFEELEGPV,			//Fee Leg Present value
		qCMPDEFLEGPV,			//Default Leg Present value 
		qCMPDURATION,			//Duration 
		qCMPACCRUED,			//Standard Accrued
		qCMPPREMIUM,			//(JLA) Premium for option products
		qCMPFWDSPREAD,			//(JLA)	FwdSpread for option products
		qCMPFWDDURATION,		//(CC)	FlatDuration for option products
		qCMPFEES,				//Fees Pricing		
		qCMPCORRELUP,			//for tranches
		qCMPCORRELDOWN,
		qCMPAVGCORRDEF,			//for CDO2			
		qCMPEL,					//for expected loss
		qCMPNPV_UPFRONT,		//For Npv with upfront (calibration only)
		//	----------------------
		//	Add here  new measures
		qCMPFLATCORR,
		qCMP_OPT_AN_DELTA,		//for option only : analytic delta from B&S
		qCMP_OPT_AN_VEGA,		//for option only : analytic vega from B&S
		qCMP_OPT_AN_GAMMA,		//for option only : analytic gamma from B&S
		//	----------------------
		//	don't change this item (should be the last)
		qCMPLAST
}qCMPMETH;

//	
//		This is the enumerate for PRICING MEASURES VECTORS, 
//		such as flows etc 
// 
//		Pricer Classes handles PriceVector() function 
//		and internal cache. 
//
typedef enum ENUM_VECTMETH
{
	qSOME,				// sample measure, not used
	qINDEXBASE,			// base vector for index curves
	qMARKETBASE,
	qTRUEMARKETBASE,
	qTRUEMARKETSPREAD,
	//	----------------------
	//	don't change this item (should be the last)
	qVECTLAST
} qVECTMETH ;

typedef enum ENUM_TERM_STRUCTURE {
        qNoTermStructure=0,
		qTermStructure=1,
		qTermStructureR=2
} qTERM_STRUCTURE;  

//		Calibration method for DefaultCurves
//
typedef enum ENUM_DEFCURVE_CALIB_ALGO 
{
	qDEFCURVE_DICHO,
	qDEFCURVE_NEWTON,
	qDEFCURVE_BRENT
} qDEFCURVE_CALIB_ALGO;  

typedef enum ENUM_INDEX_CMPT_METHOD {
        qAVERAGE,
		qHOMOTHETIE

} qINDEX_CMPT_METHOD;  

/*
typedef enum ENUM_CAL_METHOD {
        qDICHOTOMY,
		qNEWTON
} qINDEX_CAL_METHOD;  
*/

typedef enum ENUM_SECURITY_TYPE {
        qRUNNING,	//Standart Credit Default Swap
		qFUNDED,	//Funded Credit Default Swap
		qINFINE,	//Infine Credit Default Swap
		qZEROCOUPON,//Zero Coupon Credit Default Swap
		qCM_TRANCHE,//Constant Matutity CDO
		qCDS_INDEX,//Constant Matutity Credit Default Swap Index
		qCM_CDS,   //Constant Matutity Credit Default Swap
		qBASKETasMEZ //pour Structuration
} qSecurity_TYPE;  

typedef enum ENUM_CREDIT_LEG_TYPE {
		qNone_Leg,					//Unuse Leg
        qRunning_Leg,				//Standart Premium Leg
		qFunded_Leg,				//Funded Premium Leg
		qInfine_Leg,				//Infine Premium Leg
		qZeroCoupon_Leg,			//ZeroCoupon Premium Leg
		qConstant_Maturity_Leg,		//Constant Matrutity Premium Leg
		qStandart_Recovery_Leg,		//Standart Defaultable Leg
		qForward_IRrates_Leg,		//Premium Leg with coupons forward	
		qCMCDS_Leg,					//Constant Matrutity CDS Leg
		qSwapLeg,					//IR Swapleg
		qInflation_Leg,				//ZeroCoupon Premium Leg
		qCorridor_Leg				//Corridor Leg
} qCredit_Leg_Type;  

typedef enum ENUM_CREDIT_LEG_STYLE {
		qStyle_None_Leg,	//Unuse Leg
        qStyle_Recovery_Leg,//Defaultable Leg
		qStyle_Premium_Leg	//Premium Leg
} qCredit_Leg_Style;  

typedef enum COPULA_TYPE {
		qNO_COPULA,		// avoid it
        qGAUSSIAN,
		qSTUDENT
} qCopula_TYPE;  

typedef enum PAYMENT_PREMIUM_LEG {
		qCONTINUE_TO_MATURITY,		//Not risky premium leg
		qACCRUED_NOT_SETTLED,		//Accrued not settled on default (period end risky notional)
		qACCRUED_SETTLED,			//Accrued settled on default (average period risky notional)
		qCOMPLETE_CURRENT_PERIOD	//Full coupon on the default period (period start risky notional)
} qPAYMENT_PREMIUM_LEG;


typedef enum INTEGRATOR_CHOICE {
	qGAUSS_LEGENDRE,
	qGAUSS_HERMITE,
	qTRAPEZE
} qIntegratorChoice;

typedef enum CORREL_BY_STRIKE {
		qStrike_NONE,
		qStrike_UP,
		qStrike_LOW
} qCorrel_By_Strike;

//JLA 
typedef enum UNDERLYING_MATURITY_STYLE 
{
	qConstantMaturity,
	qResidualMaturity
} qUnderlying_Maturity_Style; 

//JLA 
typedef enum KO_STYLE
{
	qKO,
	qNKO
} qKoStyle; 

//JLA 
typedef enum ACCELERATION_STYLE
{
	qACC,
	qNACC
} qAccelerationStyle; 

typedef enum RESCALING_TYPE
{
	qRescal_Std,
	qRescal_Eqty,
	qRescal_Eloss,
	qRescal_Std_Maturity,
	qRescal_Eqty_Maturity,
	qRescal_Eloss_Maturity,
	qRescal_Eqty_NormByEL_Maturity,
	qRescal_Eqty_Mixed_Maturity,
	qRescal_Eqty_Mixed_NormByEL_Maturity,
	qRescal_Eqty_Digital_Maturity,
	qRescal_ING_Maturity
} qRescalType; 

typedef enum ENUM_TYPEGENODATES {
        qACC_START_DATE,
		qACC_END_DATE,
		qPAY_DATE,
		qOBS_START_DATE,
		qOBS_END_DATE,
} qTYPEGENODATES; 

//Types used for Interpolation  
/**
typedef enum ENUM_INTERPOL_TYPE{
        qINTERPOL_CONSTANT, 
		qINTERPOL_LINEAR
} qINTERPOL_TYPE; 
**/ 
typedef enum TWO_FACTORS_CORRELATION_TYPE{
	TFCT_FULL, 
	TFCT_SAME_INTER_DIFF_INTRA,
	TFCT_SAME_INTER_SAME_INTRA
} qTWO_FACTORS_CORRELATION_TYPE;

//Types used for distribution type 
typedef enum ENUM_DISTRIB_TYPE{
        qDISTRIB_STD, 
		qDISTRIB_COLLFWD, 
		qDISTRIB_STD_TSR,
		qDISTRIB_COLLFWD_TSR,
		qDISTRIB_VN,
		qDISTRIB_APPROX,
		qDISTRIB_STD_NEW
} qDISTRIB_TYPE; 


//Types used for Optimisation  
typedef enum ENUM_OPTIMIZE_TYPE{
		qDICHOTOMY,
		qNEWTON,
        qOPT_NAG_NOLIN_BC, 
		qOPT_NAG_NOLIN_BC_BFGS,
		qOPT_LEVMAR,
		qOPT_GLOBAL_UNCON,
		qOPT_GANSO,
		qOPT_PSWARM
} qOPTIMIZE_TYPE; 

typedef enum DEF_RAN_GEN{
	q_NAG , 
	q_RAN1 , 
	q_RAN2,
	q_DEF,
	q_RANMAR,
	q_RNG_STR,
	q_KISS, //  Last Random uniform

	// first  random norm 
	q_INV_CUM_NORM_ACKLAM=10,
	q_INV_CUM_NORM_MORO=11

} qRAN_GEN;


class ICM_EnumsCnv
{
public:
	template<class T> static void cnv(const std::string&name,T&ref)
	{
		std::string myList; 
		bool ret; ICM_EnumsCnv::cnv(name,ref,ret,myList); 
		if (!ret) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_EnumsCnv::cnv: can't convert "<<name<<" : valid values are "<< myList); 
	}
	template<class T> static void toString(T ref,std::string& name) 
	{
		const ICM_Traits<T>::list_t& theList = ICM_Traits<T>::getList(); 
		ICM_Traits<T>::list_t::const_iterator it = theList.begin(); 
		while (it!=theList.end()) 
		{
			if (it->second==ref) { name=it->first; return ; }
			++it; 
		}
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_EnumsCnv::toString: can't find enum value "<<ref ); 
	}
	template<class T> static std::string toString(T ref) 
	{
		std::string name; 
		ICM_EnumsCnv::toString(ref,name); 
		return name; 
	}
	static void cnv(const std::string& item,qAccelerationStyle&out,bool&isOK,std::string&valueList); 
	static void cnv(const std::string& item,qKoStyle&out,bool&isOK,std::string&valueList); 
	static void cnv(const std::string& item,qUnderlying_Maturity_Style&out,bool&isOK,std::string&valueList); 
	static void cnv(const std::string& item,qCDS_ADJ&out,bool&isOK,std::string&valueList);
	static void cnv(const std::string& item,qCMPMETH&out,bool&isOK,std::string&valueList); 
// 	static void cnv(const std::string& item,qINTERPOL_TYPE&out,bool&isOK,std::string&valueList); 
	static void cnv(const std::string& item,qSENSITIVITY_TYPE&out,bool&isOK,std::string&valueList); 
	static void cnv(const std::string& item,qPAYMENT_PREMIUM_LEG&out,bool&isOK,std::string&valueList); 
	static void cnv(const std::string& item,qTWO_FACTORS_CORRELATION_TYPE&out,bool&isOk,std::string&valueList); 
	static void cnv(const std::string& item,qVECTMETH&out,bool&isOk,std::string&valueList); 
	static void cnv(const std::string& item,qDEFCURVE_CALIB_ALGO&out,bool&isOk,std::string&valueList);
	// 
	static void cnv(const int& ,qTERM_STRUCTURE& ,bool& );
	static void cnv(const int& ,qRescalType&,bool&);
	static void toString(qAccelerationStyle,std::string&); 
	static void toString(qKoStyle,std::string&); 
	static void toString(qUnderlying_Maturity_Style,std::string&); 
	static void toString(qCDS_ADJ,std::string&);
	static void toString(qCMPMETH,std::string&);
	static void toString(qRescalType,std::string&);
	static void toInt(qRescalType,int&);
	static void toString(const qRAN_GEN& type ,std::string& RandomName );
	static void cnv(const std::string& name,qRAN_GEN& ref,bool& ret,std::string& list);
	// 
}; 
template <class T> class ICM_Traits
{
public:
	typedef std::map<std::string,T> list_t ; 
	static const list_t& getList(); 
} ;

class ORDER_CHAR
{
	public :
		char* id ;
		double tau ;
		double m_tol;
	public :
		ORDER_CHAR(char* i, double d, double tol=0.05):id(i), tau(d),m_tol(tol){}
		bool operator < (const ORDER_CHAR & rhs) const {return ((tau < rhs.tau) && (fabs(tau - rhs.tau)>m_tol));}
} ;

class ORDER_DOUBLE_WITH_ORDER
{
	public :
		double m_value ;
		double m_tol;
	public :
		ORDER_DOUBLE_WITH_ORDER(double d, double tol=DB_TOL):m_value(d),m_tol(tol){}
		bool operator < (const ORDER_DOUBLE_WITH_ORDER & rhs) const {return ((m_value < rhs.m_value) && (fabs(m_value - rhs.m_value)>m_tol));}
		ORDER_DOUBLE_WITH_ORDER(const ORDER_DOUBLE_WITH_ORDER& value):m_value(value.m_value),m_tol(value.m_tol){}
} ;

template <class T> class ORDER_INDEX_DBL
{
	public :
		double	m_index;
		T		m_value;
		double	m_tol;

	public :
		ORDER_INDEX_DBL(double indx, double value, double tol=DB_TOL):m_index(indx),m_value(value),m_tol(tol){}
		bool operator < (const ORDER_INDEX_DBL & rhs) const {return ((m_index < rhs.m_index) && (fabs(m_index - rhs.m_index)>m_tol));}
} ;

//Types used for Fwd 
typedef enum ENUM_COUPON_TYPE{
        qCPN_FIXED, 
		qCPN_FIXED_AND_FWD, 
		qCPN_FWD,
		qCPN_LAST_FWD
} qCOUPON_TYPE; 


//Types used for Calibration of Index Correlation  
typedef enum ENUM_CALIBRATION_INDEX_TYPE{
        qCAL_BASE_CORRELATION, 
		qCAL_BASE_CORRELATION_TS,
		qCAL_BASE_CORRELATION_LINEAR_MATURITY,
		qCAL_PWC_CORREL,
		qCAL_PWL_CORREL,
		qCAL_BASE_CORRELATION_TSR
} qCAL_INDEX_CORR_TYPE; 

//Types used for Riksy Periods 
typedef enum ENUM_RISK_TYPE{
		qNoRisk, 
		qRiskToDate,
		qStdRisk
} qRISK_TYPE; 

typedef enum CORRELATION_FIT_TYPE{
	qCORREL_1F_SINGLE_VALUE,
	qCORREL_1F_BETAS,
	qCORREL_2F_SAME_INTER_SAME_INTRA,
	qCORREL_2F_SAME_INTER_DIFF_INTRA,
	qCORREL_2F_DIFF_INTER_DIFF_INTRA
} qCORRELATION_FIT_TYPE;


typedef enum BS_CONVEXITY_TYPE{
	qBS_CONVEXITY_STD,
	qBS_CONVEXITY_LRM		// Linear Rate Model
} qBS_CONVEXITY_TYPE;

// evoltion of an option on a default underlying
typedef enum DEF_MAT_OPT{
	q_KO = 1, //1
	q_KO_NO_ACC =0	, //0
	q_NO_KO_NO_ACC =-1 //-1
} qDEF_MAT;




#endif  //_ICM_ENUMS_H_
