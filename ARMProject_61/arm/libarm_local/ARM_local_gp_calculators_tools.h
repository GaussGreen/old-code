#ifndef ARMLOCAL_GP_CALCULATORS_TOOLS_H
#define ARMLOCAL_GP_CALCULATORS_TOOLS_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include "ARM_result.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <gpbase/gpvector.h>
#include <gpbase/gpmatrix.h>

using ARM::ARM_GramFctorArg;

class ARM_ReferenceValue;

///////////////////////////////////////////////
/// Function to get the association object <-> XL name
/// Used in Set and Get Calculator fonctions
///////////////////////////////////////////////
/// General types to Get from a Calculator
#define GC_ACCESS_SECURITY     "SECURITY"
#define GC_ACCESS_MODEL        "MODEL"
#define GC_ACCESS_CALIB        "CALIBMETHOD"
#define GC_ACCESS_PRECALIB     "PRECALIBMETHOD"
#define GC_ACCESS_MDM          "MDM"

/// General types of sub product to Get from a Calculator
#define GC_ACCESS_CAP          "CAP"
#define GC_ACCESS_FLOOR        "FLOOR"
#define GC_ACCESS_FUNDING      "FUNDING"
#define GC_ACCESS_STDLEG       "STDLEG"
#define GC_ACCESS_RFLEG        "RFLEG"
#define GC_ACCESS_STDSWAP      "STDSWAP"
#define GC_ACCESS_RFSWAP       "RFSWAP"
#define GC_ACCESS_BERMUDA      "BERMUDA"

/// Specific types to Get from a Calculator
#define GC_ACCESS_OSW_PORT              "OSWPORTFOLIO"
#define GC_ACCESS_OSW_SECOND_PORT       "OSWSECONDPORTFOLIO"
#define GC_ACCESS_SHORT_TERM_PORT       "STMPORTFOLIO"
#define GC_ACCESS_STD_PORT				"STDPORTFOLIO"
#define GC_ACCESS_EXO_SWAPTION          "EXOSWAPTION"
#define GC_ACCESS_EXOSWAP			    "SWAP"
#define GC_ACCESS_CF_PORT               "CFPORTFOLIO"
#define GC_ACCESS_DOM_OSW_PORT          "DOMOSWPORTFOLIO"
#define GC_ACCESS_FOR_OSW_PORT          "FOROSWPORTFOLIO"
#define GC_ACCESS_FX_PORT               "FXPORTFOLIO"
#define GC_ACCESS_FLOORED_FX_PORT       "FLOOREDFXPORTFOLIO"
#define GC_ACCESS_CAPPED_FX_PORT        "CAPPEDFXPORTFOLIO"
#define GC_ACCESS_REDEMPTION_FX_PORT    "REDEMPTIONFXPORTFOLIO"
#define GC_ACCESS_EXTRA_FX_PORT         "EXTRAFXPORTFOLIO"
#define GC_ACCESS_EXTRA_FX_CALIB        "EXTRAFXCALIBMETHOD"
#define GC_ACCESS_POWER_REVERSE			"POWERREVERSE"
#define GC_ACCESS_DFBSMODEL				"DFBSMODEL"
#define GC_ACCESS_OSW_BSMODEL			"OSWBSMODEL"
#define GC_ACCESS_QGM_SKEWPARAM			"SKEWCALIBPARAM"
#define GC_ACCESS_SO_PORT				"SOPORTFOLIO"
#define GC_ACCESS_SO_SECOND_PORT		"SOSECONDPORTFOLIO"
#define GC_ACCESS_SO_CSTMANAGER         "CSOCSTMANAGER"
#define GC_ACCESS_CMSLONG_PORT			"CMSLNGPORTFOLIO"
#define GC_ACCESS_CMSSHORT_PORT			"CMSSHORTPORTFOLIO"
#define GC_ACCESS_OPT_PORT              "OPTIONPORTFOLIO"
#define GC_ACCESS_LOCAL_PORT              "LOCALPORTFOLIO"
#define GC_ACCESS_CSO_PORT				"CSOPORTFOLIO"
#define GC_ACCESS_CALIB_FLAG			"CCSOCALIBFLAG"

extern long ARMLOCAL_DateStripCombiner_Create(
         const VECTOR<long>&    C_dateStripIds,
		 const string&			funcToMerge,
         ARM_result&			result, 
         long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_DateStripCombiner_GetData(
         long dateStripCombinerId,
	     long dataType,
		 size_t dateStripNb,
		 VECTOR<double>& Data,
		 ARM_result& result );

extern long ARMLOCAL_DateStripCombiner_GetMergeData(
	long dateStripCombinerId,
	VECTOR<double>& Data,
    ARM_result& result );

extern long ARMLOCAL_ARM_GetCcyFromGenCalculator(long calcId,
        const CCString& ccyType,
        ARM_result& result,
        long objId = -1);

extern long ARMLOCAL_Calculator_Initialize(
        const long& gcId,
        const long& mktMangerId,
        const CCString& toCalSigma,
        const CCString& toCalMrs,
		const CCString& strikeTypeToCalMrs,
        const CCString& toAdjKcap,
        const CCString& toAdjKfloor,
        const CCString& modelType,
        const CCString& toCalKskew,
        const double&    kShift,
		const long& frontier,
        ARM_result&	result, 
        long        objId );

extern string GCGetTypeToClass(const string& typeToGet, long gcId);

extern long ARMLOCAL_GC_Get(
        long gcId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_GC_Set(
        long calculatorId,
        long dataId,
		const vector<string>& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_GC_GetPricingData(
		long gcId,
		const string& key,
		ARM_GramFctorArg& argResult,
		ARM_result& result );

extern long ARMLOCAL_CallOnMepiVanillaArgCreate(
		const string& CurveName,
		const string& EquityName,
		double startDate,
		double endDate,
		long resetFreq,
		double riskFactor,
		double strike,
		double maxBorrow,
		double protectionCurveStart,
		double protectionCurveEnd,
		double startingPortfolio,
		double startingCash,
		double minInvested,
		double leverageCost,
		double cashSpread,
		double fees,
		double alreadyAsianed,
		long asianDatesNb,
        ARM_result&  result,
        long objId);

int ARMLOCAL_IsWarning_OnObject( long C_ObjIdId,
	ARM_result&  result );

long ARMLOCAL_GetWarning_OnObject( const long& C_ObjIdId,
	ARM_result&  result,
	long objId = -1
	);


/****************************************************************
	General Swaption Converter to Variable Notional Swaption
*****************************************************************/
long ARMLOCAL_ConvertToVarNotionalSwaption(long yieldCurveId,long swaptionId,
										   ARM_result& result,long objId);

/****************************************************************
	Basis converter
*****************************************************************/

extern long ARMLOCAL_BasisConverter(
	const double& asOfDate,
	const string& domCcyStr,
	const string& fgnCcyStr,
	const long& domDateStripId,
	const long& forDateStripId,
	const long& fundDateStripId,
	const string& domDayCount,
	const string& domFreq,
	const string& fgnDayCount,
	const string& fgnFreq,
	const long& domZcId,
	const long& fngZcId,
	const long& domDiscZcId,
	const long& fgnDiscZcId,
	const long& forexId,
	const long& domNotinalCrv,
	const long& forNotionalCrv,
	const long& forSpreadCrv,
	VECTOR<double>& Margin,
	ARM_result& result);


#endif

//-----------------------------------------------------------------------------
/*---- End of file ----*/
