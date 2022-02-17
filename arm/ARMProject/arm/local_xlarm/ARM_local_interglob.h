#ifndef _ARM_LOCAL_INTERGLOB_H
#define _ARM_LOCAL_INTERGLOB_H


#include <libCCtools++\CCString_STL.h>

#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_init.h>


#include <ARM\libarm\ARM_message.h>

#include <ARM\libarm\ARM_result.h>

#include "ICMKernel\glob\icm_enums.h" 

#define ARM_DEFAULT_ERR			-990  

#ifdef STL_WIN32
#include <map>
#include <vector>
#define MAP		std::map
#else
#include <map.h>
#include <vector.h>
#define MAP		map
#endif	// STL_WIN32

/* mimic the class AUTO_MODE from modes.h   */
struct INTERFACE_AUTO_MODE
{
    enum _type
    {
        AUTO_MODE_NONE      = 1,
        AUTO_MODE_IRG       = AUTO_MODE_NONE << 1,          // calibrate on swaptions
        AUTO_MODE_IDX       = AUTO_MODE_IRG << 1,           // calibrate on caps
        AUTO_MODE_SWOPT_IRG = AUTO_MODE_IDX <<1,            // calibrate on swaptions, optimize on caps
        AUTO_MODE_IRG_SWOPT = AUTO_MODE_SWOPT_IRG << 1,     // calibrate on caps, optimize on swaptions
        AUTO_MODE_DUAL      = AUTO_MODE_IRG_SWOPT << 1,     // calibrate swaptions and captions, 2 vols curves
        AUTO_MODE_AUTOMEAN  = AUTO_MODE_DUAL << 1,          // optimize the mean reversion rate
        AUTO_MODE_BASISDISC = AUTO_MODE_AUTOMEAN << 1,		// basis discounting
		AUTO_MODE_CAPDIAG   = AUTO_MODE_BASISDISC << 1,
		AUTO_MODE_DIAG		= AUTO_MODE_CAPDIAG << 1
    };
};



typedef MAP<const char*, ARM_message*, ltstr> ARM_OBJECTS_LIST_t;
extern ARM_OBJECTS_LIST_t* ARM_OBJECTS_LIST;

typedef MAP<const char*, ARM_message*, ltstr> ARM_ERRORS_LIST_t;
extern ARM_ERRORS_LIST_t* ARM_ERRORS_LIST;



extern ARMLOCAL_Init* armlocal_init;

#define XLLOCALARM_INIFILE_NAME		"C:\\Program Files\\Arm\\local_xlarm.ini"
#define DEFAULT_CURRENCY			"EUR"
#define DATA_FOLDER			        "P:\\Samba\\arm-import"


#define ARM_OK					0
#define ARM_KO					-1





#define ARM_NULL_OBJECT_ID		-1

/// The default size for error message is 1000
// 
//	JLA: we don't want handle message size limitation.
//
// #define ARM_RESULT()					result.setRetCode(ARM_KO); \
// 										char* var = new char[1000]; \
// 										x.GetErrorMessage(var); \
// 										msg = var; \
// 										result.setMsg(msg); \
// 										long retCodeTmp = result.getRetCode (); \
// 										delete var;  \
// 										return (retCodeTmp); 
// 

#define ARM_RESULT()			{ result.setRetCode(ARM_KO); \
								std::string var = x.GetErrorString(); \
								msg = var.c_str(); \
								result.setMsg(msg); \
								long retCodeTmp = result.getRetCode ();  \
								return (retCodeTmp); }

//	This macro allows to fill an ARM_Result with a stream syntax
//	usage : 
//	ICMRESULT( result, "Some message, i="<<i<<" and string="<<myArg" ) ;
#ifndef	ICMRESULT
	#define ICMRESULT(argICMRESULT,msgICMRESULT) \
	{ \
		std::stringstream _sstr; _sstr<<msgICMRESULT<<std::endl; \
		argICMRESULT.setMsgString(_sstr.str()); \
	}
#endif	// ICMRESULT


extern ARM_message* BuildErrorMessage (const CCString& errValue, bool storeInLog = false);

extern CCString GetCurStringCellCoordinates ();
extern CCString GetLastCurCellEnvValue ();
extern CCString GetLastCurActiveCellEnvValue ();

extern CCString LocalMakeObjectId (long objId, const CCString& objClass);
std::string	LocalGetObjectLabel(long objId,const std::string& objClass); 
extern CCString LocalGetStringObjectClass (const CCString& stringObjectId);
extern CCString LocalGetStringObjectClassAndError(const CCString& stringObjectId);
extern long LocalGetNumObjectId (const CCString& stringObjectId);

void LocalSetCurCellEnvValue (const CCString& curClass, long objId);
void LocalSetCurCellEnvValue (const CCString& curClass, long objId,const CCString& caller);

extern void SetCurCellErrValue (const CCString& errValue);
extern void FreeCurCellErr ();
extern long FreeCurCellContent ();
extern long FreeCurCellContentWithoutFreeingARMObject();
extern long FreeCurCellContent (const CCString & envVal);
extern CCString GetEnvVar (const CCString& cellCoord);
extern CCString GetErrValue (const CCString& cellCoord);
extern void PrintErrorsList ();

extern long LOCALARM_PersistentListsInit ();
extern long LOCALARM_PersistentListsDelete ();
extern long LOCALARM_PersistentListsClear ();
extern void LOCALARM_IniFileRead(void);
extern long LOCALARM_DeconnexionEToolkit ();

extern long ARM_ConvPriceYield(const CCString& aPriceYield, ARM_result& result);
extern long ARM_ConvIrIndName (const CCString& aIrIndName, ARM_result& result);
extern long ARM_ConvTmIxName (const CCString& aTmIxName, ARM_result& result);
extern long ARM_ConvFrequency (const CCString& aFrequency, ARM_result& result);

extern long ARM_ConvCMIndName (const CCString& aIrType, ARM_result& result);
extern long ARM_ConvCompMeth (const CCString& aCompMeth, ARM_result& result);
extern long ARM_ConvCallOrPut (const CCString& aCorP, ARM_result& result);
extern long ARM_ConvMktType (const CCString& aMkt, ARM_result& result);
extern long ARM_ConvInterpMethod (const CCString& aIntMeth, ARM_result& result);
/*MICHAEL*/
extern long ARM_ConvWhatIsInterp (const CCString& aIntMeth, ARM_result& result);
/*MICHAEL*/


/// common function for interpolation
extern long ARM_ConvCPInterpMethod( const CCString& InfMethod, ARM_result& result, long Flag );
extern long ARM_ConvCPIMonthlyInterpMethod( const CCString& InfMethod, ARM_result& result);
extern long ARM_ConvCPIDailyInterpMethod( const CCString& InfMethod, ARM_result& result);
extern long ARM_ConvCPIExtrapolMethod( const CCString& InfMethod, ARM_result& result);
extern long ARM_ConvInfSwapType (const CCString& input, ARM_result& result);
extern long ARM_ConvCapOrFloor (const CCString& aCorF, ARM_result& result);
extern long ARM_ConvExerciseType (const CCString& aOptType, ARM_result& result);
extern long ARM_ConvRecOrPay (const CCString& aRorP, ARM_result& result);
extern long ARM_ConvLongOrShort (const CCString& aLorS, ARM_result& result);
extern long ARM_ConvFwdRule (const CCString& aFwdRule, ARM_result& result);
extern long ARM_ConvRule (const CCString& rule, ARM_result& result);
extern long ARM_ConvObjectClass (const CCString& idUnderlying, ARM_result& result);
extern long ARM_ConvSvtyParam (const CCString& aParam, ARM_result& result);
extern long ARM_ConvCurvSvtyParam (const CCString& aParam, ARM_result& result);
extern long ARM_ConvUpDownDouble (const CCString& aUpDownDouble, ARM_result& result);
extern long ARM_ConvInOut (const CCString& aInOut, ARM_result& result);
extern long ARM_ConvStdBoost(const CCString& aInOut, ARM_result& result);
extern long ARM_ConvVolType (const CCString& aVolType, ARM_result& result);
extern long ARM_ConvForwardYieldMethod (const CCString& aForwardYieldMeth, ARM_result& result);
extern long ARM_ConvYieldOrVol (const CCString& yieldOrVol, ARM_result& result);
extern long ARM_ConvCalcMod (const CCString& calcMod, ARM_result& result);
extern long IsCall (const CCString& aString);
extern long IsPut (const CCString& aString);
extern long IsCap (const CCString& aString);
extern long IsFloor (const CCString& aString);
extern long IsReceive (const CCString& aString);
extern long IsPay (const CCString& aString);
extern long IsBermudan (const CCString& aString);
extern long StrikeCode (const CCString& aString, ARM_result& result);
extern long ARM_StrikeCode (const CCString& aString, ARM_result& result);
extern long ARM_GetDimType( const CCString& aRule, ARM_result& result);
extern long ARM_ConvParam (const CCString& param, ARM_result& result);
extern long retrieveArmFile (const CCString& serverFileName, const CCString& clientFileName);
extern long extractCurveFromFile (const CCString& curveFileName, VECTOR<CCString>& matu, VECTOR<double>& rate);
extern long ARM_ConvIrIndNameToFreq (const CCString& aIrIndName, ARM_result& result);
extern long ARM_ConvIrIndIdToFreqId (long IrIndId, ARM_result& result);
extern long ARM_ConvAutoMode (const CCString& aAutoMode, ARM_result& result);
extern long ARM_ConvAutoMode2 (const CCString& aAutoMode, ARM_result& result);
extern long ARM_ConvFineMode (const CCString& aFineMode, ARM_result& result);
extern long ARM_ConvMcMethod (const CCString& aMcMethod, ARM_result& result);
extern long ARM_ConvShapeType (const CCString& aShapeType, ARM_result& result);
extern long ARM_ConvCalculationMethod (const CCString& aCalcMeth, ARM_result& result);
extern long ARM_ConvAmortMethod (const CCString& aAmortMethod, ARM_result& result);
extern long ARM_ConvModelType (const CCString& aModelType, ARM_result& result);
extern long ARM_ConvVolType2 (const CCString& aVolType, ARM_result& result);
extern long ARM_ConvSmileNoSmile (CCString* aSmileNotSmile, ARM_result& result);
extern long ARM_ConvYesOrNo (const CCString& YesOrNo, ARM_result& result);
bool ARM_ConvYesOrNo (const std::string& item); 
extern long ARM_ConvSmiledModelFlag (const CCString& SmiledModelType, ARM_result& result);
extern long ARM_NotionalType(const CCString& aRule, ARM_result& result);
extern long ARM_ConvTypeValues (const CCString& aValuesType, ARM_result& result);
extern long ARM_ConvTypeDates (const CCString& aDatesType, ARM_result& result);
extern long ARM_ConvSigRhoNu (const CCString& param, ARM_result& result);
extern long ARM_ConvPricerType (const CCString& param, ARM_result& result);
extern long ARM_ConvBootstrap (const CCString& param, ARM_result& result);
extern long ARM_NotionalType( const CCString& aRule, ARM_result& result);
extern long ARM_ConvSmiledModelFlag (const CCString& SmiledModelType, ARM_result& result);
extern CCString ARM_ConvTypeDeal(const CCString& typeDeal, ARM_result& result);
extern long ARM_Convhedge(const CCString& typeHedge, ARM_result& result);
extern long ARM_ConvDateStripDataType(const CCString& aValuesType, ARM_result& result);
long ARM_ConvDatePowerReverseDataType(const CCString& aValuesType, ARM_result& result);
extern long ARM_ConvCorrPfType(const CCString& corrPfType, ARM_result& result);
extern long ARM_ConvSeasonalityMode(const CCString& seasonMode, ARM_result& result);
extern long ARM_ConvCopulaPricerType (const CCString& PricerType, ARM_result& result);
extern long ARM_ConvMaturityCapMode (const CCString& Mode, ARM_result& result);
extern long ARM_ConvMaturityCapCalibrationMode (const CCString& MaturityCapCalibrationMode, ARM_result& result);
extern long ARM_ConvSummitFormulae(const CCString& param, ARM_result& result);
// extern long ARM_ConvInterpolType(const CCString& Type, ARM_result& result);


// Fonctions de convertions pour le credit ---------------------------------------
extern long ARM_ConvStringToBool (const CCString& aBoolean, ARM_result& result);
extern long ARM_ConvCalibMeth (const CCString& Type, ARM_result& result);
extern long ARM_ConvCBOptionType (const CCString& aOptionType, ARM_result& result);
extern long ARM_ConvCBOptionStrikeType (const CCString& aStrikeType, ARM_result& result);
extern long ARM_ConvCBOptionAccruedOnEx (const CCString& aAccrued, ARM_result& result);
extern long ARM_ConvCBOptionBarrierType (const CCString& aBarrierType, ARM_result& result);
extern long ARM_ConvCBOptionBarrierStrikeType (const CCString& aStrikeType, ARM_result& result);
// extern long ARM_ConvAccOnDef (const CCString& aRorP, ARM_result& result);
extern long ARM_ConvPaysOnDef (const CCString& aRorP, ARM_result& result);
extern long ARM_CheckDate (const double& date, ARM_result& result);
//extern long ICM_ConvCurvType (const CCString& aOptType, ARM_result& result);
extern long ICM_ConvCptType (const CCString& option, ARM_result& result);
extern long ICM_ConvGreekType (const CCString& option, ARM_result& result);
extern long ARM_ConvTypeGenDates (const CCString& aDatesType, ARM_result& result);
extern CCString ARM_ConvTypeModel(const CCString& typeModel, ARM_result& result);
extern long ARM_ConvTypePayLag (const CCString& aDatesType,ARM_result& result);
long ARM_ConvCreditPricerType (const CCString& PricerType, ARM_result& result);
long ARM_ConvCreditPricerType (const std::string&name); 
extern long ARM_ConvCreditMezzType (const CCString& MezzType, ARM_result& result);
qCDS_ADJ ARM_ConvAdjCalCDS (const CCString& aDatesType, ARM_result& result);
// CCString ARM_ConvAdjCalCDSToString (qCDS_ADJ AdjCDS, ARM_result& result);
extern long ARM_ConvLegType (const CCString& Type, ARM_result& result);
extern long ARM_ConvCalibType (const CCString& Type, ARM_result& result);
extern long ARM_ConvUnderlyingMatuStyle (const CCString& MaturityCapMode, ARM_result& result) ;
extern long ARM_ConvTypeUpOrLow (const CCString& UpOrLow,ARM_result& result);
// extern long ARM_Conv_Sectorial_Correlation (const CCString& Type,ARM_result& result);

// End Fonctions de convertions pour le credit -----------------------------------

extern long LocalExtractCurveFromFileMO (const CCString& curveFileName, std::vector<CCString>& matu, std::vector<double>& rate, long adjOrNotId);
extern long LocalExtractCurveFromFileBS (const CCString& curveFileName, std::vector<CCString>& matu, std::vector<double>& rate);
extern long LocalExtractVolFromFile (const CCString& FileName, std::vector<CCString>& yearTerm, std::vector<CCString>& tenor, std::vector<double>& vol);
extern long LocalExtractSmileFromFile (const CCString& FileName, std::vector<CCString>& yearTerm, std::vector<double>& strike, std::vector<double>& vol);
extern long ExtractVectorDoubleFromFile (const CCString& FileName, std::vector<double>& vect, long& vecSize);
extern long ExtractVectorDateFromFile (const CCString& FileName, std::vector<CCString>& vect, long& vecSize);
extern long LocalExtractInfoFromFilePRCS(const CCString& FileName, std::vector<CCString>& listeString, std::vector<double>& listeDouble, int& nbCol);

extern long ARM_ConvSummitManual (const CCString& aSummitManual, ARM_result& result);

extern long ARM_ConvReplicMode( const CCString& replicMode, ARM_result& result);
extern long ARM_ConvStopMode( const CCString& stopMode, ARM_result& result);
extern long ARM_ConvIndexMethod(const CCString& IndexMethod, ARM_result& result);
extern long ARM_ConvKO_Or_NoKO (const CCString& param, ARM_result& result);
extern long ARM_ConvKOStyle (const CCString& param, ARM_result& result);
extern long ARM_ConvAccStyle (const CCString& param, ARM_result& result);


extern long ARM_ConvINFDigitalPayoffType (const CCString& payOffType, ARM_result& result);

#endif	// _ARM_LOCAL_INTERGLOB_H

// EOF %M%
