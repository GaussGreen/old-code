#ifndef ARM_LOCAL_PARSEXML_H
#define ARM_LOCAL_PARSEXML_H

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml_common.h>

#include <set>
#include <gpmodels\enummodel.h>
using ARM::ARM_PricingModelType ;

class ARM_ZeroCurve;
class ARM_ZeroLInterpol;
class ARM_VolLInterpol;
class ARM_VolCube;
class ARM_Date;
class ARM_Matrix;
class ARM_Vector;
class ARM_Swaption;
class ARM_Option;
class ARM_PowerReverse;
class ARM_Security;
class ARM_SwapLeg;
class ARM_ReferenceValue;
class ARM_FlexAccretSwaption;
class ARM_OptionPortfolio;
class ARM_CapFloor;
class ARM_SpreadOption;
class ARM_StdPortfolio;
class ARM_VolCurve;
class ARM_StripOption;
class ARM_FxSpreadStripOption;
class ARM_BasisCurve;
class ARM_CorridorLeg;
class ARM_CorridorDblCondition;

namespace ARM{
	class ARM_CRFCalculator;
	class ARM_MaturityCapCalculator;
	class ARM_TARNCalculator;
	class ARM_BermudaSwaptionCalculator;
	class ARM_CallableSnowBallCalculator;
	class ARM_CaptionCalculator;
	class ARM_FXVanillaCalculator;
	class ARM_SeasonalityManager;
	class ARM_ResetManager;
	class ARM_InfCurv;
};


#define SUMMITSTRING "ARMNT_CALL/ASSETINFO/"
#define ETKSTRING "Response/"

// position ds ARM_CONFIG_VECTOR pour CRF
#define CRF_OSW_CFG_VECT_IDX   7
#define CRF_CAP_CFG_VECT_IDX   6
#define CRF_FLOOR_CFG_VECT_IDX 5


// Flow Structure necessary for retrieving a Flow attributes
typedef struct _SummitFlow
{
    ARM_Date fStartDate;
    ARM_Date fEndDate;
    ARM_Date fFixingDate;
	ARM_Date fPayDate;

	int      fIntDays;
    double   fRate;
	double   fFxRate;
	double   fForward;
	double   fSpread;
	double   fDecompRate;
	double   fNotional;
    double   fInterimInterest;
	double   fFlows;
	double   fStrike;
	string   fType;//char  fType[20];
	int      fDays;
	string	 fCcy;// char  fCcy[20];
	double   fZeroRate;
	double   fDiscFactor;
	double   fPV;
	double   fAICRate;
}
ARM_SummitFlow;

// For Retrieving A Summit Flow

extern void RetrieveSummitFlow(const CCString& tradeId, const CCString& tradeType,
						       const CCString& curveId,
						       const ARM_Date& inDate, 
						       ARM_SummitFlow* outFlow);

extern double convPlotInYearTerm(const CCString& plot, ARM_Date AsOf, ARM_Date Settle, char* currency);


extern int ARMLOCAL_ParseXMLForInterpolator(const char* chaineXML);

extern ARM_ZeroLInterpol* ARMLOCAL_ParseXMLForZC(const char* chaineXML,
												 ARM_Date aSdate,
												 const char* sCcy,
												 int interp);

extern ARM_ZeroLInterpol* ARMLOCAL_ParseXMLForMYAndZC(const char* chaineXML,
													  ARM_Date aSdate,
													  const char* sCcy,
													  long adjOrNotId = 1,
													  long rawId = K_DEFAULT_CURVEMOD,
													  long swapFrqId = K_DEF_FREQ);

extern ARM_ZeroLInterpol* ARMLOCAL_CreateZC(const CCString& index,
											const CCString& currency,
											const CCString& cvName,
											ARM_Date aSdate,
											const CCString& raw = "DEFAULT",
											long adjOrNotId = 1,
											long swapFrqId = K_DEF_FREQ);

extern ARM_BasisCurve* ARMLOCAL_CreateZCSpreaded(const CCString& index,
												 const CCString& currency,
												 const CCString& cvName,
												 ARM_Date aSdate,
												 long adjOrNotId = 1,
												 const CCString& raw = "DEFAULT",
												 long swapFrqId = K_DEF_FREQ,
												 long mmFrqId = K_DEF_FREQ,
												 long interpId = K_LINEAR,
												 ARM_ZeroCurve* ZeroCurve = NULL);

extern ARM::ARM_InfCurv* ARMLOCAL_ParseXMLForInfZC(const char* chaineXML,
												   double aSdate,
												   const CCString& index);

extern long ARMLOCAL_ParseXMLForMY(const char* chaineXML,
								   long adjOrNotId,
								   VECTOR<CCString>* matu,
								   VECTOR<double>* yield);

void ARMLOCAL_ParseXMLForMYCalypso(const std::string& chaineXML,
								   bool doAdjust,
								   std::vector<std::string>& matu,
								   std::vector<double>& yield);

extern VECTOR<CCString> ARMLOCAL_GetListTenorsFromXML(const char* chaineXML,
													  VECTOR<CCString>* listRequetes,
													  CCString NomParam);

extern VECTOR<CCString> ARMLOCAL_GetListStrikesFromXML(const char* chaineXML,
													   VECTOR<CCString>* listRequetes);


extern long ARMLOCAL_GetDatesAndDatasForVol(const char* chaineXML,
											ARM_Date asof,
											long index,
											ARM_Matrix* &Volatility,
											ARM_Vector* &vYearTerms,
											const char* currency,
											long nbStrikes = 0,
											VECTOR<CCString>* listYearTerms = NULL);

extern long ARMLOCAL_GetStringAndDatasForVol(const char* chaineXML,
											 ARM_Date asof,
											 long index,
											 ARM_Matrix* &Volatility,
											 VECTOR<CCString> &yearTerms,
											 long nbStrikes = 0);

extern long ARMLOCAL_GetDatesAndDatasForCurve(const char* chaineXML,
											  const CCString& currency,
											  ARM_Date asof,
											  long index,
											  ARM_Matrix* &Volatility,
											  ARM_Vector* &vYearTerms,
											  long nbStrikes = 0,
											  long isSmile = 1,
											  VECTOR<CCString>* listYearTerms = NULL,
											  const CCString& impOrHist = "");

extern CCString ARMLOCAL_GetCcyCalFromXML(const char* chaineXML);

extern CCString ARMLOCAL_ParseXMLForGetDateLastWarm();

extern ARM_Date ARMLOCAL_GetFwdDate(const char* chaineXML);

extern double ARMLOCAL_ParseFxCurve(const char* chaineXML);

extern long ARMLOCAL_GetDatesAndDatasForFxVol(const char* xmlResponse,
											  ARM_Date asOfDate,
											  ARM_Matrix* &Volatility,
											  ARM_Vector* &vYearTerms);

extern ARM_Object* ARMLOCAL_ParseExotic(const char* chaineXML,
										const ARM_Date& date,
										const char* filter,
										CCString& structureId,
										CCString& bookName,
										CCString& custId,
										CCString& dealId,
										VECTOR <string>& listAssetId,
										int aResetFreqForCorridorOptim = 365);

extern ARM_Object* ARMLOCAL_ParseGlobalExotic( const char* chaineXML,
											   const ARM_Date& date,
											   const char* filter, 
											   CCString& structureId,
											   CCString& bookName,
											   CCString& custId,
											   CCString& dealId, 
											   VECTOR<string>& listAssetId,
											   int aResetFreqForCorridorOptim = 365);

extern ARM_Object* ARMLOCAL_ParseCorridorSpreadOption(const char* chaineXML,
													  const ARM_Date& date, 
													  CCString& bookName,
													  const char* filter,
													  CCString& custId,
													  CCString& dealId,
													  VECTOR <string>& listAssetId,
													  long isEtk,
													  long isSwap);

extern ARM_Object* ARMLOCAL_ParseMemorySO(const char* chaineXML,
										  const ARM_Date& date,
										  CCString& bookName,
										  CCString& custId,
										  CCString& dealId,
										  long isEtk = 1);

extern ARM_Object* ARMLOCAL_ParseMemoryCap(const char* chaineXML,
										   const ARM_Date& date,
										   CCString& bookName,
										   CCString& custId,
										   CCString& dealId,
										   long isEtk = 1);

extern ARM_Object* ARMLOCAL_ParseObject(const char* chaineXML,
										  const char* typeDeal,
										  const ARM_Date& date,
										  const char* filter,
										  CCString& bookName,
										  CCString& structureId,
										  CCString& custId,
										  CCString& dealId);

extern ARM_Swaption* ARMLOCAL_ParseSwaption(const char* chaineXML,
											CCString& bookName,
											CCString& structureId,
											CCString& custId,
											CCString& dealId,
											ARM_ReferenceValue& rPremium,
											long isEtk = 1);

extern ARM_Option* ARMLOCAL_ParseFxOption(const char* chaineXML,
										  CCString& bookName,
										  CCString& structureId,
										  ARM_Date& primeDate,
										  double& prime,
										  CCString& primeCcy,
										  CCString& custId,
										  CCString& dealId);

extern ARM_PowerReverse* ARMLOCAL_ParseGenPRCS(const char* chaineXML,
											   const ARM_Date& date,
											   CCString& bookName,
											   CCString& structureId,
											   CCString& custId,
											   CCString& dealId);

extern ARM_PowerReverse* ARMLOCAL_ParsePRCS(const char* chaineXML,
											const ARM_Date& date,
											CCString& bookName,
											CCString& structureId,
											CCString& custId,
											CCString& dealId);

extern ARM_PowerReverse* ARMLOCAL_ParseMONOPRCS(const char* chaineXML,
												const ARM_Date& date,
												CCString& bookName,
												CCString& structureId,

												CCString& custId,
												CCString& dealId);

extern ARM_PowerReverse* ARMLOCAL_ParseUnderlyingPRCS ( const char* chaineXML,
														const ARM_Date& date,
														CCString& bookName,
														CCString& structureId,
														CCString& custId,
														CCString& dealId);

extern ARM_FlexAccretSwaption* ARMLOCAL_ParseFlexAcc(const char* chaineXML,
													 const ARM_Date& date,
													 CCString& bookName,
													 CCString& custId,
													 CCString& dealId,
													 long isEtk = 1);

extern ARM::ARM_CRFCalculator* ARMLOCAL_ParseCRF(const char* chaineXML,
												 const ARM_Date& date,
												 CCString& bookName,
												 CCString& custId,
												 CCString& dealId,
												 long isEtk = 1,
                                                 double FxSpot = -1);

extern ARM_OptionPortfolio* ARMLOCAL_ParseCRA(const char* chaineXML,
											  const ARM_Date& date,
											  CCString& bookName,
											  CCString& custId,
											  CCString& dealId,
											  ARM_Vector* CraPricing,
											  long isEtk = 1);

extern ARM_CapFloor* ARMLOCAL_ParseCap(const char* chaineXML, 
									   const ARM_Date& date,
									   CCString& bookName, 
									   CCString& custId, 
									   CCString& dealId,
									   ARM_ReferenceValue& Premium,
									   long isEtk = 1);

extern ARM_SpreadOption* ARMLOCAL_ParseSpreadOption(const char* chaineXML, 
												   const ARM_Date& date,
												   CCString& bookName, 
												   CCString& custId, 
												   CCString& dealId,
												   ARM_ReferenceValue& Premium,
												   long isEtk = 1,
												   int aResetFreqForCorridorOptim = 365);

extern ARM::ARM_MaturityCapCalculator* ARMLOCAL_ParseMaturityCap(const char* chaineXML, 
															 const ARM_Date& date, 
															 CCString& bookName, 
															 CCString& custId, 
															 CCString& dealId, 
															 ARM_ReferenceValue& Premium, 
															 long isEtk=1);

extern ARM::ARM_TARNCalculator* ARMLOCAL_ParseTarn(	const char* chaineXML, 
													const ARM_Date& date, 
													CCString& bookName, 
													CCString& custId, 
													CCString& dealId, 
													ARM_ReferenceValue& Premium, 
													long isEtk=1 );

extern ARM::ARM_BermudaSwaptionCalculator* ARMLOCAL_ParseBermudaSwaption(	const char* chaineXML, 
																			const ARM_Date& date, 	
																			CCString& bookName, 
																			CCString& custId, 
																			CCString& dealId, 
																			ARM_ReferenceValue& Premium, 
																			long isEtk=1 );

extern ARM::ARM_CaptionCalculator* ARMLOCAL_ParseCaption(const char* chaineXML, 
														 const ARM_Date& date, 	
														 CCString& bookName, 
														 CCString& custId, 
														 CCString& dealId, 
														 ARM_ReferenceValue& Premium, 
														 long isEtk=1 );

extern ARM_FxSpreadStripOption* ARMLOCAL_ParseFxStrip(	const char* chaineXML, 
														const ARM_Date& date, 	
														CCString& bookName,
														CCString& structureId,
														CCString& custId, 
														CCString& dealId,
														long isEtk=1 );

extern ARM::ARM_FXVanillaCalculator* ARMLOCAL_ParseFxStripCalculator(const char* chaineXML, 
																const ARM_Date& date, 	
																CCString& bookName,
																CCString& structureId,
																CCString& custId, 
																CCString& dealId,
																long isEtk = 1 );

extern ARM_Swap* ARMLOCAL_ParseSwap(const char* chaineXML,
									const ARM_Date& date, 
									CCString& bookName,
									CCString& custId, 
									CCString& dealId,
									VECTOR<string>& listAssetId,
									long isEtk = 1);

extern double ARMLOCAL_GetFxCurve(const CCString& ccy1,
								  const CCString& ccy2,
								  ARM_Date asOf,
								  double amount,
								  const CCString& cvName);

extern ARM_StdPortfolio*  ARMLOCAL_ParseCorridorDblCondition(const char* chaineXML,
										   const ARM_Date& date,
										   CCString& bookName,
										   CCString& custId,
										   CCString& dealId,
										   long isEtk = 1);

ARM_VolLInterpol* etoolkit_GetVolATMFromSummit(const CCString& index,
											   const CCString& currency,
											   const CCString& cvName,
											   ARM_Date date,
											   const CCString& vtype,
											   const CCString& impOrHist = "IRFWDVOL");

ARM_VolLInterpol* etoolkit_GetSmileFromSummit(const CCString& index,
											  const CCString& currency,
											  const CCString& cvName,
											  ARM_Date date,
											  const CCString& vtype,
											  const CCString& matuIndex,
											  int smileType,
											  const CCString& impOrHist = "");

ARM_VolCube* etoolkit_GetVolCubeFromSummit(const CCString& index,
										   const CCString& currency,
										   const CCString& cvName,
										   ARM_Date date,
										   const CCString& vtype,
										   VECTOR<CCString>& tenors,
										   bool vCorrelCube = false);

long etoolkit_GetInitialVolFromSummit (const CCString& index,
									   const CCString& currency,
									   const CCString& cvName,
									   double date,
									   const CCString& vtype,
									   const CCString& matuIndex,
									   VECTOR<CCString> *maturities,
									   VECTOR<CCString> *tenors,
									   VECTOR<double> *vol);

ARM_VolLInterpol* etoolkit_GetFXVolATMFromSummit(const CCString& ccy1,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const CCString& impOrHist="FXVOL");

ARM_VolLInterpol* etoolkit_GetFXCorrelFromSummit(const CCString& ccy1,
												 const CCString& index,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const VECTOR<CCString>& tenors);

ARM_VolLInterpol* etoolkit_GetCorrelFromSummit(const CCString& ccy1,
											   const CCString& index1,
											   const CCString& ccy2,
											   const CCString& index2,
											   ARM_Date date,
											   const CCString& cvName);

double ARMLOCAL_ParseXMLForMeanRev(const CCString& xmlResponse);

double etoolkit_getXMLMeanRevFromSummit(const CCString& C_ccy,
										const CCString& C_index,
										const CCString& C_cvname,
										const ARM_Date& Date,
										const CCString& C_NumFactor="3F");

double ARMLOCAL_ParseXMLForCutOff(const CCString& xmlResponse);

double etoolkit_getXMLCutOffFromSummit(const CCString& C_ccy,
									   const CCString& C_index,
									   const CCString& C_cvname,
									   const CCString& C_NumFactor,
									   const ARM_Date& Date);

double ARMLOCAL_ParseXMLForFixing(const CCString& xmlResponse);

double etoolkit_getXMLFixingFromSummit(const CCString& C_source,
									   const CCString& C_index,
									   const CCString& C_tenor,
									   const CCString& C_ccy,
									   const ARM_Date& Date);

ARM_Vector* ARMLOCAL_ParseXMLForQmodParam(const CCString& xmlResponse);

ARM_Vector* etoolkit_getXMLQmodParamFromSummit(const CCString& C_ccy,
											   const CCString& C_index,
											   const CCString& C_cvname,
											   const ARM_Date& Date);

ARM_VolLInterpol* ARMLOCAL_ParseXMLForQFXParam(const CCString& xmlResponse, const ARM_Date& Date);

ARM_VolLInterpol* etoolkit_getXMLQFXParamFromSummit(const CCString& C_ccy,
												    const CCString& C_index,
												    const CCString& C_cvname,
												    const ARM_Date& Date);

ARM_VolLInterpol* etoolkit_GetXMLFxSmileFromSummit(const CCString& ccy1,
												   const CCString& ccy2,
												   const CCString& cvName,
												   ARM_Date date,
												   const char* CallOrPut = "C",
												   int IncOrDec = K_INCREASING);

ARM_VolLInterpol* etoolkit_GetPivotsAndInterpForFxVolCurve(const CCString& ccy1,
														   const CCString& index,
														   const CCString& cvName,
														   ARM_Date date);

extern void ARMLOCAL_GetListFxStrikesFromXML(const char* chaineXML,
											 VECTOR<CCString>* listRequetes);

long etoolkit_GetInitialFXVolFromSummit (const CCString& ccy1,
										 const CCString& ccy2,
										 const CCString& cvName,
										 const ARM_Date& date,
										 const CCString& impOrHist,
										 VECTOR<CCString>* vMaturities,
										 ARM_Vector* vVolATM);

long etoolkit_GetInitialFXSmileFromSummit (const CCString& ccy1,
										   const CCString& ccy2,
										   const CCString& cvName,
										   const ARM_Date& date,
										   const CCString& impOrHist,
										   ARM_Matrix* mSmile,
										   ARM_Vector* vTenors);

ARM_VolLInterpol* etoolkit_CreateFXVolFromSummit(const CCString& ccy1,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const CCString& impOrHist,
												 const CCString& VolType);

ARM_VolCurve* ARMLOCAL_CreateNewFXVolFromSummitWithCurves(const CCString& ccy1,
															const CCString& ccy2,
															ARM_Date& Date,
															const CCString& cvName,
															ARM_ZeroCurve* domZc,
															ARM_ZeroCurve* forZc,
															double fxSpot,
															long WhatIsInterpolated,
															long correctSplineWithLinear,
															long isATM);

ARM_VolCurve* ARMLOCAL_CreateNewFXVolFromSummitWithForwards(const CCString& ccy1,
															  const CCString& ccy2,
															  ARM_Date& Date,
															  const CCString& cvName,
															  const VECTOR<double>& Forwards,
															  long WhatIsInterpolated,
															  long correctSplineWithLinear,
															  long isATM);

ARM::ARM_SeasonalityManager* ARMLOCAL_ParseXMLForSeasonMgr(const CCString& Index,
														   const CCString& Ccy,
														   const CCString& CvName,
														   const ARM_Date& asof,
														   long modeId);

ARM::ARM_ResetManager* ARMLOCAL_ParseXMLForResetMgr(const ARM_Date& asof,
													const CCString& Index,
													const CCString& Source,
													const CCString& Ccy,
													long isInflatIndexId,
													const CCString& Term);

ARM_ReferenceValue* ARMLOCAL_ParseXMLCRFMeanRevParam(const CCString& xmlResponse, 
													  const ARM_Date& Date, 
													  const char* ccy, 
													  const char* index);

VECTOR<CCString> ARMLOCAL_ParseListTenors(const CCString& index,
										  const CCString& ccy,
										  const CCString& cvname,
										  ARM_Date date,
										  const CCString& cvtype);

ARM_ReferenceValue* etoolkit_getXMLModelParamFromSummit(const ARM_Date& date,
													  const CCString& model,
													  const CCString& type,
													  const CCString& factorName,
													  const CCString& ccy,
													  const CCString& index,
													  const CCString& cvName,
													  long calcModId = K_LININTERPOL_REF);

ARM_ReferenceValue* ARMLOCAL_ParseXMLForModelParam(const CCString& xmlResponse,
												   long calcModId);

extern ARM_StripOption* ARMLOCAL_ParseFxOptionStrip(const char* chaineXML, const ARM_Date& date, 
													CCString& bookName, CCString& structureId, 
													CCString& custId, CCString& dealId, 
													long isEtk=1);

void MakeCorridorSpreadOption(ARM_CorridorLeg* corridorLeg,
							 vector<ARM_Security*>& ListAssets,
							 vector<double>& vWeights,
							 vector<double>& vMktPrices,
							 vector<double>& vVegas,
							 const char* id);

//CALYPSO Parser
ARM_Object* ARMLOCAL_ParseCalypsoObject(string xmlContent, const string type, const string modelType,const ARM_Date& date);

extern ARM::ARM_CRFCalculator* ARMLOCAL_ParseCalypsoCRF(string chaineXML,
												 const ARM_Date& date,
												 const string modeType,
												 string& bookName,
												 string& custId,
												 string& dealId);
extern ARM::ARM_CRFCalculator* ARMLOCAL_ParseCalypsoCRF(string chaineXML,
									 const ARM_Date& asOfDate,
									 ARM_PricingModelType::ModelType	modelType,
									 string& bookName,
									 string& custId,
									 string& dealId);

//extern ARM::ARM_BermudaSwaptionCalculator* ARMLOCAL_ParseCalypsoBermudaSwaption(string chaineXML,
extern ARM_Object* ARMLOCAL_ParseCalypsoBermudaSwaption(string chaineXML,
																			 const ARM_Date& date,
																			 const string modeType,
																			 string& bookName,
																			 string& custId,
																			 string& dealId);

extern ARM::ARM_BermudaSwaptionCalculator* ARMLOCAL_ParseCalypsoBermudaSwaption(string chaineXML,
									 const ARM_Date& asOfDate,
									 ARM_PricingModelType::ModelType	modelType,
									 string& bookName,
									 string& custId,
									 string& dealId);

//long GetVolFromCalypso();


// Calypso vol function

extern ARM_VolCurve* GetVolSmileFromCalypso(ARM_Date asOfDate,ARM_Currency* armCurrency, vector<double> valueVector ,
						vector<string> expirySetVector,vector<double> strikeSetVector, double tenorCount, double tenorRank = 0);

extern ARM_VolCurve* GetVolSmileFromCalypso(ARM_Date asOfDate, string xmlInput);
extern ARM_VolCube* GetVolCubeFromCalypso(ARM_Date asOfDate,string cvName, string xmlInput, string suffix="");
extern ARM_VolLInterpol* GetVolSurfaceFromCalypso(ARM_Date asOfDate, string xmlInput);

ARM_ReferenceValue* ARMLOCAL_ParseFixingFromInstrument(const CCString& xmlResponse);

#endif 
