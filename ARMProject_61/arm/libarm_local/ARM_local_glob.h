#ifndef ARM_LOCAL_GLOB_H
#define ARM_LOCAL_GLOB_H

/// forward declaration
//class ARM_Matrix;
class ARM_AbstractMarketClass;


#include <ARM\libarm\ARM_result.h>
#include <glob\dates.h>
#include <util\refvalue.h>
#include <deque>

namespace ARM
{
	template <typename T> class ARM_GP_T_Matrix;
	typedef ARM_GP_T_Matrix<double> ARM_GP_Matrix;
}


#define GETDEFAULTVALUE -1111
#define GETDEFAULTVALUESTR "DEFAULT"

extern int DBL_INT(double x);

void Local_XLDATE2ARMDATE (double xldate,
						   char* myArmDate);
void Local_XLDATE2ARMDATE (double xldate,ARM_Date&armDate) ;
ARM_Date ConvertToARMDATE (double xldate);

double Local_ARMDATE2XLDATE (const CCString& myArmDate);


ARM::ARM_GP_Matrix * CreateARM_GP_MatrixFromVECTOR(const VECTOR<VECTOR<double> >& param);

ARM_Matrix* CreateARMMatrixFromVECTOR(const VECTOR<double>& param, int nbrows, int nbcolumns);
ARM_Vector* CreateARMMatuVectorFromStrVECTOR(const VECTOR<CCString>& param);

extern long ARMLOCAL_GetCurrency(long ObjId, ARM_result& result);
extern long ARMLOCAL_ARM_GetDefaultCurrency (ARM_result& result);
extern long ARMLOCAL_ARM_SetDefaultCurrency (const CCString& isoCCy, ARM_result& result);

extern long ARMLOCAL_ARM_GetFRFCurrency (ARM_result& result);

extern long ARMLOCAL_bsflexible (double F, double V,double B, 
								 double K,double CallPutFlag,
						  ARM_result& result);

extern long ARMLOCAL_ARM_Price (long secId,
								long modId,
								ARM_result& result);

extern long ARMLOCAL_FreeObject (long secId,
								 ARM_result& result);

extern long ARMLOCAL_FreeAllObjects (ARM_result& result);

extern long ARMLOCAL_NextBusinessDay (double date,
									  const CCString& cal,
									  long days,
									  ARM_result& result);

extern long ARMLOCAL_IsBusinessDay (double date,
									const CCString& isoccyname,
									ARM_result& result);

extern long ARMLOCAL_FutDelivery (const CCString& fut,
								  const CCString& currency,
								  ARM_result& result);

extern long ARMLOCAL_ARM_Accrued (long secId,
								  double fwdDate,
								  long modId,
								  ARM_result& result);

extern long ARMLOCAL_IMPLIEDVOL (long instId,
								 long modelId,
								 double price,
								 bool LnVol,
								 ARM_result& result);

extern long ARMLOCAL_ARM_View (long instId,
							   ARM_result& result,
							   bool aXMLResult = false);

extern long ARMLOCAL_SetNotional (long secId,
								  long rId,
								  double percentRemainder,
								  ARM_result& result,
								  int interpolDates=K_PAY_DATES);

extern long ARMLOCAL_GetExpiry (long secId,
								ARM_result& result);

extern long ARMLOCAL_Sensitivity (long secId,
								  long modId,
								  long paramId,
								  ARM_result& result);

extern long ARMLOCAL_GetFRMShortRateVols(long modId,
										 VECTOR<double>* matu,
										 VECTOR<double>* rate,
										 ARM_result& result);

extern long ARMLOCAL_XCccyAdjustment(long startDate,
									 long endDate,
									 long payFreq,
									 const CCString& domCcy,
									 long forIndexType,
									 const CCString& forCcy,
									 long spreadsId,
									 long zcDomId,
									 long discDomId,
									 long zcForId,
									 long discForId,
									 double FX,
									 long couponId, 
									 long domDc,
									 long forDc,
									 ARM_result& result,
									 long objId = -1);

extern long ARMLOCAL_ADJUSTTOBUSDATE (double date,
									  const CCString& currency,
									  long ruleId,
									  ARM_result& result);

extern long ARMLOCAL_FwdPrice (long secId,
							   long modId,
							   double fwdDate,
							   ARM_result& result);

extern long ARMLOCAL_CvSensitivity (long secId,
									long modId,
									long paramId,
									ARM_result& result);

extern long ARMLOCAL_ARM_BetweenDates (long date1,
									   long date2,
									   long daycountId,
									   long isYearFrac,
									   ARM_result& result);

extern long ARMLOCAL_ARM_CountBusinessDays (long date1,
											long date2,
											CCString calendar,
											ARM_result& result);

extern long ARMLOCAL_ARM_ADDYEARS  (long date,
									long nb,
									long ruleId,
									const CCString& Ccy,
									ARM_result& result);

extern long ARMLOCAL_ARM_ADDMONTHS  (long date,
									 long nb,
									 long ruleId,
									 const CCString& Ccy,
									 ARM_result& result);

extern long ARMLOCAL_FxConvert (const CCString& ccy1,
								const CCString& ccy2,
								double asOfDate,
								double amount,
								const CCString& cvname,
								ARM_result& result);

extern long ARMLOCAL_ARM_ADDPERIOD  (double date,
									 long freq,
									 const CCString& ccy,
									 long nbPeriods,
									 long adjRuleId,
									 long goToEndOfMonth,
									 ARM_result& result);

extern long ARMLOCAL_ARM_Cover(long secId,
							   long modId,
							   ARM_result& result);

extern long ARMLOCAL_ARM_Price_OptUnder (long secId,
										 long modId,
										 ARM_result& result);

extern long ARMLOCAL_ParallelShift (long crvId,
									double value,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_ClonedAndSetNotional (long secId,
										   long rId,
										   double percentRemainder,
										   ARM_result& result,
										   int interpolDates=K_PAY_DATES,
										   long objId = -1);

extern long ARMLOCAL_ClonedAndSet (long secId,
								   long valtosetType,
								   double valtoset,
								   const CCString typetoset,
								   ARM_result& result,
								   long objId = -1);

extern long ARMLOCAL_Clone (long objectId,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_DoPastReset (long secId,
								  VECTOR<CCString>& resetMgrIds,
								  long AsOf,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_FIXRATES (double swapLegId,
							   VECTOR<double>& rate,
							   ARM_result& result,
							   long objId = -1);

extern long ARMLOCAL_ARM_DisplayScheduleDates (long instId,
											   long datesType,
											   long recId,
											   long viewInitExchId,
											   ARM_result& result);
extern long ARMLOCAL_ARM_DisplayScheduleValues (long instId,
												long valuesType,
												long recId,
												long modelId,
												ARM_result& result);

long ARMLOCAL_ARM_DisplayReplicPortfolio(long instId,
										 long WeightOrStrike,
                                         long PayoffOrSensi,
                                         long recId,
                                         long modelId,
                                         VECTOR<double>& DataResult,
										 ARM_result& result);

extern long ARMLOCAL_INTERPOL (const VECTOR<double>& vecX,
							   const VECTOR<double>& vecY,
							   double X,
							   long interpId,
							   ARM_result& result);

extern long ARMLOCAL_TRIANGULARINTERPOL (const VECTOR<double>& vecX,
										 const VECTOR<double>& vecY,
										 const VECTOR<double>& matZ,
										 double X,
										 double Y,
										 ARM_result& result);

extern long ARMLOCAL_KImp (long secId,
						   long modId,
						   double price,
						   long param,
						   ARM_result &result);

extern long ARMLOCAL_ARM_DisplayZC(long zcId,
								   ARM_result &result);

extern long ARMLOCAL_BSSpot (long secId,
							 long modId,
							 double date,
							 ARM_result& result);


extern long ARMLOCAL_ARM_GetMeanRevFromSummit (const CCString& C_ccy,
											   const CCString& C_index,
											   const CCString& C_cvname,
											   double C_date,
											   const CCString& C_NumFactor,
											   ARM_result& result);

extern long ARMLOCAL_ARM_GetCutOffFromSummit (const CCString& C_ccy,
											  const CCString& C_index,
											  const CCString& C_cvname,
											  const CCString& C_NumFactor,
											  double C_date,
											  ARM_result& result);

extern long ARMLOCAL_HEDGE (long secId,
							long hedgeId,
							ARM_result& result);


extern long ARMLOCAL_ARM_GetFxSmileFromSummit (const CCString& C_ccy1,
											  const CCString& C_ccy2,
											  const CCString& C_cvname,
											  double C_date,
											  ARM_result& result,
											  long objId = -1);

extern long ARMLOCAL_GETINSTRUMENTFROMSUMMIT (const CCString& idSummit,
											  CCString& typeId,
											  double asOf,
											  const CCString& filter,
											  ARM_result& result,
											  long objId = -1);

extern long ARMLOCAL_GETINSTRUMENTFROMCALYPSO (CCString& idCalypso,
											   CCString& typeId,
											   CCString& modelType,
											   double asOf,
											   ARM_result& result,
											   long objId = -1);

extern long ARMLOCAL_GETINFOFROMPRCS (long prcsId,
									  const CCString& datatype,
									  ARM_result& result);

extern long ARMLOCAL_GETOBJINFOFROMPRCS (long prcsId,
										 const CCString& datatype,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_INITCRF(long crfId,
							 long zcCpnId,
							 long swVolId,
							 long capVolId,
							 long rhoCapId,
							 long nuCapId,
						     long betaCapId,
							 long rhoSwaptId,
							 long nuSwaptId,
							 long betaSwaptId,
							 long zcFundId,
						     long zcCpnBasisId,
						     long zcFundBasisId,
						     double fxSpot,
							 long isUpdate,	
							 CCString& C_modelType,
							 long meanReversionId,
                             long skewRecalibFlag,
							 long inSABRSigmaOrAlpha,
                             ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_INITMATCAP(long matcapId,
							    long zcId,
							    long capVolId,
							    long rhoCapId,
							    long nuCapId,
							    long betaCapId,
								long nbpas,
                                long inSABRSigmaOrAlpha,
								ARM_result& result,
							    long objId = -1);

extern long ARMLOCAL_INITTARN(long tarnId,
					          long zcCpnId,
					          long swVolId,
					          long capVolId,
					          long rhoCapId,
					          long nuCapId,
					          long betaCapId,
					          long rhoSwoptId,
					          long nuSwoptId,
					          long betaSwoptId,
					          long zcFundId,
					          long zcCpnBasisId,
					          long zcFundBasisId,
					          double fxSpot,
					          CCString& inModelType,
                              double betaCorrel,
                              double hump,
                              long   SABRSigmaOrAlpha,
                              ARM_result& result,
					          long objId = -1);

extern long ARMLOCAL_INITCSB(long csbId,
							 long zcCpnId,
							 long swVolId,
							 long capVolId,
							 long rhoCapId,
							 long nuCapId,
							 long betaCapId,
							 long rhoSwoptId,
							 long nuSwoptId,
							 long betaSwoptId,
							 double hump,
							 double betaCorrel,
							 double reCorrel,
							 long inSABRSigmaOrAlpha,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_INITBERMUDASWAPTION (long bsId,
										  long mktDataManagerId,
										  vector<string> mdmKeys,
										  vector<int>* controlVariates,
										  vector<double>* controlPrices,
										  vector<string> modelParams,
										  bool mrsCalibFlag,
										  bool betaCalibFlag,
										  int numMethodType,
										  int amcIter,
										  int mcIter,
										  int maxBucketSize,
										  string genType,
										  int treeSteps,
										  vector<int> portfolioMode,
										  bool boundaryFlag,
										  bool approxMarginFlag,
										  bool freezeBetasFlag,
										  int modelType,
										  bool calculateProbaFlag,
										  ARM_result& result,
										  long objId = -1);

extern long ARMLOCAL_INITSWAPTIONBERMUDA (  long bsId,
											vector<int>* controlVariates,
											vector<double>* controlPrices,
											vector<string> modelParams,
											bool mrsCalibFlag,
											bool atmCalibFlag,
											int numMethodType,
											int amcIter,
											int mcIter,
											int maxBucketSize,
											string genType,
											int treeSteps,
											vector<int> portfolioMode,
											bool boundaryFlag,
											bool approxMarginFlag,
											bool freezeBetasFlag,
											int modelType,
											bool calculateProbaFlag,
											long zcCpnId,
											long swoptVcId,
											long capVcId,
											long capRoId,
											long capNuId,
											long capBetaId,
											long swoptRoId,
											long swoptNuId,
											long swoptBetaId,
											long normalModelId,
											long inSABRSigmaOrAlpha,
											ARM_result& result,
											long objId = -1 );

extern long ARMLOCAL_INITCAPTION (long captionId,
								  long mktDataManagerId,
								  vector<string> mdmKeys,
								  int nbFactors,
								  int SFRMVolType,
								  string SwoptCalibrationMode,
								  string BetaCalibrationMode,
								  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_INITCSO (long csoId,
							  long zcId,
							  long capVolId,
							  long swoptVolId,
							  long inSABRSigmaOrAlpha,
							  long rhoCapId,
							  long nuCapId,
							  long betaCapId,
							  long rhoSwoptId,
							  long nuSwoptId,
							  long betaSwoptId,
							  long flatVolId,
							  long convAdjustVolId,
							  long convAdjustManagerId,
							  long correlDiagCapId,
							  long correlDiagSwoptId,
							  long correlCorrId,
							  long mrsId,
							  long correlId,
							  long volRatioId,
							  long mrsSpreadId,
							  vector<double> modelParams,
							  vector<string> calibParams,
							  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
							  long forexId,
							  long fundZcId,
							  long domBasisZcId,
							  long fundBasisZcId,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_INITCSO (long csoId,
							  long mktDataManagerId, 
							  // vector<string> mdmKeys, : Now irrelevant
							  vector<double> modelParams,
							  vector<string> calibParams,
							  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_ARM_GetFixing (const CCString& source,
									const CCString& index,
									const CCString& tenor,
									const CCString& ccy,
									double asof,
									ARM_result& result);

extern long ARMLOCAL_SetSecurityData(long secId,
							         const VECTOR<double>& data,
					                 ARM_result& result,
					                 long objId);

extern long ARMLOCAL_GetSecurityData(long secId,
							         const CCString& data,
					                 ARM_result& result);

extern long ARMLOCAL_GetLastDateWarm(ARM_result& result);

extern long ARMLOCAL_SecurityFlows(	vector<CCString>& C_labels,
									vector<double>& C_values,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_INITCRASPREAD(long ccsoId,
								  long zcId,
								  long capVolId,
								  long rhoCapId,
								  long nuCapId,
								  long betaCapId,
								  long swoptVolId,
								  long rhoSwoptId,
								  long nuSwoptId,
								  long betaSwoptId,
								  long sigmaOrAlpha,
								  long convAdjustVolCapId,
								  long convAdjustVolSwoptId,
								  long convAdjustType,
								  long correlCorrId,
								  long correlDiagCapId,
								  long correlDiagSwoptId,
								  double mrs,
								  double volRatio,
								  double mrsSpread,
								  double correl,
								  vector<string> modelParams,
								  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
								  const vector<string>& localCalibFlags,
								  ARM_result& result,
								  long objId = -1);


extern long ARMLOCAL_INITGLOBALCAP(long calculatorId,
							       long zcId,
							       long capVolId,
							       long rhoId,
							       long nuId,
							       long betaId,
							       double mrs,
							       ARM_Vector* calibParams,
							       int nbSteps,
							       vector<string> randGenerator,
							       int samplerType,
							       std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ C_productsToPrice,
							       ARM_result& C_result,
							       long objId = -1);

extern long ARMLOCAL_INITTARNFX(long calculatorId,
								vector<CCString> MCParams,
								int factorNb,
								int timeStepNb,
								int spaceStepNb,
								double stdDevNb,
								int skipPDE,
								int rescalling,
								string modelType,
								int smileFlag,
								int mixCalib,
								int oneFactorFlag,
								string correlType,
								vector<CCString> zeroCurvesId,
								vector<CCString> basisCurvesId,
								vector<CCString> forexId,
								vector<CCString> ATMVolId, //for swopt BSGen
								vector<CCString> fxVolId, //for BS fx models
								vector<CCString> mixtureParamsId, //for mixture fx models
								vector<CCString> mrsParamsId,
								vector<CCString> QParamsId,
								long correlMatrixId,
								ARM_result& result,
								long objId = -1);

/*--- For setting Summit keys to a curve (ZC or Vol) object ---*/
/*
void SetMktExternalCharacteristics(ARM_AbstractMarketClass* mktObj,
                                   const CCString& index,
							       const CCString& currency,
							       const CCString& cvName,
                                   const string& type);


*/

extern long ARMLOCAL_ViewCell(
	const long&				ObjectId,
    ARM_result&				result, 
	long					objId = -1);

extern long ARMLOCAL_MatrixVectorViewer(
	const long&			ObjectId,
	long&				nbRows,
    long&				nbCols,
	vector<double>&     values,
	ARM_result&			result);

extern long ARMLOCAL_GetFixingsFromInstrument(const CCString& idSummit,
											  ARM_result& result,
											  long objId = -1);

#endif	// ARM_LOCAL_GLOB_H

/*------------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
