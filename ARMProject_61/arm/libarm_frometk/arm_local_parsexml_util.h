#ifndef ARM_LOCAL_PARSEXML_UTIL_H
#define ARM_LOCAL_PARSEXML_UTIL_H

//#include <ARM\libarm_local\firstToBeIncluded.h>
//#include <libCCTools++/CCString.h>
class CCString;

class ARM_ReferenceValue;
class ARM_ExerciseStyle;
class ARM_Currency;
class ARM_Vector;
class ARM_Date;
class ARM_SwapLeg;

//// Pour FromRefValueToARM_Curve
#include "gpbase/curvetypedef.h"
#include "gpbase/curve.h"
#include "XMLTools.h"

using ARM::ARM_Curve;
using ARM::ARM_Interpolator;
using ARM::ARM_GP_Vector;

#ifndef XML_DEFINE
#define XML_DEFINE
#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;
#endif


ARM_Date GetDateFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName);
double GetDoubleFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName, bool throwExcept = false);
int GetIntFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName);
CCString GetCCStringFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName);

// Accesseur au Niveau du noeud l'ENV
CCString GetCustumerId(MSXML2::IXMLDOMNode* node);
CCString GetDealId(MSXML2::IXMLDOMNode* node);
CCString GetBook(MSXML2::IXMLDOMNode* node);

// Accesseur au Niveau du noeud OPTION
int GetPorC(MSXML2::IXMLDOMNode* node);
ARM_Currency* GetOptionCcy(MSXML2::IXMLDOMNode* node);


ARM_ReferenceValue*  GetCRASpread(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetCRABarrierDown(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetCRABarrierUp(MSXML2::IXMLDOMNode* node);
int GetCRAFirstPeriodRank(MSXML2::IXMLDOMNode* node);

///// BLOB FUNCTIONS /////
// Generic
int		 GetBLOBNum		(MSXML2::IXMLDOMNode* node, int n=1); // n : index of Num field (ex : n=3 <=> Num3)
int		 GetBLOBTiming	(MSXML2::IXMLDOMNode* node, int n=1); // n : index of timing field (ex : n=3 <=> Timing3)
int		 GetBLOBBasis	(MSXML2::IXMLDOMNode* node, int n=1); // n : index of Basis field (ex : n=3 <=> Basis3)
int		 GetBLOBResetGap(MSXML2::IXMLDOMNode* node, int n=1); // n : index of Reset Gap field (ex : n=2 <=> ResetGap2)
int		 GetBLOBIdx		(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int n=1); // n : index of Index field (ex : n=3 <=> Index3)
double	 GetBLOBAmount	(MSXML2::IXMLDOMNode* node, int n=1); // n : index of Amount field (ex : n=3 <=> Amount3)
double	 GetBLOBRate	(MSXML2::IXMLDOMNode* node, int n=1); // n : index of Rate field (ex : n=3 <=> Rate3)
char*	 GetBLOBIdxTerm	(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int n=1); // n : index of Index field (ex : n=3 <=> Index3)
char*	 GetBLOBResetCal(MSXML2::IXMLDOMNode* node, int n=1); // n : index of Reset Calendar field (ex : n=2 <=> ResetCal2)
char*	 GetBLOBString  (MSXML2::IXMLDOMNode* node, int n=1); // n : index of String field (ex : n=2 <=> String2)
ARM_Date GetBLOBDate	(MSXML2::IXMLDOMNode* node, int n=1); // n : index of Date field (ex : n=3 <=> Date3)
ARM_ReferenceValue* GetBLOBSchedAmount(MSXML2::IXMLDOMNode* node, int schedIdx=1, int amountIdx=1); // default: field Amount1 from Schedule1
ARM_ReferenceValue* GetBLOBSchedRate(MSXML2::IXMLDOMNode* node,const char* dateName, int schedIdx=1, int rateIdx=1); // default: field Rate1 from Schedule1
ARM_ReferenceValue* GetBLOBSchedRateMulti(MSXML2::IXMLDOMNode* node,const char* dateName, int schedIdx=1, int rateIdx=1,double factorMulti=1); // default: field Rate1 from Schedule1
ARM_ReferenceValue* GetBLOBSchedRateMultiWithLastDate(MSXML2::IXMLDOMNode* node,const char* dateName,const ARM_Date& date, int schedIdx=1, int rateIdx=1, double factorMulti=1); // default: field Rate1 from Schedule1



ARM_Date GetOptionExpiry(MSXML2::IXMLDOMNode* node);
ARM_Date GetOptionStartExpiry(MSXML2::IXMLDOMNode* node);
ARM_Vector* GetFees(MSXML2::IXMLDOMNode* node);
ARM_Vector* GetBermudDates(MSXML2::IXMLDOMNode* node);
ARM_Vector* GetExerciseDates(MSXML2::IXMLDOMNode* node);
int GetExerciceStyleType(MSXML2::IXMLDOMNode* node);
ARM_ExerciseStyle * GetExerciceStyle(MSXML2::IXMLDOMNode* node);
double GetStrike(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetStrikes(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetExerFees(MSXML2::IXMLDOMNode* node);

// soecific to CRF because it has its own input structure
int GetCRFDayCount(MSXML2::IXMLDOMNode* node);
int GetCRFIndexDayCount(MSXML2::IXMLDOMNode* node);
int GetCRFTiming(MSXML2::IXMLDOMNode* node);
int GetCRFResetGap(MSXML2::IXMLDOMNode* node);
void GetCRFIndexTerm(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, char* indexTerm);
void GetCRFResetCal(MSXML2::IXMLDOMNode* node, char* resetCal);
ARM_Date GetCRFFixEndDate(MSXML2::IXMLDOMNode* node);

ARM_ReferenceValue*  GetCRFLeverage(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetCRFOneLeverage(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetCRFCpnMin(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetCRFOneCpnMin(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetCRFCpnMax(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetCRFOneCpnMax(MSXML2::IXMLDOMNode* node);

// Accesseur au Niveau du noeud ASSET
int GetStubRule(MSXML2::IXMLDOMNode* node);
int GetAssetDayCount(MSXML2::IXMLDOMNode* node);
int GetInterestDayCount(MSXML2::IXMLDOMNode* node);
ARM_Date GetStartDate(MSXML2::IXMLDOMNode* node);
ARM_Date GetEndDate(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetNotionalConst(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetNotional(MSXML2::IXMLDOMNode* node,bool isKernel = true);
ARM_ReferenceValue*  GetNotionalCust(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue*  GetNotionalWithPayDates(MSXML2::IXMLDOMNode* node, bool isCust);
ARM_ReferenceValue*  GetFundNotionalWithPayDates(MSXML2::IXMLDOMNode* node);

//Noeud ASSET: Tarn calculator
string GetTarnNodeName(MSXML2::IXMLDOMNode* node, string nodeName);
ARM_GP_Vector* GetTarnCustomDates(MSXML2::IXMLDOMNode* node);
ARM_GP_Vector* GetPastFixings(MSXML2::IXMLDOMNode* node, ARM_Date asOfDate);
ARM_GP_Vector* GetPastStart(MSXML2::IXMLDOMNode* node, ARM_Date asOfDate);

int GetCompoundingFreq(IXMLDOMNode* node);
ARM_Currency* GetCcy(MSXML2::IXMLDOMNode* node);
int GetYieldFreq(MSXML2::IXMLDOMNode* node);
int GetResetFreq(MSXML2::IXMLDOMNode* node);
int GetPayFreq(MSXML2::IXMLDOMNode* node);
int GetResetGap(MSXML2::IXMLDOMNode* node);
int GetPayResetGap(MSXML2::IXMLDOMNode* node);
int GetPayGap(MSXML2::IXMLDOMNode* node);
int GetResetTiming(MSXML2::IXMLDOMNode* node);
int GetRngProductType(MSXML2::IXMLDOMNode* node);
int GetPayTiming(MSXML2::IXMLDOMNode* node);
int GetPorS(MSXML2::IXMLDOMNode* node);
int GetLDPricingMeth(MSXML2::IXMLDOMNode* node);
int GetIndexType(MSXML2::IXMLDOMNode* node);
double GetSpread(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetSpreadVariable(MSXML2::IXMLDOMNode* node);
void GetPayCalendar(MSXML2::IXMLDOMNode* node, char* payCal);
void GetResetCalendar(MSXML2::IXMLDOMNode* node, char* resetCal);
int GetIntRule(MSXML2::IXMLDOMNode* node);
int GetPayIntRule(MSXML2::IXMLDOMNode* node);
int GetDecompFreq(MSXML2::IXMLDOMNode* node);
int GetNotionalExchangeFlag(MSXML2::IXMLDOMNode* node);
double GetMeanRev(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetCRFStrike(MSXML2::IXMLDOMNode* node);
int GetCapFloor(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetCapStrikes(MSXML2::IXMLDOMNode* node);
double GetCapStrike(MSXML2::IXMLDOMNode* node);
void GetResetRollDate(MSXML2::IXMLDOMNode* node, char* resetRollDate);
void GetPayRollDate(MSXML2::IXMLDOMNode* node, char* payRollDate);
int GetDecompFlag(MSXML2::IXMLDOMNode* node);
void GetStubDate1(MSXML2::IXMLDOMNode* node, char* stubdate);
void GetStubDate2(MSXML2::IXMLDOMNode* node, char* stubdate);
ARM_Date GetFstEffDate(MSXML2::IXMLDOMNode* node,ARM_Date defaultDate);
void GetCustom(MSXML2::IXMLDOMNode* node, char* cust);
int GetAmortFreq(MSXML2::IXMLDOMNode* node);
int GetFwdRule(MSXML2::IXMLDOMNode* node);
int GetPayIndexFwdRule(MSXML2::IXMLDOMNode* node);
int GetCompFreq(MSXML2::IXMLDOMNode* node);
int GetCompMode(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetPayFixedRate(MSXML2::IXMLDOMNode* node);
void GetPayIndexResetCal(MSXML2::IXMLDOMNode* node, char* resetCal);

//PL : value set in "pay roll" and "reset roll" on asset main screen
int GetScheduleRule(MSXML2::IXMLDOMNode* node);

// Specific Data
ARM_Date GetStartDateSpecific(MSXML2::IXMLDOMNode* node);
double GetStrikeSpecific(MSXML2::IXMLDOMNode* node);
double GetNotionalSpecific(MSXML2::IXMLDOMNode* node);

// Accesseur sur CONFIG_VECTOR (0 for the first...)
bool GetConfigVectorFlag(MSXML2::IXMLDOMNode* node, int pos);

int ParseCapCashFlows(MSXML2::IXMLDOMNode* node, 
						ARM_Vector * flowStartDates, 
						ARM_Vector * flowEndDates,
						ARM_Vector * paymentDates, 
						ARM_Vector * resetDates, 
						ARM_Vector * intDays, 
						ARM_Vector * strike,
						ARM_Vector * fwd,
						ARM_ReferenceValue * spread, 
						ARM_ReferenceValue * Notional);
ARM_ReferenceValue* GetPrimes(MSXML2::IXMLDOMNode* node);
ARM_Vector * GetFixing(ARM_Vector * fwd, ARM_Vector * fixingDates, ARM_Date date);
void GetFormula(MSXML2::IXMLDOMNode* node, char* formula);
ARM_ReferenceValue* GetRefValueFromFormula(MSXML2::IXMLDOMNode* node, int nthArg);
void GetFormulaInMemorySO(MSXML2::IXMLDOMNode* node, char* formula);
void GetNthArgInFormula(MSXML2::IXMLDOMNode* node, int n, char* Arg, int mode = 0);

int GetFormulaIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int mode = 0);
int GetSpreadIndexType1(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int isCMT = 0);
int GetSpreadIndexType2(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy);

/// In formula "FORM1 ( FORM2 (Arg1) )"
/// returns index type represented by FORM2
double GetSubFormulaIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy);

double GetSpreadPayOff(MSXML2::IXMLDOMNode* node);
double GetSpreadSlopeFlag(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetVariableSpreadPayOff(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetLowBarrier(MSXML2::IXMLDOMNode* node, int resetGap, char* ccy);
ARM_ReferenceValue* GetUpBarrier(MSXML2::IXMLDOMNode* node, int resetGap, char* ccy);
ARM_ReferenceValue* GetBoostedRate(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetCorridorPastFixings(MSXML2::IXMLDOMNode* node, int indexType, const ARM_Date& asOfDate);
ARM_ReferenceValue* GetCorridorStepUpSpread(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetSwaplegStepUpSpread(MSXML2::IXMLDOMNode* node,bool adjustStartDate);
ARM_ReferenceValue* GetStepUpFixCoupon(MSXML2::IXMLDOMNode* node);

int GetCorridorIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int refOrPaid);

int GetSpreadPayoffLibor(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy);
double GetWeight1(MSXML2::IXMLDOMNode* node);
double GetWeight2(MSXML2::IXMLDOMNode* node);

ARM_Vector* GetSpreadFirstFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date);
ARM_Vector* GetSpreadSecondFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date);
ARM_ReferenceValue* GetRefValSpreadFirstFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date);
ARM_ReferenceValue* GetRefValSpreadSecondFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date);

void GetPayIdxPastFixings(IXMLDOMNode* node, const ARM_Date& asOf, int indexType, ARM_ReferenceValue* payIdxPastFixings);
ARM_SwapLeg* GetPayLeg(MSXML2::IXMLDOMNode* node, int indexType, ARM_Currency* discountCcy, int basisDaycount, int payFreq, const ARM_Date& endDateNA, ARM_Date asOf = ARM_DEFAULT_DATE);
ARM_IRIndex* GetPayIndex(MSXML2::IXMLDOMNode* node, int indexType, ARM_Currency* IndexCcy, int basisDaycount, int payFreq);

void GetEcheancier(MSXML2::IXMLDOMNode* node, int resetFReq, int payFreq, ARM_Vector** startDates, ARM_Vector** endDates, ARM_Vector** resetDates, ARM_Vector** paymentDates, ARM_Date asOf = ARM_DEFAULT_DATE);
void GetSpreadCustomFixings(MSXML2::IXMLDOMNode* node, ARM_Vector** startDates1, ARM_Vector** resetDates1, ARM_Vector** fixings1, ARM_Vector** startDates2, ARM_Vector** resetDates2, ARM_Vector** fixings2);

/// For swap legs
void GetSwapCustomSchedule(MSXML2::IXMLDOMNode* node, ARM_Vector** startDates, ARM_Vector** endDates, ARM_Vector** resetDates, ARM_Vector** paymentDates);

// Convertion ARM_ReferenceValue -> ARM_Curve
ARM_Curve* FromRefValueToARM_Curve(const ARM_ReferenceValue* refVal, const ARM_Date& date, ARM_Interpolator<double,double>* interpolator=NULL);

void GetFormulaForCredit(MSXML2::IXMLDOMNode* node, char* formula);

void ParseSpreadFollowTxt(MSXML2::IXMLDOMNode* node, 
						  double& spread,
						  string& ccy,
						  string& Term,
						  string& IndexName,
						  bool& spreadonly);

// For caption data
int GetCaptionIndex(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy);
int GetCaptionNotifDays(MSXML2::IXMLDOMNode* node);
int GetCaptionUnderPorS(MSXML2::IXMLDOMNode* node);
ARM_Date GetCaptionLockout(MSXML2::IXMLDOMNode* node);
ARM_Date GetCaptionFinalExpiry(MSXML2::IXMLDOMNode* node);
ARM_Vector* GetCaptionExerciseDates(MSXML2::IXMLDOMNode* node, 
									ARM_Date AsOf, 
									ARM_INDEX_TYPE indexType,
									int resetTiming,
									int payTiming,
									ARM_Currency* ccy,
									int notifDays,
									char* resetCal,
									char* payCal);
int GetCaptionCpnDayCount(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetCaptionCpnSpread(MSXML2::IXMLDOMNode* node);
int GetCaptionFundDayCount(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetCaptionFundSpread(MSXML2::IXMLDOMNode* node);
ARM_ReferenceValue* GetFxFixings(MSXML2::IXMLDOMNode* node);

// For Global Cap
int GetGlobalCapIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy);

void GetGlobalCapInfo(	MSXML2::IXMLDOMNode* node, 
						ARM_ReferenceValue*& spreads,
						ARM_ReferenceValue*& fixedRates,
						ARM_ReferenceValue*& barriers,
						ARM_SwapLeg* leg,
						double& numerator,
						double& denominator,
						double& nbIter,
						ARM_Date asOf);

// Restrikable

int		GetIfResrtikableType(MSXML2::IXMLDOMNode* node, int& payIndexType, int& refIndexType);
double	GetRestrikableAlpha(MSXML2::IXMLDOMNode* node);
double	GetRestrikableBeta(MSXML2::IXMLDOMNode* node);
double	GetRestrikableRange(MSXML2::IXMLDOMNode* node);
int		GetRestrikablePayIndex(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy);
int		GetRestrikableRefIndex(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy);


// Parsing Calypso XSPSwap 8.3 function

string calypso2ARMCalendar(string& cal);
int GetCalypsoPorS(MSXML2::IXMLDOMNodePtr& node);
void parseFormula(string formula,double& strike,double& leverage,string& index);

// Function for Calypso 8004
ARM_ReferenceValue* createRefValue(vector<double> dates, vector<double> values, int calcMeth, int extrapolMeth, bool checkIfConstant = true, bool useDateEvenIfConstant = false);
	


#endif 
