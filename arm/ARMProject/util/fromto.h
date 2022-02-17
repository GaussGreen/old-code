/*----------------------------------------------------------------------------*/
#ifndef _FROMTO_H
#define _FROMTO_H


#include <ctype.h>

//#include "zerocurv.h"
//#include "volint.h"
#include "refvalue.h"
#include "gpbase/datestrip.h"
#include "currency.h"

class ARM_Swap;
class ARM_SwapLeg;
class ARM_VolCurve;
using ARM::ARM_DateStrip;


#define ARM_DEFAULT_ERR			-990  




/*--- Convert absolute Vega to relative Vega(BS) ---*/
extern double FromAbsoluteVegaToRelativeVega(char* expiry,
                                             char* matu,
                                             ARM_VolCurve* vol,
                                             ARM_ZeroCurve* zc,
                                             double absoluteVega);                                             
extern ARM_Swap* MakeGenSwap(ARM_Date& curDate,
                             double expiry,
                             double matu,
                             ARM::ARM_Currency* ccy,
							 long daycountId = -1);

extern double ApproximatedForwardRate(double Mat,
                                      double Tenor,
                                      ARM_ZeroCurve* zc,
                                      int isSwapRate,
									  long daycountId = -1);

extern void FromRRSTRtoVol(ARM_Matrix* RR, ARM_Matrix* STR, ARM_Vector* pivotVol);


extern double FromStrMatuToDouble(const char* matu, ARM_Date* asOfDate = NULL);

extern double FromRateToRate(double yield, double term,
                             int fromCompMeth, int toCompMeth);

extern double FromRateToRate(double yield, 
                             ARM_Date& startDate, ARM_Date& endDate,
                             int fromRateType, int toRateType);

extern double FromRateToZero(double yield, double term,
                             int yieldCompMeth);

extern double FromZeroToRate(double zero, double term,
                             int yieldCompMeth);

extern int FromFrequencyToLiborType(int freq);

extern int FromFrequencyToPiborType(int freq);

extern int FromFrequencyToEuriborType(int freq);

extern int FromFrequencyToXiborType(int freq, char *ccy);

extern int FromLiborTypeToFrequency(int index);

extern int FromReuterCodeToMonth(char code);

extern void FromFutureCodeToDates(char* Terms, ARM_Date& dateContrat,
                                  ARM_Date& dateUnder);

extern void FromTOYCodeToDates(char* terms, ARM_Date& dateContrat,
                               ARM_Date& dateUnder);

ARM_PRODUCT_TYPE StringToProductType(char* str);
 
ARM_GEN_MATU StringToGenMatu(char* str);
 
double StringMatuToYearTerm(char* str);

/// the opposite function
string YearTermToStringMatu(double d);

string ConvertYearTermToStringMatu(double d);
 
double GenMatuToYearTerm(ARM_GEN_MATU gm);

vector< string> splitStringIntoPieces( const string& s, const string& delimiter );

extern int FromStringToIndexType(const char* term, const char* index);

extern void GetSummitCurveIndexFromCurrency(char* ccy, char* index);
std::string GetSummitCurveIndexFromCurrency(const std::string& ccy); 



extern int FromSummitDaycountToARMDaycount(const char* daycount);

extern int FromSummitFreqToARMFreq(const char* freq);

extern int FromSummitGapToARMGap(const char* dayGap);

extern double ConvertCouponAsSummit(double coupon, int Freq,
                                    int DayCount, int newFreq,
                                    int newDayCount);

extern double CalcFwdFXSpot(ARM_Date& AsofDate,
                            double Spot,
                            ARM_Date& aFwdDate,
                            ARM_ZeroCurve* NumDiscountCurve,
                            ARM_ZeroCurve* UndDiscountCurve);

extern void GetFrqAndNbFromString(const char* period,
                                  int& Nb,
                                  long& freqId);

// functions moved from spreadleg.cpp
extern int IsFixedIndex(ARM_INDEX_TYPE idx);
extern int IsLiborIndex(ARM_INDEX_TYPE idx);
extern int IsCMSIndex(ARM_INDEX_TYPE idx);
extern int IsCMTIndex(ARM_INDEX_TYPE idx);

// functions moved from interglob.h
extern long ARM_ConvPayResetRule (const char* aRule);
extern long ARM_ConvIrType (const char* aIrType);
extern long ARM_ConvIrTypeWithoutException (const char*  aIrType);
extern long ARM_ConvDayCount (const char* aDayCount);
extern long ARM_ConvDecompFrequency (const char* aFrequency);
extern long ARM_ConvIntRule (const char* aRule);
	   long ARM_ConvStartAdjRule(const string& aRule);
extern long ARM_ConvStubRule (const char* aRule);
extern long ARM_NotionalExchange(const char* aRule);
extern long ARM_ConvPricCorridorLD (const char* param);
extern long ARM_ConvWhatIsInterp (const char* aIntMeth);
extern long ARM_ConvFwdRule (const char* aFwdRule);
extern long ARM_ConvCvMethod (const char* aCvMeth);
extern char* ARM_FromIndexTypeToStringIndex(const ARM_INDEX_TYPE aIndexType);
extern long ARM_ConvMonth (const char* aMonth);

// function moved from arm_local_glob.h
extern ARM_Vector* CreateARMVectorFromVECTOR(const std::vector<double>& param, int newSize = -1);

extern char* ARM_STRDUP(char* inStr);

extern ARM_ReferenceValue* ARM_FromStartRefToResetRef(ARM_ReferenceValue* valstartRef,
													  ARM_SwapLeg* refLeg);
extern ARM_ReferenceValue* ARM_FromStartDSToResetDS(ARM_ReferenceValue* valstartRef,
													  const ARM_DateStrip& DS);

// Calypso function
extern long calypso2ARMStubRule (const char* aRule);
extern int calypso2ARMDaycount(const char* daycount);
extern int calypso2ARMFreq(const char* freq);
extern int calypso2ARMPayOrRec(int payRec);
extern int calypso2ARMTiming(const char* calypsoTiming) ;
extern int calypso2ARMInterestRule(string calypsoRule);
extern string calypso2ARMNotionalExchange(bool initExchange, bool finalExchange);
extern int calypso2ARMIndexType(string indexName,string indexTenor); 
extern int calypso2ARMIndexType(string calypsoIndex);
extern int calypso2ARMDateRoll (string dateRollStr);


extern double FromIndexTypeToTerm(ARM_INDEX_TYPE indexType);

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
