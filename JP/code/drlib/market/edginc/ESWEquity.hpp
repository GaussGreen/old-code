//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ESWEquity.hpp
//
//   Description   Equity stream in ESW
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ESW_EQUITY_HPP
#define EDG_ESW_EQUITY_HPP

#include "edginc/FXAsset.hpp"
#include "edginc/ESWAverage.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/OutputRequestUtil.hpp"


DRLIB_BEGIN_NAMESPACE

///////////////////////////////
class MARKET_DLL ESWEquity : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultESWEquity(){
        return new ESWEquity();
    }

    enum accrualType {NONE=0,
                      LIBOR_WITHOUT_SPREAD=1,
                      LIBOR_WITH_SPREAD=2,
                      FIXED=3};

    void populateEq();

    /** needed to run example */
    ESWEquity();

    string              ccyTreatment;
    CAssetSP                asset; // $unregistered

    /** return an cashflow array consisting of equity paydates (includes avg out paydates) */
    CashFlowArraySP populateCF(CashFlowArraySP& notional, 
                               bool             calNotional, 
                               const string&   ccyISOCode,
                               OutputRequestUtil::KnownCashFlows& knownCFs, 
                               DateTimeArray& payment_dates);


    void init(CAssetWrapper assetWrapper,YieldCurveWrapper yc,const DateTime& valueDate,
              const DateTime& endDate,bool isCallable,bool isFixedNotional, 
              DayCountConventionSP payDCC,
              string swapTypeIn);

    void curtailEq(const DateTime& endDate);
    void curtailAvOut(const DateTime& endDate);
        
    void Validate(DateTimeArray &liborAccrueDates, bool useExplicitLinking, bool isFloating,
                  const DateTime& callDate,bool isCallable, const DateTime& valueDate);

    double priceEq(YieldCurveWrapper discount, const Control* control, CResults* results,
        OutputRequestUtil::KnownCashFlows& knownCFs, DateTimeArray& payment_dates);

        
    /** calculate num of shares for div leg for a given date list 
        taking into account of contract size change (with new fixing levels) */
    void calcNumOfShares(const DateTimeArray& exDates, DoubleArray& numOfShares) const;
        
    /** calculate notional size for libor leg for a given date list 
        taking into account of contract size change (with new fixing levels) */
    void calcNotional(const DateTimeArray& liborAccrueDates, /* (I) */
                      const CIntArray& liborPeriod, /* (I) */
                      bool explicitLink, /* (I) */
                      const CBoolArray      isCompounding, /* (I) */
                      DoubleArray& notional /* (M) */) const;

    void setFixing(const DateTime& valDate, 
                   const DateTime& rollDate, 
                   CAssetWrapper assetWrapper,
				   const Theta* shift);
        
    DateTimeArray       payDates;

    
	string          LegName;    
    int             accrueEqType;
    bool            isReceiving;
    CIntArray       periodIndex;
    DateTimeArray   eqRefixDates;
        
    DoubleArray             size;
    DoubleArray             eqStartLevel;
    DoubleArray             eqEndLevel;
    DoubleArray             eqAccrueDF;
    DoubleArray             fxFixing;
    ESWAvInSP               avIn;
    ESWAvOutSP              avOut;

    //for Brazil Swap
    double spreadEq;



    // unregistered
    DateTime            valueDate; // get it from instrument $unregistered
    bool                isFixedNotional; // get it from instrument $unregistered
    DateTime            EqCallSettleDate; // $unregistered
    DoubleArray         eqGrowthFactor; // $unregistered
    DayCountConventionSP payDCC; // $unregistered


    void getInLevelAndNumShares(int iPayment, double& inLevel, double& numShares, 
                                const DateTime& refDate = DateTime(0,0), bool soFar = false) const;
    int getOutLevel(int iPayment, double numShares) const;
    int getDivEqIndex(const DateTime &exDate) const;
    int getLiborEqIndex(const DateTime &liborAccrueDate, int liborPeriod, bool explicitLink) const;
    static int NumHistoricalDates(const DateTimeArray &dates, const DateTime &currentDate, bool inclusive = true);
	//double getFracSharesInAvgPeriod(const DateTime &exDate) const;
	double getAvInNotional(const DateTime& valueDate, int iAveSample) const;	
    bool isCcyStruck() const;
    static void removeDuplicates(CashFlowArray & cf);
    void validateHistoricalSamples(const DateTime& valueDate);

private:
    // to help in KNOWN_CASHFLOWS computation
    vector<bool> eqStartIsKnown; // $unregistered
    vector<bool> eqEndIsKnown; // $unregistered

    string swapType;   //STANDARD: standard swap, BRZ: Brazil swap $unregistered

    //spreadEq = 0, return 1
    //spreadEq <> 0, return (1+ spreadEq)^yearFrac
    //used for Brazil Equity Swap
    double getScalingFactorFromSpread(int iRefixDate) const;
    void validateBrazilSwap() const;

};

typedef smartPtr<ESWEquity> ESWEquitySP;
typedef smartConstPtr<ESWEquity> ESWEquityConstSP;

DRLIB_END_NAMESPACE
#endif
