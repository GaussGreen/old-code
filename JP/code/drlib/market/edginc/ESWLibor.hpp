//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ESWLibor.hpp
//
//   Description   libor leg for equity swap.
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ESW_LIBOR_HPP
#define EDG_ESW_LIBOR_HPP

#include "edginc/DayCountConvention.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/FloatRate.hpp"
#include "edginc/PayStream.hpp"
#include "edginc/ESWEquity.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ESWLibor : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultESWLibor(){
        return new ESWLibor();
    }

//    // flag for ave
//    enum ESWAveEnum {NO_AVE, AVE_IN, AVE_OUT} ESWAveType;
    string ESWAveType;

    ESWLibor();
    /** validate Libor leg */
    virtual void Validate(ESWEquityConstSP eq,const DateTime& valDate, 
                          const DateTime& startDate, const DateTime& endDate);

    /** initialisation */
    void init(ESWEquitySP eq, const DateTime& valDate,
              const DateTime& endDate, 
              const DateTime& eqSettleCallDate,
              bool isCallable, 
              YieldCurveSP discount,
              string swapTypeIn);


    void curtailLibor(const DateTime& endDate,const DateTime& eqSettleCallDate);

    /** calculation of accrual factors from settle to payment date for equity (avOut)
        leg. The factors are saved in eq->eqAccrueDF for regular EAP, and in 
        avOut->discountFactors for avOut EAP. */
    void eqAccrualFactors(ESWEquitySP eq, const DateTime& valueDate, const YieldCurveSP &discCurve);

    double calcStubFactor(const DateTime& startDate,
                                      const DateTime& endDate,
                                      bool  addSpread) const;

    /** calculation of accrual factor for two specified dates. */
    double calcAccrueFactor(const DateTime& valueDate,
                                  const DateTime& accruStartDate,
                                  const DateTime& accrueEndDate,
                                  const DoubleArray& rates, // this is libor rates (possibly with spread)
                                  const double histDF,
                                  bool  addSpread) const;

    // stairs or linear interp spread
    double InterpSpread(const DateTime& aDate, bool isLinear = true) const;

    void populateFloatRate(CDoubleArray& ratesInUse, bool addspreads) const;

    void setFixing(const DateTime& valDate, 
				   const DateTime& rollDate, 
				   const ESWAvInSP& avIn, 
				   const YieldCurveSP& discount);

    /** return a floater definition */
    FloatRateSP getFloater(const string& type, const DateTimeArray &fixingDates,
                                 const string& frontType, const string& backType) const;
	
    /** return a floater with isMulti set to false */
    FloatRateSP	getFloater(const string& type) const;

    /** return a payDCC */
    DayCountConventionSP&	getpayDCC();
    DayCountConventionSP&	getrateDCC();
    //HolidayConstSP&         getHolidays();

    HolidayConstSP getHolidays() const;

    /** price libor leg */
    double priceIR(Control* control, CResults* results, OutputRequestUtil::KnownCashFlows& knownCFs, 
                               DateTimeArray& payment_dates);

	/** find equity accrue start date that implicitly corresponds to accrue date. */
	// I.e. find the accrue start A satisfying eq acc start <= lib acc < eq acc end
	int findImplicitMatchingIndex(ESWEquityConstSP & eq, int accIdxMatch) const;
        
    string                  LegName;
    bool                    isReceiving;
    HolidayWrapper          hol;
    bool                    isFloating;
    DayCountConventionSP    payDCC;
    DayCountConventionSP    rateDCC;
    string                  rateType;
    bool                    compoundFlat;
    bool                    useExplicitLinking;
    string                  frontStubType;
    string                  backStubType;
    bool                    needSpread;

    CIntArray       periodIndex;
    DateTimeArray   refixDates;
    DateTimeArray   accrueDates;
    DateTimeArray   payDates;
    CDoubleArray        spreads;
    CDoubleArray        liborFixings;
    CBoolArray      isCompounding;
    YieldCurveWrapper       couponCurveWrapper;

    //only used for fraction of floator for Brazil swap, optional input
    CDoubleArray        weights;
    string swapType;   //STANDARD: standard swap, BRZ: Brazil swap $unregistered

   
    BadDayConventionConstSP  badDayConvention;
    
    // not instrument data
    DoubleArray     Notionals; // $unregistered

    void GetMarket(const IModel* model, const CMarketDataSP market);    
	void curtailBack(int Idx,
					  DateTimeArray &refixDates,
					  DateTimeArray &payDates,
					  DateTimeArray &accrueDates,
					  DoubleArray   &RefixLevels,
					  DoubleArray   &spreads,
					  DoubleArray   &Notionals);


    void computeGenAccrualFactors(int accrueType, 
                                    const DateTimeArray & xAccrueStartDates,
                                    const DateTimeArray & xAccrueEndDates,
                                    const DoubleArray & xAccrueDiscFactors,
                                    const DateTime& valueDate, 
                                    const CAssetSP & asset,
                                    YieldCurveSP discCurve,
                                    DoubleArray &growthFactor);

	void validateHistoricalSamples(const DateTime& valDate) const;

	void setCouponCurve(const YieldCurveSP& discount);


private:
    // unregistered
    YieldCurveSP    couponCurve; // $unregistered
    DateTime        valueDate; // $unregistered

    // array to handle ave-out, one stream for each ave-out 
    PayStreamArray   LiborPayArr; // $unregistered
    
    // array to handle each av in before consolidation
    PayStreamArray   AvInPayArr; // $unregistered

    // floater definition

    FloatRateSP RegularFloater; // $unregistered

    const IModel* Model; // $unregistered
    CMarketDataSP Market; // $unregistered

    // trim libor array
    void CreateAvOutLiborStream(int i, const DateTimeArray& avDates, 
                          const ESWEquitySP& eq, const DateTime& valueDate,               
                          DateTimeArray& irefixDates, DateTimeArray& iaccrueDates, DateTimeArray& ipayDates,
                          DoubleArray& iNotionals, DoubleArray& iRefixLevels,
                          DoubleArray& iSpreads, DoubleArray& iWeights, BoolArray& iComp);

	bool findMatchingIndex(ESWEquityConstSP & eq,int periodIndexToMatch,int & matchedPdIdx) const;
	
	/** returns the index of the first false element of the compounding array */
	int	firstNoCompound(const BoolArray & compnd) const;
    
	/** removes first numToRemove elements of all libor data.  This is used for averaging in, for both consolidated and non consolidated case. */
    void TrimLibor(const ESWEquitySP&	eq,
                    DoubleArray&		Notionals,
                    DateTimeArray&		refixDates,
                    DateTimeArray&		accrueDates,
                    DateTimeArray&		payDates,
                    DoubleArray&		liborFixings,
                    DoubleArray&		spreads,
                    DoubleArray&		weights,
                    BoolArray&			isCompounding,
                    int					numToRemove) const;
    
    int getNCLibAccAfterIn(const DateTime& lastAvInDate) const;

    double tabulatePayArray(
                                  const PayStreamArray&    payStreamArray,
                                  Control*                 control, 
                                  CResults*                results,
                                  CashFlowArraySP&         irCF) const;

    /** loops through each pay stream and appends its knownCF'S and payment dates */
    void addKnownCFsAndPaymtDts(const PayStreamArray&    payStreamArray,
                                      OutputRequestUtil::KnownCashFlows& knownCFs, 
                                      DateTimeArray& payment_dates) const;

    double computeAccruedInterest(const PayStreamArray&    payStreamArray) const;
	int getLibAccAfterSample(const DateTime& inAccStartDate) const;
	void validateExplicitLinking(ESWEquityConstSP eq) const;

	void validateBrazilSwap() const;

};


typedef smartPtr<ESWLibor> ESWLiborSP;

typedef smartConstPtr<ESWLibor> ESWLiborConstSP;

DRLIB_END_NAMESPACE
#endif
