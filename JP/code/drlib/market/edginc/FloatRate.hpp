/*****************************************************************************
 *
 *    Group       : Equity Derivatives Research
 *
 *    Description : Defines a floating rate
 *
 *    Author      : Stephen Hope
 *
 *    Date        : 23 May 2001
 *
 *
 *******************************************************************************/

#ifndef _EDG_FLOATRATE_H
#define _EDG_FLOATRATE_H

#include "edginc/DateTime.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/CompoundBasis.hpp" 


DRLIB_BEGIN_NAMESPACE

/** Implementation of structures and functions to
    support the generation of payments based upon
    interest rates and notionals. Compounding is
    supported. */
class MARKET_DLL FloatRate: public CObject {
public:
    static CClassConstSP const TYPE;
    friend class FloatRateHelper;

/** Constructor */
    FloatRate(const YieldCurve* coupCcy,
              const Holiday* holidays,
              const DayCountConvention* rateDayCountConv,
              const BadDayConvention* rateBadDayConv,
              const MaturityPeriod* resetInterval,
              const MaturityPeriod* matInterval,
              const bool isMultiPeriod,
              const MaturityPeriodArray* = NULL,
              const DateTimeArray* = NULL,
              const bool inArrears = false,
              const int inBasis = CompoundBasis::SIMPLE);


    /** overrides default */
    virtual void validatePop2Object(); 

    /** Returns the floatRate Ccy */
    string getCcy()const;

	/** Returns the floatRate YC Name */
	string getYCName() const;

    /** Get the yield curve from the market data cache */
    void getMarketData(const IModel* model, const CMarketDataSP market);

    /** Calculate the money market forward interest rates on the given 
        dates on or after today's date */
    void calculateRates(const DateTime& today,             /* (I) ignore dates before this day */
                        CashFlowArray&  refixList,
                        DateTimeArray&  payDates)const;   /* (M) calc rates on these dates */

    /** Calculate the money market forward interest rate on the given date */ 
    double calculateRate(const DateTime& date)const;

	/** Calculate the money market forward interest rate between two dates */ 
    double fwd(const DateTime& loDate, const DateTime& hiDate, const DayCountConvention* dcc)const;

    /** calculate the maturity date given the start date and an interval */
    DateTime findMatDate(const DateTime& indexStartDate)const;   /* (M) When to start from */

    // this is here to allow ESW libor object access and use this class
    // or alternatively add constructors, get/set methods, make some field public in here
    friend class ESWLibor;

    // return adjust for in Arrears (convexity/delay adjustment).
    bool isinArrears() const;
   
private:
    FloatRate();
    FloatRate(const FloatRate &rhs);
    FloatRate& operator=(const FloatRate& rhs);


    // interpolate refix dates to compute correct MaturityPeriod for theDate.
    // only should be called if isMultiPeriods == TRUE
    void InterpRefixDates(const DateTime &theDate, MaturityPeriodSP&) const;

    // find MatInterval
	MaturityPeriodSP findMatInterval(const DateTime &theDate) const;

    YieldCurveWrapper           coupCcy;
    MaturityPeriodSP            resetInterval;
    MaturityPeriodSP            matInterval;
    HolidayWrapper              holidays;
    BadDayConventionConstSP     rateBadDayConv;
    DayCountConventionConstSP   rateDayCountConv;
    bool			isMultiPeriod;
    MaturityPeriodArraySP	MultiPeriods;		  
    DateTimeArraySP refixDates; // to interpolate to compute period interval

    // is apply the convexity / delay adjustment (default = false)
    bool    inArrears;


    //allow basis other than CompoundBasis::SIMPLE
    int basis;

};

typedef smartConstPtr<FloatRate> FloatRateConstSP;
typedef smartPtr<FloatRate> FloatRateSP;



DRLIB_END_NAMESPACE

#endif


