//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RangeAccrue.hpp
//
//   Description : 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_RangeACCRUE_HPP
#define EDR_RangeACCRUE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/IAggregate.hpp"
#include "edginc/MonteCarlo.hpp"

DRLIB_BEGIN_NAMESPACE
class RangeAccrue;

class PRODUCTS_DLL RangeAccrueMaker : public CObject {
public:
    static CClassConstSP const TYPE;
    friend class RangeAccrueMakerHelper;
    friend class RangeAccrue;
    
    virtual RangeAccrue* getDiscounted(const DateTime baseDate, const YieldCurve* discount);

    DateTimeArray getMonitorDates();
    
    void validatePop2Object();

    // for reflection
    RangeAccrueMaker(): CObject(TYPE){};
    virtual ~RangeAccrueMaker() {};

private:
    // fields
    bool          isInside;             // Inside or not
    DateTimeArray monitoringDates;      // list of range monitoring date.
    DateTimeArray paymentDates;         // list of payment date.
    DoubleArray   amounts;              // list of amounts
    DoubleArray   lowLevels;            // lower level or range
    DoubleArray   highLevels;           // higher level or range
};

typedef smartPtr<RangeAccrueMaker> RangeAccrueMakerSP;

class PRODUCTS_DLL RangeAccrue {
public:
    RangeAccrue(const DateTime    baseDate, 
                const YieldCurve* discount,
                RangeAccrueMaker* maker); 

    void scalePvAmounts(const double scaling);

    double getValue(const IPathGenerator*  pathGen,
                    const int endStep,
                    IAggregateSP& assetBasketForRange,
                    SimpleDoubleArray& assetCompsForRange);
    
    double getValue(const DoubleArray& levels, const int iStep);

    double getValue(const double level, const int iStep, const bool isPVed = true);

    DateTimeArray getMonitorDates(){
        return maker->monitoringDates;
    };

    //return a copy of payment dates
    DateTimeArray getPayDates(){
        return maker->paymentDates;
    };

    void makeKnownCashFlow(const IPathGenerator*  pathGen,
                    const int endStep,
                    const bool hasFuture,
                    IAggregateSP& assetBasketForRange,
                    SimpleDoubleArray& assetCompsForRange);

    // return a copy of knownCFL
    CashFlowArraySP getKnownCashFlows(){
        CashFlowArraySP cfl(new CashFlowArray(0));
        for (int i=0; i<knownCFL.size();i++)
            cfl->push_back(knownCFL[i]);
        return cfl;
    };

    void setIsMonitorStep(const DateTimeArray& timeStepDate);

private:
    RangeAccrueMaker* maker;
    DoubleArray       pvAmounts;        // PV value of amounts
    IntArray          isMonitorStep;    // is monitor date?  
                                        // The array size is time steps of engine (MC simDates or Tree TimeSteps)
    CashFlowArray     knownCFL;         // knownCFL
    DateTime          valueDate;        //valuationDate of instruments

    double            remainedValue;    // known coupon, but not paid yet.  Need to add to the value.

};
typedef refCountPtr<RangeAccrue> RangeAccrueSP;

DRLIB_END_NAMESPACE

#endif

