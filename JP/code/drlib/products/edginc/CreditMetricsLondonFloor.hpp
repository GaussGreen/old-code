//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditMetricsLondonFloor.hpp
//
//   Description : Credit Metrics with London floor adjustment Algorithm
//
//   Date        : April 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_CREDIT_METRICS_LONDON_FLOOR_HPP
#define EDR_CREDIT_METRICS_LONDON_FLOOR_HPP

#include "edginc/CreditMetricsModel.hpp"

DRLIB_BEGIN_NAMESPACE
class ITrancheLossCalculatorLegacy;
/** 
 * CCM model with London Floor
 * */
class PRODUCTS_DLL CreditMetricsLondonFloor : public CreditMetricsModel {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~CreditMetricsLondonFloor();
    
    /** Overridden to apply a 'London Floor' adjustment */
    virtual ITrancheLossCalculator* createLossCalculator(
        const DateTimeArray&           timeline,    /* (I) */
        CreditTrancheLossConfigConstSP tranche,     /* (I) */
        CounterPartyCreditConstSP      cpty) const; /* (I) */

    // Not sure that the London Floor adjustment makes sense for recovered
    // notional, so do not override createRecoveredNotionalCalculator

    /** Invoked by the containing model after the instrument data is fetched
        ie after CInstrument::GetMarket is invoked. Overrides
        CreditMetricsModel in order to get london floor curve */
    virtual void postFetchMarketData(IModel*            model,
                                    MarketDataConstSP market);

    /** Called immediately after object constructed */
    void validatePop2Object();

protected:
    /** Only build instances of that class using reflection */
    CreditMetricsLondonFloor(const CClassConstSP& clazz);
    
private:
    /** Default constructor */
    static IObject* defaultConstructor();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** For CreditMetricsLondonFloor::IIntoProduct */
    static void loadIntoProduct(CClassSP& clazz);

protected:    
    // ------
    // FIELDS
    // ------

    /** London floor curve */
    ICDSParSpreadsWrapper londonFloor; 
    
    /** London floor, equity strike */
    double londonFloorEqStrike;
    
    /** London floor, senior strike */
    double londonFloorSeniorStrike;
};


DRLIB_END_NAMESPACE

#endif

