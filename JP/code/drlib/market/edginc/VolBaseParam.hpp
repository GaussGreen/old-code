//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBaseParam.hpp
//
//   Description : A parameterised CVolBase which supports VolRequestRaw
//
//   Date        : 17 May 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLBASEPARAM_HPP
#define EDR_VOLBASEPARAM_HPP

#include "edginc/VolBaseParamSurface.hpp"

DRLIB_BEGIN_NAMESPACE

/** An abstract parameterised CVolBase which supports VolRequestRaw. Unclear why
    we are deriving from CVolBaseParamSurface (looks suspect) */
class MARKET_DLL VolBaseParam: public CVolBaseParamSurface,
                    public virtual IVolProcessed {
public:
    static CClassConstSP const TYPE;

    virtual ~VolBaseParam();

    const TimeMetric& getTimeMetric() const;
    virtual TimeMetricConstSP GetTimeMetric() const; // unfortunate - to tidy up

    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const;

    /** Retrieves time metric from parent's vol surface */
    virtual void getMarket(const IModel*     model, 
                           const MarketData* market,
                           const string&     name);

    /** Calls above method with name=getName() */
    virtual void getMarket(const IModel* model, 
                           const MarketData* market);

    IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      asset) const;

    IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

protected:
    VolBaseParam(CClassConstSP clazz);
    VolBaseParam(CClassConstSP clazz, const string& name);
    TimeMetricSP      timeMetric; // transient
private:
    static void load(CClassSP& clazz);
    
};
typedef smartPtr<VolBaseParam> VolBaseParamSP;

DRLIB_END_NAMESPACE
#endif
