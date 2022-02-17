//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : ModelFilter.hpp
//
//   Contains    : ModelFilter, PriceCounter, ExposureHighlighter,
//                 IRVegaPointwiseExposureHighlighter
//
//   Description : ModelFilter is an abstract model that doesn't actually price.
//                 How useful is that?  Handy for when you want to find out things
//                 like how may pricings are there (PriceCounter) or what market
//                 objects am I sensitive to (ExposureHighlighter and
//                 IRVegaPointwiseExposureHighlighter) without a full-blown pricing run.
//----------------------------------------------------------------------------
#ifndef MODEL_FILTER_H
#define MODEL_FILTER_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/Model.hpp"
#include "edginc/IRVegaPointwise.hpp"
#include "edginc/LastSensDate.hpp"

DRLIB_BEGIN_NAMESPACE

class SensControl;

class RISKMGR_DLL ModelFilter: public CModel,
                           virtual public LastProductSensDate {
public:
    static CClassConstSP const TYPE;

    // make me one - wraps a genuine model to drive market data selection etc.
    ModelFilter(IModel* model);

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results) = 0;

    // these all tediously delegate down to the "real" model
    virtual CResults* Run(CInstrument* instrument, 
                          CControl*    control);

    virtual CResultsArraySP RunMulti(IInstrumentCollectionSP instruments, 
                                     CControl* control);

    /** calculate prices and store results */
    virtual void PriceMulti(IInstrumentCollectionSP instruments, 
                            CControl* control, 
                            CResultsArraySP results);
 
    virtual MarketObjectSP getMarketObject(
        const MarketData* market,
        const string&     name,
        CClassConstSP     type,
        const string&     domesticYCName) const;

    virtual void setDomesticYCName (string discountYieldCurveName) const;

    virtual void getInstrumentAndModelMarket (const MarketData*  market,
                                              CInstrument* inst);

    virtual void getInstrumentsAndModelMarket (MarketDataConstSP market,
                                               IInstrumentCollectionSP insts);

    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const;

    virtual void getComponentMarketData(const MarketData*    market,
                                        MarketObjectSP       mo) const;

    virtual void getComponentMarketData(const MarketData* market,
                                        MarketObjectSP    mo,
                                        const string&     domesticYCName) const;

    virtual CorrelationBaseSP getCorrelation(const string&     corrName,
                                             CClassConstSP     type1,
                                             CClassConstSP     type2,
                                             CClassConstSP     corrType,
                                             const MarketData* market) const;
 
    virtual MarketObjectSP modifyMarketData(
        const MarketData*     market,
        const CClassConstSP&  clazz,     // what type was originally requested
        const MarketObjectSP& mo) const; /* what GetMarket returned or what was
                                          * "inline" already */

    virtual void getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments);

    virtual const string& getDomesticYCName() const;

    virtual DateTime getValueDate() const;

    virtual SensControl* AlterControl(
        const SensControl* currSensControl) const;

    virtual void flush();

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    // override clone so that count is accurate
    virtual IObject* clone() const;

    // when to stop tweaking (for product provided end date)
    virtual DateTime endDate(const CInstrument* instrument,
                             const Sensitivity* sensControl) const;

protected:
    ModelFilter(CClassConstSP clazz);
    ModelFilter(IModel* model, CClassConstSP clazz);
    ModelFilter(const ModelFilter &rhs, CClassConstSP clazz);
    
    IModelSP realModel;

private:
    friend class ModelFilterHelper;

    ModelFilter& operator=(const ModelFilter& rhs);

    // For use by clone method (see) to answer a shallow copy of the appropriate
    // derived type when invoked through a base class (i.e. ModelFilter) pointer
    // or reference.
    virtual ModelFilter* shallowCopy() const = 0;
};

class RISKMGR_DLL PriceCounter : public ModelFilter {
public:
    static CClassConstSP const TYPE;

    PriceCounter(IModel* model);

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    // how often was this model called for pricing?
    int pricings() const;

protected:
    PriceCounter(CClassConstSP clazz);
    PriceCounter(IModel* model, CClassConstSP clazz);
    PriceCounter(const PriceCounter &rhs, CClassConstSP clazz);

private:
    friend class PriceCounterHelper;

    virtual PriceCounter* shallowCopy() const;

    int count;
};

class RISKMGR_DLL ExposureHighlighter: public ModelFilter {
public:
    static CClassConstSP const TYPE;

    // make me one - wraps a genuine model to drive market data selection etc.
    ExposureHighlighter(IModel* model);

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    virtual CResultsArraySP RunMulti(IInstrumentCollectionSP instruments, 
        CControl* control);

protected:
    ExposureHighlighter(CClassConstSP clazz);
    ExposureHighlighter(IModel* model, CClassConstSP clazz);
    ExposureHighlighter(const ExposureHighlighter &rhs, CClassConstSP clazz);

private:
    friend class ExposureHighlighterHelper;

    virtual ExposureHighlighter* shallowCopy() const;

    double dummyPrice;
};

// A custom ExposureHighlighter that implements IRVegaPointwise::ISensitivePoints
// Thus, it is used to report IRVol exposure (for models that use IRVol)
class RISKMGR_DLL IRVegaPointwiseExposureHighlighter: public ExposureHighlighter,
    public virtual IRVegaPointwise::ISensitivePoints {
public:
    static CClassConstSP const TYPE;

    IRVegaPointwiseExposureHighlighter(IModel* model);

    /** Essentially relies on instrument implementing ISensitiveIRVolPoints.
    If not returns null. */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const;

protected:
    IRVegaPointwiseExposureHighlighter(CClassConstSP clazz);
    IRVegaPointwiseExposureHighlighter(IModel* model, CClassConstSP clazz);
    IRVegaPointwiseExposureHighlighter(const IRVegaPointwiseExposureHighlighter &rhs,
        CClassConstSP clazz);

private:
    friend class IRVegaPointwiseExposureHighlighterHelper;

    virtual IRVegaPointwiseExposureHighlighter* shallowCopy() const;
};

DRLIB_END_NAMESPACE
#endif
