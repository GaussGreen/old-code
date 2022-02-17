//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ClosedFormCDSPS.hpp
//
//   Description : Closed form model for taking CDS par spreads as input
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : August 30, 2001
//
//
//----------------------------------------------------------------------------

#ifndef CLOSEDFORMCDSPS_HPP
#define CLOSEDFORMCDSPS_HPP
#include "edginc/Model.hpp"
#include "edginc/IRVegaPointwise.hpp"
#include "edginc/QuantoCDSAlgorithm.hpp"
#include "edginc/IRGridPointCache.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

/** Credit default swaps can be priced closed form in a few different ways. 
    One way is using CDS par spreads as market input. Another is with credit 
    DRs' asset-based closed-form model. The methods would take different 
    market data and so each is getting its own model object. This one is
    for the CDS par spreads.
*/

class PRODUCTS_DLL ClosedFormCDSPS: public CModel,
                       public virtual IHasForwardRatePricer,
                       public virtual QuantoCDSParSpreads::IAlgorithmBuilder,
                       public virtual IRVegaPointwise::ISensitivePoints{
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormCDSPSHelper;

    /** Simple constructor */
    ClosedFormCDSPS();

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    /** overridden to [reference] copy QuantoCDSAlgorithmSP */
    virtual IObject* clone() const;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(ClosedFormCDSPS* model,
                           Control*    control, 
                           CResults*   results) const = 0;
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct{
    public:
        friend class ClosedFormCDSPSHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormCDSPS* model) const = 0;
    };

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    
    /** Returns an instance of IAlgorithm. Typically the quantoCDSParSpreads
        parameter would be ignored but is there in case you want to switch
        the algorithm dependent upon some property of the quanto'd curve */
    virtual QuantoCDSParSpreads::IAlgorithmSP cdsQuantoAlgorithm(
        const QuantoCDSParSpreads* quantoCDSParSpreads) const;

    /** Uses quanto calibration parameters to identify relevant points */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const;

    /** sets debug state to specified value. If true, then any calls to
        cdsQuantoAlgorithm will return an object that will cache debug data.
        This can be retrieved via getDebugInfo() */
    virtual void selectDebugState(bool switchOn);
    /** Returns debug info - may be null eg if the algorithm has not been 
        used */
    virtual IObjectSP getDebugInfo() const;

    /** For quanto, need to adjust instEndDate to take into account IR spot vol
        calibration which looks at coupons of swap going beyond instEndDate */
    DateTime endDate(const Sensitivity* sensControl,
                     const CInstrument* inst,
                     const DateTime&    instEndDate) const;

    /** Whether to enable RiskMapping when computing sensitivities; returns
        riskMappingIrrelevant.  See IModel::wantsRiskMapping(). */
    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    //------------------------------
    // IHasForwardRatePricer methods
    //------------------------------

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;

    //------------------------------
    // CModel methods
    //------------------------------
    /** invoked after instruments has got its market data.  Allow model
     ** to get extra data */
    virtual void getMarket(const MarketData *market,
                           IInstrumentCollectionSP instruments);


protected:
    ClosedFormCDSPS(CClassConstSP clazz);

private:
    ClosedFormCDSPS(const ClosedFormCDSPS &rhs);
    ClosedFormCDSPS& operator=(const ClosedFormCDSPS& rhs);
    void createQuantoAlgorithm(bool debugOn) const;
    //// fields ////
    string quantoCalibrationStyle;  // eg CMS - for quanto
    string quantoCalibrationMaturity; // eg 10Y - for quanto
    IForwardRatePricerSP forwardRateModel;    /* used for determining fee leg cashflows */

    mutable QuantoCDSAlgorithmSP quantoAlgorithm; // transient, not registered $unregistered
    mutable IRGridPointCacheSP   irGridPtsCache;  // transient, not registered $unregistered
};

DRLIB_END_NAMESPACE
#endif
