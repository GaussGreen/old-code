//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ClosedFormCDSBasket.hpp
//
//   Description : Closed form model for pricing a basket of CDSs (plus some 
//                 basis adjustment, if required)
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef CLOSEDFORMCDSBASKET_HPP
#define CLOSEDFORMCDSBASKET_HPP

#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

/** Credit Index Swaps (CISs) can be priced closed form in a few different 
 * ways. One way is pricing the basket of (index basis adjusted) CDSs on the
 * underlying names in the index. */
class PRODUCTS_DLL ClosedFormCDSBasket: public CModel,
                                        virtual public IHasForwardRatePricer {
public:
    static CClassConstSP const TYPE;

    ClosedFormCDSBasket();

    /** Create a MarketDataFetcher which will be used to retrieve market 
     * data etc */
    virtual MarketDataFetcherSP createMDF() const;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct {
    public:
        virtual void price(ClosedFormCDSBasket* model,
                           Control*    control, 
                           CResults*   results) const = 0;
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormCDSBasket* model) const = 0;
    };

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument* instrument, 
                       CControl*    control, 
                       CResults*    results);

    /** Whether to enable RiskMapping when computing sensitivities; returns
        riskMappingIrrelevant.  See IModel::wantsRiskMapping(). */
    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** Let the model determine which is the latest date it shows 
     * sensitivity to */
    DateTime endDate(const Sensitivity* sensControl,
                     const CInstrument* inst,
                     const DateTime&    instEndDate) const;

    //------------------------------
    // IHasForwardRatePricer methods
    //------------------------------

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;
    
private:
    // For reflection
    static void load(CClassSP& clazz);
    static IObject* defaultClosedFormCDSBasket();

    ClosedFormCDSBasket(const ClosedFormCDSBasket &rhs);
    ClosedFormCDSBasket& operator=(const ClosedFormCDSBasket& rhs);

    // Fields
    bool                 calculateIndexBasis; /* If false, use index basis provided. Otherwise
                                               * compute it. Default: true */
    bool                 useIndexBasis;       /* If false, the index basis will not be used, ie,
                                               * produce a "dummy" index basis. Default: true */
    IForwardRatePricerSP forwardRateModel;    /* used for determining fee leg cashflows */
};

DRLIB_END_NAMESPACE
#endif
