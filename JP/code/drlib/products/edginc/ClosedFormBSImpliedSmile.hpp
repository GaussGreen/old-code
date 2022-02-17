//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Filename    : ClosedFormBSImplied.hpp
//
//   Description : Closed-form Black-Scholes implied vol smile 
//                 Model
//
//   Author      : Charles Morcom
//
//   Date        : January 23, 2006
//
//
//----------------------------------------------------------------------------

#ifndef QR_CLOSEDFORMBSIMPLIEDSMILE_HPP
#define QR_CLOSEDFORMBSIMPLIEDSMILE_HPP

#include "edginc/Model.hpp"
#include "edginc/Results.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ClosedFormBSImpliedSmile)
FORWARD_DECLARE(CInstrument)
FORWARD_DECLARE(Control)
FORWARD_DECLARE(MarketDataFetcher)

/** Pricing and tweaks are with respect to the BS implied vol smile at fixed strikes. */
  class PRODUCTS_DLL ClosedFormBSImpliedSmile: public CModel,
					       public virtual IHasForwardRatePricer {
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormBSImpliedSmileHelper;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(ClosedFormBSImpliedSmile* model,
                           Control*        control, 
                           CResults*       results) const = 0;
        virtual ~IProduct() {};
     
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class ClosedFormBSImpliedSmileHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormBSImpliedSmile* model) const = 0;
    };

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    virtual MarketDataFetcherSP createMDF() const;

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    ClosedFormBSImpliedSmile();
    ClosedFormBSImpliedSmile(CClassConstSP clazz);
    virtual ~ClosedFormBSImpliedSmile();

    /** Key method providing access to the pricer */
    virtual IForwardRatePricerSP getForwardRatePricer() const;

private:

    IForwardRatePricerSP forwardRateModel;    /* used for determining fee leg cashflows */
    //class MarketDataFetcherBSSmile;
    //ClosedFormBSImpliedSmile(const ClosedFormBSImpliedSmile &rhs);
    //ClosedFormBSImpliedSmile& operator=(const ClosedFormBSImpliedSmile& rhs);
};

DRLIB_END_NAMESPACE
#endif
